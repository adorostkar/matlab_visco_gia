"""Time-stepping solver for the incompressible-limit Maxwell GIA model.

Model (2-D plane strain, flat earth, small perturbations about a hydrostatic
pre-stressed state):

    div(sigma) + rho0 g grad(u_y) = 0,   sigma = -p I + tau
    div(u) + p/kappa = 0
    d tau/dt + tau/t_M = 2 mu d e/dt,    e = deviatoric strain,  t_M = eta/mu

The gravity ("pre-stress advection") term is kept in its weak form
``int rho0 g u_y div(v)`` which contains the restoring buoyancy force of the
deflected free surface.  The Maxwell law is integrated exactly for a strain
rate that is constant over each step:

    tau_{n+1} = a tau_n + 2 mu b (e_{n+1} - e_n),
    a = exp(-dt/t_M),  b = (1 - a) t_M/dt,

so every step is one linear solve with a constant matrix (``mu_eff = b mu``).
Stresses are non-dimensionalised by a reference shear modulus for conditioning.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import fem
from .mesh import Mesh
from .model import SECONDS_PER_YEAR, Config


def _graded(n, length, ratio):
    """n intervals summing to ``length``; each is ``ratio`` times the previous."""
    if n < 1:
        return np.zeros(0)
    w = ratio ** np.arange(n)
    return length * w / w.sum()


def build_mesh(cfg: Config) -> Mesh:
    li = cfg.ice.half_width
    if not 0 < li < cfg.L:
        raise ValueError("ice half_width must lie in (0, L)")
    dx = np.concatenate(
        [_graded(cfg.nx_ice, li, 1.0), _graded(cfg.nx_far, cfg.L - li, cfg.x_ratio)]
    )
    xv = np.concatenate([[0.0], np.cumsum(dx)])
    dy, layer_of = [], []
    # layers are listed from the surface downwards; the mesh runs bottom-up
    for k, layer in enumerate(cfg.layers):
        d = _graded(cfg.ny_layer, layer.thickness, cfg.y_ratio)  # surface-first
        dy.append(d)
        layer_of.append(np.full(cfg.ny_layer, k))
    dy_top_down = np.concatenate(dy)
    yv = -np.concatenate([[0.0], np.cumsum(dy_top_down)])[::-1]
    mesh = Mesh(xv, yv)
    mesh.layer_of_row = np.concatenate(layer_of)[::-1]
    return mesh.refine(cfg.refine) if cfg.refine else mesh


@dataclass
class Result:
    mesh: Mesh
    times: np.ndarray  # s
    surface_x: np.ndarray  # surface node x coordinates (m)
    surface_ux: np.ndarray  # (nt, nsurf) m
    surface_uy: np.ndarray  # (nt, nsurf) m
    u: np.ndarray  # final displacement [ux; uy]
    p: np.ndarray  # final pressure (Pa)
    snapshots: List[dict] = field(default_factory=list)
    # final element-centre stress (Pa): columns sxx, syy, sxy (sigma = -p I + tau)
    stress: np.ndarray = None


def simulate(cfg: Config, mesh: Mesh = None, callback=None, snapshot_times=()) -> Result:
    mesh = mesh or build_mesh(cfg)
    nD, nP = mesh.nD, mesh.nP

    mus = np.array([l.mu for l in cfg.layers])
    kaps = np.array([l.kappa for l in cfg.layers])
    etas = np.array([l.eta for l in cfg.layers])
    row_layer = np.repeat(mesh.layer_of_row, mesh.nx)
    mu = mus[row_layer]
    kappa = kaps[row_layer]
    eta = etas[row_layer]
    mu_ref = mus.max()

    dt = cfg.dt
    with np.errstate(divide="ignore", over="ignore"):
        tM = eta / mu
        r = dt / tM
        a = np.exp(-r)
        b = np.where(r > 1e-12, -np.expm1(-r) / np.where(r > 1e-12, r, 1.0), 1.0)
    rho_g = cfg.rho_earth * cfg.g / mu_ref
    rho_adv = rho_g if (cfg.gravity and cfg.advection) else 0.0

    # boundary conditions: symmetry at x=0, roller at x=L, fixed bottom
    fixed = np.unique(
        np.concatenate(
            [
                mesh.left(),
                mesh.right(),
                mesh.bottom(),
                nD + mesh.bottom(),
            ]
        )
    )
    free = np.setdiff1d(np.arange(2 * nD + nP), fixed)

    def factor(mu_eff):
        S = fem.assemble_system(mesh, mu_eff / mu_ref, kappa / mu_ref, rho_adv)
        if cfg.gravity:
            S = S + fem.winkler_matrix(mesh, rho_g)
        Sf = S[free][:, free].tocsc()
        if cfg.solver == "direct":
            return spla.splu(Sf).solve
        if cfg.solver != "gmres":
            raise ValueError("solver must be 'direct' or 'gmres'")
        # symmetric diagonal equilibration: displacement and pressure rows
        # differ by many orders of magnitude when lengths are in metres
        d = 1.0 / np.sqrt(np.abs(Sf.diagonal()))
        D = sp.diags(d)
        Ss = (D @ Sf @ D).tocsc()
        ilu = spla.spilu(Ss, drop_tol=1e-8, fill_factor=50, diag_pivot_thresh=1.0)
        M = spla.LinearOperator(Ss.shape, ilu.solve)

        def solve_it(b):
            if not b.any():
                return np.zeros_like(b)
            y, info = spla.gmres(Ss, d * b, M=M, rtol=cfg.gmres_tol, restart=100,
                                 maxiter=50)
            if info != 0:
                raise RuntimeError(f"GMRES did not converge (info={info})")
            return d * y

        return solve_it

    lu_el = factor(mu)
    lu_visc = factor(mu * b)

    # surface load at the Gauss points of the top-row elements (units of mu_ref)
    gx = fem.top_gauss_x(mesh)
    if cfg.ice.profile is None:
        shape = (gx < cfg.ice.half_width).astype(float)
    else:
        shape = np.asarray(cfg.ice.profile(gx), float)

    def load_vec(t):
        ty = -cfg.ice.density * cfg.g * cfg.ice.height(t) / mu_ref * shape
        f = np.zeros(2 * nD + nP)
        f[: 2 * nD] = fem.top_traction_rhs(mesh, ty)
        return f

    def solve(lu, rhs):
        x = np.zeros(2 * nD + nP)
        x[free] = lu(rhs[free])
        return x

    # elastic response at t = 0
    x = solve(lu_el, load_vec(0.0))
    u = x[: 2 * nD]
    e = fem.strain_dev(mesh, u)
    tau = 2 * (mu / mu_ref)[:, None, None] * e

    surf = mesh.top()
    nsteps = int(round(cfg.t_end / dt))
    times = np.arange(nsteps + 1) * dt
    sux = np.zeros((nsteps + 1, len(surf)))
    suy = np.zeros_like(sux)
    sux[0], suy[0] = u[surf], u[nD + surf]
    snaps = []
    snap_idx = {int(round(s / dt)) for s in snapshot_times}

    def stress_of(x, tau):
        """Element-centre stress sigma = -p I + tau (Pa)."""
        pe = x[2 * nD :][mesh.conn1].mean(axis=1)
        t = tau.mean(axis=1)
        return mu_ref * np.stack(
            [t[:, 0] - pe, t[:, 1] - pe, t[:, 2]], axis=1
        )

    def emit(k, x, tau):
        if k in snap_idx:
            snaps.append(
                dict(time=times[k], ux=x[:nD].copy(), uy=x[nD : 2 * nD].copy(),
                     p=x[2 * nD :] * mu_ref, stress=stress_of(x, tau))
            )
        if callback:
            callback(k, times[k], x)

    emit(0, x, tau)
    mub = (mu * b / mu_ref)[:, None, None]
    af = a[:, None, None]
    for k in range(1, nsteps + 1):
        hist = 2 * mub * e - af * tau
        rhs = load_vec(times[k])
        rhs[: 2 * nD] += fem.history_rhs(mesh, hist)
        x = solve(lu_visc, rhs)
        u = x[: 2 * nD]
        e_new = fem.strain_dev(mesh, u)
        tau = af * tau + 2 * mub * (e_new - e)
        e = e_new
        sux[k], suy[k] = u[surf], u[nD + surf]
        emit(k, x, tau)

    return Result(
        mesh=mesh,
        times=times,
        surface_x=mesh.node_x[surf],
        surface_ux=sux,
        surface_uy=suy,
        u=u,
        p=x[2 * nD :] * mu_ref,
        snapshots=snaps,
        stress=stress_of(x, tau),
    )
