"""Finite element kernels (vectorised over elements) for Q2-Q1 elements."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .mesh import Mesh


def _shape1d(x, order):
    x = np.asarray(x)
    if order == 2:
        N = np.stack([x * (x - 1) / 2, 1 - x**2, x * (x + 1) / 2], axis=-1)
        dN = np.stack([x - 0.5, -2 * x, x + 0.5], axis=-1)
    else:
        N = np.stack([(1 - x) / 2, (1 + x) / 2], axis=-1)
        dN = np.stack([-0.5 * np.ones_like(x), 0.5 * np.ones_like(x)], axis=-1)
    return N, dN


_gp, _gw = np.polynomial.legendre.leggauss(3)


def _tensor(order):
    """Values and reference derivatives at the 3x3 Gauss points.

    Returns arrays of shape (9 gauss points, nodes): N, dN/dxi, dN/deta.
    Gauss point g = q*3+p -> (xi_p, eta_q); local node a = j*n+i.
    """
    N1, dN1 = _shape1d(_gp, order)
    N = np.einsum("qj,pi->qpji", N1, N1).reshape(9, -1)
    Nxi = np.einsum("qj,pi->qpji", N1, dN1).reshape(9, -1)
    Neta = np.einsum("qj,pi->qpji", dN1, N1).reshape(9, -1)
    return N, Nxi, Neta


N2, N2xi, N2eta = _tensor(2)
N1, _, _ = _tensor(1)
W2 = np.einsum("q,p->qp", _gw, _gw).ravel()  # reference weights, g = q*3+p
# Gauss point coordinates in the reference element, g = q*3+p
GXI = np.tile(_gp, 3)
GETA = np.repeat(_gp, 3)


def _grads(mesh: Mesh):
    sx = (2.0 / mesh.hx)[:, None, None]
    sy = (2.0 / mesh.hy)[:, None, None]
    Gx = N2xi[None] * sx
    Gy = N2eta[None] * sy
    W = W2[None] * (mesh.hx * mesh.hy / 4.0)[:, None]
    return Gx, Gy, W


def _coo(mesh, rows_c, cols_c, vals, shape):
    r = np.broadcast_to(rows_c[:, :, None], vals.shape)
    c = np.broadcast_to(cols_c[:, None, :], vals.shape)
    return sp.coo_matrix((vals.ravel(), (r.ravel(), c.ravel())), shape=shape)


def assemble_system(mesh: Mesh, mu, kappa, rho_g):
    """Saddle-point matrix for unknowns [ux; uy; p] (per-element mu, kappa).

        sum_e int 2 mu (eps_dev(u):eps(v)) + rho g u_y div v - p div v = f.v
        int q (div u + p/kappa) = 0

    ``rho_g`` is the pre-stress advection coefficient (rho0*g, 0 to disable).
    """
    nD, nP = mesh.nD, mesh.nP
    Gx, Gy, W = _grads(mesh)
    mu = np.asarray(mu)[:, None, None]
    Wm = W[:, :, None] * mu[:, :1, :1]  # (ne, g, 1) weight times mu

    def form(A, B):  # int A_a B_b weighted, (ne, a, b); first index = test
        return np.einsum("egi,egj->eij", A * W[:, :, None], B)

    def mform(A, B):
        return np.einsum("egi,egj->eij", A * Wm, B)

    Kxx = mform(Gx, Gx) * (4 / 3) + mform(Gy, Gy)
    Kyy = mform(Gy, Gy) * (4 / 3) + mform(Gx, Gx)
    Kxy = 2 * (-mform(Gx, Gy) / 3 + 0.5 * mform(Gy, Gx))
    Kyx = np.swapaxes(Kxy, 1, 2)

    Nq = np.broadcast_to(N2[None], Gx.shape)
    Axy = rho_g * form(Gx, Nq)  # test v_x, trial u_y
    Ayy = rho_g * form(Gy, Nq)

    Np = np.broadcast_to(N1[None], (mesh.ne, 9, 4))
    B1 = form(Gx, Np)  # (test q2 node, trial q1 node) int phi_x psi
    B2 = form(Gy, Np)
    M = form(Np, Np)
    Cp = M / np.asarray(kappa)[:, None, None]

    c2, c1 = mesh.conn2, mesh.conn1
    n = 2 * nD + nP
    blocks = [
        (c2, c2, Kxx, 0, 0),
        (c2, c2, Kxy + Axy, 0, nD),
        (c2, c2, Kyx, nD, 0),
        (c2, c2, Kyy + Ayy, nD, nD),
        (c2, c1, -B1, 0, 2 * nD),
        (c2, c1, -B2, nD, 2 * nD),
        (c1, c2, np.swapaxes(B1, 1, 2), 2 * nD, 0),
        (c1, c2, np.swapaxes(B2, 1, 2), 2 * nD, nD),
        (c1, c1, Cp, 2 * nD, 2 * nD),
    ]
    mats = []
    for rc, cc, v, ro, co in blocks:
        mats.append(_coo(mesh, rc + ro, cc + co, v, (n, n)))
    return sum(mats[1:], mats[0]).tocsr()


def strain_dev(mesh: Mesh, u):
    """In-plane components (exx, eyy, exy) of the 3-D deviatoric strain at
    the Gauss points, shape (ne, 9, 3); ``u = [ux; uy]``."""
    Gx, Gy, _ = _grads(mesh)
    ux = u[: mesh.nD][mesh.conn2]
    uy = u[mesh.nD : 2 * mesh.nD][mesh.conn2]
    uxx = np.einsum("ega,ea->eg", Gx, ux)
    uxy_ = np.einsum("ega,ea->eg", Gy, ux)
    uyx = np.einsum("ega,ea->eg", Gx, uy)
    uyy = np.einsum("ega,ea->eg", Gy, uy)
    div = uxx + uyy
    return np.stack([uxx - div / 3, uyy - div / 3, 0.5 * (uxy_ + uyx)], axis=-1)


def history_rhs(mesh: Mesh, s):
    """int s : eps(v) for in-plane stress-like field s (ne, 9, 3)."""
    Gx, Gy, W = _grads(mesh)
    fx = np.einsum("eg,ega,eg->ea", W, Gx, s[..., 0]) + np.einsum(
        "eg,ega,eg->ea", W, Gy, s[..., 2]
    )
    fy = np.einsum("eg,ega,eg->ea", W, Gy, s[..., 1]) + np.einsum(
        "eg,ega,eg->ea", W, Gx, s[..., 2]
    )
    out = np.zeros(2 * mesh.nD)
    np.add.at(out, mesh.conn2.ravel(), fx.ravel())
    np.add.at(out, mesh.nD + mesh.conn2.ravel(), fy.ravel())
    return out


def top_traction_rhs(mesh: Mesh, ty):
    """Vertical traction on the top surface; ``ty`` has shape (nx, 3) with the
    traction values at the 3 Gauss points of every top-row element."""
    f = np.zeros(2 * mesh.nD)
    top_nodes = mesh.top()
    N, _ = _shape1d(_gp, 2)  # (3 gauss, 3 nodes)
    for e in range(mesh.nx):
        ids = top_nodes[2 * e : 2 * e + 3]
        f[mesh.nD + ids] += (mesh.hx[e] / 2) * (_gw * ty[e]) @ N
    return f


def top_gauss_x(mesh: Mesh):
    """x coordinates of the 3 Gauss points of each top-row element, (nx, 3)."""
    return mesh.xc[: mesh.nx, None] + 0.5 * mesh.hx[: mesh.nx, None] * _gp[None]


def winkler_matrix(mesh: Mesh, rho_g):
    """Buoyancy restoring force of the deflected free surface,
    rho g int u_y v_y ds on the top boundary (full system size)."""
    n = 2 * mesh.nD + mesh.nP
    M1 = np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]]) / 30.0
    top = mesh.top()
    ids = np.stack([top[2 * e : 2 * e + 3] for e in range(mesh.nx)]) + mesh.nD
    vals = rho_g * mesh.hx[: mesh.nx, None, None] * M1[None]
    return _coo(mesh, ids, ids, vals, (n, n)).tocsr()
