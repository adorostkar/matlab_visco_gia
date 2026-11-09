"""Verification against the analytic response of an incompressible Maxwell
half-space with surface buoyancy under a sinusoidal load q cos(kx):

    u_y(t) = -q/(rho g) [1 - 2k mu/(rho g + 2k mu) exp(-t/tau)],
    tau    = eta/mu + 2 k eta/(rho g)
"""
import numpy as np
import pytest

from gia import Config, Ice, Layer, SECONDS_PER_YEAR, simulate


def run(nx=16, ny=16, dt_years=500.0, t_end_years=40000.0, y_ratio=1.0):
    lam = 2.0e6
    k = 2 * np.pi / lam
    eta, E, nu = 1.0e21, 2.6e11, 0.4999
    layer = Layer(thickness=lam, E=E, nu=nu, eta=eta)
    q = 1.0e6  # Pa (small enough to be linear anyway)
    ice = Ice(half_width=lam / 4, thickness=q / (917.0 * 9.81), density=917.0,
              profile=lambda x: np.cos(k * x))
    cfg = Config(L=lam / 2, layers=[layer], ice=ice, nx_ice=nx // 2,
                 nx_far=nx - nx // 2, x_ratio=1.0, ny_layer=ny, y_ratio=y_ratio,
                 dt=dt_years * SECONDS_PER_YEAR,
                 t_end=t_end_years * SECONDS_PER_YEAR)
    res = simulate(cfg)
    rho_g = cfg.rho_earth * cfg.g
    mu = layer.mu
    tau = eta / mu + 2 * k * eta / rho_g
    t = res.times
    exact = -q / rho_g * (1 - 2 * k * mu / (rho_g + 2 * k * mu) * np.exp(-t / tau))
    return res.surface_uy[:, 0], exact, tau, q / rho_g


def test_sinusoidal_load_matches_analytic_relaxation():
    num, exact, tau, uinf = run()
    # elastic response and the full relaxation curve
    assert abs(num[0] - exact[0]) < 0.03 * abs(exact[0])
    assert np.max(np.abs(num - exact)) < 0.03 * uinf
    # reaches isostatic equilibrium -q/(rho g)
    assert abs(num[-1] + uinf) < 0.02 * uinf


def test_viscous_relaxation_is_monotone_and_converges():
    n1, e1, *_ = run(nx=8, ny=8)
    n2, e2, *_ = run(nx=16, ny=16)
    err1 = np.max(np.abs(n1 - e1))
    err2 = np.max(np.abs(n2 - e2))
    assert err2 < err1
    assert np.all(np.diff(n2) <= 1e-9)  # sinks monotonically
