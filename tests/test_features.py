import numpy as np

from gia import Config, SECONDS_PER_YEAR, build_mesh, presets, simulate


def test_elastic_uniform_load_exact():
    cfg = presets.elastic_uniform_load()
    res = simulate(cfg)
    layer = cfg.layers[0]
    lam = layer.E * layer.nu / ((1 + layer.nu) * (1 - 2 * layer.nu))
    exact = -1.0e6 * layer.thickness / (lam + 2 * layer.mu)
    assert np.allclose(res.surface_uy[0], exact, rtol=1e-8)
    assert np.allclose(res.stress[:, 1], -1.0e6, rtol=1e-6)  # syy = -q


def test_refine_and_gmres_agree_with_direct():
    base = dict(nx_ice=4, nx_far=4, ny_layer=2, dt=1000 * SECONDS_PER_YEAR,
                t_end=5000 * SECONDS_PER_YEAR)
    m0 = build_mesh(Config(**base))
    assert build_mesh(Config(refine=1, **base)).ne == 4 * m0.ne
    d = simulate(Config(**base)).surface_uy
    g = simulate(Config(solver="gmres", **base)).surface_uy
    assert np.allclose(d, g, rtol=1e-6, atol=1e-6)


def test_glaciation_cycle_rebounds():
    cfg = presets.glaciation_cycle(t_lgm_yr=4000, t_eog_yr=8000, run_yr=40000,
                                   dt_yr=500, nx_ice=4, nx_far=8, ny_layer=4)
    res = simulate(cfg)
    uy = res.surface_uy[:, 0]
    assert uy[0] == 0 and uy.min() < -100
    assert abs(uy[-1]) < 0.1 * abs(uy.min())  # nearly fully rebounded
