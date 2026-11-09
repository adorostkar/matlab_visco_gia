"""Ready-made configurations (ports of the MATLAB test problems)."""
from __future__ import annotations

import numpy as np

from .model import SECONDS_PER_YEAR, Config, Ice, Layer


def standard(emagn: float = 1.0, **kw) -> Config:
    """Default GIA test (Visco_parameters.m): two layers, 2 km boxcar ice sheet.

    ``emagn`` scales the Young's modulus of the lower layer (the E jump between
    the two subdomains in the MATLAB code)."""
    layers = [Layer(2.0e6), Layer(2.0e6, E=4.0e11 * emagn)]
    return Config(layers=layers, **kw)


def glaciation_cycle(
    t_beg_yr=0.0, t_lgm_yr=90000.0, t_eog_yr=100000.0, run_yr=120000.0,
    dt_yr=1000.0, **kw,
) -> Config:
    """Ice grows linearly to the glacial maximum, melts linearly, then the
    earth relaxes (MATLAB ``height_ice`` test case 11)."""
    s = SECONDS_PER_YEAR
    ice = Ice(history=Ice.glaciation(t_beg_yr * s, t_lgm_yr * s, t_eog_yr * s))
    return Config(ice=ice, t_end=run_yr * s, dt=dt_yr * s, **kw)


def elastic_uniform_load(nu=0.3, E=1.0e11, depth=1.0e6, width=1.0e6, q=1.0e6,
                         **kw) -> Config:
    """Pure-elasticity benchmark (MATLAB test problems 0/8/9): uniform surface
    pressure ``q`` on a laterally constrained elastic column, no gravity.
    Exact surface displacement: ``-q*depth/(lambda + 2 mu)``."""
    layer = Layer(depth, E=E, nu=nu, eta=float("inf"))
    ice = Ice(half_width=width / 2, thickness=q / (917.0 * 9.81),
              profile=lambda x: np.ones_like(x))
    return Config(L=width, layers=[layer], ice=ice, gravity=False,
                  advection=False, nx_ice=4, nx_far=4, x_ratio=1.0,
                  ny_layer=4, t_end=0.0, dt=SECONDS_PER_YEAR, **kw)
