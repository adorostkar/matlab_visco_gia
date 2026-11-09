# Glacial Isostatic Adjustment

Python finite-element model of glacial isostatic adjustment (GIA) of a
layered Maxwell viscoelastic flat Earth (2-D plane strain) loaded by an ice sheet.
The original MATLAB implementation is archived in `matlab_archive/`.

## Model

- Momentum: `div(sigma) + rho0 g grad(u_y) = 0`, `sigma = -p I + tau`
- Compressibility: `div(u) + p/kappa = 0` (nearly incompressible, `nu = 0.4999`)
- Maxwell rheology for the deviatoric stress: `dtau/dt + tau/t_M = 2 mu de/dt`, `t_M = eta/mu`
  (integrated exactly per step: `tau_{n+1} = a tau_n + 2 mu b (e_{n+1} - e_n)`,
  `a = exp(-dt/t_M)`, `b = (1-a) t_M/dt`)
- Gravity: bulk pre-stress advection plus the buoyancy (Winkler) restoring force of
  the deflected free surface. Without the latter the surface never reaches isostatic equilibrium.
- Ice load: vertical traction `-rho_ice g h(t)` on the top surface for `x <= half_width`.
- Q2-Q1 Taylor-Hood elements on a structured, gradable rectangular mesh; symmetry
  at `x=0`, roller at `x=L`, fixed bottom; one sparse LU factorisation reused for all steps.

## Usage

```
pip install -r requirements.txt
python run_gia.py            # default two-layer experiment, writes out/
python -m pytest tests       # verification
```

```python
from gia import Config, Layer, Ice, SECONDS_PER_YEAR, simulate
cfg = Config(layers=[Layer(1e5, eta=float("inf")), Layer(3.9e6, eta=1e21)],
             t_end=1e4 * SECONDS_PER_YEAR, dt=100 * SECONDS_PER_YEAR)
res = simulate(cfg)       # res.surface_uy[time, x]
```

## Verification

`tests/test_analytic.py` compares against the analytic response of an
incompressible Maxwell half-space to a sinusoidal load:
`u_y(t) = -q/(rho g) [1 - 2k mu/(rho g + 2k mu) exp(-t/tau)]`,
`tau = eta/mu + 2 k eta/(rho g)` (error < 0.5%). The default experiment also
reaches the isostatic depression `rho_ice h / rho_earth` (about 556 m for 2 km of ice).
