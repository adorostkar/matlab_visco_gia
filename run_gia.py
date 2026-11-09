"""Run the default GIA experiment: a 2 km thick, 1000 km half-width ice
boxcar placed on a two-layer Maxwell earth at t = 0.

    python run_gia.py [--years 20000] [--dt 250] [--out out] [--no-plot]
"""
import argparse
import os

import numpy as np

from gia import SECONDS_PER_YEAR, presets, simulate
from gia.plotting import plot_field, plot_mesh, plot_surface


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--years", type=float, default=20000.0)
    ap.add_argument("--dt", type=float, default=250.0, help="time step (years)")
    ap.add_argument("--out", default="out")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--no-advection", action="store_true",
                    help="drop the bulk pre-stress advection term")
    ap.add_argument("--preset", choices=["standard", "glaciation"],
                    default="standard", help="standard: ice on from t=0; "
                    "glaciation: growth, melt and relaxation")
    ap.add_argument("--refine", type=int, default=0, help="uniform refinements")
    ap.add_argument("--emagn", type=float, default=1.0,
                    help="E jump of the lower layer (standard preset)")
    ap.add_argument("--solver", choices=["direct", "gmres"], default="direct")
    ap.add_argument("--snapshots", type=float, nargs="*", default=[],
                    help="years at which to save full fields (npz)")
    args = ap.parse_args()

    common = dict(refine=args.refine, solver=args.solver,
                  advection=not args.no_advection)
    if args.preset == "glaciation":
        cfg = presets.glaciation_cycle(run_yr=args.years, dt_yr=args.dt, **common)
    else:
        cfg = presets.standard(args.emagn, t_end=args.years * SECONDS_PER_YEAR,
                               dt=args.dt * SECONDS_PER_YEAR, **common)
    res = simulate(cfg, snapshot_times=[y * SECONDS_PER_YEAR for y in args.snapshots])

    os.makedirs(args.out, exist_ok=True)
    np.savez(os.path.join(args.out, "gia_surface.npz"),
             times_years=res.times / SECONDS_PER_YEAR, x=res.surface_x,
             ux=res.surface_ux, uy=res.surface_uy)
    for snap in res.snapshots:
        np.savez(os.path.join(args.out, f"data-y{snap['time'] / SECONDS_PER_YEAR:.0f}.npz"),
                 **snap)
    yrs = res.times / SECONDS_PER_YEAR
    step = max(1, len(yrs) // 10)
    for k in range(0, len(yrs), step):
        print(f"year {yrs[k]:8.0f}: centre uy = {res.surface_uy[k, 0]:9.2f} m, "
              f"max |ux| = {np.abs(res.surface_ux[k]).max():8.2f} m")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        plot_surface(res, ax[0], scale=1e3, xlim=(0, 4000))
        ax[0].set_xlabel("x (km)")
        ax[1].plot(yrs, res.surface_uy[:, 0])
        ax[1].set_xlabel("time (yr)"); ax[1].set_ylabel("uy under ice centre (m)")
        plot_field(res.mesh, res.u[res.mesh.nD:], ax[2], "final uy (m)")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "gia_surface.png"), dpi=120)
        fig, a = plt.subplots(figsize=(8, 3))
        plot_mesh(res.mesh, a)
        fig.savefig(os.path.join(args.out, "gia_mesh.png"), dpi=120)


if __name__ == "__main__":
    main()
