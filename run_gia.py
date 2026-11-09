"""Run the default GIA experiment: a 2 km thick, 1000 km half-width ice
boxcar placed on a two-layer Maxwell earth at t = 0.

    python run_gia.py [--years 20000] [--dt 250] [--out out] [--no-plot]
"""
import argparse
import os

import numpy as np

from gia import Config, SECONDS_PER_YEAR, simulate


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--years", type=float, default=20000.0)
    ap.add_argument("--dt", type=float, default=250.0, help="time step (years)")
    ap.add_argument("--out", default="out")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--no-advection", action="store_true",
                    help="drop the bulk pre-stress advection term")
    args = ap.parse_args()

    cfg = Config(t_end=args.years * SECONDS_PER_YEAR,
                 dt=args.dt * SECONDS_PER_YEAR, advection=not args.no_advection)
    res = simulate(cfg)

    os.makedirs(args.out, exist_ok=True)
    np.savez(os.path.join(args.out, "gia_surface.npz"),
             times_years=res.times / SECONDS_PER_YEAR, x=res.surface_x,
             ux=res.surface_ux, uy=res.surface_uy)
    yrs = res.times / SECONDS_PER_YEAR
    step = max(1, len(yrs) // 10)
    for k in range(0, len(yrs), step):
        print(f"year {yrs[k]:8.0f}: centre uy = {res.surface_uy[k, 0]:9.2f} m, "
              f"max |ux| = {np.abs(res.surface_ux[k]).max():8.2f} m")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for k in np.linspace(0, len(yrs) - 1, 6).astype(int):
            ax[0].plot(res.surface_x / 1e3, res.surface_uy[k],
                       label=f"{yrs[k]:.0f} yr")
        ax[0].set_xlim(0, 4000)
        ax[0].set_xlabel("x (km)"); ax[0].set_ylabel("vertical displacement (m)")
        ax[0].legend()
        ax[1].plot(yrs, res.surface_uy[:, 0])
        ax[1].set_xlabel("time (yr)"); ax[1].set_ylabel("uy under ice centre (m)")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "gia_surface.png"), dpi=120)


if __name__ == "__main__":
    main()
