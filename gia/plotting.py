"""Visualisation helpers (port of Bvisual_mesh and the MATLAB figures)."""
from __future__ import annotations

import numpy as np


def _plt():
    import matplotlib.pyplot as plt

    return plt


def plot_mesh(mesh, ax=None, nodes=False):
    """Draw element edges (and optionally the Q2 nodes)."""
    ax = ax or _plt().gca()
    for x in mesh.xv:
        ax.plot([x, x], [mesh.yv[0], mesh.yv[-1]], "k-", lw=0.5)
    for y in mesh.yv:
        ax.plot([mesh.xv[0], mesh.xv[-1]], [y, y], "k-", lw=0.5)
    if nodes:
        ax.plot(mesh.node_x, mesh.node_y, "r.", ms=2)
    ax.set_aspect("equal")
    return ax


def plot_surface(res, ax=None, scale=1.0, indices=None, xlim=None):
    """Surface vertical displacement profiles at selected time indices."""
    ax = ax or _plt().gca()
    idx = indices if indices is not None else np.linspace(
        0, len(res.times) - 1, 6).astype(int)
    yrs = res.times / (365.25 * 24 * 3600)
    for k in idx:
        ax.plot(res.surface_x / scale, res.surface_uy[k], label=f"{yrs[k]:.0f} yr")
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_xlabel("x"); ax.set_ylabel("vertical displacement (m)")
    ax.legend()
    return ax


def plot_field(mesh, values, ax=None, title=None):
    """Filled plot of a per-node (Q2) or per-element field."""
    plt = _plt()
    ax = ax or plt.gca()
    v = np.asarray(values)
    if v.size == mesh.nD:
        m = ax.pcolormesh(mesh.x2, mesh.y2, v.reshape(mesh.n2y, mesh.n2x),
                          shading="gouraud")
    else:
        m = ax.pcolormesh(mesh.xv, mesh.yv, v.reshape(mesh.ny, mesh.nx))
    plt.colorbar(m, ax=ax)
    if title:
        ax.set_title(title)
    return ax
