"""Structured Q2-Q1 Taylor-Hood mesh on a rectangle."""
from __future__ import annotations

import numpy as np


class Mesh:
    """Rectilinear mesh. ``xv``/``yv`` are the (ascending) vertex coordinates.

    Q2 (velocity/displacement) node (I, J), 0<=I<=2nx, 0<=J<=2ny, has id
    ``J*(2nx+1)+I``; Q1 (pressure) vertex (i, j) has id ``j*(nx+1)+i``.
    Element e = ey*nx + ex. Local node a = j*3+i (Q2) or j*2+i (Q1).
    """

    def __init__(self, xv, yv):
        self.xv = np.asarray(xv, float)
        self.yv = np.asarray(yv, float)
        self.nx = len(self.xv) - 1
        self.ny = len(self.yv) - 1
        nx, ny = self.nx, self.ny
        self.ne = nx * ny
        self.hx = np.tile(np.diff(self.xv), ny)
        self.hy = np.repeat(np.diff(self.yv), nx)
        self.xc = np.tile(0.5 * (self.xv[:-1] + self.xv[1:]), ny)
        self.yc = np.repeat(0.5 * (self.yv[:-1] + self.yv[1:]), nx)

        # Q2 node coordinates
        x2 = np.empty(2 * nx + 1)
        x2[0::2] = self.xv
        x2[1::2] = 0.5 * (self.xv[:-1] + self.xv[1:])
        y2 = np.empty(2 * ny + 1)
        y2[0::2] = self.yv
        y2[1::2] = 0.5 * (self.yv[:-1] + self.yv[1:])
        self.x2, self.y2 = x2, y2
        self.n2x, self.n2y = 2 * nx + 1, 2 * ny + 1
        self.nD = self.n2x * self.n2y
        self.nP = (nx + 1) * (ny + 1)
        X, Y = np.meshgrid(x2, y2)
        self.node_x = X.ravel()
        self.node_y = Y.ravel()

        ey, ex = np.divmod(np.arange(self.ne), nx)
        j, i = np.divmod(np.arange(9), 3)
        self.conn2 = (2 * ey[:, None] + j) * self.n2x + (2 * ex[:, None] + i)
        j, i = np.divmod(np.arange(4), 2)
        self.conn1 = (ey[:, None] + j) * (nx + 1) + (ex[:, None] + i)

    # boundary node sets (Q2 ids)
    def left(self):
        return np.arange(self.n2y) * self.n2x

    def right(self):
        return self.left() + self.n2x - 1

    def bottom(self):
        return np.arange(self.n2x)

    def top(self):
        return (self.n2y - 1) * self.n2x + np.arange(self.n2x)

    def refine(self, times: int = 1) -> "Mesh":
        """Uniform refinement (each element split in four); layer bookkeeping
        is carried over."""
        xv, yv = self.xv, self.yv
        layer = getattr(self, "layer_of_row", None)
        for _ in range(times):
            xv = np.sort(np.concatenate([xv, 0.5 * (xv[:-1] + xv[1:])]))
            yv = np.sort(np.concatenate([yv, 0.5 * (yv[:-1] + yv[1:])]))
            if layer is not None:
                layer = np.repeat(layer, 2)
        m = Mesh(xv, yv)
        if layer is not None:
            m.layer_of_row = layer
        return m

    def element_at(self, x, y):
        """Index of the element containing the point (x, y)."""
        i = np.clip(np.searchsorted(self.xv, x, side="right") - 1, 0, self.nx - 1)
        j = np.clip(np.searchsorted(self.yv, y, side="right") - 1, 0, self.ny - 1)
        return j * self.nx + i
