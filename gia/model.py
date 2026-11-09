"""Physical model description: earth layers, ice load, numerical settings."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional

SECONDS_PER_YEAR = 365.25 * 24 * 3600.0


@dataclass
class Layer:
    """Maxwell viscoelastic layer. ``eta=inf`` gives a purely elastic layer."""

    thickness: float  # m
    E: float = 4.0e11  # Young's modulus, Pa
    nu: float = 0.4999  # Poisson ratio
    eta: float = 1.0e21  # viscosity, Pa s

    @property
    def mu(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def kappa(self) -> float:
        return self.E / (3.0 * (1.0 - 2.0 * self.nu))

    @property
    def maxwell_time(self) -> float:
        return self.eta / self.mu


@dataclass
class Ice:
    """Boxcar ice sheet occupying 0 <= x <= half_width (symmetry axis at x=0).

    ``history`` maps time in seconds to a thickness factor in [0, 1]; the
    default is an ice sheet that is fully present from t=0.
    """

    half_width: float = 1.0e6  # m
    thickness: float = 2.0e3  # m
    density: float = 917.0  # kg/m^3
    history: Optional[Callable[[float], float]] = None
    # optional load shape f(x) in [0, 1] over the whole surface (replaces the
    # boxcar; used e.g. for sinusoidal verification loads)
    profile: Optional[Callable] = None

    def height(self, t: float) -> float:
        f = 1.0 if self.history is None else float(self.history(t))
        return self.thickness * f


@dataclass
class Config:
    L: float = 5.0e7  # horizontal extent (m); x in [0, L]
    layers: List[Layer] = field(
        default_factory=lambda: [Layer(2.0e6), Layer(2.0e6)]
    )
    rho_earth: float = 3300.0  # kg/m^3
    g: float = 9.81  # m/s^2
    ice: Ice = field(default_factory=Ice)
    gravity: bool = True  # buoyancy restoring force of the free surface
    advection: bool = True  # pre-stress advection term in the bulk
    # discretisation
    nx_ice: int = 8  # elements under the ice
    nx_far: int = 24  # elements outside the ice
    x_ratio: float = 1.25  # geometric growth of far-field elements
    ny_layer: int = 8  # elements per layer
    y_ratio: float = 1.0  # geometric growth towards depth (>1 refines surface)
    # time stepping
    t_end: float = 20000.0 * SECONDS_PER_YEAR
    dt: float = 100.0 * SECONDS_PER_YEAR

    @property
    def depth(self) -> float:
        return sum(l.thickness for l in self.layers)
