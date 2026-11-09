"""Terrain wind exposure and the evaporative demand it drives.

Three layers, in the order a calculation uses them:

* **climatology** -- speed and direction conventions, directional constancy,
  the neutral log-profile roughness inversion, a seasonal x diurnal harmonic
  model and Weibull gust shape.  What the wind does at a point, over years.
* **exposure** -- how terrain redistributes that wind across a block.  Two
  schemes, deliberately both: `sx_raster`, the upwind shelter angle over a
  short fetch (Winstral), and `micromet_weight`, the Liston and Elder slope-
  plus-curvature weighting, which unlike Sx is two-sided and so can represent
  speed-up over a crest as well as sheltering in a lee.
* **demand** -- FAO-56 reference ET0, the aerodynamic term that radiation
  alone cannot supply.

Requires `pyterrainsol <https://github.com/jborevitz/pyterrainsol>`_ for net
radiation on the inclined, shaded surface.  The dependency is one way: wind
needs the solar chain, the solar chain never needs wind.
"""
from .wind import (speed_direction, components, directional_constancy,
                   roughness_length, scale_speed, harmonic_design, fit_harmonic,
                   harmonic_summary, weibull_fit, exceedance)
from .exposure import (upwind_profile, horizon_angle, sx, sx_rose, sx_raster,
                       slope_aspect, curvature, micromet_weight)
from .demand import vapour_pressure_deficit, wind_at_2m, penman_monteith

__all__ = ["speed_direction","components","directional_constancy",
           "roughness_length","scale_speed","harmonic_design","fit_harmonic",
           "harmonic_summary","weibull_fit","exceedance",
           "upwind_profile","horizon_angle","sx","sx_rose","sx_raster",
           "slope_aspect","curvature","micromet_weight",
           "vapour_pressure_deficit","wind_at_2m","penman_monteith"]
__version__ = "0.1.0"
