"""Wind-aware evaporative demand: FAO-56 reference evapotranspiration.

The aerodynamic half of the demand calculation.  Energy, net radiation and the
thermodynamic helpers come from
`pyterrainsol <https://github.com/jborevitz/pyterrainsol>`_; this module adds
the term that radiation alone cannot supply, proportional to ``u2 x VPD``.

The dependency runs one way only.  Penman-Monteith needs net radiation on the
inclined, shaded surface, so it needs the solar chain; the solar chain never
needs wind.  That is why these are layered packages rather than peers.

Use :func:`penman_monteith` when wind matters -- which in this landscape means
winter: net radiation is near zero in June but the air keeps moving, and the
aerodynamic share of ET0 peaks near 61%.  Use
:func:`pyterrainsol.demand.priestley_taylor` when the signal you resolve is
radiation, or when you want a wind-free baseline to difference against.

Formulae follow FAO Irrigation and Drainage Paper 56 (Allen et al. 1998).
"""

import numpy as np

from pyterrainsol.demand import (saturation_vapour_pressure, svp_slope,
                                 psychrometric_constant)


def vapour_pressure_deficit(t_min_c, t_max_c, vapour_pressure):
    """FAO-56 vapour pressure deficit (kPa) from daily extremes and actual ea.

    FAO-56 defines es as the MEAN of the saturation pressures at Tmin and Tmax,
    not the saturation pressure at mean temperature -- the curve is convex, so
    the two differ by several per cent and always in the same direction.

    This matters for SILO users in particular: SILO publishes a ``vp_deficit``
    variable that is NOT this quantity.  At 845 Avonside SILO's own value
    averages 0.836 kPa against 0.614 kPa here, and substituting it inflates
    annual ET0 by about 14%.  Compute the deficit from tmin/tmax/vp instead.
    """
    es = 0.5 * (saturation_vapour_pressure(t_min_c)
                + saturation_vapour_pressure(t_max_c))
    return np.maximum(es - np.asarray(vapour_pressure, float), 0.0)


def wind_at_2m(speed, height_m):
    """Scale a wind speed from `height_m` to the FAO-56 reference 2 m.

    The FAO-56 logarithmic correction for the reference grass surface,
    4.87 / ln(67.8 z - 5.42).  Reanalysis and most weather stations report at
    10 m, where the factor is 0.748.
    """
    z = np.asarray(height_m, float)
    return np.asarray(speed, float) * 4.87 / np.log(67.8 * z - 5.42)


def penman_monteith(net_rad, t_c, elevation_m, wind_2m, vpd,
                    soil_heat=0.0, split=False):
    """FAO-56 reference evapotranspiration ET0, mm per day.

    The companion to :func:`priestley_taylor`.  Both answer "how much water can
    this atmosphere take", but they differ in what they are allowed to know:

    * ``priestley_taylor`` sees energy only.  Advection is absorbed into a
      fixed ``alpha``, so it cannot respond to a windy day at all.  Use it when
      the signal you resolve is radiation, and when you want a wind-free
      baseline to difference against.
    * ``penman_monteith`` adds the aerodynamic term, proportional to
      ``u2 x VPD``.  Use it when wind matters -- which at 845 Avonside means
      midwinter, where the aerodynamic share of ET0 peaks near 61% because net
      radiation is near zero while the air keeps moving.

    Parameters are daily means: `net_rad` and `soil_heat` in MJ m-2 d-1, `t_c`
    in degrees C, `wind_2m` in m/s at 2 m (see :func:`wind_at_2m`), `vpd` in
    kPa (see :func:`vapour_pressure_deficit`).

    With ``split=True`` returns ``(radiation_term, aerodynamic_term)`` instead
    of their sum, which is how the 60/40 annual partition in the README was
    measured.  Negative values are clipped to zero only on the total, so a
    genuinely negative midwinter radiation term stays visible in the split.
    """
    delta = svp_slope(t_c)
    gamma = psychrometric_constant(elevation_m)
    u2 = np.asarray(wind_2m, float)
    available = np.asarray(net_rad, float) - np.asarray(soil_heat, float)
    denom = delta + gamma * (1.0 + 0.34 * u2)
    rad = 0.408 * delta * available / denom
    aero = gamma * (900.0 / (np.asarray(t_c, float) + 273.0)) * u2 \
        * np.asarray(vpd, float) / denom
    if split:
        return rad, aero
    return np.maximum(rad + aero, 0.0)
