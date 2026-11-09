"""Penman-Monteith, the wind-aware companion to Priestley-Taylor.

Also pins the two unit traps that cost real debugging time: FAO-56 defines es
as the mean of the saturation pressures rather than the saturation pressure of
the mean, and SILO publishes vapour pressure in hPa where FAO-56 wants kPa.
"""

import inspect

import numpy as np
import pytest

from pyterrainwind.demand import (penman_monteith, vapour_pressure_deficit,
                                  wind_at_2m)
from pyterrainsol.demand import priestley_taylor, saturation_vapour_pressure


def test_wind_at_2m_matches_fao_table():
    # FAO-56 Table 4: the 10 m conversion factor is 0.748.
    assert wind_at_2m(1.0, 10.0) == pytest.approx(0.748, abs=0.001)
    # 2 m is the reference height, so the factor is 1 to within rounding
    # (4.87/ln(67.8*2-5.42) = 1.00022, not exactly unity).
    assert wind_at_2m(3.3, 2.0) == pytest.approx(3.3, rel=1e-3)


def test_vpd_uses_mean_of_saturation_not_saturation_of_mean():
    # The SVP curve is convex, so es(Tmin),es(Tmax) averaged always exceeds
    # es of the mean temperature.  A build that uses the latter under-reads.
    tmin, tmax = 2.0, 24.0
    vpd = vapour_pressure_deficit(tmin, tmax, 0.0)
    naive = saturation_vapour_pressure(0.5 * (tmin + tmax))
    assert vpd > naive
    # Never negative, even when the air is supersaturated in the record.
    assert vapour_pressure_deficit(5.0, 6.0, 3.0) == 0.0


def test_penman_monteith_assembles_the_fao56_equation():
    """Check the ET0 assembly against the equation written out by hand.

    Delta and gamma are typed as literals here rather than taken from the
    module, so this exercises the assembly of FAO-56 eq. 6 independently of
    svp_slope and psychrometric_constant (which have their own tests).
    """
    T, Z, u2, Rn, G, vpd = 30.2, 2.0, 2.0, 14.33, 0.14, 1.15
    delta, gamma = 0.2458, 0.06735          # at T = 30.2 C and 2 m altitude
    expected = (0.408 * delta * (Rn - G)
                + gamma * (900.0 / (T + 273.0)) * u2 * vpd) \
        / (delta + gamma * (1.0 + 0.34 * u2))
    assert expected == pytest.approx(5.246, abs=0.002)   # pins the arithmetic
    assert penman_monteith(Rn, T, Z, u2, vpd, soil_heat=G) \
        == pytest.approx(expected, rel=1e-3)


def test_penman_monteith_split_sums_to_total():
    rad, aero = penman_monteith(12.0, 18.0, 861.0, 2.4, 0.9, split=True)
    total = penman_monteith(12.0, 18.0, 861.0, 2.4, 0.9)
    assert rad + aero == pytest.approx(total, rel=1e-12)
    assert aero > 0.0


def test_aerodynamic_term_vanishes_in_still_air_but_radiation_does_not():
    # With u2 = 0 Penman-Monteith keeps only the radiation term, which is the
    # structural difference from Priestley-Taylor: PT has no u at all.
    rad, aero = penman_monteith(12.0, 18.0, 861.0, 0.0, 0.9, split=True)
    assert aero == 0.0
    assert rad > 0.0


def test_wind_raises_demand_monotonically():
    u = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    et = penman_monteith(12.0, 18.0, 861.0, u, 0.9)
    assert np.all(np.diff(et) > 0)


def test_priestley_taylor_ignores_wind_by_construction():
    # Guards the "two options" contract: PT must stay a wind-free baseline, so
    # anyone adding a u argument to it has made an error.
    assert "wind" not in inspect.signature(priestley_taylor).parameters


def test_silo_vapour_pressure_unit_trap():
    """SILO publishes vp in hPa; feeding it straight in zeroes the aero term.

    This is a real failure mode: es - ea goes negative, the clip sends VPD to
    zero, and ET0 silently loses its entire aerodynamic term rather than
    raising.  The test documents the symptom so it is recognised next time.
    """
    tmin, tmax, vp_hpa = 4.8, 18.4, 9.93
    wrong = vapour_pressure_deficit(tmin, tmax, vp_hpa)
    right = vapour_pressure_deficit(tmin, tmax, vp_hpa / 10.0)
    assert wrong == 0.0
    assert right == pytest.approx(0.495, abs=0.01)
