# pyterrainwind — methods

Conventions, traps and limitations. Findings from the block this was developed on
are in the [README](../README.md). The solar chain this package sits on has its own
methods page in
[`pyterrainsol`](https://github.com/jborevitz/pyterrainsol/blob/main/docs/methods.md).

## Why wind is a separate package

Penman-Monteith needs net radiation on the inclined, shaded surface, so it needs the
solar chain. The solar chain never needs wind. The dependency runs one way, so the
packages are layered rather than peers, and the radiation code exists once rather
than twice.

The practical reason is that wind is an extra position to take. Aspect-driven
radiation is predictable geometry and is necessary at any site with relief; wind
exposure requires a transfer function from terrain to speed that is scheme-dependent
and, below, only partly constrained. A collaborator modelling light should not have
to install or argue with that.

## Wind conventions

Penman-Monteith's aerodynamic term scales as **u × VPD**, so a wind speed without
humidity multiplies an unknown and contributes nothing. Both are required.

Three conventions, each of which has cost someone a day:

- **Direction is where wind comes FROM**, degrees clockwise from north. A westerly
  (blowing toward the east, `u > 0`) reads 270.
- **Never average directions.** Average the components. Use `directional_constancy`
  — vector-mean over scalar-mean — to find out whether a mean direction means
  anything before quoting one.
- **Reanalysis timestamps are usually UTC.** Shift to local time before any diurnal
  analysis, or the afternoon maximum lands at dawn.

### Fitting a wind climatology, and how to be wrong at it

- **Do not calibrate on one year.** At a synoptically dominated site the best year in
  a twenty-year record can have half again the explanatory power of the pooled fit,
  and the largest annual amplitude. Fit the pool; use single years to measure spread.
- **Keep the model small.** Annual and diurnal sine/cosine pairs plus their
  interaction — nine terms — is usually the ceiling. The interaction is what lets the
  diurnal amplitude vary through the year and is often worth as much as the annual
  cycle itself. Higher orders fit weather.
- **Expect a low R².** Where fronts rather than the sun set the wind, a climatology
  explains a small fraction of variance by construction. That is a property of the
  place, not a failure of the fit.
- **Report amplitudes and phases, not coefficients**, and check which is which. A
  tensor product can be laid out seasonal-major or diurnal-major, and reading the
  wrong pair silently swaps the annual and daily cycles.
- **Weibull `k` is dimensionless and scale-free.** Doubling every speed moves `c` and
  leaves `k` alone, so a low `k` is never an artefact of a low background — it means
  a genuinely wider relative spread, typically a mixture of calm and synoptic
  regimes. Rank sites or seasons by mean speed and you can get their extremes
  backwards, since exceedance is `exp(-(x/c)^k)` and depends sharply on `k`.
- **Fit roughness on near-neutral steps only.** A stable nocturnal boundary layer
  steepens the two-height profile and inverts to a spuriously large z0. And if the
  two speeds come from a reanalysis you recover the *model's* effective roughness,
  a property of its land-surface configuration rather than a measurement of the
  landscape.

### Plotting a direction climatology

Use a **cyclic** colormap so 359° and 1° are the same colour and north carries no
seam, and a compass ring rather than a linear bar — a bar still implies the ends are
far apart. Note that `twilight_shifted` is near-white at its midpoint, so masked
cells must take a tone that appears nowhere in the ramp or they read as southerlies.
Scale opacity by directional constancy rather than by speed: fading what is
directionless is more honest than fading what is merely slow, and a hard speed cut
erases real light breezes.


## Humidity and vapour pressure, where the units bite

- **VPD is the mean of the saturation pressures, not the saturation pressure of the
  mean.** The curve is convex, so the two differ by several per cent and always in
  the same direction.
- **SILO publishes a `vp_deficit` that is not the FAO-56 deficit.** At the Avonside
  cell SILO's own value averages 0.836 kPa against 0.614 kPa computed the FAO way,
  and substituting it inflates annual ET0 by about 14%.
- **SILO publishes vapour pressure in hectopascals; FAO-56 wants kilopascals.** This
  one fails silently and is worth knowing by sight: with `ea` ten times too large,
  `es − ea` goes negative, a clip sends VPD to zero, and the whole aerodynamic term
  disappears without an error. ET0 drops to the radiation term alone and still looks
  like a plausible number. The symptom is a radiation/aerodynamic split that comes
  back 100/0. Pinned by `test_silo_vapour_pressure_unit_trap`.
- **Wind speed must be fitted as a scalar.** Fitting the u and v components and
  taking the magnitude of the fitted pair collapses the annual cycle, because
  opposing directions cancel in a vector mean. Components are right for the
  direction arrow; the scalar daily mean is right for the aerodynamic term.
- **Keep one clear-sky reference throughout.** Two clearness ratios coexist in this
  chain — observed shortwave over a flat unobstructed modelled clear sky, which
  scales the per-cell field, and observed over the FAO-56 `Rso`, which drives the
  longwave cloudiness term. They are different numbers. Mixing them biases any
  comparison between a real year and a climatology.

## Terrain exposure: two schemes, and what each can say

`sx_raster` implements the Winstral upwind shelter angle; `micromet_weight`
implements the Liston and Elder slope-plus-curvature weighting. Both are provided
because they answer different questions, and the difference is not cosmetic:

- **Sx is one-sided by construction.** It measures sheltering by upwind terrain and
  has no mechanism for speed-up above the free stream; its zero means "no upwind
  obstruction", not "accelerated". Converting it to a multiplier that exceeds 1
  asserts a speed-up the index never measured — that two-sidedness comes from
  re-centring on the field median, which is a normalisation choice, not physics.
- **MicroMet is genuinely two-sided.** Slope in the wind direction plus curvature
  gives acceleration over a convex crest and deceleration in a concave hollow, which
  is the mechanism a resident describes on a ridge.
- **Mean-preserving normalisation.** A reanalysis supplies the grid-cell areal
  average, so dividing by the field mean before applying the multiplier is the
  defensible way to redistribute that wind rather than change it. Note the raw mean
  is generally not 1.
- **Sx is dominated by the nearest terrain**, so near-field sampling detail matters
  more than fetch length. Two traps, both now commented in `exposure.py` and pinned
  by tests: an affine transform maps a coordinate to a cell **corner** while array
  interpolation treats index 0 as the cell **centre** (subtract half a cell, or every
  profile is displaced); and the ray must start at the **adjacent cell**, not at the
  step size.
- **Curvature at the array edge.** Fetching neighbours with `np.roll` wraps around
  the edge. Where a mask touches the window boundary with no margin this silently
  contaminates cells with extreme values. Use `np.pad(mode="edge")`.
- **Aspect azimuth ordering.** The downslope azimuth is `atan2(-dzdx, -dzdy)` — east
  component first, north second. Transposing the arguments rotates aspect by exactly
  90°, which zeroes any slope-in-wind-direction term computed from it and leaves a
  result that looks plausible rather than wrong.
- **Synthetic-terrain tests** must place the observation point at a cell centre, or
  bilinear interpolation halves the sampled ridge height. On an exact tilted plane
  every point along a ray subtends the same angle, so a "controlling distance"
  assertion is meaningless there; use a flat plain with one isolated ridge.

## The transfer function is the weak step

Nothing in either scheme converts geometry to metres per second from first
principles at this scale. Sx is in degrees of arc; MicroMet returns a dimensionless
weight. Jackson and Hunt's linear theory does give a crest speed-up of 2H/L, but it
is valid only to about H/L = 0.3, and steep ground exceeds that — the flow separates
and the theory over-predicts. Treat the chosen multiplier range as a **declared
scenario** and say so in any result that depends on it. The way to close it is two
anemometers, in the most and least exposed cells, and a measured ratio.

## What this does not do

- **No mountain-scale wind model.** A 4.4 km reanalysis such as BARRA-C2 already
  resolves valley channelling and delivers wind at the site pixel. Re-deriving that
  here would re-derive, badly, what the reanalysis did. This package reads regional
  terrain to *explain* the wind, never to *predict* it. What a reanalysis cannot
  contain is the last 500 m, and that is what the exposure fields are for.
- **No stability correction.** The roughness inversion assumes a neutral log
  profile. A stable nocturnal layer steepens the two-height ratio and inverts to a
  spuriously large z0.
- **No canopy.** Shelterbelts and tree cover change exposure substantially and are
  not represented; this is bare-earth terrain only.
- **No water limit.** This is demand, not actual evapotranspiration. What the
  surface actually gives up is a soil-water question and belongs in another package.
