# IM2 Hillslope Hydrology

## Overview

The IM2 hillslope hydrology scheme (`use_IM2_hillslope_hydrology`) enables
subgrid lateral water redistribution among topographic units (topounits)
within a land gridcell. Rather than treating each topounit's surface and
subsurface runoff as immediately leaving the gridcell, the scheme routes a
fraction of those fluxes downslope to the next lower topounit, where they
are received as an additional water input to the soil surface. This provides
a simple, computationally inexpensive representation of hillslope-scale
hydrologic connectivity that does not require solving a full lateral flow
equation.

The flag is set via the namelist parameter `use_IM2_hillslope_hydrology` in
the `elm_inparm` group and defaults to `.false.`.

## Conceptual Framework

Each land gridcell contains one or more topounits ordered by elevation.
When `use_IM2_hillslope_hydrology = .true.`, the scheme:

1. **Routes lateral outflow downslope.** At the end of each timestep a
   fixed fraction $f_{\downarrow}$ of three column-level runoff fluxes is
   redirected to the next-lower (downhill) topounit rather than exiting the
   gridcell:
   - surface runoff ($q_{\text{surf}}$),
   - perched water table drainage ($q_{\text{drain,perched}}$), and
   - surface water runoff ($q_{\text{h2osfc,surf}}$).

2. **Stores received water at the topounit level.** The redirected water is
   accumulated in the topounit-level state variable `from_uphill` (kg m⁻²).
   Area weighting is applied so that mass is conserved when the uphill and
   downhill topounits cover different fractions of the gridcell.

3. **Distributes received water to columns.** Each timestep, a fixed
   fraction $f_{\uparrow}$ of the `from_uphill` store is partitioned among
   the soil/crop/pervious-road columns of the receiving topounit in
   proportion to each column's topounit weight, and added to the net water
   input at the top of the soil profile ($q_{\text{top_soil}}$).

The two fractions are hard-coded constants defined in `elm_varcon.F90`:

| Constant          | Variable name      | Value | Description |
| :---------------- | :----------------- | :---: | :---------- |
| $f_{\downarrow}$  | `frac_to_downhill` | 1.0   | Fraction of column runoff fluxes transferred to the downhill topounit |
| $f_{\uparrow}$    | `frac_from_uphill` | 0.5   | Fraction of the `from_uphill` state delivered to columns per timestep |

With $f_{\downarrow} = 1.0$, all of the selected runoff fluxes are
redirected downhill. With $f_{\uparrow} = 0.5$, the stored water is
delivered with an exponential decay (a linear reservoir), smoothing the
input over multiple timesteps.

## Topounit Connectivity

Topounits within a gridcell are connected according to their elevation.
The downhill neighbor of each topounit is determined at initialization
(`initGridCellsMod.F90`) by finding the nearest lower-elevation topounit
on the same gridcell:

```
For each topounit t1:
    find t2 on the same gridcell such that
        elevation(t2) < elevation(t1)  and
        elevation(t2) is the maximum among all such t2
    downhill_ti(t1) = t2   (or -1 if no lower topounit exists)
```

The lowest topounit on any gridcell has `downhill_ti = -1` and sends no
water downhill; its runoff exits the gridcell normally.

## Water Routing Algorithm

### Step 1 — Downhill transfer (in `HydrologyDrainage`)

For every column $c$ on topounit $t$ that has a downhill neighbor
$t' = \text{downhill_ti}(t)$:

$$
\begin{aligned}
q_{\text{to_downhill}}(c)
&= f_{\downarrow}\,\max\left(0,\,q_{\text{surf}}(c)\right) \\
&\quad + f_{\downarrow}\,\max\left(0,\,q_{\text{drain,perched}}(c)\right) \\
&\quad + f_{\downarrow}\,\max\left(0,\,q_{\text{h2osfc,surf}}(c)\right)
\end{aligned}
$$

The corresponding column fluxes are reduced by the same amount:

$$
q_{\text{surf}}(c) \leftarrow q_{\text{surf}}(c) - f_{\downarrow}\,\max\left(0,\,q_{\text{surf}}(c)\right)
$$

and similarly for $q_{\text{drain,perched}}$ and $q_{\text{h2osfc,surf}}$.

Only positive (outgoing) fluxes are transferred; negative values (which
would represent water returning from the stream) are left unchanged.

The redirected water is added to the `from_uphill` state of the downhill
topounit, with an area-weighting factor that converts from the uphill
column's area to the downhill topounit's area:

$$
\Delta S_{t'} = q_{\text{to_downhill}}(c)\;\Delta t
\;\frac{w_c^{\text{gcell}}}{w_{t'}^{\text{gcell}}}
$$

where $w_c^{\text{gcell}}$ is the fractional weight of column $c$ relative
to the gridcell and $w_{t'}^{\text{gcell}}$ is the fractional weight of the
downhill topounit relative to the gridcell. This ensures mass conservation
when topounits differ in areal extent.

The cumulative update over all columns on all uphill topounits gives:

$$
S_{t'}^{t} = S_{t'}^{\,t-1} + \sum_{c \in \text{uphill}} \Delta S_{t'}(c)
$$

### Step 2 — Uphill input to columns (in `SurfaceRunoff`)

At the beginning of the next call to `SurfaceRunoff`, the water stored in
`from_uphill` is distributed to the columns of its topounit.

First, the sum of column weights on the topounit (for the hydrology filter,
which includes soil, crop, and pervious-road columns) is accumulated:

$$
W_t = \sum_{c \in \text{hydrology filter on } t} w_c^{\text{topounit}}
$$

Then, the per-column flux is:

$$
q_{\text{from_uphill}}(c)
= \frac{w_c^{\text{topounit}}}{W_t}
\;\frac{f_{\uparrow}\;S_t}{\Delta t}
$$

where $S_t$ is the current `from_uphill` state (kg m⁻²) and $\Delta t$ is
the model timestep. This flux is added to the net top-of-soil water input:

$$
q_{\text{top_soil}}(c) \leftarrow q_{\text{top_soil}}(c) + q_{\text{from_uphill}}(c)
$$

### Step 3 — State update (in `BalanceCheck`)

After the fluxes have been applied the `from_uphill` store is reduced by
the delivered amount:

$$
S_t \leftarrow \max\left(0,\;S_t - q_{\text{from_uphill}}(c)\,\Delta t\right)
$$

A floor of $10^{-20}$ kg m⁻² is applied to prevent spurious negative
values from floating-point round-off during recession:

$$
\text{if } S_t < 10^{-20} \text{ then } S_t = 0
$$

## Water Balance

The lateral fluxes enter the column water budget explicitly. In the balance
check the column error is:

$$
\begin{aligned}
\epsilon(c) &= \Delta W(c) \\
&\quad - \left(P(c) + q_{\text{flood}}(c) + q_{\text{from_uphill}}(c) + q_{\text{irrig}}(c)\right)\Delta t \\
&\quad + \left(E(c) + q_{\text{surf}}(c) + q_{\text{h2osfc,surf}}(c) + q_{\text{to_downhill}}(c)\right. \\
&\qquad \left. +\; q_{\text{drain}}(c) + q_{\text{drain,perched}}(c) + q_{\text{snwcp,ice}}(c)\right)\Delta t
\end{aligned}
$$

where $P$ is precipitation, $E$ is total evapotranspiration, and $\Delta W$
is the change in column total water storage. A non-zero $\epsilon$ indicates
a mass conservation error.

## Variables Summary

| Symbol | Code variable | Units | Description |
| :----- | :------------ | :---: | :---------- |
| $S_t$ | `top_ws%from_uphill(t)` | kg m⁻² | Water stored at topounit from uphill transfer |
| $q_{\text{to_downhill}}$ | `col_wf%qflx_to_downhill(c)` | mm s⁻¹ | Column flux sent to downhill topounit |
| $q_{\text{from_uphill}}$ | `col_wf%qflx_from_uphill(c)` | mm s⁻¹ | Column flux received from uphill topounit |
| $w_c^{\text{topounit}}$ | `col_pp%wttopounit(c)` | — | Column weight relative to topounit |
| $w_c^{\text{gcell}}$ | `col_pp%wtgcell(c)` | — | Column weight relative to gridcell |
| $w_t^{\text{gcell}}$ | `top_pp%wtgcell(t)` | — | Topounit weight relative to gridcell |
| $W_t$ | `top_pp%uphill_wt(t)` | — | Sum of column weights on topounit (hydrology filter) |
| — | `top_pp%downhill_ti(t)` | — | Topounit index of downhill neighbor (−1 if none) |
| $f_{\downarrow}$ | `frac_to_downhill` | — | Fraction of runoff sent downhill (= 1.0) |
| $f_{\uparrow}$ | `frac_from_uphill` | — | Fraction of stored water delivered per timestep (= 0.5) |

## Code Structure

| File | Subroutine | Role |
| :--- | :--------- | :--- |
| `src/main/elm_varctl.F90` | — | Declares `use_IM2_hillslope_hydrology` flag |
| `src/main/elm_varcon.F90` | — | Defines `frac_to_downhill` and `frac_from_uphill` constants |
| `src/main/controlMod.F90` | `control_init` | Reads flag from namelist; broadcasts via MPI |
| `src/main/initGridCellsMod.F90` | `set_topounit` | Sets `downhill_ti` for each topounit |
| `src/biogeophys/SoilHydrologyMod.F90` | `SurfaceRunoff` | Computes `qflx_from_uphill`; adds to `qflx_top_soil` |
| `src/biogeophys/HydrologyDrainageMod.F90` | `HydrologyDrainage` | Computes `qflx_to_downhill`; updates `from_uphill` state |
| `src/biogeophys/BalanceCheckMod.F90` | `BalanceCheck` | Decrements `from_uphill` state; checks water balance |
