# Results

View profiles, trends, and transverse temperature results after running the simulation.

## Overview

The Results tab displays simulation output after execution. It provides interactive visualization of flow variables along the pipeline and over time.

## Running a Simulation

1. Configure all input tabs (Fluids through General)
2. Click **Run** in the sidebar
3. Results appear automatically in this tab upon completion

## Result Types

### Profiles (Spatial)

Spatial distribution of variables along the pipe at a given time:

- **Pressure** [kgf/cm²]
- **Temperature** [°C]
- **Liquid holdup** [-]
- **Velocity** (liquid, gas, mixture) [m/s]
- **Flow pattern** [-]
- **Density** [kg/m³]
- **Measured depth** [m]

Branch schema output blocks related to profiles:

- `productionProfile`
- `serviceProfile`
- `crossProductionProfile`
- `crossServiceProfile`

### Trends (Temporal)

Time evolution of variables at specific locations:

- Pressure at outlet/inlet
- Flow rates
- Temperature evolution
- Holdup transients

Branch schema output blocks related to trends:

- `productionTrend`
- `serviceTrend`
- `crossProductionTrend`
- `crossServiceTrend`

### Transverse Temperature

Radial temperature distribution through the pipe wall layers at selected positions.

This is configured in branch schema through `crossProductionProfile` and `crossServiceProfile`.

## Interacting with Results

- **Select result set** — Choose between different output datasets
- **X-axis** — Typically measured depth for profiles, time for trends
- **Y-axis** — Select one or more variables to plot
- **Time slider** — For transient profiles, select the time instant to display
- **Animation** — Animate profiles over time

## Export

Results can be exported as:

- CSV files (tabular data)
- JSON (raw output)
- Images (plot screenshots)

## Troubleshooting

| Issue | Possible Cause |
|-------|---------------|
| No results shown | Simulation has not been run yet |
| Empty result set | Simulation failed — check terminal output |
| Unexpected values | Verify boundary conditions and fluid properties |

!!! tip
    Use the steady-state solution as a sanity check before running long transient simulations.
