# Time

Set simulation time, initial conditions, and time step controls for transient simulations.

## Overview

The Time tab configures temporal aspects of the simulation: whether it runs as steady-state or transient, the total simulation duration, time-stepping strategy, and initial conditions.

## Basic Settings

| Parameter | Description |
|-----------|-------------|
| **Transient simulation** | `true` for transient; `false` for steady-state only |
| **Final time [s]** | Total simulation duration |
| **Initial fluid ID (prod)** | ID of the fluid initially filling the production system |
| **Initial condition** | How the simulation is initialized |

### Initial Condition Options

| Code | Description |
|------|-------------|
| 0 | User-defined (explicit initial state) |
| 1 | Steady-state solution (compute steady-state first, then go transient) |
| 2 | Snapshot file (`.snp`) — resume from a saved state |
| 3 | Gas-lift discharge procedure |

When using initial condition `3` (gas-lift discharge), branch schema also supports:

- `initialConfig.fluidSalinity`
- `initialConfig.gasLineInterfaceLength`
- `initialConfig.prodLineInterfaceLength`
- `initialConfig.dischargeControl`
- `initialConfig.dischargeParameters` (pressure, flow, latency, and temperature controls)

## Maximum Time Increments

Define time breakpoints and corresponding maximum time steps:

| Column | Unit | Description |
|--------|------|-------------|
| **tempos** | s | Time breakpoint |
| **dtmax** | s | Maximum allowed Δt after this breakpoint |

The solver uses adaptive time-stepping but never exceeds the specified `dtmax`. Smaller steps give better accuracy but longer run times.

**Example:**

| tempos [s] | dtmax [s] |
|------------|-----------|
| 0 | 0.1 |
| 10 | 1.0 |
| 100 | 5.0 |

This allows small steps during initial transients and larger steps once the system stabilizes.

## Snapshot Recording

Specify simulation times at which to save the full system state to a `.snp` file. These snapshots can later be used as initial conditions for subsequent simulations.

## Segregation Control

Time-dependent control of phase segregation in the pipe:

| Column | Description |
|--------|-------------|
| **tempoSegrega [s]** | Time at which segregation setting changes |
| **segrega [0/1]** | `0` = no segregation; `1` = segregation enabled |

## JSON Structure

```json
{
  "initialConfig": {
    "transient": true,
    "initialCondition": 1,
    "initialFluidId": 0
  },
  "time": {
    "finalTime": 3600,
    "times": [0, 10, 100],
    "maxDT": [0.1, 1.0, 5.0],
    "saveSnapshot": [1800, 3600],
    "segregationTime": [0, 600],
    "segregation": [0, 1]
  }
}
```

!!! tip
    Use `initialCondition: 1` (steady-state first) for most production simulations. This ensures the transient starts from a physically consistent state.

!!! note
    Branch schema also includes `initialConfig.snapshotFile` for `initialCondition: 2`, and `initialConfig.steadyStateOrder` for steady-state numerical order selection.
