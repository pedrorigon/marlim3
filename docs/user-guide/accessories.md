# Accessories

Configure flow sources, valves, pumps, and pig launchers along the pipeline.

## Overview

Accessories are devices or boundary elements placed at specific measured depths along the production pipe. They include flow sources (inlets), valves, artificial lift equipment, and pigging devices.

## Sources

### IPR (Inflow Performance Relationship)

Defines reservoir inflow using a productivity index or IPR curve:

- Reservoir pressure and temperature
- Productivity index (PI)
- IPR model type (linear, Vogel, etc.)
- Position along the pipe (measured depth)

### Liquid Source

Constant or time-varying liquid injection at a point along the pipe.

Branch schema highlights:

| Field | Description |
|-------|-------------|
| **prodFluidId** | Production-fluid ID injected by this source |
| **measuredLength** | Position in production line [m] |
| **time** | Event times [s] |
| **temperature** | Source temperature profile [degC] |
| **beta** | Complementary-fluid fraction profile |
| **liquidFlowRate** | Standard liquid flow profile [sm3/d] |

### Mass Source

Mass flow rate injection at a specified position.

Branch schema highlights:

| Field | Description |
|-------|-------------|
| **thermType** | `0`: equilibrium gas from fluid model, `1`: explicit gas mass flow |
| **totalMassFlowRate** | Total mass-flow profile [kg/s] |
| **complementaryMassFlowRate** | Complementary-fluid mass-flow profile [kg/s] |
| **gasMassFlow** | Gas mass-flow profile [kg/s], valid when `thermType=1` |

### Gas Source

Gas injection at a specified position (e.g., for gas-lift).

Branch schema highlights:

| Field | Description |
|-------|-------------|
| **dry** | `true`: dry gas (uses gasFluid); `false`: rich gas via production fluid |
| **prodFluidId** | Required for rich-gas source (`dry=false`) |
| **gasFlowRate** | Gas flow profile [sm3/d] |
| **temperature** | Source temperature profile [degC] |

### Pressure Source

Fixed pressure boundary at a point along the pipe.

### Radial Pore

Distributed inflow from the formation (radial pore model).

In branch schema this object is `porousRadialSource`.

### 2D Pore

Two-dimensional pore inflow model for complex reservoir coupling.

In branch schema this object is `porous2DSource`.

## Valves

### DHSV / Generic Valve

Downhole safety valve or generic restriction with:

- Opening/closing schedule
- Cv (flow coefficient)
- Position along pipe

Branch schema highlights:

| Field | Description |
|-------|-------------|
| **cvCurve** | `0`: opening as area ratio; `1`: opening via valve-stem displacement curve |
| **opening** | Time-dependent opening profile |
| **cd** | Discharge coefficient (Sachdeva model) |
| **x1, cv1** | Stem-displacement vs Cv curve (when `cvCurve=1`) |

### Gas Lift Valves (VGL)

Gas-lift valve configuration:

- Injection orifice diameter
- Opening/closing pressures
- Position on production and service lines

### Master Valve (Production)

Surface production master valve.

### BCS Valve

Subsurface controlled valve.

## Pumps

### ESP (Electrical Submersible Pump)

- Performance curves (head vs. flow)
- Speed (RPM) schedule
- Pump efficiency
- Position along pipe

### Volumetric Pump

Positive-displacement pump with fixed flow rate.

### Delta Pressure

Simplified pump model as a fixed pressure differential.

## Pig

Pigging device with:

- Launch time and position
- Pig friction/resistance characteristics
- Detection logic

## JSON Structure

```json
{
  "liquidSource": [...],
  "massSource": [...],
  "gasSource": [...],
  "pressureSource": [...],
  "porousRadialSource": [...],
  "porous2DSource": [...],
  "valve": [...],
  "gasLiftSource": [...],
  "esp": [...],
  "volumetricPump": [...],
  "pressureDrop": [...],
  "masterValve": { ... },
  "masterValve2": { ... },
  "pig": [...]
}
```

!!! warning
    Accessory positions must fall within the pipe discretization range. Position `0` is the pipe inlet.
