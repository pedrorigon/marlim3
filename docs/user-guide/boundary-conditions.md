# Boundary Conditions

Set pressure and flow boundary conditions at system inlet and outlet.

## Overview

Boundary conditions define the system constraints at the pipe extremities. Different combinations are used depending on whether the system is a production well, injection well, or flowline segment.

In `schema_branch.json`, if no inlet BC object is provided (`pressureCondition` or `flowRatePressureCondition`), the inlet is treated as closed and inflow must come from sources (IPR, mass source, liquid source, gas source, etc.).

## Downstream Pressure (Separator)

The most common BC for production systems — sets a pressure at the downstream end (separator/platform):

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Active** | — | Enable/disable this BC |
| **Pressure schedule** | kgf/cm² | Time-varying separator pressure |

The pressure can be constant or vary with time for transient simulations (e.g., slug catcher pressure oscillations).

## Downstream Flow

A downstream flow or IPR boundary condition is configured by adding the corresponding accessory (mass source, liquid source, gas source, or IPR) at a small but non-zero measured depth in the **Accessories** tab.

!!! info
    No additional configuration is needed here — once the accessory is placed, it automatically acts as a boundary condition.

## Upstream Pressure

Pressure condition at the pipe upstream end:

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Active** | — | Enable/disable this BC |
| **Pressure** | kgf/cm² | Upstream pressure |
| **Temperature** | °C | Upstream fluid temperature |
| **Fluid quality** | — | Gas mass fraction at inlet |
| **Beta ratio** | — | Water cut ratio |

## Upstream Flow + Pressure

Fully determines the system from the upstream side (steady-state only):

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Pressure** | kgf/cm² | Upstream pressure |
| **Temperature** | °C | Upstream fluid temperature |
| **Mass flow rate** | kg/s | Total mass flow rate |
| **Beta ratio** | — | Water cut ratio |

## Gas Injection (Service Line)

Boundary condition for gas-lift injection on the service line:

| Parameter | Unit | Description |
|-----------|------|-------------|
| **BC type** | — | `0` — Injection pressure; `1` — Injection flow rate |
| **Temperature** | °C | Gas injection temperature |
| **Gas flow rate** | sm³/d | Standard gas flow rate |
| **Injection pressure** | kgf/cm² | Gas injection pressure |
| **Guess initial flow rate** | — | Use first element as initial guess (helps convergence) |

## Injection Well

For injection well simulations:

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Fluid type** | — | `0` User fluid, `1` Water, `2` CO₂ flash, `3` CO₂ compositional |
| **Salinity** | g/(kg water) | Required when fluid type = 1 |
| **PVTSim file** | — | `.tab`/`.ctm` for CO₂ models |
| **BC combination** | — | Specifies which pair of variables is fixed |
| **Injection temperature** | °C | Surface injection temperature |
| **Flow rate** | sm³/d | Injection flow rate |
| **Injection pressure** | kgf/cm² | Surface injection pressure |
| **Bottom-hole pressure** | kgf/cm² | Downhole pressure constraint |

### BC Combinations

| Code | Description |
|------|-------------|
| 0 | Flow + IPR |
| 1 | Injection pressure + IPR |
| 2 | Bottom pressure + IPR |
| 3 | Injection + bottom pressure |
| 4 | Flow + injection pressure |
| 5 | Flow + bottom pressure |

## Choke Boundary Objects

Branch schema also defines dedicated choke objects for boundary restrictions:

| Object | Description |
|--------|-------------|
| **surfaceChoke** | Production-line choke at surface/outlet with opening and coefficient controls |
| **injectionChoke** | Injection-line choke object for injector workflows |

## JSON Structure

```json
{
  "separator": {
    "active": true,
    "time": [0],
    "pressure": [50.0]
  },
  "initialConfig": {
    "pressureCondition": {
      "active": true,
      "time": [0],
      "pressure": [200.0],
      "temperature": [80.0]
    }
  },
  "gasInj": {
    "active": true,
    "bcType": 0,
    "time": [0],
    "temperature": [40.0],
    "gasFlowRate": [100000.0],
    "injectionPressures": [150.0]
  }
}
```
