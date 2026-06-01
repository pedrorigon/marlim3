# Fluids

Configure production fluids, PVT models, and complementary/gas fluid properties.

## Overview

The Fluids tab defines the thermodynamic fluid model used in the simulation. Marlim3 supports multiple PVT modeling approaches that can be combined depending on the system complexity.

## Fluid Model

| Parameter | Description |
|-----------|-------------|
| **Compositional model** | Uses compositional models (only when flash table is disabled) |
| **Flash table model** | Uses PVT flash table instead of black oil or compositional |
| **Dynamic table model** | Post-simulation table from compositional model, reusing BO P/T mapping (steady-state only). Can significantly speed up network simulations |
| **Latent heat from PVTSim** | Uses PVTSim table interpolation for enthalpy (requires PVTSim file) |
| **Fluid type** | `0` — Liquid dominated (black oil); `1` — Gas dominated (title ratio) |
| **Cp model** | `0` — Black oil; `1` — PVTSim table interpolation |
| **JTL model** | `0` — Standard; `1` — d(ρ_l)/dT from PVTSim |
| **PVTSim file** | `.tab` or `.ctm` file with PVT data |
| **Table P** | Compressibility table |
| **Table G** | Gas compressibility table |

## Fluid Definitions

Each fluid entry represents a complete set of properties for a given fluid. You can define multiple fluids (e.g., different wells producing different oils).

In branch/network schema files, this list is named `productionFluid` (instead of `fluid`).

### Key Properties

- **API gravity** — Oil API density
- **Gas specific gravity** — Gas relative density (air = 1)
- **Water cut** — Water fraction
- **GOR** — Gas-Oil Ratio [sm³/sm³]
- **BSW** — Basic Sediment & Water [%]
- **Bubble point** — Saturation pressure [kgf/cm²]
- **Salinity** — Water salinity [g/kg]

### Extended Production-Fluid Options (Branch Schema)

`schema_branch.json` includes additional production-fluid controls, especially for viscosity and compositional behavior:

| Parameter | Description |
|-----------|-------------|
| **emulsionType** | Emulsion viscosity model selector |
| **srBpModel** | RS correlation model |
| **deadOilModel** | Dead-oil viscosity model |
| **liveOilModel** | Live-oil viscosity model |
| **undersaturatedOilModel** | Undersaturated-oil viscosity model |
| **blackOilViscModel** | In flash-table mode, choose table viscosity or black-oil correlations |
| **blackOilWaterModel** | In flash-table mode, choose Joule-Thomson treatment |
| **userMolarFraction** | If true, read composition from input instead of .ctm |
| **molarFraction** | User-defined pseudocomponent molar fractions |
| **userGORComp** | Correct compositional fractions to match user GOR |

## Complementary Fluid

A secondary fluid used for specific scenarios (e.g., injection fluid, kill fluid). Configured with its own density, viscosity, and thermal properties.

In branch schema this object is `complementaryFluid`, with `complementaryFluidType`:

| Type | Meaning |
|------|---------|
| **0** | Standard complementary fluid (full property set required) |
| **1** | Water (salinity-based setup) |
| **2** | Similar to type 0, with internal friction-reducer handling |

## Gas Fluid

Properties specific to the gas phase when using gas-dominated systems or gas-lift operations:

- Gas density at standard conditions
- Gas viscosity correlations
- Gas thermal properties

In branch schema, `gasFluid` also includes:

| Parameter | Description |
|-----------|-------------|
| **criticalCorrelation** | Gas pseudo-critical P/T model selector |
| **useFlashTable** | Use flash table for service-line gas properties |
| **userMolarFraction** | Use user gas composition rather than .ctm |
| **molarFraction** | Gas pseudocomponent molar fractions |

## Compressibility Table Configuration

When `initialConfig.pressureTable` or `initialConfig.gasTable` is enabled, branch schema defines table ranges under `compTable`:

- `numPoints`
- `minPressure`, `maxPressure`
- `minTemperature`, `maxTemperature`

## JSON Structure

```json
{
  "initialConfig": {
    "compositionalFluidModel": false,
    "flashTableFluidModel": false,
    "dynamicTableModel": false,
    "latentHeat": false,
    "fluidType": 0,
    "cpModel": 0,
    "jtlModel": 0,
    "pvtFile": "example.tab"
  },
  "productionFluid": [
    {
      "id": 0,
      "active": true,
      "api": 28.0,
      "gasDensity": 0.65,
      "gor": 150.0,
      "bsw": 0.0,
      "srBpModel": 0,
      "deadOilModel": 3
    }
  ],
  "complementaryFluid": {
    "active": false,
    "complementaryFluidType": 0
  },
  "gasFluid": {
    "active": false,
    "gasDensity": 0.65,
    "criticalCorrelation": 1
  },
  "compTable": {
    "active": false,
    "numPoints": 40,
    "minPressure": 1.0,
    "maxPressure": 600.0,
    "minTemperature": 4.0,
    "maxTemperature": 180.0
  }
}
```

!!! tip
    For compositional simulations, ensure you have uploaded the corresponding `.tab` or `.ctm` PVTSim file in the sidebar before selecting it.
