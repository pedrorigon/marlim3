# Fluids

This section explains how Marlim3 represents real fluids and phase behavior.

## Why Fluid Modeling Matters

In multiphase simulation, nearly every governing equation depends on fluid properties: density, viscosity, enthalpy, compressibility, and phase split. Inaccurate or inconsistent fluid data propagates into pressure and temperature predictions and can dominate overall model error.

## Modeling Families

Marlim3 supports three thermodynamic model families, selected via `initialConfig`:

| Model family | Key flags | Typical use | Main tradeoff |
|--------------|-----------|-------------|---------------|
| Black-oil | `flashTableFluidModel = false`, `compositionalFluidModel = false` | Fast engineering studies, screening | Simplified composition physics |
| Compositional | `compositionalFluidModel = true` | Rich gas, condensate, complex PVT | Higher computational cost |
| Flash-table | `flashTableFluidModel = true` | Precomputed lookup from PVTSim `.tab`/`.ctm` files | Accurate within table range; depends on table quality |

!!! note
    In transient simulations, black-oil is the lightest option. Consider `pressureTable = true` to pre-build compressibility tables and reduce runtime further.

## Global Fluid Configuration (`initialConfig`)

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `flashTableFluidModel` | bool | `false` | Use PVT table lookup as primary property source |
| `compositionalFluidModel` | bool | `false` | Use compositional treatment (only if flash-table is off) |
| `dynamicTableModel` | bool | `false` | Build local table from compositional run after initial BO pass (steady-state networks) |
| `fluidType` | int | `0` | `0`: liquid-dominated (black-oil mass-transfer logic); `1`: gas-dominated (quality-based mass-transfer) |
| `cpModel` | int | `0` | Specific-heat model: `0` = black-oil; `1` = PVT table interpolation |
| `jtlModel` | int | `0` | `1`: use d(ρ_l)/dT from PVTSim table even in black-oil mode |
| `pvtFile` | string | — | Name of `.tab` or `.ctm` PVT file |
| `pressureTable` | bool | `false` | Pre-build compressibility table for production fluids |
| `gasTable` | bool | `false` | Pre-build compressibility table for service-line gas |
| `latentHeat` | bool | `false` | Use PVT table enthalpy (requires `pvtFile`) |
| `freeGasDensityCorrectionBO` | bool | `false` | Correct in-situ gas density for pressure-dependent lighter-hydrocarbon release |

## Production Fluid Definition (`productionFluid`)

Each element in the `productionFluid` array defines one fluid that can be referenced by sources and initial conditions via its `id`.

### Minimum Required Fields

| Field | Unit | Description |
|-------|------|-------------|
| `id` | — | Integer identifier (referenced by sources, `initialFluidId`, etc.) |
| `api` | °API | Oil API gravity (black-oil models) |
| `gasDensity` | — | Gas relative density (gas/air at standard conditions) |
| `gor` | Sm³/Sm³ | Gas-oil ratio |
| `bsw` | m³/m³ | Basic Sediment and Water fraction (default `0`) |

### Viscosity Model Selectors

Black-oil viscosity is built from three layered correlations:

| Field | Default | Options |
|-------|---------|---------|
| `deadOilModel` | `3` | `0`: ASTM; `1`: Beggs & Robinson; `2`: Modified B&R; `3`: Glaso; `4`: Kartoatmodjo-Schmidt; `5`: Petrosky-Farshad; `6`: Beal; `7`: user temperature-viscosity table |
| `liveOilModel` | `0` | `0`: Beggs & Robinson; `1`: Kartoatmodjo-Schmidt; `2`: Petrosky-Farshad |
| `undersaturatedOilModel` | `0` | `0`: Vasquez & Beggs; `1`: Kartoatmodjo-Schmidt; `2`: Petrosky-Farshad; `3`: Beal; `4`: Khan |

When `deadOilModel = 0` (ASTM), provide `temp1`, `visc1`, `temp2`, `visc2`.
When `deadOilModel = 7` (table), provide `deadOilTemp` and `deadOilVisc` arrays.

### Solution-Gas Ratio (RS) Correlations

| `srBpModel` | Correlation |
|-------------|-------------|
| `0` | Vazquez & Beggs (default) |
| `1` | Lasater |
| `2` | Standing |
| `3` | Glaso |
| `4` | Livia Fulchignoni |

### Emulsion Options

| `emulsionType` | Model | Notes |
|----------------|-------|-------|
| `0` | Water-quality-weighted mixture | Default |
| `1` | Weak Woelflin | — |
| `2` | Medium Woelflin | — |
| `3` | Strong Woelflin | — |
| `4` | Exponential | Requires `emulsionCoefA`, `emulsionCoefB` |
| `5` | Pal-Rhodes | Requires `phi100` |
| `6` | User BSW-multiplier table | Requires `bswVec`, `emulVec` |
| `7` | Oil viscosity below saturation BSW | — |

### Critical-Properties Correlation

| `criticalCorrelation` | Name | Notes |
|----------------------|------|-------|
| `0` | Standard | — |
| `1` | Brown et al | Better for CO₂-rich fluids |
| `2` | Piper et al | Better for CO₂-rich fluids |

### Compositional Options

| Field | Description |
|-------|-------------|
| `userMolarFraction` | If `true`, uses `molarFraction` array; if `false`, reads from `.ctm` file |
| `molarFraction` | Pseudocomponent molar fractions (same order as in `.ctm` file) |
| `userGORComp` | If `true`, adjusts molar fractions to match the provided `gor` value |

## Complementary Fluid (`complementaryFluid`)

A third liquid phase (besides hydrocarbons and water) that travels with the stream without dissolving in gas or slipping relative to the mixture. Examples include:

- Glycol or ethanol cushions (hydrate prevention before master valve opening)
- Chemical inhibitors
- Completion fluids (brine)

| `complementaryFluidType` | Behavior |
|--------------------------|----------|
| `0` | Generic user-defined fluid (requires density, compressibility, viscosity, etc.) |
| `1` | Water-based (only `salinity` is required) |
| `2` | Generic fluid with internal friction-reducer treatment |

Key properties for type `0`:

| Field | Unit |
|-------|------|
| `density` | kg/m³ |
| `compressibility` | 1/Pa |
| `thermalExpansivity` | 1/K |
| `surfaceTension` | N/m |
| `specificHeat` | J/(kg·K) |
| `conductivity` | W/(m·K) |
| `temp1`, `visc1`, `temp2`, `visc2` | °C, cP (ASTM viscosity) |

## Gas Fluid (`gasFluid`)

Used for dry-gas definitions in service/injection contexts (when `initialConfig.gasLine = true`):

| Field | Description |
|-------|-------------|
| `gasDensity` | Gas relative density (gas/air at ambient conditions) |
| `CO2Fraction` | CO₂ molar fraction (default `0`) |
| `criticalCorrelation` | `1`: Brown et al; `2`: Piper et al |
| `useFlashTable` | Use flash table for service-line gas properties |
| `userMolarFraction` | Provide gas composition explicitly |
| `molarFraction` | Pseudocomponent molar fractions |

## Precomputed Property Grids (`compTable`)

When `pressureTable` or `gasTable` is enabled, `compTable` defines the P-T grid for interpolation:

| Field | Unit | Purpose |
|-------|------|---------|
| `numPoints` | — | Number of grid points in each dimension |
| `minPressure` | kgf/cm² | Lower pressure bound |
| `maxPressure` | kgf/cm² | Upper pressure bound |
| `minTemperature` | °C | Lower temperature bound |
| `maxTemperature` | °C | Upper temperature bound |

!!! warning
    Choose ranges that fully cover the expected operating envelope. Extrapolation outside the table is unreliable.

## Example JSON

```json
{
  "initialConfig": {
    "flashTableFluidModel": false,
    "compositionalFluidModel": false,
    "fluidType": 0,
    "cpModel": 0,
    "jtlModel": 0,
    "pressureTable": true,
    "pvtFile": "PVT_CASE.ctm"
  },
  "productionFluid": [
    {
      "id": 0,
      "active": true,
      "api": 28.0,
      "gasDensity": 0.65,
      "gor": 150.0,
      "bsw": 0.10,
      "deadOilModel": 3,
      "liveOilModel": 0,
      "criticalCorrelation": 1,
      "emulsionType": 0
    }
  ],
  "compTable": {
    "active": true,
    "numPoints": 50,
    "minPressure": 1.0,
    "maxPressure": 500.0,
    "minTemperature": 4.0,
    "maxTemperature": 120.0
  }
}
```

!!! tip
    Start with one fluid and validate a simple steady-state case before adding advanced options (compositional corrections, emulsions, multiple fluids) incrementally.
