# Materials

Define materials for pipe wall layers (steel, insulation, concrete, etc.).

## Overview

Materials define the thermal properties of each layer in a pipe cross-section. They are referenced by ID in the Cross Sections tab when building the pipe wall structure.

## Material Types

| Type | Description |
|------|-------------|
| **0 — Solid** | Standard solid material (steel, concrete, insulation) |
| **1 — Fluid (user)** | User-defined fluid layer (e.g., annular fluid) |
| **2 — Water** | Seawater (predefined properties with temperature dependence) |
| **3 — Air** | Air layer (predefined properties) |

## Properties

| Property | Unit | Description |
|----------|------|-------------|
| **Conductivity (k)** | W/(m·°C) | Thermal conductivity |
| **Specific Heat (Cp)** | J/(kg·°C) | Specific heat capacity |
| **Density (ρ)** | kg/m³ | Material density |

## Typical Values

| Material | k [W/(m·°C)] | Cp [J/(kg·°C)] | ρ [kg/m³] |
|----------|--------------|-----------------|-----------|
| Carbon steel | 50 | 500 | 7800 |
| Stainless steel | 15 | 500 | 8000 |
| Concrete | 1.5 | 880 | 2300 |
| PU foam insulation | 0.03 | 1500 | 60 |
| Polypropylene | 0.22 | 1800 | 900 |

## JSON Structure

```json
{
  "material": [
    {
      "id": 0,
      "active": true,
      "label": "Steel",
      "type": 0,
      "conductivity": 50.0,
      "specificHeat": 500.0,
      "rho": 7800.0
    }
  ]
}
```

!!! note
    Material IDs are referenced by the cross-section layer definitions. Ensure material IDs match when configuring pipe layers.
