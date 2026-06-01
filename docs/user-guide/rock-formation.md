# Rock Formation

Define rock lithology properties for heat exchange with the surrounding formation.

## Overview

The Rock Formation tab specifies the thermal properties of the geological formation surrounding the pipeline or wellbore. These properties control the heat exchange between the flowing fluid and the environment.

## Production Time

| Parameter | Unit | Description |
|-----------|------|-------------|
| **Production time** | days | Duration of production prior to the simulation. Determines the thermal penetration radius into the formation |

!!! info
    The production time is critical for transient heat transfer. A longer production time means the formation has been heated (or cooled) over a larger radius, affecting the current heat transfer rate.

## Rock Properties

Each rock definition provides:

| Property | Unit | Description |
|----------|------|-------------|
| **Conductivity (k)** | W/(m·°C) | Thermal conductivity of the rock |
| **Specific Heat (Cp)** | J/(kg·°C) | Specific heat capacity |
| **Density (ρ)** | kg/m³ | Rock density |

## Typical Values

| Rock Type | k [W/(m·°C)] | Cp [J/(kg·°C)] | ρ [kg/m³] |
|-----------|--------------|-----------------|-----------|
| Sandstone | 2.0–4.0 | 800–900 | 2200–2600 |
| Shale | 1.5–2.5 | 800–1000 | 2300–2700 |
| Limestone | 2.5–3.5 | 800–900 | 2500–2700 |
| Salt | 5.0–6.0 | 850 | 2200 |

## JSON Structure

```json
{
  "initialConfig": {
    "formation": {
      "productionTime": 365,
      "properties": [
        {
          "id": 0,
          "label": "Sandstone",
          "conductivity": 2.5,
          "specificHeat": 850.0,
          "density": 2500.0
        }
      ]
    }
  }
}
```

!!! note
    Rock formation IDs are referenced by pipe segments to define which formation surrounds each section of the pipeline.
