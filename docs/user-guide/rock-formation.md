# Rock Formation

Formation properties control long-term thermal interactions between the wellbore/pipeline and surrounding geological media.

## Why Formation Matters

Transient thermal response extends beyond the pipe wall. Heat diffuses into surrounding rock, which acts as a large thermal reservoir. This determines:

- **Cooldown rate** during shutdowns (formation supplies or absorbs heat)
- **Warmup behavior** during restarts
- **Steady-state temperature profiles** in deep wells (geothermal gradient coupling)

The formation model in Marlim3 uses an analytical/numerical approach based on prior production time to estimate the thermally disturbed radius around the wellbore.

## Configuration

Formation properties are defined inside `initialConfig.formation`:

```json
"initialConfig": {
  "formation": {
    "productionTime": 365,
    "properties": [...]
  }
}
```

Each pipe segment references a formation entry via `formationId`.

## Key Parameters

| Field | Unit | Meaning | Impact |
|-------|------|---------|--------|
| `productionTime` | days | Prior continuous operating duration | Longer time → wider heated zone → milder initial cooldown rates |
| `conductivity` | W/(m·°C) | Heat propagation speed in rock | Higher → faster heat dissipation to formation |
| `specificHeat` | J/(kg·°C) | Thermal storage capacity per unit mass | Higher → more thermal inertia |
| `density` | kg/m³ | Rock mass density | Higher → more thermal inertia |

**Thermal diffusivity** $\alpha = k / (\rho \cdot C_p)$ governs how fast the thermal front advances into the formation.

## Common Rock Properties

| Rock type | `conductivity` | `specificHeat` | `density` |
|-----------|---------------|---------------|-----------|
| Sandstone | 2.0–4.0 | 800–1000 | 2200–2600 |
| Shale | 1.0–2.5 | 800–1000 | 2300–2700 |
| Limestone | 2.5–3.5 | 850–950 | 2500–2700 |
| Salt | 5.0–6.0 | 850 | 2200 |

## Modeling Guidance

- **Depth consistency:** Use lithology-appropriate properties at each depth zone. Assign different `formationId` values to pipe segments crossing different rock layers.
- **Avoid underestimation:** Unrealistically low conductivity exaggerates thermal insulation from the formation, leading to optimistic cooldown predictions.
- **Production time calibration:** For new wells, use a short `productionTime` (e.g., 1–30 days). For mature wells, use actual production history (e.g., 365+ days). This significantly affects early transient thermal behavior.
- **ID synchronization:** Each `productionPipe` segment with `environment = 0` (buried/formation) should have a valid `formationId` pointing to an entry in the `properties` array.

## Example JSON

### Single Formation

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

### Multiple Formation Zones

```json
{
  "initialConfig": {
    "formation": {
      "productionTime": 180,
      "properties": [
        {
          "id": 0,
          "label": "Reservoir sandstone",
          "conductivity": 3.0,
          "specificHeat": 900.0,
          "density": 2400.0
        },
        {
          "id": 1,
          "label": "Overburden shale",
          "conductivity": 1.5,
          "specificHeat": 950.0,
          "density": 2500.0
        }
      ]
    }
  }
}
```

!!! tip
    For cooldown sensitivity studies, vary `productionTime` and formation `conductivity` together — they interact strongly in determining early thermal transient rates.
