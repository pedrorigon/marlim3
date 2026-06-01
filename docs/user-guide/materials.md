# Materials

Materials define the thermal properties of solid layers and media used in heat-transfer calculations between the flowing fluid, pipe wall, and surrounding environment.

## Concept

Marlim3 computes radial heat conduction through concentric layers (defined in `crossSection`). Each layer references a `material` entry that provides thermal conductivity, heat capacity, and density. These properties determine:

- How quickly heat penetrates the wall (thermal diffusivity = k / (ρ·Cp))
- How much energy the wall can store (thermal inertia = ρ·Cp·V)
- The equilibrium temperature profile across layers

## Material Types

| `type` | Meaning | Required properties |
|--------|---------|---------------------|
| `0` | Solid (steel, concrete, polymer insulation, cement) | `conductivity`, `specificHeat`, `rho` |
| `1` | User-defined fluid-like layer (stagnant fluid annulus) | `conductivity`, `specificHeat`, `rho` |
| `2` | Water (properties computed internally from correlations) | None (auto-computed) |
| `3` | Air (properties computed internally from correlations) | None (auto-computed) |

!!! note
    For types `2` and `3`, internal correlations provide temperature-dependent properties. Only `id`, `type`, and optionally `label` need to be specified.

## Key Thermal Properties (Type 0 and 1)

| Property | Unit | Physical interpretation | Typical range |
|----------|------|-------------------------|---------------|
| `conductivity` | W/(m·°C) | Rate of heat flow through the material | Steel: 50; PU foam: 0.03; Cement: 0.6 |
| `specificHeat` | J/(kg·°C) | Energy to raise 1 kg by 1°C | Steel: 500; Insulation: 1500; Cement: 1000 |
| `rho` | kg/m³ | Mass density (thermal mass contribution) | Steel: 7800; PU foam: 60; Cement: 500–2000 |

**Thermal diffusivity** $\alpha = k / (\rho \cdot C_p)$ controls how fast temperature fronts propagate. High-conductivity, low-density materials equilibrate faster.

## Common Material Library

| Material | `conductivity` | `specificHeat` | `rho` | Notes |
|----------|---------------|---------------|-------|-------|
| Carbon steel | 50 | 500 | 7800 | — |
| Stainless steel | 15 | 500 | 8000 | — |
| PU insulation | 0.03 | 1500 | 60 | Sensitive to water ingress |
| Polypropylene | 0.22 | 1800 | 900 | Solid insulation |
| Cement | 0.6 | 1000 | 500–2000 | Varies with formulation |
| Concrete coating | 1.5 | 880 | 2300 | Weight coating |

## Engineering Guidance

- **Consistent units:** All values use SI convention as per schema (W/(m·°C), J/(kg·°C), kg/m³).
- **Realistic insulation:** Over-optimistic insulation conductivity (< 0.02 W/(m·°C)) strongly distorts cooldown predictions. Use manufacturer-rated wet values for subsea applications.
- **Reuse material IDs:** Share material definitions across cross sections to avoid duplication and inconsistency. A single material array serves the entire model.
- **Completion fluids:** Use `type: 2` (water) for completion-fluid annuli rather than manually specifying water properties.

## Example JSON

```json
{
  "material": [
    {
      "id": 0,
      "active": true,
      "label": "Carbon steel",
      "type": 0,
      "conductivity": 50.0,
      "specificHeat": 500.0,
      "rho": 7800.0
    },
    {
      "id": 1,
      "active": true,
      "label": "PU insulation",
      "type": 0,
      "conductivity": 0.03,
      "specificHeat": 1500.0,
      "rho": 60.0
    },
    {
      "id": 2,
      "active": true,
      "label": "Completion fluid",
      "type": 2
    },
    {
      "id": 3,
      "active": true,
      "label": "Cement",
      "type": 0,
      "conductivity": 0.6,
      "specificHeat": 1000.0,
      "rho": 500.0
    }
  ]
}
```

!!! tip
    For cooldown and restart studies, verify that insulation and cement conductivity values reflect actual field conditions (aged, wet, or damaged), not ideal lab values.
