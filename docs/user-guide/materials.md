# Materials

Materials define the thermal properties of the solid (or fluid) layers that surround the flow area. They are the bridge between geometry and temperature: they determine how fast heat moves through walls, how much energy walls can store, and ultimately how fluid temperature evolves along the system.

> **JSON key:** `material` (EN) · `material` (PT) — top-level array

---

## Thermal Conductivity

Conductivity describes how easily heat flows through the material by conduction. High conductivity means the material transmits heat quickly (e.g., steel); low conductivity means the material resists heat flow (e.g., insulation foam).

In practice, conductivity controls the **steady-state temperature gradient** across a layer: higher conductivity → smaller temperature difference across the layer.

> **JSON key:** `conductivity` (EN) · `condutividade` (PT) — unit: W/(m·°C)

## Specific Heat

Specific heat describes how much energy is needed to raise the temperature of one kilogram of the material by one degree. Materials with high specific heat absorb more energy before their temperature changes.

In practice, specific heat (together with density) controls the **thermal inertia** — how long it takes for temperature changes to propagate through the layer during transient events (cooldown, warmup, restart).

> **JSON key:** `specificHeat` (EN) · `calorEspecifico` (PT) — unit: J/(kg·°C)

## Density

Mass density contributes to thermal inertia. Combined with specific heat, it determines the volumetric heat capacity ($\rho \cdot C_p$) — the energy stored per unit volume per degree.

> **JSON key:** `rho` (EN) · `rho` (PT) — unit: kg/m³

## Thermal Diffusivity

The combination $\alpha = k / (\rho \cdot C_p)$ is the thermal diffusivity — it governs how fast temperature fronts propagate through the material. High diffusivity means rapid equilibration; low diffusivity means slow response.

---

## Material Types

Not all layers are solids with user-specified properties. Marlim3 supports four material types:

### Solid Material (Type 0)

Standard solid: steel, concrete, polymer insulation, cement. All three thermal properties (conductivity, specific heat, density) must be specified.

### User-Defined Fluid Layer (Type 1)

A stagnant or quasi-stagnant fluid layer (e.g., a fluid-filled annulus where convection is not explicitly modeled). Properties must be provided by the user.

### Water (Type 2)

An internal water model computes temperature-dependent properties automatically. Only the material type needs to be specified — no manual property input is needed. Use this for completion-fluid annuli, water-filled gaps, or seawater layers.

### Air (Type 3)

An internal air model computes properties automatically. Use this for gas-filled annuli or atmospheric gaps.

> **JSON key:** `type` (EN) · `tipo` (PT)
> Values: `0` = solid, `1` = user fluid, `2` = water, `3` = air

!!! note
    For types 2 and 3, properties are computed internally from correlations. Only `id`, `type`, and optionally `label` need to be specified.

---

## Common Material Properties

| Material | Conductivity [W/(m·°C)] | Specific Heat [J/(kg·°C)] | Density [kg/m³] |
|----------|------------------------|--------------------------|-----------------|
| Carbon steel | 50 | 500 | 7800 |
| Stainless steel | 15 | 500 | 8000 |
| PU insulation (foam) | 0.03 | 1500 | 60 |
| Polypropylene (solid) | 0.22 | 1800 | 900 |
| Cement | 0.6 | 1000 | 500–2000 |
| Concrete weight coating | 1.5 | 880 | 2300 |

---

## Practical Guidance

- **Realistic insulation values:** Over-optimistic conductivity (< 0.02 W/(m·°C)) strongly distorts cooldown predictions. Use manufacturer-rated wet/aged values for subsea applications.
- **Reuse material IDs:** Define materials once in the top-level array and reference them by ID across multiple cross sections. Avoids duplication and inconsistency.
- **Completion fluids:** Use `type: 2` (water) for water-based completion-fluid annuli rather than manually specifying properties.
- **Unit consistency:** All values follow the schema convention (SI-like: W/(m·°C), J/(kg·°C), kg/m³).

---

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
    For cooldown and restart studies, verify that insulation and cement conductivity values reflect actual field conditions (aged, wet, or damaged), not ideal laboratory values.
