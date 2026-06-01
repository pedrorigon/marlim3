# Rock Formation

The formation surrounding a wellbore or buried pipeline acts as a large thermal reservoir. Heat diffuses into (or from) the rock during production, and this exchange controls long-term thermal behavior: cooldown rates during shutdowns, warmup during restarts, and steady-state temperature profiles in deep wells.

> **JSON key:** `initialConfig.formation` (EN) · `configuracaoInicial.formacao` (PT)

---

## Prior Production Time

The formation around a producing well is not at virgin temperature — years of production have heated a radial zone around the wellbore. The **production time** estimates the radius of this thermally disturbed region. A wider heated zone means the near-wellbore formation is warmer and cooldown starts from a higher baseline.

- **Short production time** (new well, 1–30 days): Small heated zone → rapid initial cooldown.
- **Long production time** (mature well, 365+ days): Wide heated zone → more thermal inertia → slower cooldown.

> **JSON key:** `productionTime` (EN) · `tempoProducao` (PT) — unit: days

---

## Formation Thermal Properties

Each formation zone is characterized by three properties that together determine heat conduction and storage:

### Thermal Conductivity

How easily heat flows through the rock. Higher conductivity means heat dissipates faster from the wellbore into the far field, leading to more heat loss during steady production but also faster thermal equilibration.

> **JSON key:** `conductivity` (EN) · `condutividade` (PT) — unit: W/(m·°C)

### Specific Heat

How much energy the rock can absorb per unit mass per degree of temperature change. Higher specific heat means the formation stores more energy and responds more slowly to temperature changes.

> **JSON key:** `specificHeat` (EN) · `calorEspecifico` (PT) — unit: J/(kg·°C)

### Density

Rock mass density. Together with specific heat, determines the volumetric heat capacity ($\rho \cdot C_p$) — the energy stored per unit volume per degree.

> **JSON key:** `density` (EN) · `densidade` (PT) — unit: kg/m³

### Thermal Diffusivity

The combination $\alpha = k / (\rho \cdot C_p)$ is the thermal diffusivity of the formation. It governs how fast temperature fronts propagate into the rock:

- High diffusivity (e.g., salt): temperature changes penetrate deep, quickly.
- Low diffusivity (e.g., shale): temperature changes are confined near the wellbore.

---

## Multiple Formation Zones

Different lithologies along a well or pipeline are represented by multiple formation entries, each with its own ID. Pipe segments reference the appropriate formation via their formation ID, allowing property variation with depth or position.

> **JSON key:** `properties` (EN) · `propriedades` (PT) — array inside formation object
> Each element has: `id`, `conductivity`, `specificHeat`, `density`

---

## Common Rock Properties

| Rock type | Conductivity [W/(m·°C)] | Specific Heat [J/(kg·°C)] | Density [kg/m³] |
|-----------|------------------------|--------------------------|-----------------|
| Sandstone | 2.0–4.0 | 800–1000 | 2200–2600 |
| Shale | 1.0–2.5 | 800–1000 | 2300–2700 |
| Limestone | 2.5–3.5 | 850–950 | 2500–2700 |
| Salt | 5.0–6.0 | 850 | 2200 |

---

## Practical Guidance

- **Depth consistency:** Assign different formation IDs to pipe segments crossing different rock layers. A sandstone reservoir section and a shale overburden section should have distinct properties.
- **Avoid underestimation:** Unrealistically low conductivity exaggerates insulation from the formation, leading to optimistic (slow) cooldown predictions.
- **Production-time sensitivity:** For cooldown studies, production time interacts strongly with formation conductivity. Vary both together to understand sensitivity.
- **Pipeline vs. well:** For buried pipelines, formation represents the surrounding soil. For exposed subsea lines, formation is not used — set the pipe environment to water or atmosphere instead.

---

## Example: Single Formation

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

## Example: Multiple Formation Zones

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
    For cooldown sensitivity studies, vary production time and formation conductivity together — they interact strongly in determining early thermal transient rates.
