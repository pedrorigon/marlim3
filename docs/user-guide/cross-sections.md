# Cross Sections

A cross section defines the radial structure of the pipe: the flow area where fluid moves, the inner wall where friction acts, and the concentric layers through which heat is conducted between fluid and environment.

---

## Flow Area and Hydraulic Resistance

The 1D axial equations require a flow area and a friction surface to compute velocities and pressure losses. Two properties define this:

### Internal Diameter

The inner diameter sets the flow cross-sectional area ($A = \pi d^2 / 4$) and the velocity scale. All hydraulic calculations (Reynolds number, friction factor, holdup geometry) depend on this dimension.

> **JSON key:** `innerDiameter` (EN) · `diametroInterno` (PT) — unit: m

### Wall Roughness

Absolute roughness of the inner pipe surface, entering the Colebrook (or Moody) friction-factor correlation. Typical values range from 1.5×10⁻⁶ m (polished tubing) to 4.6×10⁻⁴ m (corroded old pipe).

> **JSON key:** `roughness` (EN) · `rugosidade` (PT) — unit: m

### Annular Flow Geometry

Some well configurations have production occurring in the annulus between tubing and casing rather than inside the tubing. When annular geometry is indicated, the hydraulic calculations use the annular area and equivalent diameter.

> **JSON key:** `annular` (EN) · `anular` (PT) — default `false`

---

## Concentric Layers (Thermal Structure)

Heat transfer between the flowing fluid and the external environment passes through concentric solid (or fluid) layers surrounding the flow area. Each cross section defines these layers from the inner wall outward:

```
┌─────────────── External boundary (formation, seawater, or atmosphere)
│  Layer N  ← cement, rock contact, or concrete coating
│  Layer 2  ← insulation (PU foam, polypropylene, etc.)
│  Layer 1  ← pipe wall (steel)
│  ═════════  Inner surface (roughness applies here)
│  Fluid flow area (innerDiameter)
└─────────────── Pipe center
```

Each layer is defined by its geometric extent, radial discretization, and material reference.

> **JSON key:** `layers` (EN) · `camadas` (PT) — array inside each cross-section object

### Layer Geometric Definition

Two ways to specify layer size:

- **By thickness** — the radial thickness of the layer, building outward from the previous boundary.
- **By outer diameter** — the absolute outer diameter of the layer.

> **JSON key:** `layerMeasurementType` (EN) · `tipoMedicaoCamada` (PT)
> Values: `"THICKNESS"` (EN) / `"ESPESSURA"` (PT) or `"OUTER_DIAMETER"` (EN) / `"DIAMETRO_EXTERNO"` (PT)
>
> Thickness: `thickness` (EN) · `espessura` (PT) — unit: m
> Outer diameter: `outerDiameter` (EN) · `diametroExterno` (PT) — unit: m

### Radial Discretization

Each layer is subdivided into radial nodes for the heat-conduction calculation. More nodes provide better resolution of temperature gradients within the layer, at the cost of additional computation.

**Guidelines:**

- Thin steel walls (< 15 mm): 1 node is sufficient.
- Thick insulation layers (> 30 mm): Use 2–4 nodes to capture transient thermal response during cooldown/warmup.
- Completion-fluid annuli: 2–3 nodes for thermal coupling resolution.

> **JSON key:** `discretization` (EN) · `discretizacao` (PT) — default `1`

### Material Reference

Each layer references a material entry (by integer ID) that provides the thermal properties (conductivity, specific heat, density) for heat-conduction calculations.

> **JSON key:** `materialId` (EN) · `idMaterial` (PT)

---

## Cross-Section Identity

Each cross section has an integer identifier. Pipe segments reference the cross section they use via this ID, allowing multiple pipe segments to share the same radial geometry.

> **JSON key:** `crossSection` (EN) · `secaoTransversal` (PT) — top-level array
> Each element: `id` — integer identifier

---

## Practical Guidance

- **Layer consistency:** Layers must nest concentrically. When using thickness mode, the outer radius of each layer is computed cumulatively from the inner diameter outward.
- **Roughness calibration:** If pressure-drop predictions deviate from field data, roughness is often the first parameter to calibrate (especially in old or scaled pipes).
- **Multiple cross sections:** Use different cross-section IDs for different pipe segments (e.g., one for the riser with thick insulation, another for the subsea flowline with concrete coating).
- **Annular wells:** When annular geometry is active, ensure the layer structure reflects the casing/tubing arrangement correctly.

---

## Example: Simple Pipeline (Steel + Insulation)

```json
{
  "crossSection": [
    {
      "id": 0,
      "active": true,
      "annular": false,
      "innerDiameter": 0.15,
      "roughness": 5e-5,
      "layers": [
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.012,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.050,
          "discretization": 3,
          "materialId": 1
        }
      ]
    }
  ]
}
```

## Example: Well Completion (Tubing + Annular Fluid + Casing + Cement)

```json
{
  "crossSection": [
    {
      "id": 1,
      "active": true,
      "annular": false,
      "innerDiameter": 0.10,
      "roughness": 1.83e-4,
      "layers": [
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.008,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.020,
          "discretization": 2,
          "materialId": 1
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.010,
          "discretization": 1,
          "materialId": 0
        },
        {
          "layerMeasurementType": "THICKNESS",
          "thickness": 0.040,
          "discretization": 2,
          "materialId": 2
        }
      ]
    }
  ]
}
```

Where material IDs map to: `0` = steel, `1` = completion fluid (type 2), `2` = cement.

!!! tip
    Start with 1 node per layer for steady-state validation, then increase to 2–3 in insulation/cement layers for transient cooldown studies where radial resolution matters.

!!! tip
    When in doubt about discretization, start with 1 node per layer for steady-state validation, then increase to 2–3 in insulation/cement layers for transient cooldown studies.
