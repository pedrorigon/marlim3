# JSON Schema Reference

Marlim3 uses JSON input files to define simulations. The full schema is documented here.

## Schema Files

- [Tramo Schema](../schema_tramo.json) — Single pipeline (tramo) input format
- [Branch Schema](../schema_branch.json) — Network branch definition

## Top-Level Structure

```json
{
  "system": "PROD",
  "versaoJSON": "3.0.0",
  "initialConfig": { ... },
  "productionFluid": [ ... ],
  "complementaryFluid": { ... },
  "gasFluid": { ... },
  "compTable": { ... },
  "material": [ ... ],
  "crossSection": [ ... ],
  "productionPipe": [ ... ],
  "servicePipe": [ ... ],
  "liquidSource": [ ... ],
  "massSource": [ ... ],
  "gasSource": [ ... ],
  "pressureSource": [ ... ],
  "porousRadialSource": [ ... ],
  "porous2DSource": [ ... ],
  "gasLiftSource": [ ... ],
  "ipr": [ ... ],
  "valve": [ ... ],
  "esp": [ ... ],
  "volumetricPump": [ ... ],
  "pressureDrop": [ ... ],
  "masterValve": { ... },
  "masterValve2": { ... },
  "pig": [ ... ],
  "separator": { ... },
  "gasInj": { ... },
  "surfaceChoke": { ... },
  "injectionChoke": { ... },
  "injectionWellBC": { ... },
  "productionProfile": { ... },
  "serviceProfile": { ... },
  "productionTrend": { ... },
  "serviceTrend": { ... },
  "crossProductionProfile": { ... },
  "crossServiceProfile": { ... },
  "crossProductionTrend": { ... },
  "crossServiceTrend": { ... },
  "screenConfig": { ... },
  "wax": { ... },
  "severeSlugging": { ... },
  "time": { ... }
}
```

## Tramo vs Branch Naming

Some names differ between schemas and GUI views. Common differences:

| Concept | Tramo/GUI naming | Branch schema naming |
|---------|-------------------|----------------------|
| System enum | `MULTIFASICO`, `INJETOR` | `PROD`, `INJ` |
| Fluid array | `fluid` | `productionFluid` |
| Accessory grouping | Nested convenience view | Top-level arrays/objects |
| Time vector key in several objects | Often presented as `times` in examples | `time` in schema |

## Key Sections

| Section | Description |
|---------|-------------|
| `system` | System type (`PROD` or `INJ` in branch schema) |
| `versaoJSON` | JSON format version |
| `initialConfig` | Global simulation settings, correlations, numerical parameters |
| `productionFluid` | Array of production fluid definitions |
| `complementaryFluid` | Optional third liquid phase (inhibitor/brine/etc.) |
| `gasFluid` | Dry-gas property definition for service line/sources |
| `compTable` | Limits and resolution for precomputed P-T property tables |
| `material` | Array of material thermal properties |
| `crossSection` | Array of pipe cross-section geometries |
| `productionPipe` | Array of production pipe segments |
| `servicePipe` | Array of service line segments |
| `liquidSource`, `massSource`, `gasSource`, `pressureSource`, `ipr` | Source/inflow definitions |
| `valve`, `esp`, `volumetricPump`, `pressureDrop`, `masterValve`, `masterValve2`, `pig` | Accessory/equipment definitions |
| `porousRadialSource`, `porous2DSource`, `gasLiftSource` | Specialized source models |
| `separator` | Downstream pressure BC |
| `gasInj` | Gas injection BC |
| `surfaceChoke`, `injectionChoke`, `injectionWellBC` | Choke and injection-well boundary objects |
| `productionProfile`, `serviceProfile`, `productionTrend`, `serviceTrend` | Longitudinal result-output configuration |
| `crossProductionProfile`, `crossServiceProfile`, `crossProductionTrend`, `crossServiceTrend` | Radial and cross-section output configuration |
| `screenConfig` | Plot and display options |
| `wax`, `severeSlugging` | Optional advanced physical-model blocks |
| `time` | Time-stepping configuration |

## Validation

Input files can be validated against the JSON schema:

```python
import json
import jsonschema

with open("schema_tramo.json") as f:
    schema = json.load(f)

with open("my_input.json") as f:
    data = json.load(f)

jsonschema.validate(data, schema)
```
