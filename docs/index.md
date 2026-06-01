# Marlim3 Documentation

Marlim3 is a 1D multiphase flow simulator developed by Petrobras for production and injection systems.

This documentation is organized around simulation concepts, so it can be used regardless of your workflow: direct JSON authoring, Python API scripting, command-line execution, or the Streamlit GUI.

## What Marlim3 Solves

Marlim3 simulates coupled fluid-flow and heat-transfer problems in well/pipeline systems:

- **Steady-state and transient** multiphase flow (oil + gas + water)
- **Thermal coupling** between fluid, pipe wall layers, and surrounding formation
- **Artificial lift**: gas-lift valves, ESPs, volumetric pumps
- **Flow control**: valves, chokes, check valves, pigging
- **Production and injection** boundary scenarios
- **Network modeling** with multiple branches, sources, and convergence logic
- **Compositional and black-oil** thermodynamic frameworks
- **Natural convection** in pipe cross-sections during shutdowns

## Concept Map

| Concept | What You Learn |
|---------|----------------|
| [Fluids](user-guide/fluids.md) | Thermodynamic models, black-oil vs. compositional, PVT tables |
| [Materials](user-guide/materials.md) | Thermal properties of pipe wall layers and media |
| [Cross Sections](user-guide/cross-sections.md) | Radial geometry, layering, and friction surface |
| [Rock Formation](user-guide/rock-formation.md) | Formation thermal properties and production-time effects |
| [Pipes](user-guide/pipes.md) | Axial geometry, inclination, discretization, and thermal environment |
| [Accessories](user-guide/accessories.md) | Sources, valves, pumps, ESPs, and pigging |
| [Boundary Conditions](user-guide/boundary-conditions.md) | Inlet/outlet closure: pressure, flow rate, and combinations |
| [Time](user-guide/time.md) | Initialization strategies, time-step scheduling, and snapshots |
| [General](user-guide/general.md) | Global physics, correlations, and numerical controls |
| [Results](user-guide/results.md) | Output configuration, interpretation, and diagnostics |

## Reference

| Reference | Content |
|-----------|---------|
| [JSON Schema Reference](reference/json-schema.md) | Complete input-file structure, field catalog, and units |
| [Tutorials](tutorials/index.md) | Step-by-step Jupyter notebook workflows |
| [Bilingual Support](dev-guide/translations.md) | English ↔ Portuguese key translation system |

## Suggested Reading Order

For new users building their first model:

1. **Fluids** — Define what flows through the system
2. **Materials and Cross Sections** — Define the pipe wall structure
3. **Pipes** — Define the axial geometry and mesh
4. **Boundary Conditions** — Close the mathematical problem
5. **Time and General** — Set simulation mode and numerical options
6. **Results** — Configure and interpret outputs

## Build This Documentation

```bash
uv sync --group docs
uv run mkdocs serve
```
