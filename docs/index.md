<h1 align="center">
	<img src="img/logo_marlim3.svg" alt="Marlim3 logo" width="320"/>
</h1>

# Marlim3 Documentation

Marlim3 is a 1D multiphase-flow simulator developed by Petrobras for production and injection systems in oil and gas applications.

This documentation is concept-oriented and workflow-agnostic: you can use it for JSON input authoring, Python scripting, command-line execution, or GUI-driven studies.

## Core Capabilities

Marlim3 supports steady-state and transient simulations for:

- **Production systems**
- **Injection systems** (water or gas)
- **Network systems**
	- Production networks
	- Injection networks
	- Gas-lift loops
- **Artificial lift models**
	- Gas-lift valves
	- ESP pumps
	- Volumetric pumps

## Advanced Modeling Scope

- **Thermal diffusion**: coupled 2D/3D diffusion models integrated with 1D flow
- **Natural convection**: 2D cross-section analysis for shutdown and segregation scenarios
- **Fluid frameworks**: black-oil, flash-table, and compositional models
- **Near-wellbore effects**: support for radial and 2D source-style modeling strategies
- **Operational events**: valves, chokes, pressure-drop devices, pigging, and control schedules

## How You Can Use Marlim3

- **Python package** for automation, integration, and parametric studies
- **Command-line executable** for direct simulation workflows
- **Streamlit GUI** for interactive model setup and result inspection
- **Bilingual keys** in English and Portuguese for input and API usage

## User Guide Map

| Concept | What You Learn |
|---------|----------------|
| [Fluids](user-guide/fluids.md) | Fluid-model families, PVT files, and property applicability by model |
| [Materials](user-guide/materials.md) | Thermal-property definitions for wall layers and special material types |
| [Cross Sections](user-guide/cross-sections.md) | Radial geometry, annular behavior, and layer construction |
| [Rock Formation](user-guide/rock-formation.md) | Formation thermal memory and rock property assignment |
| [Pipes](user-guide/pipes.md) | Segment geometry, discretization, ambient coupling, and thermal links |
| [Accessories](user-guide/accessories.md) | Sources, valves, pumps, pressure devices, and pigging |
| [Boundary Conditions](user-guide/boundary-conditions.md) | Inlet/outlet closure strategies and injection boundary modes |
| [Time](user-guide/time.md) | Initialization strategies, step schedule, segregation windows, snapshots |
| [General](user-guide/general.md) | Global physics switches, correlations, and advanced numerical controls |
| [Results](user-guide/results.md) | Profile/trend/cross outputs and interpretation strategy |

## Reference and Tutorials

| Section | Content |
|---------|---------|
| [JSON Schema Reference](reference/json-schema.md) | Full field catalog, object structure, and units |
| [Tutorials](tutorials/index.md) | Practical notebook-based workflows |
| [Bilingual Support](dev-guide/translations.md) | English ↔ Portuguese key mapping and usage |

## Recommended Reading Path

For a first model:

1. **Fluids**
2. **Materials** and **Cross Sections**
3. **Pipes**
4. **Boundary Conditions** and **Accessories**
5. **Time** and **General**
6. **Results**

## Build This Documentation

```bash
uv sync --group docs
uv run mkdocs serve
```
