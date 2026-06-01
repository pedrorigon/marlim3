# Marlim3 Documentation

Welcome to the **Marlim3** documentation — a 1D multiphase flow simulator developed by Petrobras for modeling production and injection systems.

## Overview

Marlim3 simulates steady-state and transient multiphase flow in pipelines, risers, and wellbores. It supports:

- **Black-oil and compositional** fluid modeling
- **Steady-state** shooting-method solvers
- **Transient** drift-flux solvers
- **Heat transfer** with radial conduction and natural convection
- **Gas-lift** and **ESP** artificial lift systems
- **Network simulation** of multi-branch systems
- **Sensitivity analysis** parameter sweeps

## User Guide

The user guide is organized to match the input tabs in the Marlim3 GUI:

| Section | Description |
|---------|-------------|
| [Fluids](user-guide/fluids.md) | Fluid property definitions (PVT, compositions, correlations) |
| [Materials](user-guide/materials.md) | Pipe and insulation material thermal properties |
| [Cross Sections](user-guide/cross-sections.md) | Pipe layer geometry and composition |
| [Rock Formation](user-guide/rock-formation.md) | Surrounding formation thermal properties |
| [Pipes](user-guide/pipes.md) | Pipeline geometry, discretization, and elevation profile |
| [Accessories](user-guide/accessories.md) | Sources, valves, pumps, and pigs |
| [Boundary Conditions](user-guide/boundary-conditions.md) | Pressure and flow boundary specifications |
| [Time](user-guide/time.md) | Transient simulation time-stepping and events |
| [General](user-guide/general.md) | Global settings, correlations, and numerical parameters |
| [Results](user-guide/results.md) | Output visualization and post-processing |

## Getting Started

```bash
# Install the package
pip install -e .

# Run the GUI
streamlit run gui/app.py

# Run a simulation from CLI
marlim3 run input.json
```

## Building the Documentation

```bash
pip install mkdocs-material
mkdocs serve     # Local preview at http://127.0.0.1:8000
mkdocs build     # Build static site to site/
```
