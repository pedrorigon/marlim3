<h1 align="center">
<img src="img/logo_marlim3.svg" alt="drawing" width="300"/>

</h1><br>

`Marlim3` is a 1D multiphase flow simulator developed by Petrobras.
  
## Core capabilities (steady-state and transient)

- **Production wells**
- **Injection wells**: Water or gas injection wells, both single-phase and multiphase
- **Networks**
    - Production networks 
    - Injection networks
    - Gas lift loops
- **Artificial Lift models**: gas lift valves, pumps

## Advanced modeling

- **Natural convection**: 2D solutions for natural convection analysis in confined spaces (single-phase or two-phase), such as pipeline cross-sections during production shutdowns
- **Compositional fluid model library**
- **Near wellbore model**: radial and 2D models to consider phenomena such as water coning
- **Thermal diffusion**: 2D and 3D coupled to the 1D flow model

## Installation

### Option 1: Install via pip

Install `Marlim3` as a Python package:

```bash
pip install marlim3
```

### Option 2: Use the executable directly

You can download the `Marlim3` executable for Linux, Windows or Mac from the Releases section on GitHub. This standalone executable allows you to run simulations directly from the terminal, without the need to install the Python package. Detailed instructions are provided below.

### Option 3: Developer setup (uv)

For development, use [uv](https://docs.astral.sh/uv/) to manage the Python environment and dependencies.

**Step 1 — Install dependencies and the Python package:**

```bash
uv sync --locked --group dev
```

This creates a `.venv` with Python 3.12+, installs all dev tools (pytest, flake8, jupyter, etc.), and installs `marlim3` in editable mode. After this step you can already use `import marlim3` in your scripts:

```bash
uv run python -c "import marlim3; print(marlim3.__version__)"
```

**Step 2 — Build and register the C++/Fortran executable** (required to run simulations):

See [Compilation](#compilation) below. After building, copy `build/Marlim3` to `marlim3/` and run:

```bash
MARLIM3_SKIP_BUILD=1 uv sync --locked
```

## Usage

### Option 1: Python Package

Use `Marlim3` as a Python library in your scripts:

```python
import marlim3

# Your simulation code here
# Example: configure and run simulations programmatically
```

For examples, refer to the tutorials available in the `docs` folder.

### Option 2: Command-Line Executable

Run `Marlim3` directly from the terminal using the compiled executable available in the Releases section on GitHub.

#### Available Commands

There are four simulation types available:

**1. Simple Production System**
```bash
./executable_name -d directory_name -i input_file
```

**2. Simple Injection System**
```bash
./executable_name -d directory_name -i input_file -s INJETOR
```

**3. Flow Network**
```bash
./executable_name -d directory_name -i input_file -s REDE
```

**4. Natural Convection in Cross-Section**
```bash
./executable_name -d directory_name -i input_file -s CONVECNAT
```

#### Command-Line Arguments

- `-d directory_name`: Output directory for simulation results
- `-i input_file`: Input file name (JSON format)
- `-s SIMULATION_TYPE`: Simulation type (INJETOR, REDE, or CONVECNAT)

#### Platform-Specific Notes

**Linux/macOS:**
```bash
./Marlim3 -d ./output -i simulation.json -s REDE
```

**Windows:**
```powershell
Marlim3.exe -d .\output -i simulation.json -s REDE
```

> **Tip:** To export results to the current working directory, use `./` (Linux/macOS) or `.\` (Windows) as the directory name.

## Compilation

Compilation is only necessary if you need to rebuild the executable from source.

### Requirements

- GCC/G++ >= 9.0
- GFortran >= 9.0
- CMake >= 3.16

### Build the executable

The project uses [CMake presets](CMakePresets.json). Available presets:

| Preset | Platform | Description |
|--------|----------|-------------|
| `gcc-release` / `gcc-debug` | Linux / macOS | GCC with partial static linking |
| `mingw-release` / `mingw-debug` | Windows | MinGW with full static linking |
| `clang-release` / `clang-debug` | Linux / macOS | Clang 20 + GFortran |

#### Linux / macOS

```bash
cmake --preset gcc-release
cmake --build --preset gcc-release -j$(nproc)
```

#### Windows (MSYS2 / MinGW64)

Ensure `g++` and `gfortran` are in your PATH (e.g., via [MSYS2](https://www.msys2.org/) with the `mingw-w64-x86_64-gcc` and `mingw-w64-x86_64-gcc-fortran` packages).

```bash
cmake --preset mingw-release
cmake --build --preset mingw-release -j%NUMBER_OF_PROCESSORS%
```

The resulting `build/Marlim3.exe` is fully statically linked and does not require external DLLs.

The compiled executable is placed at `build/Marlim3`.

To use it with the Python package, copy it to the `marlim3/` directory:

```bash
cp build/Marlim3 marlim3/
```

Then activate the package locally (skipping recompilation):

```bash
MARLIM3_SKIP_BUILD=1 uv sync --locked
```

### Run tests

```bash
uv run pytest tests/ -v
```

### Run the GUI

```bash
uv sync --group gui
uv run streamlit run gui/app.py
```

The GUI auto-detects the executable from `build/Marlim3`.

## Note

Several resources and portions of the source code are currently written in Portuguese. We plan to gradually translate all content into English.