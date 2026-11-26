"""RsPbTable: the in-situ Rs is the table interpolated between the grid nodes that bracket each cell."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
DEMO = ROOT / "demos" / "simplifiedProduction"
MODEL = "simplifiedProduction.mr3"
REL_TOL = 1e-9
pytestmark = pytest.mark.simulacao

# (numPoints, (pmin, pmax), (tmin, tmax)); the demo's profile spans about 81-90 C and 11-77 kgf/cm2.
INTERIOR = (7, (10.0, 130.0), (50.0, 110.0))        # 10 C spacing; the bottom cell sits on the 90 C node
FIRST_INTERVAL = (5, (10.0, 130.0), (80.0, 200.0))  # all cells in the first temperature interval
LAST_INTERVAL = (5, (10.0, 80.0), (0.0, 95.0))      # all cells in the last temperature interval


@pytest.fixture(scope="module")
def engine():
    configured = os.environ.get("MARLIM3_TEST_EXECUTABLE")
    suffix = ".exe" if sys.platform == "win32" else ""
    candidates = [Path(configured)] if configured else [
        ROOT / "build" / f"Marlim3{suffix}",
        ROOT / "tools" / "marlim3" / f"Marlim3{suffix}",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    pytest.skip("Compiled engine not found; set MARLIM3_TEST_EXECUTABLE")


def run_with_table(engine, directory, points, pressures, temperatures):
    """Run the demo with the RS table; return the table and (p, T, Rs) of the cells inside it."""
    directory.mkdir(parents=True)
    model = json.loads((DEMO / MODEL).read_text(encoding="utf-8"))
    # pressureTable makes the reader parse compTable
    model["initialConfig"].update({"RsPbTable": True, "pressureTable": True})
    model["compTable"] = {
        "active": True, "numPoints": points,
        "minPressure": pressures[0], "maxPressure": pressures[1],
        "minTemperature": temperatures[0], "maxTemperature": temperatures[1],
    }
    model["productionProfile"]["SR"] = True
    (directory / "model.mr3").write_text(json.dumps(model), encoding="utf-8")
    process = subprocess.run(
        [str(engine), "--input", "model.mr3", "--dir", str(directory)], cwd=directory,
        capture_output=True, text=True, errors="replace", timeout=300,
    )
    assert process.returncode == 0, process.stdout[-3000:] + process.stderr[-3000:]
    table = [
        [float(value) for value in line.split(";") if value.strip()]
        for line in (directory / "perfilRSLivia").read_text().splitlines() if line.strip()
    ]
    lines = [
        line for line in (directory / "PERFISP-0.dat").read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    head = [name.strip() for name in lines[0].split(";")]
    columns = [
        head.index(next(name for name in head if name.startswith(prefix)))
        for prefix in ("Pressure", "Temperature", "Solubility ratio")
    ]
    cells = [tuple(float(line.split(";")[k]) for k in columns) for line in lines[1:]]
    n = len(table) - 1
    inside = [
        (p, t, rs) for p, t, rs in cells
        if table[1][0] < p < table[n][0] and table[0][1] < t < table[0][n]
    ]
    return table, inside


def bracket(axis, value):
    """Index k with axis[k] <= value <= axis[k + 1]."""
    return max(k for k in range(len(axis) - 1) if axis[k] <= value)


def interpolate(table, pressure, temperature):
    """Bilinear interpolation between the bracketing rows and columns."""
    n = len(table) - 1
    temperatures = [table[0][j] for j in range(1, n + 1)]
    pressures = [table[i][0] for i in range(1, n + 1)]
    j = bracket(temperatures, temperature)
    i = bracket(pressures, pressure)
    rt = (temperature - temperatures[j]) / (temperatures[j + 1] - temperatures[j])
    rp = (pressure - pressures[i]) / (pressures[i + 1] - pressures[i])
    low = (1 - rp) * table[i + 1][j + 1] + rp * table[i + 2][j + 1]
    high = (1 - rp) * table[i + 1][j + 2] + rp * table[i + 2][j + 2]
    return (1 - rt) * low + rt * high


@pytest.mark.parametrize("grid", [INTERIOR, FIRST_INTERVAL, LAST_INTERVAL],
                         ids=["interior", "first-interval", "last-interval"])
def test_in_situ_rs_is_the_table_between_the_bracketing_columns(engine, tmp_path, grid):
    table, cells = run_with_table(engine, tmp_path / "run", *grid)
    assert len(cells) >= 10
    for p, t, rs in cells:
        assert rs == pytest.approx(interpolate(table, p, t), rel=REL_TOL), (p, t)


def test_in_situ_rs_on_a_grid_temperature_is_that_column(engine, tmp_path):
    table, cells = run_with_table(engine, tmp_path / "run", *INTERIOR)
    n = len(table) - 1
    temperatures = [table[0][j] for j in range(1, n + 1)]
    pressures = [table[i][0] for i in range(1, n + 1)]
    on_grid = [(p, t, rs) for p, t, rs in cells if t in temperatures]
    assert on_grid
    for p, t, rs in on_grid:
        j = temperatures.index(t) + 1
        i = bracket(pressures, p)
        rp = (p - pressures[i]) / (pressures[i + 1] - pressures[i])
        assert rs == pytest.approx((1 - rp) * table[i + 1][j] + rp * table[i + 2][j], rel=REL_TOL), (p, t)
