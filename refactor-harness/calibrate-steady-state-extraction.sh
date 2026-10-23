#!/usr/bin/env bash
# Check every extracted block and prove that the checker rejects numerical and binding changes.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$script_dir/steady-state-extraction-layout.py" --calibrate
