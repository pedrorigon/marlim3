#!/usr/bin/env python3
"""Build one standalone Marlim3 desktop artifact for the current platform."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIST_DIR = ROOT / "dist"
WORK_DIR = ROOT / "build" / "desktop"
SPEC_FILE = ROOT / "gui" / "packaging" / "marlim3-desktop.spec"


def run(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
    cwd: Path = ROOT,
) -> None:
    print(f"+ {' '.join(command)}", flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def normalized_machine() -> str:
    machine = platform.machine().lower()
    return {
        "amd64": "x86_64",
        "x64": "x86_64",
        "aarch64": "arm64",
    }.get(machine, machine)


def platform_config() -> tuple[str, str, str]:
    system = platform.system()
    machine = normalized_machine()
    supported = {
        ("Linux", "x86_64"): (
            "gcc-release",
            "Marlim3",
            "Marlim3-desktop-linux-x64",
        ),
        ("Windows", "x86_64"): (
            "mingw-release",
            "Marlim3.exe",
            "Marlim3-desktop-windows-x64.exe",
        ),
        ("Darwin", "arm64"): (
            "gcc-release",
            "Marlim3",
            "Marlim3-desktop-macos-arm64.dmg",
        ),
    }
    try:
        return supported[(system, machine)]
    except KeyError as exc:
        raise RuntimeError(
            f"Unsupported desktop build platform: {system} {machine}. "
            "Supported targets are Linux x64, Windows x64, and macOS ARM64."
        ) from exc


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool not found in PATH: {name}")


def build_engine(preset: str, engine_name: str) -> None:
    require_tool("cmake")
    jobs = str(os.cpu_count() or 1)
    run(["cmake", "--preset", preset])
    run(["cmake", "--build", "--preset", preset, f"-j{jobs}"])

    engine = ROOT / "marlim3" / engine_name
    if not engine.is_file():
        raise FileNotFoundError(f"Native engine was not generated: {engine}")


def sync_dependencies() -> None:
    require_tool("uv")
    env = os.environ.copy()
    env["MARLIM3_SKIP_BUILD"] = "1"
    env["MARLIM3_SKIP_EXECUTABLE_RESOLUTION"] = "1"
    run(["uv", "sync", "--locked", "--group", "gui", "--group", "desktop"], env=env)


def build_pyinstaller_bundle() -> None:
    require_tool("uv")
    shutil.rmtree(DIST_DIR, ignore_errors=True)
    DIST_DIR.mkdir(parents=True)
    run(
        [
            "uv",
            "run",
            "--group",
            "gui",
            "--group",
            "desktop",
            "pyinstaller",
            str(SPEC_FILE),
            "--noconfirm",
            "--clean",
        ]
    )


def validation_environment() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    env.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")
    return env


def validate_executable(executable: Path, smoke_test: bool) -> None:
    run([str(executable), "--print-config"])
    if smoke_test:
        run([str(executable), "--window-smoke-test"], env=validation_environment())


def package_linux(smoke_test: bool) -> Path:
    raw_executable = DIST_DIR / "Marlim3-Desktop"
    if not raw_executable.is_file():
        raise FileNotFoundError(f"PyInstaller output not found: {raw_executable}")
    validate_executable(raw_executable, smoke_test)

    artifact = DIST_DIR / "Marlim3-desktop-linux-x64"
    raw_executable.replace(artifact)
    return artifact


def package_windows(smoke_test: bool) -> Path:
    raw_executable = DIST_DIR / "Marlim3-Desktop.exe"
    if not raw_executable.is_file():
        raise FileNotFoundError(f"PyInstaller output not found: {raw_executable}")
    validate_executable(raw_executable, smoke_test)
    artifact = DIST_DIR / "Marlim3-desktop-windows-x64.exe"
    raw_executable.replace(artifact)
    return artifact


def package_macos(smoke_test: bool) -> Path:
    app_bundle = DIST_DIR / "Marlim3.app"
    executable = app_bundle / "Contents" / "MacOS" / "Marlim3"
    if not executable.is_file():
        raise FileNotFoundError(f"PyInstaller app bundle not found: {app_bundle}")
    validate_executable(executable, smoke_test)

    require_tool("hdiutil")
    dmg_root = WORK_DIR / "dmg-root"
    shutil.rmtree(dmg_root, ignore_errors=True)
    dmg_root.mkdir(parents=True)
    shutil.copytree(app_bundle, dmg_root / "Marlim3.app")
    (dmg_root / "Applications").symlink_to("/Applications")

    artifact = DIST_DIR / "Marlim3-desktop-macos-arm64.dmg"
    run(
        [
            "hdiutil",
            "create",
            "-volname",
            "Marlim3",
            "-srcfolder",
            str(dmg_root),
            "-ov",
            "-format",
            "UDZO",
            str(artifact),
        ]
    )
    shutil.rmtree(app_bundle)
    return artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-engine-build",
        action="store_true",
        help="Use the existing native engine in marlim3/.",
    )
    parser.add_argument(
        "--skip-dependency-sync",
        action="store_true",
        help="Do not run uv sync before packaging.",
    )
    parser.add_argument(
        "--skip-smoke-test",
        action="store_true",
        help="Skip the offscreen desktop window smoke test.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    preset, engine_name, expected_artifact = platform_config()
    print(
        f"Building Marlim3 desktop for {platform.system()} {normalized_machine()}",
        flush=True,
    )

    if not args.skip_engine_build:
        build_engine(preset, engine_name)
    elif not (ROOT / "marlim3" / engine_name).is_file():
        raise FileNotFoundError(
            f"--skip-engine-build requires marlim3/{engine_name} to exist"
        )

    if not args.skip_dependency_sync:
        sync_dependencies()
    build_pyinstaller_bundle()

    smoke_test = not args.skip_smoke_test
    if platform.system() == "Linux":
        artifact = package_linux(smoke_test)
    elif platform.system() == "Windows":
        artifact = package_windows(smoke_test)
    else:
        artifact = package_macos(smoke_test)

    if artifact.name != expected_artifact or not artifact.is_file():
        raise RuntimeError(f"Expected desktop artifact was not generated: {artifact}")
    print(f"Desktop artifact: {artifact}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
