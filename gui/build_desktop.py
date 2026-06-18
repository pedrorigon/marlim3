#!/usr/bin/env python3
"""Build one standalone Marlim3 desktop artifact for the current platform."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIST_DIR = ROOT / "dist"
WORK_DIR = ROOT / "build" / "desktop"
SPEC_FILE = ROOT / "gui" / "packaging" / "marlim3-desktop.spec"
LINUX_PACKAGING_DIR = ROOT / "gui" / "packaging" / "linux"
ICON_FILE = ROOT / "assets" / "icons" / "app-icon.png"
APPIMAGE_TOOL_URL = (
    "https://github.com/AppImage/AppImageKit/releases/download/continuous/"
    "appimagetool-x86_64.AppImage"
)


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
            "Marlim3-desktop-linux-x64.AppImage",
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

    app_dir = WORK_DIR / "AppDir"
    shutil.rmtree(app_dir, ignore_errors=True)
    (app_dir / "usr" / "bin").mkdir(parents=True)
    icon_dir = app_dir / "usr" / "share" / "icons" / "hicolor" / "1024x1024" / "apps"
    icon_dir.mkdir(parents=True)
    applications_dir = app_dir / "usr" / "share" / "applications"
    applications_dir.mkdir(parents=True)
    metainfo_dir = app_dir / "usr" / "share" / "metainfo"
    metainfo_dir.mkdir(parents=True)
    shutil.copy2(raw_executable, app_dir / "usr" / "bin" / "Marlim3")
    shutil.copy2(ICON_FILE, app_dir / "marlim3.png")
    shutil.copy2(ICON_FILE, app_dir / ".DirIcon")
    shutil.copy2(ICON_FILE, icon_dir / "marlim3.png")
    app_run = app_dir / "AppRun"
    shutil.copy2(LINUX_PACKAGING_DIR / "AppRun", app_run)
    app_run.chmod(0o755)
    desktop_name = "br.com.petrobras.marlim3.desktop"
    desktop_file = LINUX_PACKAGING_DIR / desktop_name
    shutil.copy2(desktop_file, app_dir / desktop_name)
    shutil.copy2(desktop_file, applications_dir / desktop_name)
    shutil.copy2(
        LINUX_PACKAGING_DIR / "br.com.petrobras.marlim3.appdata.xml",
        metainfo_dir / "br.com.petrobras.marlim3.appdata.xml",
    )

    tool = WORK_DIR / "appimagetool-x86_64.AppImage"
    tool.parent.mkdir(parents=True, exist_ok=True)
    if not tool.is_file():
        print(f"Downloading {APPIMAGE_TOOL_URL}", flush=True)
        urllib.request.urlretrieve(APPIMAGE_TOOL_URL, tool)
    tool.chmod(0o755)

    artifact = DIST_DIR / "Marlim3-desktop-linux-x64.AppImage"
    env = os.environ.copy()
    env.update({"APPIMAGE_EXTRACT_AND_RUN": "1", "ARCH": "x86_64"})
    run([str(tool), str(app_dir), str(artifact)], env=env)
    raw_executable.unlink()

    validation_env = os.environ.copy()
    validation_env["APPIMAGE_EXTRACT_AND_RUN"] = "1"
    run([str(artifact), "--print-config"], env=validation_env)
    if smoke_test:
        env = validation_environment()
        env["APPIMAGE_EXTRACT_AND_RUN"] = "1"
        run([str(artifact), "--window-smoke-test"], env=env)
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
