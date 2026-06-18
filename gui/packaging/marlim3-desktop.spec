# -*- mode: python ; coding: utf-8 -*-

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules, copy_metadata


root = Path(SPECPATH).resolve().parents[1]

icons_dir = root / "assets" / "icons"
branding_dir = root / "assets" / "branding"
windows_icon = icons_dir / "app-icon.ico"
macos_icon = icons_dir / "app-icon.icns"
runtime_icon = icons_dir / "app-icon.png"
logo = branding_dir / "logo.svg"

version_namespace = {}
exec((root / "_version.py").read_text(encoding="utf-8"), version_namespace)
app_version = version_namespace["__version__"]

required_assets = (windows_icon, macos_icon, runtime_icon, logo)
missing_assets = [str(path) for path in required_assets if not path.is_file()]
if missing_assets:
    raise FileNotFoundError(f"Missing desktop assets: {', '.join(missing_assets)}")

datas = []
binaries = []
hiddenimports = []

streamlit_datas, streamlit_binaries, streamlit_hiddenimports = collect_all("streamlit")
datas += streamlit_datas
binaries += streamlit_binaries
hiddenimports += streamlit_hiddenimports

hiddenimports += [
    "jsonschema",
    "lxml",
    "matplotlib",
    "nbformat",
    "networkx",
    "numpy",
    "openpyxl",
    "pandas",
    "plotly",
    "plotly.graph_objects",
    "PySide6.QtCore",
    "PySide6.QtGui",
    "PySide6.QtNetwork",
    "PySide6.QtSvg",
    "PySide6.QtWebChannel",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWidgets",
    "requests",
    "seaborn",
    "sklearn",
    "xmltodict",
]

for package in ("streamlit", "PySide6", "plotly", "pandas", "jsonschema", "marlim3"):
    try:
        datas += copy_metadata(package)
    except Exception:
        pass

marlim_datas = collect_data_files("marlim3")
datas += marlim_datas
hiddenimports += collect_submodules("marlim3")

datas += [
    (str(root / "gui" / "app.py"), "gui"),
    (str(root / "demos"), "demos"),
    (str(runtime_icon), "assets/icons"),
    (str(logo), "assets/branding"),
]

engine_name = "Marlim3.exe" if sys.platform.startswith("win") else "Marlim3"
engine_path = root / "marlim3" / engine_name
if not engine_path.is_file():
    raise FileNotFoundError(f"Missing native Marlim3 engine: {engine_path}")
binaries.append((str(engine_path), "marlim3"))

a = Analysis(
    [str(root / "gui" / "desktop.py")],
    pathex=[str(root)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

if sys.platform == "darwin":
    exe = EXE(
        pyz,
        a.scripts,
        [],
        exclude_binaries=True,
        name="Marlim3",
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,
        console=False,
        disable_windowed_traceback=False,
        argv_emulation=False,
        target_arch=None,
        icon=str(macos_icon),
        codesign_identity=None,
        entitlements_file=None,
    )
    coll = COLLECT(
        exe,
        a.binaries,
        a.datas,
        strip=False,
        upx=False,
        upx_exclude=[],
        name="Marlim3",
    )
    app = BUNDLE(
        coll,
        name="Marlim3.app",
        icon=str(macos_icon),
        bundle_identifier="br.com.petrobras.marlim3",
        info_plist={
            "CFBundleName": "Marlim3",
            "CFBundleDisplayName": "Marlim3",
            "CFBundleShortVersionString": app_version,
            "CFBundleVersion": app_version,
            "NSHighResolutionCapable": True,
        },
    )
else:
    exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.datas,
        [],
        name="Marlim3-Desktop",
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,
        upx_exclude=[],
        runtime_tmpdir=None,
        console=False,
        disable_windowed_traceback=False,
        argv_emulation=False,
        target_arch=None,
        icon=str(windows_icon) if sys.platform.startswith("win") else None,
        codesign_identity=None,
        entitlements_file=None,
    )
