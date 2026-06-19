# -*- mode: python ; coding: utf-8 -*-

import sys
from pathlib import Path

from PyInstaller.depend.bindepend import get_imports
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
    "networkx",
    "numpy",
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


def collect_linux_runtime_libraries():
    import PySide6

    pyside_dir = Path(PySide6.__file__).resolve().parent
    qt_dir = pyside_dir / "Qt"
    runtime_roots = (
        qt_dir / "lib" / "libQt6WebEngineCore.so.6",
        qt_dir / "plugins" / "platforms" / "libqxcb.so",
        qt_dir / "plugins" / "platforms" / "libqwayland.so",
        qt_dir / "plugins" / "platforms" / "libqoffscreen.so",
    )
    glibc_libraries = {
        "ld-linux-x86-64.so.2",
        "libanl.so.1",
        "libc.so.6",
        "libdl.so.2",
        "libm.so.6",
        "libmvec.so.1",
        "libnsl.so.1",
        "libpthread.so.0",
        "libresolv.so.2",
        "librt.so.1",
        "libutil.so.1",
    }
    pending = [path for path in runtime_roots if path.is_file()]
    processed = set()
    portable_libraries = {}
    missing_libraries = set()

    while pending:
        binary = pending.pop()
        resolved_binary = binary.resolve()
        if resolved_binary in processed:
            continue
        processed.add(resolved_binary)

        for library_name, library_path in get_imports(str(resolved_binary)):
            if library_name in glibc_libraries:
                continue
            if not library_path:
                missing_libraries.add(library_name)
                continue
            soname_path = Path(library_path)
            resolved_library = soname_path.resolve()
            if pyside_dir in resolved_library.parents:
                continue
            portable_libraries[library_name] = soname_path
            pending.append(resolved_library)

    if missing_libraries:
        names = ", ".join(sorted(missing_libraries))
        raise RuntimeError(f"Missing Linux desktop runtime libraries: {names}")

    return [
        (str(library_path), ".")
        for _, library_path in sorted(portable_libraries.items())
    ]


if sys.platform.startswith("linux"):
    binaries += collect_linux_runtime_libraries()


a = Analysis(
    [str(root / "gui" / "desktop.py")],
    pathex=[str(root)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[str(root / "gui" / "packaging" / "hooks")],
    runtime_hooks=[],
    excludes=[
        "IPython",
        "_tkinter",
        "jedi",
        "nbformat",
        "openpyxl",
        "pytest",
        "scipy",
        "sklearn",
        "streamlit.hello",
        "streamlit.testing",
        "tkinter",
    ],
    noarchive=False,
    optimize=1,
    hooksconfig={
        "matplotlib": {
            "backends": ["Agg"],
        },
    },
)

if sys.platform.startswith("linux"):
    def required_linux_qt_data(entry):
        destination = entry[0].replace("\\", "/")
        name = Path(destination).name

        if name == "qtwebengine_devtools_resources.pak":
            return False
        if "/qtwebengine_locales/" in destination:
            return name in {"en-US.pak", "pt-BR.pak"}
        if name.startswith("qtwebengine_") and name.endswith(".qm"):
            return False
        return True

    a.datas = type(a.datas)(filter(required_linux_qt_data, a.datas))

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
