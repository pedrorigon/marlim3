"""Desktop launcher for the Streamlit GUI.

The launcher starts the bundled Streamlit app on localhost and renders it in a
Qt WebEngine window. Qt and its Chromium runtime are included in desktop builds.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path


APP_NAME = "Marlim3"
APP_ID = "br.com.petrobras.marlim3"
DESKTOP_FILE_NAME = APP_ID
DEFAULT_HOST = "127.0.0.1"
MAX_PORT_ATTEMPTS = 10


def get_bundle_root() -> Path:
    """Return the root that contains gui/, demos/ and marlim3/."""
    env_root = os.environ.get("MARLIM3_DESKTOP_ROOT")
    if env_root:
        return Path(env_root).resolve()
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS).resolve()
    return Path(__file__).resolve().parents[1]


def get_gui_app_path(bundle_root: Path | None = None) -> Path:
    root = bundle_root or get_bundle_root()
    return root / "gui" / "app.py"


def executable_name() -> str:
    return "Marlim3.exe" if sys.platform.startswith("win") else "Marlim3"


def get_engine_path(bundle_root: Path | None = None) -> Path:
    root = bundle_root or get_bundle_root()
    return root / "marlim3" / executable_name()


def get_app_icon_path(bundle_root: Path | None = None) -> Path:
    root = bundle_root or get_bundle_root()
    return root / "assets" / "icons" / "app-icon.png"


def get_logo_path(bundle_root: Path | None = None) -> Path:
    root = bundle_root or get_bundle_root()
    return root / "assets" / "branding" / "logo.svg"


def configure_platform_app_id() -> None:
    """Set the native application identifier used for taskbar grouping."""
    if sys.platform.startswith("win"):
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_ID)


def find_free_port(host: str = DEFAULT_HOST) -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((host, 0))
        return int(sock.getsockname()[1])


def streamlit_flag_options(port: int, host: str = DEFAULT_HOST) -> dict[str, object]:
    return {
        "server.address": host,
        "server.port": port,
        "server.headless": True,
        "server.fileWatcherType": "none",
        "browser.gatherUsageStats": False,
        "global.developmentMode": False,
    }


def configure_environment(bundle_root: Path, engine_path: Path) -> None:
    os.environ["MARLIM3_DESKTOP_ROOT"] = str(bundle_root)
    os.environ.setdefault("MARLIM3_EXEC", str(engine_path))
    os.environ.setdefault("MARLIM3_SKIP_EXECUTABLE_RESOLUTION", "1")
    os.environ.setdefault("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")


def wait_for_server(
    url: str,
    process: subprocess.Popen | None = None,
    timeout: float = 30.0,
) -> None:
    health_url = f"{url}/_stcore/health"
    deadline = time.monotonic() + timeout
    last_error: Exception | None = None
    while time.monotonic() < deadline:
        if process is not None and process.poll() is not None:
            raise RuntimeError(
                f"Streamlit exited before becoming ready (exit code {process.returncode})"
            )
        try:
            with urllib.request.urlopen(health_url, timeout=1.0) as response:
                if response.status == 200 and response.read().decode().strip() == "ok":
                    return
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
        time.sleep(0.25)
    raise RuntimeError(f"Streamlit did not become ready at {health_url}") from last_error


def run_streamlit_server(app_path: Path, port: int, host: str = DEFAULT_HOST) -> None:
    from streamlit.web import bootstrap

    options = streamlit_flag_options(port=port, host=host)
    bootstrap.load_config_options(options)
    bootstrap.run(
        str(app_path),
        False,
        [],
        options,
    )


def open_desktop_window(url: str, close_after_ms: int | None = None) -> None:
    from PySide6.QtCore import QTimer, QUrl
    from PySide6.QtGui import QIcon
    from PySide6.QtWebEngineWidgets import QWebEngineView
    from PySide6.QtWidgets import QApplication, QMainWindow

    configure_platform_app_id()

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationDisplayName(APP_NAME)
    app.setApplicationName(APP_NAME)
    app.setDesktopFileName(DESKTOP_FILE_NAME)
    app.setQuitOnLastWindowClosed(True)

    icon_path = get_app_icon_path()
    app_icon = QIcon(str(icon_path)) if icon_path.exists() else QIcon()
    if app_icon.isNull():
        raise RuntimeError(f"Invalid desktop icon: {icon_path}")
    app.setWindowIcon(app_icon)

    window = QMainWindow()
    window.setWindowTitle(APP_NAME)
    window.setMinimumSize(960, 640)
    window.setWindowIcon(app_icon)

    browser = QWebEngineView(window)
    browser.setUrl(QUrl(url))
    window.setCentralWidget(browser)
    window.showMaximized()

    if close_after_ms is not None:
        QTimer.singleShot(close_after_ms, window.close)

    app.exec()
    browser.close()
    browser.deleteLater()


def build_runtime_config(port: int | None = None, host: str = DEFAULT_HOST) -> dict[str, str | int]:
    bundle_root = get_bundle_root()
    chosen_port = port if port is not None else 0
    return {
        "bundle_root": str(bundle_root),
        "gui_app": str(get_gui_app_path(bundle_root)),
        "engine": str(get_engine_path(bundle_root)),
        "app_icon": str(get_app_icon_path(bundle_root)),
        "logo": str(get_logo_path(bundle_root)),
        "host": host,
        "port": chosen_port,
        "url": f"http://{host}:{chosen_port}",
    }


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch the Marlim3 desktop GUI.")
    parser.add_argument("--host", default=DEFAULT_HOST)
    parser.add_argument("--port", type=int, default=None)
    parser.add_argument("--print-config", action="store_true", help="Print resolved runtime paths and exit.")
    parser.add_argument("--serve", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--parent-pid", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--window-smoke-test", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def server_command(host: str, port: int) -> list[str]:
    args = [
        "--serve",
        "--host",
        host,
        "--port",
        str(port),
        "--parent-pid",
        str(os.getpid()),
    ]
    if getattr(sys, "frozen", False):
        return [sys.executable, *args]
    return [sys.executable, "-m", "gui.desktop", *args]


def stop_server(process: subprocess.Popen) -> None:
    if sys.platform.startswith("win"):
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True,
            check=False,
        )
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=5)
        return

    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


def monitor_parent(parent_pid: int, interval: float = 0.5) -> None:
    while True:
        try:
            os.kill(parent_pid, 0)
        except (OSError, ProcessLookupError):
            os._exit(0)
        time.sleep(interval)


def start_server(
    host: str,
    requested_port: int | None = None,
) -> tuple[subprocess.Popen, str]:
    attempts = 1 if requested_port is not None else MAX_PORT_ATTEMPTS
    last_error: Exception | None = None

    for _ in range(attempts):
        port = requested_port if requested_port is not None else find_free_port(host)
        url = f"http://{host}:{port}"
        popen_options = (
            {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
            if sys.platform.startswith("win")
            else {"start_new_session": True}
        )
        process = subprocess.Popen(server_command(host, port), **popen_options)
        try:
            wait_for_server(url, process=process)
            return process, url
        except RuntimeError as exc:
            last_error = exc
            stop_server(process)

    raise RuntimeError(
        f"Could not start Marlim3 on a free local port after {attempts} attempt(s)"
    ) from last_error


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    config = build_runtime_config(port=args.port, host=args.host)

    if args.print_config:
        print(json.dumps(config, indent=2, sort_keys=True))
        return 0

    bundle_root = Path(str(config["bundle_root"]))
    app_path = Path(str(config["gui_app"]))
    engine_path = Path(str(config["engine"]))
    required_files = {
        "Streamlit app": app_path,
        "Marlim3 engine": engine_path,
        "desktop icon": Path(str(config["app_icon"])),
        "Marlim3 logo": Path(str(config["logo"])),
    }
    for description, path in required_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"{description} not found: {path}")

    configure_environment(bundle_root, engine_path)

    if args.serve:
        if args.port is None:
            raise ValueError("--serve requires --port")
        if args.parent_pid is not None:
            threading.Thread(
                target=monitor_parent,
                args=(args.parent_pid,),
                daemon=True,
            ).start()
        run_streamlit_server(app_path, args.port, args.host)
        return 0

    process, url = start_server(args.host, requested_port=args.port)
    try:
        close_after_ms = 3000 if args.window_smoke_test else None
        open_desktop_window(url, close_after_ms=close_after_ms)
    finally:
        stop_server(process)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
