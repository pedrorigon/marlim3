import os
import sys
from pathlib import Path
from unittest.mock import Mock

os.environ.setdefault("MARLIM3_SKIP_EXECUTABLE_RESOLUTION", "1")

from gui import desktop


def test_desktop_paths_resolve_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("MARLIM3_DESKTOP_ROOT", str(tmp_path))

    assert desktop.get_bundle_root() == tmp_path.resolve()
    assert desktop.get_gui_app_path() == tmp_path / "gui" / "app.py"
    assert desktop.get_engine_path().parent == tmp_path / "marlim3"
    assert desktop.get_app_icon_path() == tmp_path / "assets" / "icons" / "app-icon.png"
    assert desktop.get_logo_path() == tmp_path / "assets" / "branding" / "logo.svg"


def test_desktop_executable_name_matches_platform():
    expected = "Marlim3.exe" if sys.platform.startswith("win") else "Marlim3"

    assert desktop.executable_name() == expected
    assert desktop.DESKTOP_FILE_NAME == desktop.APP_ID


def test_configure_environment_sets_desktop_defaults(monkeypatch, tmp_path):
    monkeypatch.delenv("MARLIM3_EXEC", raising=False)
    monkeypatch.delenv("MARLIM3_SKIP_EXECUTABLE_RESOLUTION", raising=False)

    engine = tmp_path / "marlim3" / desktop.executable_name()
    desktop.configure_environment(tmp_path, engine)

    assert os.environ["MARLIM3_DESKTOP_ROOT"] == str(tmp_path)
    assert os.environ["MARLIM3_EXEC"] == str(engine)
    assert os.environ["MARLIM3_SKIP_EXECUTABLE_RESOLUTION"] == "1"


def test_streamlit_flags_are_headless_and_localhost():
    flags = desktop.streamlit_flag_options(port=8765)

    assert flags["server.address"] == "127.0.0.1"
    assert flags["server.port"] == 8765
    assert flags["server.headless"] is True
    assert flags["server.fileWatcherType"] == "none"


def test_build_runtime_config_uses_expected_url(monkeypatch, tmp_path):
    monkeypatch.setenv("MARLIM3_DESKTOP_ROOT", str(tmp_path))

    config = desktop.build_runtime_config(port=9012)

    assert config["url"] == "http://127.0.0.1:9012"
    assert Path(config["gui_app"]) == tmp_path / "gui" / "app.py"
    assert Path(config["app_icon"]) == tmp_path / "assets" / "icons" / "app-icon.png"
    assert Path(config["logo"]) == tmp_path / "assets" / "branding" / "logo.svg"


def test_build_runtime_config_uses_automatic_port_by_default(monkeypatch, tmp_path):
    monkeypatch.setenv("MARLIM3_DESKTOP_ROOT", str(tmp_path))

    config = desktop.build_runtime_config()

    assert config["port"] == 0
    assert config["url"] == "http://127.0.0.1:0"


def test_wait_for_server_checks_streamlit_health_endpoint(monkeypatch):
    response = Mock()
    response.status = 200
    response.read.return_value = b"ok"
    response.__enter__ = Mock(return_value=response)
    response.__exit__ = Mock(return_value=False)
    urlopen = Mock(return_value=response)
    monkeypatch.setattr(desktop.urllib.request, "urlopen", urlopen)

    desktop.wait_for_server("http://127.0.0.1:9012", timeout=0.1)

    urlopen.assert_called_once_with(
        "http://127.0.0.1:9012/_stcore/health",
        timeout=1.0,
    )


def test_server_command_relaunches_desktop_module(monkeypatch):
    monkeypatch.setattr(desktop.sys, "frozen", False, raising=False)
    monkeypatch.setattr(desktop.os, "getpid", lambda: 4567)

    command = desktop.server_command("127.0.0.1", 9012)

    assert command == [
        sys.executable,
        "-m",
        "gui.desktop",
        "--serve",
        "--host",
        "127.0.0.1",
        "--port",
        "9012",
        "--parent-pid",
        "4567",
    ]


def test_run_streamlit_server_loads_port_options(monkeypatch, tmp_path):
    from streamlit.web import bootstrap

    load_config_options = Mock()
    run = Mock()
    monkeypatch.setattr(bootstrap, "load_config_options", load_config_options)
    monkeypatch.setattr(bootstrap, "run", run)
    app_path = tmp_path / "app.py"

    desktop.run_streamlit_server(app_path, 9012)

    options = desktop.streamlit_flag_options(9012)
    load_config_options.assert_called_once_with(options)
    run.assert_called_once_with(str(app_path), False, [], options)


def test_main_stops_server_when_window_closes(monkeypatch, tmp_path):
    app_path = tmp_path / "gui" / "app.py"
    engine_path = tmp_path / "marlim3" / desktop.executable_name()
    icon_path = tmp_path / "assets" / "icons" / "app-icon.png"
    logo_path = tmp_path / "assets" / "branding" / "logo.svg"
    app_path.parent.mkdir(parents=True)
    engine_path.parent.mkdir(parents=True)
    icon_path.parent.mkdir(parents=True)
    logo_path.parent.mkdir(parents=True)
    app_path.touch()
    engine_path.touch()
    icon_path.touch()
    logo_path.touch()

    process = Mock()
    monkeypatch.setattr(
        desktop,
        "build_runtime_config",
        lambda **kwargs: {
            "bundle_root": str(tmp_path),
            "gui_app": str(app_path),
            "engine": str(engine_path),
            "app_icon": str(icon_path),
            "logo": str(logo_path),
            "host": "127.0.0.1",
            "port": 0,
            "url": "http://127.0.0.1:0",
        },
    )
    monkeypatch.setattr(
        desktop,
        "start_server",
        Mock(return_value=(process, "http://127.0.0.1:9012")),
    )
    open_window = Mock()
    stop_server = Mock()
    monkeypatch.setattr(desktop, "open_desktop_window", open_window)
    monkeypatch.setattr(desktop, "stop_server", stop_server)

    assert desktop.main([]) == 0
    open_window.assert_called_once_with(
        "http://127.0.0.1:9012",
        close_after_ms=None,
    )
    stop_server.assert_called_once_with(process)


def test_main_rejects_missing_desktop_assets(monkeypatch, tmp_path):
    app_path = tmp_path / "gui" / "app.py"
    engine_path = tmp_path / "marlim3" / desktop.executable_name()
    app_path.parent.mkdir(parents=True)
    engine_path.parent.mkdir(parents=True)
    app_path.touch()
    engine_path.touch()

    start_server = Mock()
    monkeypatch.setattr(
        desktop,
        "build_runtime_config",
        lambda **kwargs: {
            "bundle_root": str(tmp_path),
            "gui_app": str(app_path),
            "engine": str(engine_path),
            "app_icon": str(tmp_path / "assets" / "icons" / "app-icon.png"),
            "logo": str(tmp_path / "assets" / "branding" / "logo.svg"),
            "host": "127.0.0.1",
            "port": 0,
            "url": "http://127.0.0.1:0",
        },
    )
    monkeypatch.setattr(desktop, "start_server", start_server)

    try:
        desktop.main([])
    except FileNotFoundError as exc:
        assert "desktop icon not found" in str(exc)
    else:
        raise AssertionError("Expected missing desktop asset failure")

    start_server.assert_not_called()


def test_main_stops_server_when_window_fails(monkeypatch, tmp_path):
    app_path = tmp_path / "gui" / "app.py"
    engine_path = tmp_path / "marlim3" / desktop.executable_name()
    icon_path = tmp_path / "assets" / "icons" / "app-icon.png"
    logo_path = tmp_path / "assets" / "branding" / "logo.svg"
    app_path.parent.mkdir(parents=True)
    engine_path.parent.mkdir(parents=True)
    icon_path.parent.mkdir(parents=True)
    logo_path.parent.mkdir(parents=True)
    app_path.touch()
    engine_path.touch()
    icon_path.touch()
    logo_path.touch()

    process = Mock()
    monkeypatch.setattr(
        desktop,
        "build_runtime_config",
        lambda **kwargs: {
            "bundle_root": str(tmp_path),
            "gui_app": str(app_path),
            "engine": str(engine_path),
            "app_icon": str(icon_path),
            "logo": str(logo_path),
            "host": "127.0.0.1",
            "port": 0,
            "url": "http://127.0.0.1:0",
        },
    )
    monkeypatch.setattr(
        desktop,
        "start_server",
        Mock(return_value=(process, "http://127.0.0.1:9012")),
    )
    monkeypatch.setattr(
        desktop,
        "open_desktop_window",
        Mock(side_effect=RuntimeError("window failed")),
    )
    stop_server = Mock()
    monkeypatch.setattr(desktop, "stop_server", stop_server)

    try:
        desktop.main([])
    except RuntimeError as exc:
        assert str(exc) == "window failed"
    else:
        raise AssertionError("Expected window failure")

    stop_server.assert_called_once_with(process)
