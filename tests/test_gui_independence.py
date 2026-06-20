from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_streamlit_gui_does_not_import_desktop_wrapper():
    source = (ROOT / "gui" / "app.py").read_text(encoding="utf-8")

    assert "import marlim3_desktop" not in source
    assert "from marlim3_desktop" not in source
    assert "sys._MEIPASS" not in source
    assert "MARLIM3_DESKTOP_ROOT" not in source


def test_desktop_wrapper_consumes_gui_as_a_data_file():
    spec = (
        ROOT / "marlim3_desktop" / "packaging" / "marlim3-desktop.spec"
    ).read_text(encoding="utf-8")

    assert 'root / "gui" / "app.py"' in spec
    assert 'desktop_root / "__main__.py"' in spec


def test_relevant_simulator_launches_do_not_use_shell():
    sources = [
        ROOT / "gui" / "app.py",
        ROOT / "marlim3" / "_tramo" / "_branch.py",
        ROOT / "marlim3" / "_tramo" / "_tramo.py",
        ROOT / "marlim3" / "_rede" / "_rede.py",
    ]

    for source_path in sources:
        assert "shell=True" not in source_path.read_text(encoding="utf-8-sig")
