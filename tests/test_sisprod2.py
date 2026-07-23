"""
test_sisprod2.py — Integração da nova arquitetura do SisProd na suíte de testes.

Este teste compila o harness ``tests/sisprod2_selftest.cpp`` junto com a
implementação da nova arquitetura ``src/core/SisProd.cpp`` (usando o cabeçalho
``src/include/SisProd2.h``) e executa a sua suíte interna, que cobre:

  - o root finder genérico (Brent) contra a referência independente (bisseção);
  - a estratégia de mistura de correntes (StreamMixing);
  - o solver de regime permanente contra a solução hidrostática analítica;
  - o determinismo e o round trip do solver;
  - o lote paralelo com OpenMP igualando o resultado sequencial;
  - a comparação automática entre implementações.

Também verifica a organização dos arquivos: o ``SisProd`` legado preservado em
``SisProd_old.cpp`` (referência para os testes de integração) e a nova
arquitetura em ``SisProd.cpp`` + ``SisProd2.h``.

Execução:
    pytest tests/test_sisprod2.py -v
"""

import os
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEW_IMPL = os.path.join(ROOT, "src", "core", "SisProd.cpp")
NEW_IMPL_R03 = os.path.join(ROOT, "src", "core", "SisProd_r03.cpp")
NEW_IMPL_R04 = os.path.join(ROOT, "src", "core", "SisProd_r04.cpp")
NEW_IMPL_R05 = os.path.join(ROOT, "src", "core", "SisProd_r05.cpp")
NEW_IMPL_R06 = os.path.join(ROOT, "src", "core", "SisProd_r06.cpp")
NEW_HEADER = os.path.join(ROOT, "src", "include", "SisProd2.h")
SELFTEST = os.path.join(ROOT, "tests", "sisprod2_selftest.cpp")
INCLUDE_DIR = os.path.join(ROOT, "src", "include")
LEGACY_SRC = os.path.join(ROOT, "src", "core", "SisProd_old.cpp")
LEGACY_HEADER = os.path.join(ROOT, "src", "include", "SisProd.h")


def _compiler():
    for candidate in ("g++", "c++", "clang++"):
        if shutil.which(candidate):
            return candidate
    return None


skip_sem_compilador = pytest.mark.skipif(
    _compiler() is None,
    reason="Nenhum compilador C++ (g++/c++/clang++) encontrado",
)


def test_legacy_sisprod_preservado():
    """O SisProd legado deve permanecer intacto em SisProd_old.cpp."""
    assert os.path.isfile(LEGACY_SRC), "SisProd_old.cpp legado não encontrado"
    assert os.path.isfile(LEGACY_HEADER), "SisProd.h legado não encontrado"
    with open(LEGACY_SRC, "r", encoding="utf-8", errors="replace") as handle:
        linhas = sum(1 for _ in handle)
    assert linhas > 20000, "SisProd_old.cpp legado parece ter sido alterado"


def test_novo_modulo_existe():
    """A nova arquitetura deve estar em SisProd.cpp + SisProd2.h."""
    assert os.path.isfile(NEW_IMPL), "SisProd.cpp (nova arquitetura) não encontrado"
    assert os.path.isfile(NEW_HEADER), "SisProd2.h não encontrado"
    assert os.path.isfile(SELFTEST), "harness sisprod2_selftest.cpp não encontrado"
    assert os.path.isfile(NEW_IMPL_R03), "SisProd_r03.cpp (R03 C0/Ud) não encontrado"
    assert os.path.isfile(NEW_IMPL_R04), "SisProd_r04.cpp (R04 C0/Ud) não encontrado"
    assert os.path.isfile(NEW_IMPL_R05), "SisProd_r05.cpp (R05 C0/Ud) não encontrado"
    assert os.path.isfile(NEW_IMPL_R06), "SisProd_r06.cpp (R06/R10 Colebrook+GLV) não encontrado"


@skip_sem_compilador
def test_suite_da_nova_arquitetura(tmp_path):
    """Compila e executa a suíte interna da nova arquitetura (deve passar 100%)."""
    binario = os.path.join(str(tmp_path), "sisprod2_tests")
    compilar = [
        _compiler(),
        "-std=c++11",
        "-O2",
        "-Wall",
        "-fopenmp",
        f"-I{INCLUDE_DIR}",
        SELFTEST,
        NEW_IMPL,
        NEW_IMPL_R03,
        NEW_IMPL_R04,
        NEW_IMPL_R05,
        NEW_IMPL_R06,
        "-o",
        binario,
    ]
    proc = subprocess.run(compilar, capture_output=True, text=True)
    assert proc.returncode == 0, f"Falha ao compilar a nova arquitetura:\n{proc.stderr}"

    executar = subprocess.run([binario], capture_output=True, text=True)
    assert executar.returncode == 0, (
        "A suíte da nova arquitetura falhou:\n" + executar.stdout + executar.stderr
    )
    assert "ALL TESTS PASSED" in executar.stdout
    assert "[FAIL]" not in executar.stdout, executar.stdout


@skip_sem_compilador
def test_flag_selecao_implementacao(tmp_path):
    """A flag MARLIM_USE_NEW_SISPROD deve alternar a implementação reportada."""

    def build_and_report(defines, label):
        binario = os.path.join(str(tmp_path), "flag_" + label)
        comando = [_compiler(), "-std=c++11", "-O2", "-fopenmp", f"-I{INCLUDE_DIR}"]
        comando += list(defines)
        comando += [SELFTEST, NEW_IMPL, NEW_IMPL_R03, NEW_IMPL_R04, NEW_IMPL_R05, NEW_IMPL_R06, "-o", binario]
        proc = subprocess.run(comando, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        run = subprocess.run([binario], capture_output=True, text=True)
        return run.stdout

    legacy_out = build_and_report([], "legacy")
    new_out = build_and_report(["-DMARLIM_USE_NEW_SISPROD"], "new")
    assert "active implementation: legacy" in legacy_out
    assert "active implementation: new" in new_out
