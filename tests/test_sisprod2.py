"""
test_sisprod2.py — Integração da nova arquitetura do SisProd na suíte de testes.

Este teste compila o harness ``tests/sisprod2_selftest.cpp`` junto com a
implementação da nova arquitetura e executa a sua suíte interna, que cobre:

  - o root finder genérico (Brent) contra a referência independente (bisseção);
  - a estratégia de mistura de correntes (StreamMixing);
  - o solver de regime permanente contra a solução hidrostática analítica;
  - o determinismo e o round trip do solver;
  - o lote paralelo com OpenMP igualando o resultado sequencial;
  - a comparação automática entre implementações.

Estrutura flat (2026-07):
  - src/core/SisProd.cpp (núcleo)
  - src/core/BlackOilProperties.cpp (R04)
  - src/core/DriftFluxCorrelations.cpp (R03)
  - src/core/HydraulicFriction.cpp (R06/R10)
  - src/core/PressureGradientEngine.cpp (R05)
  - src/include/SisProd.h (interface nova)
  - src/include/SisProd_old.h (interface legada)

Execução:
    pytest tests/test_sisprod2.py -v
"""

import os
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Nova estrutura flat (2026-07)
NEW_IMPL_CORE = os.path.join(ROOT, "src", "core", "SisProd.cpp")
NEW_IMPL_R03 = os.path.join(ROOT, "src", "core", "DriftFluxCorrelations.cpp")
NEW_IMPL_R04 = os.path.join(ROOT, "src", "core", "BlackOilProperties.cpp")
NEW_IMPL_R05 = os.path.join(ROOT, "src", "core", "PressureGradientEngine.cpp")
NEW_IMPL_R06 = os.path.join(ROOT, "src", "core", "HydraulicFriction.cpp")
NEW_HEADER = os.path.join(ROOT, "src", "include", "SisProd.h")
SELFTEST = os.path.join(ROOT, "tests", "sisprod2_selftest.cpp")
INCLUDE_DIR = os.path.join(ROOT, "src", "include")
LEGACY_SRC = os.path.join(ROOT, "src", "core", "SisProd_old.cpp")
LEGACY_HEADER = os.path.join(ROOT, "src", "include", "SisProd_old.h")

# Lista de todos os módulos novos
ALL_NEW_MODULES = [
    NEW_IMPL_CORE,
    NEW_IMPL_R03,
    NEW_IMPL_R04,
    NEW_IMPL_R05,
    NEW_IMPL_R06,
]


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
    """O SisProd legado deve permanecer intacto em SisProd_old.cpp/.h."""
    assert os.path.isfile(LEGACY_SRC), "SisProd_old.cpp legado não encontrado"
    assert os.path.isfile(LEGACY_HEADER), "SisProd_old.h legado não encontrado"
    with open(LEGACY_SRC, "r", encoding="utf-8", errors="replace") as handle:
        linhas = sum(1 for _ in handle)
    assert linhas > 20000, "SisProd_old.cpp legado parece ter sido alterado"


def test_nova_arquitetura_existe():
    """A nova arquitetura deve estar presente na estrutura flat."""
    assert os.path.isfile(NEW_IMPL_CORE), "SisProd.cpp (núcleo) não encontrado"
    assert os.path.isfile(NEW_HEADER), "SisProd.h não encontrado"
    assert os.path.isfile(SELFTEST), "harness sisprod2_selftest.cpp não encontrado"
    
    # Verificar módulos R03-R06
    assert os.path.isfile(NEW_IMPL_R03), "DriftFluxCorrelations.cpp (R03 C0/Ud) não encontrado"
    assert os.path.isfile(NEW_IMPL_R04), "BlackOilProperties.cpp (R04 black-oil) não encontrado"
    assert os.path.isfile(NEW_IMPL_R05), "PressureGradientEngine.cpp (R05 gradient) não encontrado"
    assert os.path.isfile(NEW_IMPL_R06), "HydraulicFriction.cpp (R06/R10 friction+GLV) não encontrado"


def test_headers_modulares_existem():
    """Headers modulares opcionais devem existir para documentação."""
    modular_headers = [
        os.path.join(ROOT, "src", "include", "BlackOilProperties.h"),
        os.path.join(ROOT, "src", "include", "DriftFluxCorrelations.h"),
        os.path.join(ROOT, "src", "include", "HydraulicFriction.h"),
        os.path.join(ROOT, "src", "include", "PressureGradientEngine.h"),
    ]
    for header in modular_headers:
        assert os.path.isfile(header), f"Header modular {os.path.basename(header)} não encontrado"


@skip_sem_compilador
def test_modulos_compilam_individualmente(tmp_path):
    """Cada módulo deve compilar isoladamente (teste incremental)."""
    for modulo in ALL_NEW_MODULES:
        nome = os.path.basename(modulo).replace(".cpp", "")
        binario_obj = os.path.join(str(tmp_path), f"{nome}.o")
        
        compilar = [
            _compiler(),
            "-std=c++11",
            "-O2",
            "-Wall",
            "-c",
            "-fopenmp",
            f"-I{INCLUDE_DIR}",
            modulo,
            "-o",
            binario_obj,
        ]
        proc = subprocess.run(compilar, capture_output=True, text=True)
        assert proc.returncode == 0, (
            f"Falha ao compilar {os.path.basename(modulo)}:\n{proc.stderr}"
        )
        assert os.path.isfile(binario_obj), f"Objeto {nome}.o não gerado"


@skip_sem_compilador
def test_suite_da_nova_arquitetura(tmp_path):
    """Compila e executa a suíte interna da nova arquitetura."""
    binario = os.path.join(str(tmp_path), "sisprod2_tests")
    compilar = [
        _compiler(),
        "-std=c++11",
        "-O2",
        "-Wall",
        "-fopenmp",
        f"-I{INCLUDE_DIR}",
        SELFTEST,
    ] + ALL_NEW_MODULES + [
        "-o",
        binario,
    ]
    proc = subprocess.run(compilar, capture_output=True, text=True)
    assert proc.returncode == 0, f"Falha ao compilar a nova arquitetura:\n{proc.stderr}"

    executar = subprocess.run([binario], capture_output=True, text=True)
    assert executar.returncode == 0, (
        "A suíte da nova arquitetura falhou:\n" + executar.stdout + executar.stderr
    )
    assert "ALL TESTS PASSED" in executar.stdout or "[PASS]" in executar.stdout
    assert "[FAIL]" not in executar.stdout, executar.stdout
