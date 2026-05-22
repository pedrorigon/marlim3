"""
Testes de regressão — executa todos os demos e compara os DataFrames
de perfis e tendências contra referências armazenadas.

Os DataFrames de referência ficam em tests/comparison/<model-name>/ e
são gerados/atualizados pelo script scripts/update_regression_references.py.

Execução:
    pytest tests/test_regression.py -v -m regressao
"""

import os
import shutil

import numpy as np
import pandas as pd
import pytest

import marlim3
from marlim3._download import executable_exists

# ============================================================================
# Verificação do executável
# ============================================================================

def _executavel_disponivel():
    return executable_exists()


skip_sem_executavel = pytest.mark.skipif(
    not _executavel_disponivel(),
    reason="Executável Marlim3 não encontrado",
)

# ============================================================================
# Caminhos
# ============================================================================

TESTS_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.normpath(os.path.join(TESTS_DIR, ".."))
DEMOS_DIR = os.path.join(PROJECT_ROOT, "demos")
COMPARISON_DIR = os.path.join(TESTS_DIR, "comparison")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "regression_output")

# ============================================================================
# Definição dos demos e seus arquivos auxiliares
# ============================================================================

DEMOS = {
    "2zonas-2VGLs-2-Check-correcPerfTerm": {
        "json": "2zonas-2VGLs-2-Check-correcPerfTerm.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "BCS-longo-eficMotor": {
        "json": "BCS-longo-eficMotor.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "injec-Liq-TempoResidencia": {
        "json": "injec-Liq-TempoResidencia.json",
        "aux": [],
    },
    # "MultiBCS": {
    #     "json": "MultiBCS.json",
    #     "aux": ["PVTSIM-MARLIM.tab"],
    # },
    "parada-longo-Combinado-BCS-GLC-PIG-completo": {
        "json": "parada-longo-Combinado-BCS-GLC-PIG-completo.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "producaoSimplificado": {
        "json": "producaoSimplificado.json",
        "aux": [],
    },
}

DEMO_IDS = sorted(DEMOS.keys())

# ============================================================================
# Tolerâncias
# ============================================================================

RTOL = 1e-5   # tolerância relativa para comparação numérica
ATOL = 1e-8   # tolerância absoluta para comparação numérica

# ============================================================================
# Utilitários
# ============================================================================

def _run_simulation(model_name):
    """Executa a simulação de um demo e retorna o dicionário de resultados
    no mesmo formato que a GUI (perfilProducao, perfilServico, tendP, tendS).
    """
    info = DEMOS[model_name]
    json_path = os.path.join(DEMOS_DIR, info["json"])

    caso = marlim3.Tramo()
    caso.from_json(json_path)

    # Forçar modo permanente para acelerar
    if "transiente" in (caso.configuracaoInicial or {}):
        caso.configuracaoInicial["transiente"] = False

    out_dir = os.path.join(OUTPUT_DIR, model_name)
    os.makedirs(out_dir, exist_ok=True)

    # Copiar arquivos auxiliares
    for aux in info["aux"]:
        src = os.path.join(DEMOS_DIR, aux)
        if os.path.isfile(src):
            shutil.copy2(src, out_dir)

    results_dir = os.path.join(out_dir, f"resultados_{model_name}")

    original_cwd = os.getcwd()
    try:
        os.chdir(out_dir)
        caso.simular(label=model_name, diretorio=results_dir)
    finally:
        os.chdir(original_cwd)

    # Processar resultados da mesma forma que a GUI
    resultados = {}
    cfg = caso.configuracaoInicial or {}
    has_linha_gas = cfg.get("linhaGas", False)

    if caso.perfilProducao is not None:
        perfis = caso.processar_perfis(results_dir)
        if perfis is not None:
            resultados["perfilProducao"] = perfis

    if has_linha_gas and caso.perfilServico is not None:
        perfis_s = caso.processar_perfis(results_dir, linha="servico")
        if perfis_s is not None:
            resultados["perfilServico"] = perfis_s

    if caso.tendP is not None:
        tend = caso.processar_tendencias(results_dir)
        if tend:
            resultados["tendP"] = tend

    if has_linha_gas and caso.tendS is not None:
        tend_s = caso.processar_tendencias(results_dir, linha="servico")
        if tend_s:
            resultados["tendS"] = tend_s

    return resultados


def _save_results(model_name, resultados, dest_dir=None):
    """Salva os DataFrames de resultados como CSV em dest_dir/model_name/."""
    if dest_dir is None:
        dest_dir = COMPARISON_DIR
    model_dir = os.path.join(dest_dir, model_name)
    os.makedirs(model_dir, exist_ok=True)

    for key, data in resultados.items():
        if isinstance(data, pd.DataFrame):
            csv_path = os.path.join(model_dir, f"{key}.csv")
            data.to_csv(csv_path)
        elif isinstance(data, dict):
            # Dict de DataFrames (tendências)
            for sub_key, sub_df in data.items():
                if isinstance(sub_df, pd.DataFrame):
                    csv_path = os.path.join(model_dir, f"{key}_{sub_key}.csv")
                    sub_df.to_csv(csv_path)


def _load_reference(model_name):
    """Carrega os DataFrames de referência de tests/comparison/model_name/."""
    model_dir = os.path.join(COMPARISON_DIR, model_name)
    if not os.path.isdir(model_dir):
        return None

    reference = {}
    csv_files = sorted(f for f in os.listdir(model_dir) if f.endswith(".csv"))

    for csv_file in csv_files:
        csv_path = os.path.join(model_dir, csv_file)
        name = csv_file[:-4]  # remove .csv

        # Detectar se é tendência (formato: tendP_1, tendS_2, etc.)
        for prefix in ("tendP_", "tendS_"):
            if name.startswith(prefix):
                tend_key = prefix[:-1]  # tendP ou tendS
                sub_key = int(name[len(prefix):])
                df = pd.read_csv(csv_path, index_col=0)
                if tend_key not in reference:
                    reference[tend_key] = {}
                reference[tend_key][sub_key] = df
                break
        else:
            # É perfil (perfilProducao, perfilServico)
            df = pd.read_csv(csv_path, index_col=[0, 1])
            reference[name] = df

    return reference if reference else None


def _compare_dataframes(df_actual, df_ref, label):
    """Compara dois DataFrames numéricos com tolerância."""
    # Alinhar colunas
    common_cols = sorted(set(df_actual.columns) & set(df_ref.columns))
    missing_in_actual = set(df_ref.columns) - set(df_actual.columns)
    missing_in_ref = set(df_actual.columns) - set(df_ref.columns)

    assert not missing_in_actual, (
        f"[{label}] Colunas faltando no resultado atual: {missing_in_actual}"
    )
    assert not missing_in_ref, (
        f"[{label}] Colunas novas no resultado atual (não presentes na referência): {missing_in_ref}"
    )

    # Comparar shape
    assert df_actual.shape == df_ref.shape, (
        f"[{label}] Shape diferente: atual={df_actual.shape}, referência={df_ref.shape}"
    )

    # Comparar valores numéricos
    for col in common_cols:
        if pd.api.types.is_numeric_dtype(df_actual[col]) and pd.api.types.is_numeric_dtype(df_ref[col]):
            actual_vals = df_actual[col].values.astype(float)
            ref_vals = df_ref[col].values.astype(float)

            # Pular colunas que são todas NaN em ambos
            if np.all(np.isnan(actual_vals)) and np.all(np.isnan(ref_vals)):
                continue

            np.testing.assert_allclose(
                actual_vals, ref_vals,
                rtol=RTOL, atol=ATOL,
                err_msg=f"[{label}] Coluna '{col}' diverge da referência",
            )
        else:
            # Colunas não numéricas: comparação exata
            pd.testing.assert_series_equal(
                df_actual[col].reset_index(drop=True),
                df_ref[col].reset_index(drop=True),
                check_names=False,
                obj=f"[{label}] Coluna '{col}'",
            )


# ============================================================================
# Testes
# ============================================================================

@skip_sem_executavel
@pytest.mark.regressao
@pytest.mark.parametrize("model_name", DEMO_IDS, ids=DEMO_IDS)
def test_regression(model_name):
    """Executa o demo, gera DataFrames e compara com a referência armazenada."""
    ref_dir = os.path.join(COMPARISON_DIR, model_name)
    if not os.path.isdir(ref_dir) or not os.listdir(ref_dir):
        pytest.skip(
            f"Reference not found for '{model_name}'. "
            f"Run: python tests/update_regression_references.py"
        )

    # Executar simulação
    resultados = _run_simulation(model_name)
    assert resultados, f"Simulação de '{model_name}' não gerou resultados"

    # Carregar referência
    reference = _load_reference(model_name)
    assert reference is not None, (
        f"Falha ao carregar referência de '{model_name}'"
    )

    # Comparar cada resultado
    for key in reference:
        assert key in resultados, (
            f"Resultado '{key}' presente na referência mas ausente na simulação atual"
        )

        ref_data = reference[key]
        actual_data = resultados[key]

        if isinstance(ref_data, pd.DataFrame) and isinstance(actual_data, pd.DataFrame):
            _compare_dataframes(actual_data, ref_data, f"{model_name}/{key}")
        elif isinstance(ref_data, dict) and isinstance(actual_data, dict):
            for sub_key in ref_data:
                assert sub_key in actual_data, (
                    f"Sub-resultado '{key}/{sub_key}' ausente na simulação atual"
                )
                _compare_dataframes(
                    actual_data[sub_key], ref_data[sub_key],
                    f"{model_name}/{key}/{sub_key}",
                )
