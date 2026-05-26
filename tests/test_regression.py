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
        "json": "pt-br/2zonas-2VGLs-2-Check-correcPerfTerm.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "BCS-longo-eficMotor": {
        "json": "pt-br/BCS-longo-eficMotor.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "injec-Liq-TempoResidencia": {
        "json": "pt-br/injec-Liq-TempoResidencia.json",
        "aux": [],
    },
    # "MultiBCS": {
    #     "json": "pt-br/MultiBCS.json",
    #     "aux": ["PVTSIM-MARLIM.tab"],
    # },
    "parada-longo-Combinado-BCS-GLC-PIG-completo": {
        "json": "pt-br/parada-longo-Combinado-BCS-GLC-PIG-completo.json",
        "aux": ["PVTSIM-MARLIM.tab"],
    },
    "producaoSimplificado": {
        "json": "pt-br/producaoSimplificado.json",
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
    no mesmo formato que a GUI (productionProfile, serviceProfile, productionTrend, serviceTrend).
    """
    info = DEMOS[model_name]
    json_path = os.path.join(DEMOS_DIR, info["json"])

    caso = marlim3.Branch()
    caso.from_json(json_path)

    # Forçar modo permanente para acelerar
    if "transient" in (caso.initialConfig or {}):
        caso.initialConfig["transient"] = False

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
        caso.simulate(label=model_name, directory=results_dir)
    finally:
        os.chdir(original_cwd)

    # Processar resultados da mesma forma que a GUI
    resultados = {}
    cfg = caso.initialConfig or {}
    has_linha_gas = cfg.get("gasLine", False)

    if caso.productionProfile is not None:
        perfis = caso._process_profiles(results_dir)
        if perfis is not None:
            resultados["productionProfile"] = perfis

    if has_linha_gas and caso.serviceProfile is not None:
        perfis_s = caso._process_profiles(results_dir, line="service")
        if perfis_s is not None:
            resultados["serviceProfile"] = perfis_s

    if caso.productionTrend is not None:
        tend = caso._process_trends(results_dir)
        if tend:
            resultados["productionTrend"] = tend

    if has_linha_gas and caso.serviceTrend is not None:
        tend_s = caso._process_trends(results_dir, line="service")
        if tend_s:
            resultados["serviceTrend"] = tend_s

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

    profile_key_map = {
        "perfilProducao": "productionProfile",
        "perfilServico": "serviceProfile",
        "productionProfile": "productionProfile",
        "serviceProfile": "serviceProfile",
    }
    trend_prefix_map = {
        "tendP_": "productionTrend",
        "tendS_": "serviceTrend",
        "productionTrend_": "productionTrend",
        "serviceTrend_": "serviceTrend",
    }

    reference = {}
    csv_files = sorted(f for f in os.listdir(model_dir) if f.endswith(".csv"))

    for csv_file in csv_files:
        csv_path = os.path.join(model_dir, csv_file)
        name = csv_file[:-4]  # remove .csv

        # Detectar se é tendência (aceita nomes antigos e novos)
        for prefix, tend_key in trend_prefix_map.items():
            if name.startswith(prefix):
                sub_key = int(name[len(prefix):])
                df = pd.read_csv(csv_path, index_col=0)
                if tend_key not in reference:
                    reference[tend_key] = {}
                reference[tend_key][sub_key] = df
                break
        else:
            # É perfil (aceita nomes antigos e novos)
            df = pd.read_csv(csv_path, index_col=[0, 1])
            mapped_name = profile_key_map.get(name, name)
            reference[mapped_name] = df

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
