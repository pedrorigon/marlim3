"""
Testes unitários do pacote marlim3.

Baseados no tutorial comparacao_horizontal_vertical.ipynb,
verificam a montagem de modelos (Tramo), serialização JSON e
a classe de cenários (Cenarios).
"""

import copy
import json
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

import marlim3


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def fluido():
    """Fluido black-oil padrão."""
    return {
        "id": 0,
        "api": 32,
        "gor": 100,
        "gasDensity": 0.7,
        "bsw": 0.0,
    }


@pytest.fixture
def material_aco():
    """Material aço."""
    return {
        "id": 0,
        "type": 0,
        "conductivity": 58,
        "specificHeat": 480,
        "rho": 7850,
    }


@pytest.fixture
def secao_transversal():
    """Seção transversal com uma camada de aço."""
    camada = {
        "materialId": 0,
        "layerMeasurementType": "THICKNESS",
        "thickness": 0.0254,
    }
    return {
        "id": 0,
        "innerDiameter": 10 * 0.0254,
        "roughness": 0.183e-3,
        "layer": [camada],
    }


@pytest.fixture
def duto_horizontal():
    """Duto horizontal de 2500 m, 20 células."""
    ncel = 20
    length_total = 2500
    linha = {
        "numCells": ncel,
        "length": length_total / ncel,
    }
    condicoes_ambiente = {
        "measuredPosition": [0, 1],
        "ambientTemp": [40, 20],
        "ambientVel": [0.5, 0.5],
    }
    return {
        "id": 0,
        "crossSectionId": 0,
        "ambientEnvironment": 2,
        "angle": 0,
        "discretization": [linha],
        "initialAndAmbientConditions": condicoes_ambiente,
    }


@pytest.fixture
def fonte_liquido():
    """Fonte de líquido a montante – 1500 sm3/d."""
    return {
        "id": 0,
        "prodFluidId": 0,
        "measuredLength": 0.1,
        "time": [0],
        "liquidFlowRate": [1500],
        "temperature": [40],
    }


@pytest.fixture
def separador():
    """Condição de contorno de pressão a jusante."""
    return {
        "time": [0],
        "pressure": [2],
    }


@pytest.fixture
def vars_saida():
    """Variáveis de saída para perfis."""
    nomes = ["pressure", "temperature", "holdup", "flowPattern", "frictionPressureGradient", "hydrostaticPressure"]
    return {"time": [0]} | {var: True for var in nomes}


@pytest.fixture
def caso_horizontal(fluido, material_aco, secao_transversal,
                    duto_horizontal, fonte_liquido, separador, vars_saida):
    """Caso base horizontal completo (sem simulação)."""
    caso = marlim3.Branch()
    caso.productionFluid = [fluido]
    caso.material = [material_aco]
    caso.crossSection = [secao_transversal]
    caso.productionPipe = [duto_horizontal]
    caso.liquidSource = [fonte_liquido]
    caso.separator = separador
    caso.productionProfile = vars_saida
    return caso


@pytest.fixture
def caso_vertical(caso_horizontal):
    """Caso base vertical (cópia do horizontal com ângulo π/2)."""
    caso = copy.deepcopy(caso_horizontal)
    caso.productionPipe[0]["angle"] = np.pi / 2
    return caso


# ============================================================================
# Testes de construção do Tramo
# ============================================================================

class TestTramoConstruction:
    """Testes de instanciação e atributos do Tramo."""

    def test_tramo_default(self):
        """Tramo vazio deve ter atributos padrão."""
        t = marlim3.Branch()
        assert t.system == "MULTIFASICO"
        assert t.jsonVersion == "1.3.9"
        assert t.productionFluid == []
        assert t.material == []
        assert t.crossSection == []
        assert t.productionPipe == []
        assert t.liquidSource == []
        assert t.separator is None
        assert t.productionProfile == {}
        assert isinstance(t.resultados, dict)

    def test_tramo_default_lists_are_independent(self):
        """Listas padrão não devem ser compartilhadas entre instâncias."""
        t1 = marlim3.Branch()
        t2 = marlim3.Branch()
        t1.productionFluid.append({"test": 1})
        assert t2.productionFluid == []

    def test_tramo_with_kwargs(self, fluido, material_aco):
        """Tramo pode ser criado passando kwargs no construtor."""
        t = marlim3.Branch(
            productionFluid=[fluido],
            material=[material_aco],
        )
        assert len(t.productionFluid) == 1
        assert t.productionFluid[0]["api"] == 32
        assert len(t.material) == 1
        assert t.material[0]["conductivity"] == 58

    def test_system_injetor(self):
        """Tramo pode ser criado com system INJETOR."""
        t = marlim3.Branch(system="INJETOR")
        assert t.system == "INJETOR"


# ============================================================================
# Testes de montagem do modelo (padrão do notebook)
# ============================================================================

class TestModelAssembly:
    """Testes que reproduzem a montagem do modelo do notebook."""

    def test_fluido_assignment(self, caso_horizontal, fluido):
        assert caso_horizontal.productionFluid == [fluido]
        assert caso_horizontal.productionFluid[0]["api"] == 32
        assert caso_horizontal.productionFluid[0]["gor"] == 100
        assert caso_horizontal.productionFluid[0]["gasDensity"] == 0.7
        assert caso_horizontal.productionFluid[0]["bsw"] == 0.0

    def test_material_assignment(self, caso_horizontal, material_aco):
        assert caso_horizontal.material == [material_aco]

    def test_secao_transversal(self, caso_horizontal):
        sec = caso_horizontal.crossSection[0]
        assert sec["innerDiameter"] == pytest.approx(10 * 0.0254)
        assert sec["roughness"] == pytest.approx(0.183e-3)
        assert len(sec["layer"]) == 1
        assert sec["layer"][0]["thickness"] == pytest.approx(0.0254)

    def test_duto_horizontal(self, caso_horizontal):
        duto = caso_horizontal.productionPipe[0]
        assert duto["angle"] == 0
        assert duto["ambientEnvironment"] == 2
        assert duto["discretization"][0]["numCells"] == 20
        assert duto["discretization"][0]["length"] == pytest.approx(125.0)

    def test_condicoes_ambiente(self, caso_horizontal):
        cond = caso_horizontal.productionPipe[0]["initialAndAmbientConditions"]
        assert cond["measuredPosition"] == [0, 1]
        assert cond["ambientTemp"] == [40, 20]
        assert cond["ambientVel"] == [0.5, 0.5]

    def test_fonte_liquido(self, caso_horizontal):
        fl = caso_horizontal.liquidSource[0]
        assert fl["liquidFlowRate"] == [1500]
        assert fl["temperature"] == [40]
        assert fl["measuredLength"] == pytest.approx(0.1)

    def test_separador(self, caso_horizontal):
        assert caso_horizontal.separator["pressure"] == [2]

    def test_perfil_producao_vars(self, caso_horizontal):
        pp = caso_horizontal.productionProfile
        assert pp["time"] == [0]
        for var in ["pressure", "temperature", "holdup", "flowPattern", "frictionPressureGradient", "hydrostaticPressure"]:
            assert pp[var] is True

    def test_length_total(self, caso_horizontal):
        """Comprimento total = numCells * length_celula."""
        duto = caso_horizontal.productionPipe[0]
        disc = duto["discretization"][0]
        total = disc["numCells"] * disc["length"]
        assert total == pytest.approx(2500.0)


# ============================================================================
# Testes de criação do caso vertical (deepcopy + alteração)
# ============================================================================

class TestVerticalCase:
    """Testes da criação do caso vertical a partir do horizontal."""

    def test_angle_vertical(self, caso_vertical):
        assert caso_vertical.productionPipe[0]["angle"] == pytest.approx(np.pi / 2)

    def test_deepcopy_independence(self, caso_horizontal, caso_vertical):
        """Alterar o caso vertical não altera o horizontal."""
        assert caso_horizontal.productionPipe[0]["angle"] == 0
        assert caso_vertical.productionPipe[0]["angle"] == pytest.approx(np.pi / 2)

    def test_vertical_preserves_other_fields(self, caso_horizontal, caso_vertical):
        """O deepcopy mantém os demais campos intactos."""
        assert caso_vertical.productionFluid == caso_horizontal.productionFluid
        assert caso_vertical.material == caso_horizontal.material
        assert caso_vertical.separator == caso_horizontal.separator


# ============================================================================
# Testes de variação de parâmetros (padrão do notebook – loop de vazões)
# ============================================================================

class TestParameterVariation:
    """Testes para variação de parâmetros como no notebook."""

    def test_variacao_vazao(self, caso_horizontal):
        vazoes = [200, 2000, 4000, 6000]
        casos = []
        for vazao in vazoes:
            c = copy.deepcopy(caso_horizontal)
            c.liquidSource[0]["liquidFlowRate"] = [vazao]
            casos.append(c)

        for i, vazao in enumerate(vazoes):
            assert casos[i].liquidSource[0]["liquidFlowRate"] == [vazao]

        # Caso original não foi alterado
        assert caso_horizontal.liquidSource[0]["liquidFlowRate"] == [1500]

    def test_variacao_diametro(self, caso_horizontal):
        diametros_pol = [5, 10, 15, 20]
        casos = []
        for d in diametros_pol:
            c = copy.deepcopy(caso_horizontal)
            c.crossSection[0]["innerDiameter"] = d * 0.0254
            casos.append(c)

        for i, d in enumerate(diametros_pol):
            assert casos[i].crossSection[0]["innerDiameter"] == pytest.approx(
                d * 0.0254
            )

    def test_variacao_gor(self, caso_horizontal):
        gors = [0, 50, 200, 500, 1000, 2000]
        casos = []
        for gor in gors:
            c = copy.deepcopy(caso_horizontal)
            c.productionFluid[0]["gor"] = gor
            casos.append(c)

        for i, gor in enumerate(gors):
            assert casos[i].productionFluid[0]["gor"] == gor

    def test_variacao_api(self, caso_horizontal):
        apis = [15, 20, 25]
        for api_val in apis:
            c = copy.deepcopy(caso_horizontal)
            c.productionFluid[0]["api"] = api_val
            assert c.productionFluid[0]["api"] == api_val

    def test_adicionar_camada_isolante(self, caso_horizontal):
        """Adiciona uma camada isolante à seção transversal."""
        c = copy.deepcopy(caso_horizontal)
        camada_isolante = {
            "materialId": 1,
            "layerMeasurementType": "THICKNESS",
            "thickness": 0.05,
        }
        c.crossSection[0]["layer"].append(camada_isolante)
        assert len(c.crossSection[0]["layer"]) == 2
        # Original não foi alterado
        assert len(caso_horizontal.crossSection[0]["layer"]) == 1


# ============================================================================
# Testes de serialização JSON (to_json / from_json)
# ============================================================================

class TestJsonSerialization:
    """Testes de serialização e desserialização JSON."""

    def test_to_json_creates_file(self, caso_horizontal):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_horizontal.to_json("test_model")
                assert os.path.isfile(os.path.join(tmpdir, "test_model.json"))
            finally:
                os.chdir(original_cwd)

    def test_to_json_valid_json(self, caso_horizontal):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_horizontal.to_json("test_model")
                with open(os.path.join(tmpdir, "test_model.json"), "r") as f:
                    data = json.load(f)
                assert isinstance(data, dict)
            finally:
                os.chdir(original_cwd)

    def test_to_json_contains_expected_keys(self, caso_horizontal):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_horizontal.to_json("test_model")
                with open(os.path.join(tmpdir, "test_model.json"), "r") as f:
                    data = json.load(f)

                assert "system" in data
                assert "productionFluid" in data
                assert "material" in data
                assert "crossSection" in data
                assert "productionPipe" in data
                assert "liquidSource" in data
                assert "separator" in data
            finally:
                os.chdir(original_cwd)

    def test_roundtrip_json(self, caso_horizontal):
        """to_json → from_json deve preservar os dados."""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_horizontal.to_json("roundtrip")

                novo = marlim3.Branch()
                novo.from_json(os.path.join(tmpdir, "roundtrip.json"))

                assert novo.system == caso_horizontal.system
                assert novo.productionFluid == caso_horizontal.productionFluid
                assert novo.material == caso_horizontal.material
                assert novo.crossSection == caso_horizontal.crossSection
                assert novo.productionPipe == caso_horizontal.productionPipe
                assert novo.liquidSource == caso_horizontal.liquidSource
                assert novo.separator == caso_horizontal.separator
            finally:
                os.chdir(original_cwd)

    def test_roundtrip_vertical(self, caso_vertical):
        """Caso vertical sobrevive a roundtrip JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_vertical.to_json("roundtrip_v")
                novo = marlim3.Branch()
                novo.from_json(os.path.join(tmpdir, "roundtrip_v.json"))
                assert novo.productionPipe[0]["angle"] == pytest.approx(np.pi / 2)
            finally:
                os.chdir(original_cwd)

    def test_to_json_empty_fields(self):
        """to_json com generate_empty_fields gera campos vazios."""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                t = marlim3.Branch()
                t.to_json("empty_model", generate_empty_fields=True)
                with open(os.path.join(tmpdir, "empty_model.json"), "r") as f:
                    data = json.load(f)
                # Com generate_empty_fields=True, campos vazios devem existir
                assert "system" in data
            finally:
                os.chdir(original_cwd)

    def test_from_json_with_dict(self, caso_horizontal):
        """from_json aceita dicionário via is_string=True."""
        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                caso_horizontal.to_json("dict_test")
                with open(os.path.join(tmpdir, "dict_test.json"), "r") as f:
                    data = json.load(f)

                novo = marlim3.Branch()
                novo.from_json(data, is_string=True)
                assert novo.productionFluid == caso_horizontal.productionFluid
            finally:
                os.chdir(original_cwd)


# ============================================================================
# Testes da classe Cenarios
# ============================================================================

class TestCenarios:
    """Testes de construção e uso da classe Cenarios."""

    def test_cenarios_default(self):
        c = marlim3.Scenarios()
        assert c.casos == {}

    def test_cenarios_with_cases(self, caso_horizontal, caso_vertical):
        """Cenarios aceita dicionário de Tramos."""
        cenarios = marlim3.Scenarios({
            "horizontal": caso_horizontal,
            "vertical": caso_vertical,
        })
        assert "horizontal" in cenarios.casos
        assert "vertical" in cenarios.casos
        assert len(cenarios.casos) == 2

    def test_cenarios_access(self, caso_horizontal):
        """Acesso ao caso dentro de Cenarios retorna o Tramo."""
        cenarios = marlim3.Scenarios({"caso1": caso_horizontal})
        assert cenarios.casos["caso1"] is caso_horizontal

    def test_cenarios_multiple_vazoes(self, caso_horizontal):
        """Padrão do notebook: criar cenários com múltiplas vazões."""
        vazoes = [200, 2000, 4000, 6000]
        casos = {}
        for vazao in vazoes:
            c = copy.deepcopy(caso_horizontal)
            c.liquidSource[0]["liquidFlowRate"] = [vazao]
            casos[f"Q = {vazao} sm3/d"] = c

        cenarios = marlim3.Scenarios(casos)
        assert len(cenarios.casos) == 4
        assert "Q = 200 sm3/d" in cenarios.casos
        assert "Q = 6000 sm3/d" in cenarios.casos

        for vazao in vazoes:
            label = f"Q = {vazao} sm3/d"
            assert cenarios.casos[label].liquidSource[0]["liquidFlowRate"] == [vazao]


# ============================================================================
# Testes da classe Rede
# ============================================================================

class TestRedeConstruction:
    """Testes básicos de instanciação da classe Rede."""

    def test_rede_default(self):
        r = marlim3.Network()
        assert r.Arquivos == []
        assert r.Conexao == []


# ============================================================================
# Testes de geometria (estrutura do duto)
# ============================================================================

class TestGeometry:
    """Testes de consistência geométrica dos modelos."""

    def test_multiple_dutos(self, caso_horizontal):
        """Modelo com dois dutos concatenados."""
        c = copy.deepcopy(caso_horizontal)
        segundo_duto = copy.deepcopy(c.productionPipe[0])
        segundo_duto["id"] = 1
        segundo_duto["angle"] = np.pi / 4  # 45°
        c.productionPipe.append(segundo_duto)
        assert len(c.productionPipe) == 2
        assert c.productionPipe[0]["angle"] == 0
        assert c.productionPipe[1]["angle"] == pytest.approx(np.pi / 4)

    def test_discretization_consistency(self, caso_horizontal):
        """Verifica que a discretização é consistente."""
        duto = caso_horizontal.productionPipe[0]
        disc = duto["discretization"][0]
        assert disc["numCells"] > 0
        assert disc["length"] > 0

    def test_secao_diametro_positive(self, caso_horizontal):
        assert caso_horizontal.crossSection[0]["innerDiameter"] > 0

    def test_secao_roughness_positive(self, caso_horizontal):
        assert caso_horizontal.crossSection[0]["roughness"] > 0


# ============================================================================
# Testes de validação de dados de entrada
# ============================================================================

class TestInputValidation:
    """Testes de consistência dos dados de entrada do modelo."""

    def test_fluido_api_range(self, fluido):
        """API gravity deve estar num intervalo razoável."""
        assert 5 <= fluido["api"] <= 70

    def test_fluido_gor_nonnegative(self, fluido):
        assert fluido["gor"] >= 0

    def test_fluido_bsw_range(self, fluido):
        assert 0.0 <= fluido["bsw"] <= 1.0

    def test_fluido_gasDensity_range(self, fluido):
        assert 0.5 <= fluido["gasDensity"] <= 1.5

    def test_pressure_separador_positive(self, separador):
        assert all(p > 0 for p in separador["pressure"])

    def test_temperature_fonte_positive(self, fonte_liquido):
        assert all(t > 0 for t in fonte_liquido["temperature"])

    def test_vazao_fonte_positive(self, fonte_liquido):
        assert all(q > 0 for q in fonte_liquido["liquidFlowRate"])

    def test_material_conductivity_positive(self, material_aco):
        assert material_aco["conductivity"] > 0

    def test_material_calor_especifico_positive(self, material_aco):
        assert material_aco["specificHeat"] > 0

    def test_material_rho_positive(self, material_aco):
        assert material_aco["rho"] > 0


# ============================================================================
# Testes de from_json com arquivos demo
# ============================================================================

class TestDemoFiles:
    """Testes de carregamento dos arquivos JSON de demonstração."""

    DEMO_DIR = os.path.join(os.path.dirname(__file__), "..", "demos")

    def _demo_files(self):
        if not os.path.isdir(self.DEMO_DIR):
            pytest.skip("Diretório demos/ não encontrado")
        return [
            f for f in os.listdir(self.DEMO_DIR) if f.endswith(".json")
        ]

    def test_demo_files_exist(self):
        demos = self._demo_files()
        assert len(demos) > 0, "Nenhum arquivo demo .json encontrado"

    @pytest.mark.parametrize("demo_file", [
        "2zones-2GLVs-2-Check-SA_en.json",
        "extended-ESP-pumpEfic.json",
        "simplifiedProduction.json",
    ])
    def test_load_demo_json(self, demo_file):
        """Arquivos demo devem ser carregáveis como Tramo."""
        filepath = os.path.join(self.DEMO_DIR, demo_file)
        if not os.path.isfile(filepath):
            pytest.skip(f"{demo_file} não encontrado")
        t = marlim3.Branch()
        t.from_json(filepath)
        assert t.system is not None
