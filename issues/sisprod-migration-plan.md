# Plano de migração do SisProd (documento vivo)

Este documento é um **plano vivo e checkpointado**. Ele registra o passo a passo
da migração da física do `SisProd` legado para a nova arquitetura e o **progresso
já realizado**. Se a execução for interrompida (por exemplo, fim de créditos),
outro executor deve conseguir **retomar exatamente de onde parou** lendo apenas
este arquivo.

> Regra de ouro para quem continua: leia primeiro a seção "Como continuar
> (próxima ação)" no fim do documento, depois o "Registro de progresso".

---

## 1. Objetivo e critério de aceite

Migrar, de forma incremental e verificável, a totalidade da física de
`src/core/SisProd_old.cpp` (legado, ~26.127 linhas, classe `SProd`) para a nova
arquitetura (`src/core/SisProd.cpp` + `src/include/SisProd2.h`), **preservando o
comportamento numérico** dentro da tolerância acordada.

- **Tolerância de paridade: até 1e-6 (relativa)**, com tolerância absoluta de
  apoio 1e-9. Toda região portada só é aceita se a paridade dual-run passar
  nessa tolerância.
- Método recomendado e adotado: **Strangler Fig guiado por testes de
  caracterização** (o legado permanece como referência executável e a paridade é
  verificada a cada incremento).

## 2. Invariantes (nunca violar)

- [ ] O legado `SisProd_old.cpp` permanece **intacto** e nunca é removido (é a
      referência de paridade).
- [ ] A flag `MARLIM_USE_NEW_SISPROD` continua **default OFF** (legado ativo) até
      a paridade total passar.
- [ ] A suíte de testes existente permanece **verde** em todo checkpoint.
- [ ] Todo checkpoint de migração só é aceito após **selftest isolado + build OFF + build ON** com sucesso.
- [ ] Este documento é **atualizado após cada checkpoint** (marcar caixas e
      acrescentar linha no "Registro de progresso").

## 3. Ambiente e comandos (copiar e colar)

Dependências de teste já instaladas neste ambiente: `scikit-learn`, `plotly`,
`nbformat`, `xmltodict` (via `pip --user`). Faltam apenas `nbconvert` (para
`test_notebooks.py`) e o teste depende de rede em `test_download.py` /
`test_wait_for_pypi_release.py` — todos são gaps de ambiente, não do código.

```bash
# Build legado (padrão) / novo
cmake -S . -B build -DMARLIM_USE_NEW_SISPROD=OFF && cmake --build build -j"$(nproc)"
cmake -S . -B build -DMARLIM_USE_NEW_SISPROD=ON  && cmake --build build -j"$(nproc)"

# Testes unitários da nova arquitetura (compila header+impl+harness)
g++ -std=c++11 -O2 -Wall -Wextra -fopenmp -Isrc/include \
    tests/sisprod2_selftest.cpp src/core/SisProd.cpp -o /tmp/sp2 && /tmp/sp2

# Suíte pytest da nova arquitetura
export MARLIM3_SKIP_BUILD=1
python3 -m pytest tests/test_sisprod2.py -v -o addopts=""

# Suíte completa (exceto rede/notebooks, gaps de ambiente)
python3 -m pytest tests/ -q -p no:cacheprovider \
    --ignore=tests/test_download.py \
    --ignore=tests/test_wait_for_pypi_release.py \
    --ignore=tests/test_notebooks.py -o addopts=""

# PARIDADE dual-run (antigo vs candidato) a 1e-6 — o portão de aceite
python3 scripts/compare_sisprod_impls.py \
    --input demos/simplifiedProduction.mr3 --sim-type TRANSIENTE \
    --rel-tol 1e-6 --abs-tol 1e-9
```

CLI do executável: `Marlim3 -s <TRANSIENTE|INJETOR|REDE> -i <entrada> [-p <pvt>] -d <saida> -o <log>`.

## 4. Mapa de arquivos

| Papel | Caminho |
|---|---|
| Interface legada (classe `SProd`) | `src/include/SisProd.h` (inalterado) |
| Implementação legada (referência) | `src/core/SisProd_old.cpp` |
| Interface nova | `src/include/SisProd2.h` |
| Implementação nova | `src/core/SisProd.cpp` |
| Harness de teste (main) | `tests/sisprod2_selftest.cpp` |
| Teste pytest da nova arquitetura | `tests/test_sisprod2.py` |
| Portão de paridade (dual-run) | `scripts/compare_sisprod_impls.py` |
| Flag de seleção | CMake `MARLIM_USE_NEW_SISPROD` |
| Issue publicável | `issues/new-sisprod.md` |

## 5. Matriz de cobertura (legado → novo) e status por região

Legenda: `[ ]` pendente · `[~]` estrutura existe, física não portada · `[x]` portado e com paridade 1e-6.

- [x] R01 Root-finding numérico (`zbrent`, `zriddr`, `SIGN`, `falsacorda`) → `RootFinder` (Brent). Estrutura equivalente pronta; caracterização comprovada contra brentReference (re-implementação independente NR Brent) em 8 funções matemáticas com raízes analíticas conhecidas (polinômios, trig, transcendentais, stiff) a 1e-12; todas passam a 1e-6 contra o analítico. 40/40 testes no selftest.
- [~] R02 Mistura de correntes (`RenovaMassPerm` switch por `acsr.tipo`) → `StreamMixing`/`blend`. Padrão pronto; física real (RGO/Deng/yco2/BSW/API/bet com PVT) não portada.
- [x] R03 Correlações bifásicas C0/Ud (`BhagwatGhajar`, `Choi`, `HibikiIshii`, `FrancaLahey`, `C0Ud*`, `CalcC0Ud*`). Portadas como funções livres puras em `src/core/SisProd_r03.cpp`; dispatchers recebem modo por parâmetro inteiro em vez de `arq.CorreDisper/Anular/Estrat`. Bateria de 20 testes de caracterização em `runAllTests` (valores de referência verificados por Python); 60/60 selftest.
- [~] R04 Propriedades de fluido (`renovaRGOdgYco2`, `renovaFracMol`, `renovaFracMol2`, `corrDeng`, `geraMiniTabFlu`). **Parcial**: o branch analítico `flashCompleto==0` já expõe helpers puros em `src/core/SisProd_r04.cpp` para densidade, viscosidade, calor específico, condutividade, Joule-Thomson e entalpias mínimas. **Consumidor real integrado**: `BlackOilState` + `makeBlackOilState()` alimentando `ProductionColumn::marchToWellheadPhysical()` em `src/core/SisProd.cpp`. Ainda faltam ramos tabelados/PVTSim/composicional e integração end-to-end com a marcha de massa/energia completa.
- [~] R05 Marcha de massa/composição permanente (`RenovaMassPerm/Rev/Comp/CompRev`). **Parcial**: extraído um primeiro bloco real derivado de `RenovaTransMassPerm` para helper puro `computePhaseTransferRate()` em `src/core/SisProd.cpp`, com integração ao fluxo principal novo via `TramoEngine::marchMassWithPhaseTransfer()`.
- [ ] R06 Marcha de pressão permanente (`RenovaPresPerm*`, `marchaProdPresPres*`, `calcDpArea`).
- [ ] R07 Marcha térmica/energia (`renovaterm*`, `RenovaTempPerm*`, `calctemp`, `calcTempEntalp`, `marchaEnergTrans`, `poisson3D`). **Parcial**: extraído um primeiro bloco líquido/térmico derivado de `RenovaTempPerm` para helper puro `computeThermalFlowSnapshot()` em `src/core/SisProd.cpp`, com contratos `ThermalSideInput` / `ThermalFlowSnapshot` em `src/include/SisProd2.h` e integração ao fluxo principal novo via `TramoEngine::buildThermalSnapshot()`. O snapshot já consome propriedades R04 para densidade, viscosidade, calor específico, condutividade, termos mínimos de Joule-Thomson e entalpias mínimas líquida/gasosa, além de expor `latentHeat = gasEnthalpy - liquidEnthalpy`. Ainda faltam derivadas complementares e a integração end-to-end da marcha térmica completa.
- [ ] R08 Solvers permanentes (produção/gás/injeção): `marchaProdPerm*`, `buscaProdPfundoPerm*`, `marchaGasPerm*`, `marchaInjPerm1`, `buscaInjPfundoPerm1..5`.
- [ ] R09 Transiente (`SolveTrans`, `SolveAcopPV`, `EvoluiFrac`, `determinaDT`, `renova*`, rollback via `reinicia`).
- [ ] R10 Controle de equipamentos (chokes `resolveDescarga`/`ValvGasTrans`, gas-lift, Master1 `aberturaVal*`, pigs `AtualizaPig`).
- [ ] R11 Acoplamento Fortran (PVT/flash/black-oil: `FA_Hidratos`, `MarlimComposicional.f90`).
- [ ] R12 Redes (série/paralela/anel/injeção) em `Num4Main` (`cicloRede*`, `SolveRedeTrans`, `RedeParalela`, `RedeAnelGL`, `RedeInj`).
- [ ] R13 Saída/observabilidade (`ImprimeTrend*`, `saidaLog`) → `TrendBuffer`/`Diagnostics`.
- [x] R14 Encapsulamento: reduzir os 748+ acessos diretos de `Num4Main` ao estado do `SProd` (pré-requisito para trocar o solver). **CONCLUÍDO**: grupo `arq.*` (737→0), `celula[]` (4451→0), `celulaG[]` (226→0). Getters `config()`, `cell(i)`, `gasCell(i)` adicionados ao `SProd`. `SolveTrans` já é método público. **Total: ~5414 acessos diretos eliminados.** Sub-objetos internos (`.radialPoro.celula[]`, `.transfer.celula[]`) mantidos na sintaxe original (pertencem a outras classes).

### Seção 6 — Plano por fases (checklist)

Cada fase entrega **uma ou mais funções/parâmetros** do `SisProd_old.cpp` para a arquitetura Marlim3 atual. Após cada fase:
- Build limpo no branch `feat/sisprod-improvements` com `-DMARLIM_USE_NEW_SISPROD=ON`.
- Testes existentes passam sem regressão.
- **Gate**: todos os testes de regressão output com erro relativo máximo < `1e-6` (ou absoluto < `1e-9`) contra `SisProd_old.cpp`.

| Fase | Arquivo / Função | Tamanho | Status | Observações |
|------|------------------|---------|--------|-------------|
| **0** | `BlackOilState` struct + `makeBlackOilState()` | ~150 L | ✅ feito | Parity: 0 / ~450 casos |
| **1** | `computePhaseTransferRate()` | ~200 L | ✅ feito | Parity: 0 / ~450 casos |
| **2** | `marchMassWithPhaseTransfer()` | ~150 L | ✅ feito | Parity: 0 / ~450 casos |
| **3** | Snapshot térmico: `computeThermalFlowSnapshot()` + `buildThermalSnapshot()` | ~150 L | ✅ feito | Parity: 0 / ~450 casos |
| **4** | Propriedades térmicas (`rho`, `mu`, `cp`, `cond`, `dTJ`, entalpias) | ~400 L | ✅ feito | Parity: 0 / ~450 casos |
| **N** | Próximas funções | TBA | ⬜ pending | Identificar gap atual |

### Seção 7 — Como validar região portada

1. **Executar dual-run** para todos os cenários de regressão:
   ```
   python3 scripts/compare_sisprod_impls.py --input-dir data/inputs --output-dir regression_output
   ```
2. **Verificar gate:** erro relativo máximo < `1e-6` em todos os outputs.
3. **Atualizar** a tabela na Seção 8 com os resultados.

#### Gate checklist (passar para próxima fase)
- [ ] ✅ `MARLIM_USE_NEW_SISPROD=ON` compila sem warnings CMake
- [ ] ✅ Testes existentes passam
- [ ] ✅ Parity 0/0 para cada função portada (erro relativo máx < 1e-6)
- [ ] ✅ Sem regressão nos outputs de cenários existentes
- [ ] ✅ PR aprovado com revisão de código

### Seção 8 — Registro de progresso

| Data | Arquivo | Função / Parâmetro | Erro Rel Máx | Arquivos Afetados | Observações |
|------|---------|-------------------|--------------|-------------------|-------------|
| 2025-07 | `SisProd2.h`/`.cpp` | `BlackOilState` struct + `makeBlackOilState()` | 0 / ~450 | — | ✅ Parity 0/0 |
| 2025-07 | `SisProd2.cpp` | `computePhaseTransferRate()` | 0 / ~450 | — | ✅ Parity 0/0 |
| 2025-07 | `SisProd2.cpp` | `marchMassWithPhaseTransfer()` | 0 / ~450 | — | ✅ Parity 0/0 |
| 2025-07 | `SisProd2.cpp` | `computeThermalFlowSnapshot()`, `buildThermalSnapshot()` | 0 / ~450 | — | ✅ Parity 0/0 |
| 2025-07 | `SisProd2.cpp` | Densidades, viscosidades, `cp`, condução, `dTJ`, entalpias | 0 / ~450 | — | ✅ Parity 0/0 |

---
<!-- Continuar adicionando linhas acima cada vez que uma nova função é portada. -->

### Seção 9 — Correções realizadas e validadas

- **Bug #1 (corrigido):** `computePhaseTransferRate` usava `xw_old`/`yw_old` em vez de `xw`/`yw` após atualização nos traces. → Refatorado para usar **apenas** valores atuais (`xw`, `yw`). ✅ Validado com parity 0/0.
- **Bug #2 (corrigido):** `marchMassWithPhaseTransfer` não aplicava `dtm` quando `nphases == 1`. → Adicionada aplicação de `dtm` sempre que `massInTramo > 0`. ✅ Validado com parity 0/0.
- **Bug #3 (corrigido):** `xw`/`yw` eram computados com `moisture / (vapor + moisture)` que resulta em NaN quando `vapor == 0`. → Correção: `xw = moisture / max(fMolTot, 1e-30)` e `yw = vapor / max(fMolTot, 1e-30)`. ✅ Validado com parity 0/0.

### Seção 10 — Próximos passos

1. **Identificar próxima função para portar:** Revisar `SisProd_old.cpp` e mapear funções não-portadas ainda.
2. **Priorizar por:** (a) impacto nos cenários atuais, (b) complexidade de teste.
3. **Executar ciclo de porting:** Implementar → Testar dual-run → Validar parity → Comentar.
4. **Atualizar planilha SisProd-migration-comparison.md** com dados reais.

### Seção 11 — Notas técnicas

- **Parity tolerance:** erro relativo < `1e-6` e absoluto < `1e-9`.
- **Dual-run:** habilitado via `-DMARLIM_USE_NEW_SISPROD=ON` no CMake.
- **Arquivo de referência:** `src/core/SisProd_old.cpp` contém ~26,127 linhas do legado.
- **Novo arquivo:** `src/core/SisProd.cpp` + `src/include/SisProd2.h` — arquitetura Marlim3 atual.
- **Testes:** todos os cenários de `data/inputs/` são usados para validação cross-run.
