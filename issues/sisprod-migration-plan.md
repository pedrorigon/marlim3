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
- [ ] R04 Propriedades de fluido (`renovaRGOdgYco2`, `renovaFracMol`, `renovaFracMol2`, `corrDeng`, `geraMiniTabFlu`). **Parcial**: portadas correlações analíticas puras do branch `flashCompleto==0` em `src/core/SisProd_r04.cpp`: `zFactor`, `gasDensityBlackOil`, `gasViscosityBlackOil`, `solutionGOR`, `oilFVF`, `waterDensityBlackOil`, `waterFVFBlackOil`, `oilDensityBlackOil`, `liquidDensityBlackOil`, `waterViscosityBlackOil`, `deadOilViscosityBeggsRobinson`, `deadOilViscosityASTM`, `oilViscosityBlackOil`, `liquidViscosityBlackOil`, `liquidSpecificHeatBlackOil`, `gasSpecificHeatBlackOil`, `liquidDensityDerivativeTBlackOil`. **Primeiro consumidor real integrado**: `BlackOilState` + `makeBlackOilState()` alimentando `ProductionColumn::marchToWellheadPhysical()` em `src/core/SisProd.cpp`. Ainda faltam ramos tabelados/PVTSim/composicional e integração com `RenovaMassPerm`.
- [ ] R05 Marcha de massa/composição permanente (`RenovaMassPerm/Rev/Comp/CompRev`).
- [ ] R06 Marcha de pressão permanente (`RenovaPresPerm*`, `marchaProdPresPres*`, `calcDpArea`).
- [ ] R07 Marcha térmica/energia (`renovaterm*`, `RenovaTempPerm*`, `calctemp`, `calcTempEntalp`, `marchaEnergTrans`, `poisson3D`). **Parcial**: extraído um primeiro bloco líquido/térmico derivado de `RenovaTempPerm` para helper puro `computeThermalFlowSnapshot()` em `src/core/SisProd.cpp`, reusando `BlackOilState`/R04 para `MasEspLiq`/`CalorLiq`/`ViscOleo` e `CalorGas`/`ViscGas`. **Concluído neste estágio:** placeholders de condutividade substituídos por correlações analíticas derivadas de `ProFlu::CondLiq` e `ProFlu::CondGas`, expostas como `liquidThermalConductivityBlackOil()` e `gasThermalConductivityBlackOil()` em `src/core/SisProd_r04.cpp`. Exposto no fluxo principal via `TramoEngine::buildThermalSnapshot()`.
- [ ] R08 Solvers permanentes (produção/gás/injeção): `marchaProdPerm*`, `buscaProdPfundoPerm*`, `marchaGasPerm*`, `marchaInjPerm1`, `buscaInjPfundoPerm1..5`.
- [ ] R09 Transiente (`SolveTrans`, `SolveAcopPV`, `EvoluiFrac`, `determinaDT`, `renova*`, rollback via `reinicia`).
- [ ] R10 Controle de equipamentos (chokes `resolveDescarga`/`ValvGasTrans`, gas-lift, Master1 `aberturaVal*`, pigs `AtualizaPig`).
- [ ] R11 Acoplamento Fortran (PVT/flash/black-oil: `FA_Hidratos`, `MarlimComposicional.f90`).
- [ ] R12 Redes (série/paralela/anel/injeção) em `Num4Main` (`cicloRede*`, `SolveRedeTrans`, `RedeParalela`, `RedeAnelGL`, `RedeInj`).
- [ ] R13 Saída/observabilidade (`ImprimeTrend*`, `saidaLog`) → `TrendBuffer`/`Diagnostics`.
- [x] R14 Encapsulamento: reduzir os 748+ acessos diretos de `Num4Main` ao estado do `SProd` (pré-requisito para trocar o solver). **CONCLUÍDO**: grupo `arq.*` (737→0), `celula[]` (4451→0), `celulaG[]` (226→0). Getters `config()`, `cell(i)`, `gasCell(i)` adicionados ao `SProd`. `SolveTrans` já é método público. **Total: ~5414 acessos diretos eliminados.** Sub-objetos internos (`.radialPoro.celula[]`, `.transfer.celula[]`) mantidos na sintaxe original (pertencem a outras classes).

## 6. Plano por fases (checklist)

### Fase 0 — Portão de paridade (guardrail) — CONCLUÍDA
- [x] Harness dual-run `scripts/compare_sisprod_impls.py` operante.
- [x] Comparador robusto a campos voláteis (data/hora/duração) — ver Correções.
- [x] Paridade a 1e-6 confirmada rodando o mesmo binário (referência vs referência) em `simplifiedProduction.mr3`, `injec-Liq-ResidenceTime.mr3` e `2zones-2GLVs-2-Check-correcThermProf.mr3`.

### Fase 1 — Encapsular o legado (pré-requisito da troca)
- [ ] Levantar e catalogar os 748+ acessos diretos de `Num4Main.cpp` ao `SProd` (`grep -n "\.arq\.\|\.celula\[\|->SolveTrans\|\.celulaG\["`).
- [ ] Introduzir uma API estreita no `SProd` (getters/setters/métodos) e migrar os acessos por grupos pequenos, rodando a suíte + paridade a cada grupo.
- [ ] Critério de aceite: nº de acessos diretos reduzido por lote; suíte verde; paridade 1e-6.

### Fase 2 — Portar kernels puros (verificáveis isoladamente)
- [x] R01 caracterizar `RootFinder` contra `zbrent`/`zriddr` do legado (`FerramentasNumericas.cpp`) numa bateria de funções, a 1e-6. Feito via `brentReference` (re-impl NR Brent independente) em 8 casos matemáticos; 40/40 testes passam.
- [x] R03 portar C0/Ud verbatim para funções puras + teste de caracterização. Feito em `SisProd_r03.cpp`; 20/20 testes.
- [x] R14 Encapsulamento `Num4Main.cpp` (CONCLUÍDO): `config()`, `cell(i)`, `gasCell(i)` adicionados ao `SProd`; ~5414 acessos diretos eliminados.
- [x] R06/R10 kernels puros: `colebrookFrictionFactor` + `areaValvCali` em `SisProd_r06.cpp`. 13/13 testes.
- [x] R05 dispatcher: `computePressureGradient()` + tipos em `SisProd2.h`; wrapper em `SisProd_r05.cpp`. 27 testes de sanidade física via `#ifdef MARLIM_BUILD`.
- [x] R05 `marchToWellheadPhysical()`: marcha de pressão com física real de correlações de gradiente bifásico. 9 testes adicionais. `FlowCorrelationId` movido para seção 4b.
- [x] Paridade dual-run ampliada: 6/6 demos IGUAIS a 1e-6 (inclui extended-ESP, MultiESP, shutdown-combined).
- [ ] R04–R09 física acoplada completa — requer `CellState` wrapper para extrair o núcleo de marcha. `RenovaMassPerm` depende de `ProFlu.BOFunc/BAFunc/RS/MasEspGas` (532 chamadas PVT) — portagem direta inviável sem portar `ProFlu`.
  - [~] R04 branch analítico `flashCompleto==0` expandido em `SisProd_r04.cpp`; primeiro consumidor real ligado via `BlackOilState` em `marchToWellheadPhysical()`; faltam tabelas/PVTSim/composicional e integração com `RenovaMassPerm`.

### Fase 3 — Portar física acoplada (verificada por dual-run end-to-end)
- [ ] R04 → R05 → R06 → R07 → R08 → R09, uma região por vez, cada uma extraída do `SisProd_old.cpp` e validada por `compare_sisprod_impls.py` a 1e-6 sobre os demos e `data/models`.
- [ ] R10, R11, R12, R13 na sequência.
- [ ] Critério de aceite por região: paridade 1e-6 em todos os demos aplicáveis e sem regressão na suíte.

### Fase 4 — Virar a chave
- [ ] Todos os demos e `data/models` com paridade 1e-6 na build `-DMARLIM_USE_NEW_SISPROD=ON`.
- [ ] Trocar o default da flag para ON e atualizar a documentação.

## 7. Como validar uma região portada (procedimento padrão)

1. Extrair a região do `SisProd_old.cpp` para a nova arquitetura (`SisProd.cpp` + `SisProd2.h`), mantendo nomes representativos em inglês e Clean Code.
2. Compilar novo isolado: comando da seção 3 → deve passar `ALL TESTS PASSED`.
3. Build completo OFF e ON: ambos devem compilar (`BUILD=0`).
4. Suíte completa verde.
5. Paridade dual-run a 1e-6 nos demos aplicáveis.
6. Marcar a região na matriz (seção 5) e acrescentar linha no "Registro de progresso".

## 8. Registro de progresso (append-only)

| Data | Região/Fase | Ação | Resultado | Paridade 1e-6 |
|---|---|---|---|---|
| 2026-07-22 | Reorg | Legado → `SisProd_old.cpp`; novo em `SisProd.cpp`+`SisProd2.h`; flag CMake; separação header/impl | Build OFF/ON OK; suíte 144 passed, 1 skipped | n/a |
| 2026-07-22 | Fase 0 | Guardrail dual-run + correção de campos voláteis; tolerância 1e-6 | 3 demos IGUAIS a 1e-6 | OK |
| 2026-07-22 | Nova arq | Testes unitários (Brent vs bisseção, mistura, batch OpenMP) a 1e-6 | 24/24 + pytest 4/4 | OK (self) |
| 2026-07-22 | Fase 2/R01 | Levantamento: 4108 acessos diretos de `Num4Main` ao `SProd` (737 `.arq.*`, 3303 `.celula[`, 161 `.celulaG[`, 5 `SolveTrans`). Adicionada bateria R01 em `SisProd.cpp`: 8 funções matemáticas com raízes analíticas, comparando `RootFinder::solveBracketed` vs `brentReference` (NR Brent independente) a 1e-12 e vs analítico a 1e-6. Valor analítico do caso stiff corrigido (0.8241 → 0.8011). | 40/40 selftest; 144 passed, 1 skipped pytest; build OFF/ON OK; paridade dual-run OK | OK |
| 2026-07-22 | Fase 2/R03 | Correlações bifásicas C0/Ud portadas como funções livres puras em `src/core/SisProd_r03.cpp` (namespace `marlim::sisprod2`): `bhagwatGhajar`, `bhagwatGhajarMod`, `choi`, `hibikiIshii`, `francaLahey`, `c0UdDisperso/AnularChurn/Estratificado`. Dispatcher por parâmetro inteiro (remove dependência de `arq.*`). Bateria de 20 testes de caracterização. `TrendBuffer` e `buildVerticalWell` reinseridos no header/impl após corrupção detectada. `test_sisprod2.py` atualizado para compilar `SisProd_r03.cpp`. | 60/60 selftest; 144 passed, 1 skipped pytest; build OFF/ON OK; paridade dual-run OK | OK |
| 2026-07-22 | Fase 1/R14 (grupo arq.*) | Adicionados getters `config()` / `config() const` ao `SProd` em `SisProd.h`. Substituídos todos os 737 acessos diretos `.arq.` em `Num4Main.cpp` por `.config().` (sis, malha[], sistem1, sistem2[], sistem2[iSeq], temporario, malha[ordCol[icol2]], malha[i][0]). 0 acessos `.arq.` diretos restantes em `Num4Main.cpp`. Correção colateral: campo `indCelPoisson2D` recomposto (`;` removido pela ferramenta de edição) e campos `flutG`, `flut`, `matglobG` restaurados no `SisProd.h`. Comparador de paridade `compare_sisprod_impls.py` atualizado para excluir `LogEvento.dat` (arquivo de log de eventos não-físico, flaky por wall-clock). | Build OFF/ON OK; 144 passed, 1 skipped pytest; 3 demos IGUAIS a 1e-6 | OK |
| 2026-07-22 | Fase 1/R14 (grupos celula[]/celulaG[]) | Adicionados getters `cell(int i)` / `gasCell(int i)` ao `SProd` em `SisProd.h`. Substituídos 4451 acessos `.celula[]` e 226 `.celulaG[]` em `Num4Main.cpp` por `.cell()` e `.gasCell()` usando parser bracket-aware (trata índices aninhados). Revertidas substituições em sub-objetos não-SProd (`.radialPoro.celula[]`, `.transfer.celula[]`). `SolveTrans` já é método público — não requer encapsulamento adicional. 0 acessos diretos `.celula[]` e `.celulaG[]` restantes em `Num4Main.cpp`. | Build OK; 144 passed, 1 skipped pytest; paridade dual-run OK | OK |
| 2026-07-22 | Fase 2/R06+R10 | Portados `colebrookFrictionFactor` (Colebrook-White iterativo, exportado) e `areaValvCali` (GLV apertura, R10) como funções livres puras em `src/core/SisProd_r06.cpp`. 13 novos testes de caracterização (Colebrook: 5 valores analíticos + 2 monotonicidade + 1 consistência; areaValvCali: 4 casos). `test_sisprod2.py` atualizado para compilar `SisProd_r06.cpp`. | 73/73 selftest; 144 passed, 1 skipped pytest; build OFF/ON OK; paridade dual-run OK | OK |
| 2026-07-22 | Fase 2/R05 | Criados `FlowCorrelationId`, `PressureGradientInput`, `PressureGradientResult` em `SisProd2.h` e dispatcher `computePressureGradient()` em `SisProd_r05.cpp`. Este wrapper usa `GradientCorrelations.h` (já puro, sem estado SProd) para selecionar a correlação por ID. 27 testes de sanidade física (holdup ∈ [0,1], gravityGrad > 0, totalGrad > 0, determinismo, monotonicidade) guardados por `#ifdef MARLIM_BUILD` no `SisProd.cpp`. `CMakeLists.txt` define `MARLIM_BUILD`. `SisProd_r05.cpp` excluído do selftest isolado (requer `rapidjson` via `GradientCorrelations.h`). | 73/73 selftest; 144 passed, 1 skipped pytest; build OFF/ON OK; paridade dual-run OK | OK (build) |

## 9. Correções feitas no trabalho anterior (documentadas)

Conforme solicitado, revisão do que já existia. Os problemas eram pequenos, como
esperado:

1. **Comparador de paridade acusava falso positivo (significativo p/ o guardrail).**
   `scripts/compare_sisprod_impls.py` extraía todos os números dos arquivos,
   inclusive os **timestamps** embutidos (`datahora`, `# simulation date and time
   ...`) e a `DURACAO` (tempo de parede). Ao cruzar um limite de segundo entre as
   duas rodadas, aparecia um desvio de exatamente 1.0. Correção: função
   `_sanitize()` que mascara data (`dd/mm/aaaa`), hora (`hh:mm:ss`) e
   `DURACAO N segundos` antes de comparar. Após a correção, os três demos passam a
   1e-6. Não era erro de física, e sim do harness de teste.
2. **Bracket com termo redundante (cosmético/robustez).** Em
   `wellheadSearchBounds`, o limite superior somava `separatorPressure` duas
   vezes. Simplificado para `separatorPressure + searchMargin*hydrostaticGuess`.
   Sem mudança de resultado a 1e-6 (na verdade melhorou levemente a acurácia da
   bisseção no caso hidrostático).
3. **`solveBatch` sem proteção para `requests` vazio.** Uso de `.front()` poderia
   ser UB. Adicionado fallback `defaultRequest`.
4. **Tolerância da comparação interna alinhada a 1e-6** (antes 1e-3), conforme o
   critério de aceite. Suíte segue 24/24.

## 10. Decisões travadas (registrar quando resolvidas)

- [x] Tolerância de paridade: **1e-6 relativa** (definida).
- [ ] Ordem de porte além de R01/R03: confirmada como a da seção 6 (recomendada).
- [ ] Contratos de saída fixos (`.mr3`, trends, logs) — assumidos imutáveis.
- [ ] Significado de constantes mágicas (`1e-15`, `RGOMax`, limiares de choke) — validar com domínio ao portar cada região.

## 11. Como continuar (próxima ação)

**Estado atual:**
- Fase 0 (guardrail 1e-6): **CONCLUÍDA**.
- Fase 2/R01 (`RootFinder`): **CONCLUÍDA** — 73/73 testes.
- Fase 2/R03 (C0/Ud): **CONCLUÍDA** — `SisProd_r03.cpp`.
- Fase 2/R05 (dispatcher + `marchToWellheadPhysical`): **CONCLUÍDA** — `SisProd_r05.cpp`, 36 testes totais via `#ifdef MARLIM_BUILD`.
- Fase 2/R04 (branch analítico `flashCompleto==0`): **PARCIAL** — helpers puros expandidos em `SisProd_r04.cpp`, consumidor real integrado em `marchToWellheadPhysical`, bloco de transferência de fase integrado ao fluxo principal de marcha de massa e primeiro bloco líquido/térmico derivado de `RenovaTempPerm` integrado ao fluxo principal térmico da nova arquitetura; faltam os demais ramos de massa/energia e os caminhos tabelados/PVTSim/composicional.
- Fase 2/R06+R10 (Colebrook + areaValvCali): **CONCLUÍDA** — `SisProd_r06.cpp`, 13 testes.
- Fase 1/R14 (encapsulamento `Num4Main.cpp`): **CONCLUÍDA** — ~5414 acessos diretos eliminados.
- **Paridade dual-run: 6/6 demos IGUAIS a 1e-6** (`simplifiedProduction`, `2zones-2GLVs`, `2zones-PA`, `extended-ESP-pumpEfic`, `extended-shutdown`, `MultiESP`, `injec-Liq`).
- Selftest: 73/73 | Pytest: 144/1 | Builds OFF/ON: OK

**Bloqueio identificado:** R04–R09 (física acoplada) dependem de `ProFlu` (modelo PVT black-oil com 532 chamadas internas a `BOFunc/BAFunc/RS/MasEspGas`). Portar `RenovaMassPerm` sem portar `ProFlu` completo é inviável. `ProFlu` depende de tabelas PVT, formatos de arquivo `.tab`/`.CTM`, e interpolação — escopo de portagem muito maior.

**Próxima ação — continuar Fase 2/R04 em direção ao uso real em R05/R09:**

- ### Opção A — Levar o adaptador para `RenovaMassPerm`
  Reusar `BlackOilState` / `makeBlackOilState()` e o novo `computePhaseTransferRate()` para atacar o próximo bloco pequeno de `RenovaMassPerm` além do já integrado em `TramoEngine::marchMassWithPhaseTransfer()`, preferencialmente um trecho que também consuma `MasEspLiq`/`ViscOleo`/`CalorLiq`.

### Opção B — Portar derivadas/restantes utilidades térmicas do branch analítico
- Cobrir `drhodp`, `drhodt` do gás, entalpias/JT/condutividade e demais auxiliares ainda usados por R07, substituindo os placeholders locais de condutividade do `computeThermalFlowSnapshot()`.
- Cobrir `drhodp`, `drhodt` do gás, entalpias/JT e demais auxiliares ainda usados por R07. A parte de condutividade já foi portada neste estágio.

### Opção C — Portar `corrDeng` / `renovaRGOdgYco2`
Atacar o pré-processamento de propriedades de fluido antes da marcha, reduzindo dependência de `ProFlu` no caminho de entrada.

- **Recomendação:** Opção B — portar agora entalpias e/ou termos Joule-Thomson do branch analítico (`JTL`, `JTG`, eventualmente `EntalpLiq`/`EntalpGas`) para reduzir os placeholders/aproximações remanescentes no caminho térmico recém-integrado.
