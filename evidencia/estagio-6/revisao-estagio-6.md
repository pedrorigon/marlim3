# Revisão do Estágio 6 — o que foi feito, o que diverge, o que falta

Auditoria feita **medindo**, não relendo. Métodos: gcov sobre um build
instrumentado à parte, grafo de chamadas do módulo, comparação
assinatura a assinatura entre header e `.cpp`, e a varredura de caracterização
contra um golden anterior aos moves.

## 1. O que o estágio entregou

22 definições saíram de `SisProd.cpp` para `src/core/SisProdGasLift.cpp`
(1.320 linhas), atrás de um `GasLiftState` que dá ao domínio uma fronteira
única. `SisProd.cpp` caiu 1.162 linhas. Nove commits: header e harness
(`1808f55`), sete moves (T073–T079), os renames (T081r) e a evidência
(`177fdf6`).

## 2. Correção — o que está provado

| verificação | resultado |
|---|---|
| L2 bit-a-bit, corpus completo | **14/14** |
| baseline ainda válido nesta máquina | **VALID** (build pristino reproduz a captura) |
| regressão L3 | 5/5 |
| consumidores compilam | 4/4 |
| `tests/comparison/` intocado | sim |
| calibração do harness | 9/9 instrumentos |
| varredura gas-lift vs golden **pré-move** | **36/36 linhas idênticas** |
| header × `.cpp`: nomes de parâmetro | 22 assinaturas, **0 divergências** |
| declarações penduradas / definições órfãs | **nenhuma** |
| SC-004 (nenhuma função > 200 linhas) | maior = 147, mediana = 65 |
| FR-034 (módulo ≤ 8.000 linhas) | 1.320 |

O golden da varredura é de 01/09; o primeiro move é de 02/09. A comparação é
legítima, não circular.

## 3. A divergência que importa: o L2 só prova o que o L2 executa

Medido com gcov, acumulando os **três cenários distintos com gas-lift** do
corpus (os demais modelos são traduções pt-br dos mesmos casos):

    executadas pelo corpus ......... 12 de 22
    nunca executadas ............... 10 de 22

O caminho de *unloading* exige `configuracaoInicial/condicaoInicial == 3`
(`arq.descarga == 1`). **Nenhum dos 14 modelos define isso** — verificado nos
arquivos `.mr3`. Não é amostragem: é uma propriedade do corpus.

Das 10 não executadas, a varredura cobre 3 (`calibratedValveArea`,
`steadyInjectionPressureDrop`, `unloadingPressureCorrection`). Restam **7
funções que nenhuma camada de execução alcança**:

| linhas | função | situação |
|---:|---|---|
| 91 | `solveUnloading` | unloading vivo, sem modelo |
| 85 | `searchUnloadingInjectionPressure` | unloading vivo, sem modelo |
| 58 | `advanceInterface` | unloading vivo, sem modelo |
| 53 | `computeGasUnloadingHydrostatics` | unloading vivo, sem modelo |
| 53 | `computeUnloadingValvePressure` | unloading vivo, sem modelo |
| 15 | `advanceBufferedGasSubStep` | **morto** — sem chamador |
| 8 | `updateBufferedGasLine` | **morto** — só o anterior o chama |
| **363** | | **27% do módulo** |

Essas 363 linhas repousam **apenas** sobre identidade token a token. É uma
garantia real e foi a garantia usada em cada move — mas não é execução, e o
verde do L2 não fala por elas.

As duas funções mortas **já eram mortas no baseline**: `subtempoGasBuf` tem
uma única ocorrência em `0f3b64f`, a própria definição. A refatoração preservou
código morto fielmente; não o introduziu.

## 4. Erro encontrado e corrigido nesta revisão

`SisProdGasLift.h` afirmava que a varredura *"is the only verification that
reaches twelve of the twenty-two functions"*. Medido: a varredura alcança
**três**. O número 12 é o de funções que o corpus **executa**, não o que a
varredura cobre — as duas coisas foram trocadas.

Superestimar cobertura de verificação num comentário é o mesmo defeito que
`calibrate-all` marcando sonda morta como `dead=0`: a evidência afirma mais do
que a corrida sustenta. O comentário foi reescrito com os números medidos e diz
explicitamente o que a versão anterior alegava.

## 5. Pendências e melhorias — por ordem de valor

**(a) Sete funções sem execução (363 linhas).** A correção é estender
`gaslift-sweep.cpp` para dirigi-las sobre células sintéticas, exatamente como o
Estágio 5 fez pelo `tempDescarga` (a tabela térmica foi de 73 para 80 linhas).
O `GasLiftState` já foi desenhado para permitir isso — é literalmente a
justificativa registrada no header. **Este é o único item com risco de
corretude.**

**(b) 14 variáveis não usadas ainda no módulo** (12 `unused-variable`, 2
`set but not used`), a mesma classe que você mandou limpar no térmico. Três
delas são âncoras de depuração mortas, guardadas por `lixo5`:

    linha  530  if ((*state.globals).lixo5 > 1000) { int para; para = 0; }
    linha 1221  double fator = 0.; if (lixo5 < 1000.) fator = 0.;
    linha  697  double verifica = state.gasCells[0].pres;

O `thermal-move.py` já tem `strip_debug_anchors` name-agnostic; a normalização
existe.

**(c) T081r ficou incompleto quanto ao seu enunciado.** Os cinco critérios
(a)–(e) estão atendidos e foram verificados um a um. Mas o enunciado do FR-039
pede *"nomes em inglês significativos"*, e **33 dos 176 locais distintos
continuam em português ou opacos**: `massica`, `areamenor`, `auxpres`,
`auxpresmax`, `presDescini`, `presRev`, `posicM2`, `rholiq`, `rhomix`,
`vazmax`, `vazvalv`, `viscliq`, `viscmix`, `precorr`, `pmed`, `tmed`, `dxmed`,
`ugsmed`, entre outros. Os critérios de aceite são mais fracos que a tarefa;
marquei a tarefa como feita pelos critérios, e registro aqui a diferença.

**(d) SC-017, um conceito um nome.** `rhog` e `rhoG` coexistem no módulo
(linhas 185 e 785), em funções diferentes, ambos densidade de gás.

**(e) `updateBufferedGasLine` tem três ramos idênticos** — os três atribuem
`VGasRBuf = gasFreeTerms[3*i + 1]`. É código morto herdado; simplificar é
seguro mas sem ganho, e mudaria tokens sob controle. Deixado como está.

## 6. Estado após os ajustes (resolvido)

Todos os cinco itens foram tratados. Números medidos, não estimados:

| item | antes | depois |
|---|---|---|
| (a) funções sem execução | 6 (355 linhas) | **0** — união corpus+varredura = 22/22 |
| (b) avisos de variável não usada | 14 | **0** (módulo 1.320 → 1.299 linhas) |
| (c) locais em português/opacos | 33 de 176 | **0** (96 renames) |
| (d) `rhog` vs `rhoG` | dois nomes, um conceito | unificado; e um nome com **dois** conceitos corrigido |
| (e) `updateBufferedGasLine` | 3 ramos idênticos | colapsado, após torná-lo observável |

A tabela da varredura foi de 36 para 116 linhas, e a referência é regerada da
árvore **anterior** ao Estágio 6 — as seis funções que os moves não podiam
verificar por execução estão verificadas retroativamente. Calibração do gaslift:
17/17. L2: 14/14 após cada passo.

## 7. Veredito

O Estágio 6 **está correto no que foi verificado**, e o que foi verificado está
medido, não presumido. Não há divergência numérica: 14/14 bit-a-bit com o
baseline provado válido nesta máquina.

A ressalva honesta não é sobre o que o estágio fez, e sim sobre o alcance da
prova: **27% do módulo nunca foi executado por nenhuma camada**, e o header
alegava cobertura maior do que existe — agora corrigido. O item (a) deveria ser
resolvido antes do Estágio 7, porque o Estágio 7 move `marchaGasPerm1/2/3`,
que tocam o mesmo `celulaG[]`, e é melhor entrar nele com a rede já esticada.
