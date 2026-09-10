# T080 — propriedade de `celulaG[]` ao fim do Estágio 6

## Veredito

**O critério de saída do Estágio 6 não é alcançável pelas tarefas do Estágio 6 —
e o plano nunca o torna verdadeiro.** Isto não é um trabalho pendente; é uma
contradição entre o critério e o cronograma, e está registrada aqui em vez de
ser contornada.

## Medição

    src/core/SisProd.cpp            204 acessos
    src/core/Num4Main.cpp           278
    src/core/LerAP.cpp               61
    src/core/SisProdGasLift.cpp       2   <- nenhum é acesso

Os dois casamentos dentro do módulo são `#include "celulaGas.h"` (substring) e
uma menção em comentário. O módulo acessa o domínio exclusivamente por
`state.gasCells` — 714 usos. Do lado do módulo, o encapsulamento está completo.

## Onde estão os 204 acessos de SisProd.cpp

| acessos | função | tarefa que a reivindica | estágio |
|---:|---|---|---|
| 43 | `marchaGasPerm2` | T094 | 7 |
| 41 | `marchaGasPerm1` | T094 | 7 |
| 39 | `marchaGasPerm3` | T094 | 7 |
| 21 | `montasistema` | T100 | 9 |
| 13 | `SolveTrans` | T127 | 8+ |
| 9 | `hidroLinServ` | T096 | 7 |
| 6 | `marchaProdPerm2` | T090 | 7 |
| 6 | `marchaProdPerm1` | T089 | 7 |
| 4 | `SProd` (construtor) | T101 | 9 |
| 4 | `determinaDT` | T122 | 8+ |
| 4 | `calcDTPseudoTrans` | T097 | 7 |
| 4 | `buscaProdPfundoPerm3` | T091 | 7 |
| 4 | `buscaGasPresPerm2` | T094 | 7 |
| 2 | `renew` | — | 9 |
| 2 | `copiaSemJson` | T101 | 9 |
| 1 | `ReiniEvolFrac` | T125 | 8+ |
| 1 | `buscaGasPresPerm3` | T094 | 7 |
| **204** | 17 funções | | |

**Nenhum destes acessos pertence a uma tarefa do Estágio 6.** Todas as 22
definições que o Estágio 6 reivindicou saíram (T073–T079), e as funções que
restam tocando `celulaG[]` são propriedade declarada de estágios posteriores.

## Por que o critério nunca se torna verdadeiro

O critério de saída diz: *"`celulaG[]` tocado exclusivamente por este módulo"*.

Mas o T094 (linha 373 de tasks.md) determina, textualmente, que
`marchaGasPerm1`, `marchaGasPerm2` e `marchaGasPerm3` — 123 dos 204 acessos,
60% do total — vão para **`src/core/SisProdSteadyState.cpp`**, e que
`buscaGasPresPerm2`/`buscaGasPresPerm3` vão para
**`src/core/SisProdSteadyStateSearch.cpp`**.

Ou seja: após o Estágio 7, `celulaG[]` será tocado pelo módulo de gas-lift **e**
por dois módulos de regime permanente. O critério de exclusividade não é apenas
prematuro no Estágio 6 — ele é incompatível com o destino que o próprio plano dá
à marcha de gás.

Isso é defensável em termos de projeto: a marcha de gás em regime permanente é
uma marcha, e o corte do Estágio 7 é por *regime*, não por *domínio de dados*.
Dois cortes legítimos cruzam-se aqui. O que não se sustenta é a redação do
critério.

## O que este estágio de fato provou

- As 22 definições movidas não deixaram cálculo de linha de gás para trás.
- O módulo não alcança `celulaG[]` diretamente: só por `GasLiftState`.
- A fronteira é atravessada por um único adaptador (`gasLiftStateOf`), o que
  torna a propriedade futura verificável num só ponto.

## Recomendação

Reescrever o critério de saída do Estágio 6 para o que é verificável agora —
*"nenhum cálculo de linha de gás reivindicado por este estágio permanece em
`SisProd.cpp`; o módulo não referencia `celulaG[]` diretamente"* — e mover a
afirmação de exclusividade para depois do T094, onde ela deve ser reformulada
como propriedade **compartilhada e explícita** entre gas-lift e regime
permanente. `Num4Main.cpp` (278) e `LerAP.cpp` (61) não são endereçados por
nenhuma tarefa do plano e precisam de decisão à parte.
