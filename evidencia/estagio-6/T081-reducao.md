# T081 — redução de `src/core/SisProd.cpp`

## Medição

    antes do Estágio 6 (24190e4)   20.483 linhas
    depois do T079                 19.321 linhas
    redução medida                  1.162
    alvo declarado no T081          1.309
    diferença                         147

## A conta fecha — o alvo é que está obsoleto

    corpos efetivamente removidos ...... 1.287
    menos os stubs de delegação ........   -67
    menos o adaptador de estado ........   -58
    = redução esperada ................. 1.162
    redução medida ..................... 1.162   confere

Os três termos, medidos e não estimados:

**Corpos removidos: 1.287.** O alvo de 1.309 vem da soma das faixas de linha
listadas no próprio T081, medidas em `SisProd.cpp` **antes do Estágio 1**. Os
Estágios 1–5 alteraram esses mesmos corpos — extraíram deriva, solvers, C0/Ud,
tendências e o cálculo térmico deles. Quando o Estágio 6 chegou, as funções já
eram 22 linhas mais curtas. As faixas do T081 nunca foram reaferidas.

**Stubs: 67 linhas, 22 funções.** Não podem ser removidos. O FR-032 exige
assinatura pública idêntica e o FR-038 proíbe renomear método público de
`SProd`; `Num4Main.cpp` chama esses métodos pelo nome. Cada stub tem 3 linhas
(assinatura, delegação, `}`), exceto `areaValvCali` com 4, por ter parâmetros
demais para uma linha.

**Adaptador: 58 linhas.** É o preço do contrato do módulo:

    gasLiftStateOf (corpo, monta 36 campos) ....... 42
    GasLiftTemperatureUpdater (3 métodos) ......... 14
    declaração antecipada de gasLiftStateOf ....... 1
    #include "SisProdGasLift.h" ................... 1

O `GasLiftTemperatureUpdater` existe porque o módulo de gas-lift precisa chamar
de volta três rotinas de temperatura que pertencem ao térmico. Sem ele, a
dependência seria circular.

## Conclusão

As 147 linhas de "diferença" não são trabalho não feito. São **22 linhas** que
os estágios anteriores já haviam removido dos mesmos corpos (contadas duas
vezes pelo alvo) e **125 linhas** de stubs e adaptador que o alvo simplesmente
não previu, embora o FR-032 e o contrato do módulo os exijam.

Perseguir o número exigiria violar o FR-032. O alvo deve ser corrigido para
1.162; a redução está completa.
