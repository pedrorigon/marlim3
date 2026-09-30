#include "SisProdSteadyStateSearch.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "celula3.h"
#include "celulaGas.h"
#include "chokegas.h"
#include "variaveisGlobais1D.h"

#include <math.h>

namespace sisprod::steady {

namespace {

/// The nine rows of the steady-state march dispatch table, named after the
/// SProd methods that run each march (SProd::marchaProdPerm1 runs
/// marchProductionSteady, and so on), so that each row reads against them.
enum class SteadyMarch {
    marchaInjPerm1,
    marchaGasPerm2,
    marchaGasPerm3,
    marchaProdPerm1,
    marchaProdPerm1Rev,
    marchaProdPerm2,
    marchaProdPresPres1,
    marchaProdPresPres1Rev,
    marchaProdPresPres2,
};

/// Resolves the four selectors to one row of the table: a pure function of
/// five values.
///
/// productionChokeOpening is the array, not the value, so that the subscript
/// happens only in the branch that needs it: Ler::copia_chokeSup leaves
/// chokep.abertura null when parserie is not positive, and reading it on every
/// dispatch would turn a conditional dereference into an unconditional one.
SteadyMarch selectSteadyMarch(int injectorWell, int prod, int tipoCC, int reverseMarch,
                              const double *productionChokeOpening) {
    if (injectorWell != 0)
        return SteadyMarch::marchaInjPerm1;
    if (prod == 0)
        return tipoCC == 0 ? SteadyMarch::marchaGasPerm2 : SteadyMarch::marchaGasPerm3;
    if (prod == 1) {
        if (tipoCC != 0)
            return SteadyMarch::marchaProdPerm2;
        return reverseMarch == 0 ? SteadyMarch::marchaProdPerm1
                                 : SteadyMarch::marchaProdPerm1Rev;
    }
    if (tipoCC != 0) {
        // Both arms select the same march, and that is correct: the secondary
        // pressure-to-pressure march already reduces to the tertiary one when the
        // choke is shut. The throat area is `abertura[0] * area`, so vazmaxSachd and
        // vazmassSachd both go to zero with it, and the residual becomes the same
        // `0. - MR` that marchProductionPressureToPressureTertiary returns literally.
        // Routing the else there would swap a guarded, general march for a narrower
        // one.
        return productionChokeOpening[0] > 1e-15 ? SteadyMarch::marchaProdPresPres2
                                                 : SteadyMarch::marchaProdPresPres2;
    }
    return reverseMarch == 0 ? SteadyMarch::marchaProdPresPres1
                             : SteadyMarch::marchaProdPresPres1Rev;
}

}  // namespace

double dispatchMarch(const SteadyStateSearchState &state, double chute, int prod, int tipoCC) {
    switch (selectSteadyMarch(state.march.input.pocinjec, prod, tipoCC, state.reverseSteady,
                              state.march.input.chokep.abertura)) {
    case SteadyMarch::marchaGasPerm2:         return marchGasSteadySecondary(state.march, chute);
    case SteadyMarch::marchaGasPerm3:         return marchGasSteadyTertiary(state.march, chute);
    case SteadyMarch::marchaProdPerm1:        return marchProductionSteady(state.march, chute);
    case SteadyMarch::marchaProdPerm1Rev:     return marchReverseProductionSteady(state.march, chute);
    case SteadyMarch::marchaProdPerm2:        return marchProductionSteadySecondary(state.march, chute);
    case SteadyMarch::marchaProdPresPres1:    return marchProductionPressureToPressure(state.march, chute);
    case SteadyMarch::marchaProdPresPres1Rev: return marchReverseProductionPressureToPressure(state.march, chute);
    case SteadyMarch::marchaProdPresPres2:    return marchProductionPressureToPressureSecondary(state.march, chute);
    // Falls out of the switch rather than returning inside it: an exhaustive
    // switch over a scoped enum still trips -Wreturn-type on GCC 11, and a
    // default: label would add an unreachable path.
    case SteadyMarch::marchaInjPerm1:         break;
    }
    return marchInjectionSteady(state.march, chute);
}

double solveSteadyRoot(const SteadyStateSearchState &state, double lowerBracket, double upperBracket, int prod, int tipoCC) {
    // Read here rather than in the solver: arq is input-deck configuration, and a
    // generic root finder has no business reading it. minit gates three early
    // returns inside the solver.
    int minit=0;
    if(state.march.input.acopColAnulPermForte == 1)minit=10;
    return rootfinding::zriddr(
        lowerBracket, upperBracket,
        [&](double guess) { return dispatchMarch(state, guess, prod, tipoCC); },
        // Domain feedback, kept on the domain side: the solver composes it as
        // monitor(objective(x)).
        [&](double residual) {
            double normalized = residual / state.march.baseConvergenceMonitor;
            if (prod != 0) {
                state.march.convergenceMonitor = fabs(normalized);
            }
            return normalized;
        },
        state.reverseSteady, minit);
}

double searchGasPressureSteadySecondary(const SteadyStateSearchState &state) {
    int valveCount = state.march.input.nvalvgas;
    double pchute;
    int maisprof = 0;
    for (int i = 1; i < valveCount; i++)
        if (state.march.productionValveCellIndices[i] < state.march.productionValveCellIndices[maisprof])
            maisprof = i;                      // busca a vÃ¡lvula mais profunda
    pchute = state.march.cells[state.march.productionValveCellIndices[maisprof]].pres; // estima uma pressao para a valvula mais profunda da linha,
    // a pressao na coluna de producao para esta valvula, observar que quando este metodo e chamado,
    // a marcha na linha de producao jÃ¡ foi feita
    double deepestValveTemperatureGuess;
    if (state.march.steadyIteration == 0) {                              // para o caso da primeira iteracao do sistema linha de producao-linha de gas
        deepestValveTemperatureGuess = state.march.cells[state.march.productionValveCellIndices[maisprof]].temp; // estimativa da temperatura na VGL mais profunda
        for (int k = state.march.gasCellCount; k > 0; k--) {           // marcha da ultima celula para a primeira
            // para fazer a estimativa da pressao de njecao
            double dx = 0.5 * (state.march.gasCells[k].dx0 + state.march.gasCells[k].dxL);
            double rhog = state.march.gasCells[k].flui.MasEspGas(pchute, deepestValveTemperatureGuess);
            pchute += rhog * 9.81 * sin(state.march.gasCells[k].duto.teta) * dx / kPascalPerKgfPerCm2; // avanco apenas pela hidrostatica
        }
    } else {
        for (int k = state.march.gasCellCount; k > 0; k--)
            pchute += state.march.updaters.steadyGasPressureDrop(k); // avanco usando os resultados da ireacao anterior, usando a hidrostatica
        // e a friccao
    }

    double delmas = 10;
    delmas = marchGasSteadySecondary(state.march, pchute); // marcha feita com a pressao de injecao estimada e a vazao de
    // injecao definida, retorna a diferenca entre a soma das vazÃµes em cada VGL calculadas na marcha
    // e a vazao de injecao definida, se negativa->pchute baixo, se positiva-> pchute alto
    double pchute2 = pchute;
    int kontaiter = 0;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    if (fabs(delmas) < 1e-3)
        return pchute;
    else {
        if (delmas < 0) { // pchute baixo, nova estimativa aumentando pchute, para encontrar um outro valor
            // de delmas positivo, para iniciar a falsa corda
            negativeResidualGuess = pchute2;
            while (delmas < 0 && kontaiter < 800) {
                pchute2 *= 1.1;
                delmas = marchGasSteadySecondary(state.march, pchute2);
                if (delmas < 0 && delmas > -0.9e10)
                    negativeResidualGuess = pchute2;
                kontaiter++;
            }
            if (kontaiter >= 800) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaGasPresPerm2 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            positiveResidualGuess = pchute2;
        } else if (delmas > 0) { // pchute alto, nova estimativa diminuindo pchute, para encontrar um outro valor
            // de delmas negativo, para iniciar a falsa corda
            positiveResidualGuess = pchute2;
            while (delmas > 0 && kontaiter < 800) {
                pchute2 *= 0.9;
                delmas = marchGasSteadySecondary(state.march, pchute2);
                if (delmas > 0 && delmas < 0.9e10)
                    positiveResidualGuess = pchute2;
                kontaiter++;
            }
            if (kontaiter >= 800) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaGasPresPerm2 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            negativeResidualGuess = pchute2;
        }
        return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 0, 0); // calculo de zero de funcao
    }
}

double searchGasPressureSteadyTertiary(const SteadyStateSearchState &state) {
    int valveCount = state.march.input.nvalvgas;
    // dois chutes de pressao a montante do choke de injecao, a convergencia para este caso,
    // choke de injecao e mais dificil, portanto, se testa mais possibiliodades de hutes:
    double pchute = state.march.injectionChoke.presEstag * 0.8; // pressao a montante do choke 20% menor que a pressao a jusante
    double pchute2 = state.march.gasCells[0].pres;         // Na primeira iteracao este valor Ã© a prÃ³pria pressao a montante

    // duas marchas para os dois chutes:
    double delmas;
    delmas = marchGasSteadyTertiary(state.march, pchute); // o valor retornado Ã© a diferenca entre a soma das vazoes VGL
    // e o valor calculado da vazao no choke de injecao, se der negativo, significa que a pressao a montante
    // se positivo, pressao a montante muito grande
    double delmas2;
    delmas2 = marchGasSteadyTertiary(state.march, pchute2);

    int kontaiter = 0;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    if (fabs(delmas) < 1e-3)
        return pchute;
    else {
        if (delmas < 0) { // chute de pressao pequeno, deve ser aumentado
            negativeResidualGuess = pchute;
            while (delmas < 0 && kontaiter < 800) {
                pchute *= 1.01;
                if (pchute > state.march.injectionChoke.presEstag)
                    pchute = state.march.injectionChoke.presEstag;
                delmas = marchGasSteadyTertiary(state.march, pchute);
                if (delmas < 0 && delmas > -0.9e10)
                    negativeResidualGuess = pchute;
                kontaiter++;
            }
            positiveResidualGuess = pchute;
        } else if (delmas >= 0) { // chute de pressao grande, deve ser diminuido
            kontaiter = 0;
            positiveResidualGuess = pchute;
            while (delmas > 0 && kontaiter < 800) {
                pchute *= 0.9;
                delmas = marchGasSteadyTertiary(state.march, pchute);
                if (delmas > 0 && delmas < 0.9e10)
                    positiveResidualGuess = pchute;
                kontaiter++;
            }
            negativeResidualGuess = pchute;
        }
        // caso nÃ£o funcione com o primeiro chute
        if (delmas2 >= 0 && kontaiter >= 800) { // chute de pressao pequeno, deve ser aumentado
            kontaiter = 0;
            positiveResidualGuess = pchute2;
            while (delmas2 > 0 && kontaiter < 800) {
                pchute2 *= 0.9;
                delmas2 = marchGasSteadyTertiary(state.march, pchute2);
                if (delmas2 > 0 && delmas2 < 0.9e10)
                    positiveResidualGuess = pchute2;
                kontaiter++;
            }
            if (kontaiter >= 800) { // falha na busca de estimativa para a falsa corda
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaGasPresPerm3 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            negativeResidualGuess = pchute2;
        } else if (delmas2 <= 0 && kontaiter >= 800) { // chute de pressao grande, deve ser diminuido
            kontaiter = 0;
            negativeResidualGuess = pchute2;
            while (delmas2 < 0 && kontaiter < 800) {
                pchute2 *= 1.01;
                if (pchute2 > state.march.injectionChoke.presEstag)
                    pchute2 = state.march.injectionChoke.presEstag;
                delmas2 = marchGasSteadyTertiary(state.march, pchute2);
                if (delmas2 < 0 && delmas2 > -0.9e10)
                    negativeResidualGuess = pchute2;
                kontaiter++;
            }
            if (kontaiter >= 800) { // falha na busca de estimativa para a falsa corda
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaGasPresPerm3 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            positiveResidualGuess = pchute2;
        }
        return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 0, 1); // calculo de zero de funcao
    }
}

namespace {

/// Brackets and solves the root for the reverse search.
///
/// Unlike the forward search this is one block rather than two arms: the reverse
/// march's residual does not split on sign the same way.
double bracketReverseRoot(const SteadyStateSearchState &state, double aumenta, double reduz, double &guessLowerBound, double &positiveResidualGuess, double &negativeResidualGuess, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double chute) {
    if (marchResidual < 0.) { // caso em que pressao a montante do choke < pressao da ultima celula, calculada pela
        // marcha, isto implica em pressao de chute alta , deve-se agora buscar
        // uma pressao de chute baixa para que marchResidual seja positivo e assim iniciar o processo
        // de calculo de zero de funcao
        negativeResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor negativo
        // a cada nova busca em que marchResidual se aproxima de zero, mas ainda negativo
        // negativeResidualGuess Ã© atualizado com a Ãºltima pressao de chute
        while (marchResidual < 0) {
            if (fabs(pchute2 - pchute) / pchute < (1. - reduz) / 10. && kontaiter > 100) {
                pchute2 *= 0.5;
            }
            pchuteAux = pchute2;
            pchute2 *= reduz; // Diminuindo a pressao na busca de marchResidual>0
            if (pchute2 <= guessLowerBound)
                pchute2 = 0.5 * (pchuteAux + guessLowerBound); // guessLowerBound inicialmente
            // e zero, mas pode acontecer de baixar demais pchute2 ao ponto de marchResidual=-1e10
            //(pressao abaixo de 0.5 no meio da marcha), neste caso, guessLowerBound se torna este valor de
            // pchute2, pois, com isto, ja se sabe que nao se pode ir abaixo de guessLowerBound
            marchResidual = marchReverseProductionSteady(state.march, pchute2);
            if (marchResidual < 0 && marchResidual > -0.9e10)
                negativeResidualGuess = pchute2; // atualizando o negativeResidualGuess
            kontaiter++;
            if (kontaiter > 100) { // limite de iteracoes, falha na busca do segundo chute
                // fim da simulacao ou aviso de falha
                if ((*state.march.globals).chaverede == 0) {
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    else {
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                } else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            while (marchResidual < -0.9e10) { // a reduÃ§Ã£o de pressao foi demais e a marcha nÃ£o foi capaz
                // de ir ate o final sem que a pressao ficasse inferior a 0.5kgf/cm2
                // deve-se aumentar a estimativa de pressao baixa
                guessLowerBound = pchute2;
                pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a menor pressao
                // que leva a marchResidual<0 e a pressao baixa demais
                marchResidual = marchReverseProductionSteady(state.march, pchute2);
                if (marchResidual < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = pchute2;
                kontaiter++;
                if (kontaiter > 100) {
                    if ((*state.march.globals).chaverede == 0) {
                        if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                            NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                        else {
                            cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    } else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        positiveResidualGuess = pchute2;
    } else if (marchResidual > 0.) { // caso em que pressao a montante do choke > pressao da ultima celula,
        // calculada pela
        // marcha, isto implica em pressao de chute baixa , deve-se agora buscar
        // uma pressao de chute alta para que marchResidual seja negativo e assim iniciar o processo
        // de calculo de zero de funcao
        positiveResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor positivo
        // a cada nova busca em que marchResidual se aproxima de zero, mas ainda positivo
        // positiveResidualGuess e atualizado com a Ãºltima pressao de chute
        while (marchResidual > 0) {
            pchuteAux = pchute2;
            pchute2 *= aumenta; // Aumentando a pressao na busca de marchResidual>0
            int limpres = 0;
            marchResidual = marchReverseProductionSteady(state.march, pchute2);
            if (marchResidual > 0 && marchResidual < 0.9e10)
                positiveResidualGuess = pchute2; // atualizando o positiveResidualGuess
            kontaiter++;
            if (kontaiter > 100) { // limite de iteracoes, falha na busca do segundo chute
                // fim da simulacao ou aviso de falha
                if ((*state.march.globals).chaverede == 0) {
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    else {
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                } else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            while (marchResidual > 0.9e10) { // o incremento de pressao foi demais e a marcha nao foi capaz
                // de ir ate o final sem que a pressao ficasse maior do que a pressao estatica
                // em um eventual IPR no meio da marcha ou maior que o limitre maximo de pressao
                // da tabela PVTSim, quando for este o caso
                // deve-se diminuir a estimativa de pressao baixa
                guessLowerBound = pchute2;
                pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a maior pressao
                // que leva a marchResidual>0 e a pressao alta demais
                marchResidual = marchReverseProductionSteady(state.march, pchute2);
                if (marchResidual > 0 && marchResidual < 0.9e10)
                    positiveResidualGuess = pchute2;
                kontaiter++;
                if (kontaiter > 100) {
                    if ((*state.march.globals).chaverede == 0) {
                        if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                            NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                        else {
                            cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;

                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    } else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        negativeResidualGuess = pchute2;
    }

    return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 0); // com as duas estimativas de pressao
    // em posicoes de sinal contrario da curva, inicia-se o processo de calculo de zero
    // de funcao
}

/// Walks the guess until the reverse march stops returning a sentinel.
bool retryUntilReverseMarchCompletes(const SteadyStateSearchState &state, double pchuteAux0, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute, double chute, double &abortValue) {
    if ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter <= 100) {
        double valtemp;
        valtemp = marchReverseProductionSteady(state.march, pchuteAux);   // marcha com pchuteAux
        while (valtemp > 0.9e10 && marchResidual > 0.9e10) { // estimativa de pressao de fundo ainda alta
            pchuteAux *= 0.99;                     // reduzindo a estimativa
            valtemp = marchReverseProductionSteady(state.march, pchuteAux);
            kontaiter++; // 50 iteracoes no maximo
        }
        while (valtemp < -0.9e10 && marchResidual < -0.9e10 && kontaiter < 100) { // estimativa de pressao de fundo ainda baixa
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection || state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection || state.march.cells[0].acsr.tipo == kAccessoryMultipleSource)
                pchuteAux *= 1.1;
            else
                pchuteAux *= 1.01; // aumentando a estimativa
            // verificando se este aumento ultrapassa o limite de pressao de uma eventual IPR no
            // fundo
            int limpres = 0;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchuteAux) < -0.01 * state.march.cells[0].acsr.ipr.Pres) {

                pchuteAux = (10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchuteAux > 1.01 * state.march.cells[0].acsr.ipr.Pres)
                    pchuteAux = 1.01 * state.march.cells[0].acsr.ipr.Pres;

                limpres = 1; // indicadpor de que este avanÃ§o de pressao esta muito alto
                // outros avancos devem ser feitos em um passo menor
            }
            // verificando se este aumento de estimativa foca acima do valor maximo de pressao
            // de uma eventual tabela PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
            valtemp = marchReverseProductionSteady(state.march, pchuteAux); // nova tentativa
            if (valtemp < -0.9e10 && limpres == 1) { // continua alto e existe o indicador
                // de que deve-se usar um passo de aumento de pressao menor
                int iterpres = 0;
                while (valtemp < -0.9e10 && iterpres < 10) { // novo laco com um passo menor,
                    // maximo de 10 iteracoes
                    pchuteAux *= 1.001;
                    valtemp = marchReverseProductionSteady(state.march, pchuteAux);
                    iterpres++;
                }
                if (iterpres >= 10) {
                    // caso em que atingiu o maximo de passos de incremento de pressa em uma situacao
                    // de passo pequeno, retorna um aviso que deu problema ou finaliza a simulacao
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        // neste caso, se finaliza a simulacao, nao tem um transiente
                        // a ser feito a seguir e nem se estÃ¡ em uma rede
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    else {
                        // apresenta apenas um aviso
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                        // se for em uma iteracao de rede, apos a primeira iteracao
                        if ((*state.march.globals).iterRede > 0)
                            {
                                abortValue = -1.1e10;
                                return true;
                            }
                        // se logo apÃƒÂ³s tem uma simulacao transiente ou se esta na primeira iteracao de rede
                        else
                            {
                                abortValue = 1.1e10;
                                return true;
                            }
                    }
                }
            }
            if (kontaiter == 100) {
                pchuteAux = 0.5 * pchuteAux0;
            }
            kontaiter++;
        }
        if (valtemp > -0.9e10 && valtemp < 0.9e10) { // conseguiu fazer a marcha atÃƒÂ© o final
            marchResidual = valtemp;
            pchute = pchuteAux;
        }
    }
    return false;
}

/// Moves the guess according to which sentinel the reverse march returned.
void classifyReverseMarchSentinel(const SteadyStateSearchState &state, double marchResidual, double &pchuteAux, double perdafric, double &taux) {
    if (marchResidual < -0.9e10) { // pressao muito baixa na marcha, deve ser aumentado o valor de chute
        // faz-se uma nova estimativa, sÃ³ que agora admitindo uma hidrostÃ¡tica de Ã¡gua,
        // o que darÃ¡ uma pressao de fundpo mais alta
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = 1000 + 0 * state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            double alfa = 0.;
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchuteAux -= state.march.cells[i - 1].acsr.delp;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) > -0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = (10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux > 1.01 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;  // limite de pchuteAux
    } else if (marchResidual > 0.9e10) { // pressÃƒÂ£o em algum ponto ficou acima de alguma pressao estatica
        // deve-se diminuir o valor do chute
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            // neste caso, utiliza-se uma fraÃ§Ã£o de vazio alta para a hidrostÃ¡tica
            double alfa = 0.8;
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.2 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchuteAux -= state.march.cells[i - 1].acsr.delp;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) > -0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = (10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux > 1.01 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;
    }
}

/// Estimates the bottom-hole pressure the reverse search starts from.
void estimateInitialReverseBottomHolePressure(const SteadyStateSearchState &state, double &completionFractionGuess, double &perdafric, double &frictionFactor, double &rmis, double &j, double &taux, double &pchute, double chute) {
    if (chute < 0) {

        if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection && fabs(state.march.cells[0].acsr.injl.QLiq) > 0.) {
            // este espaco faz uma estimativa de quanto deve ser a perda de carga media
            // a partir do valor da vazao no inicio da tubulacao, isto e feito apenas
            // se existir uma fonte de liquido, celula[0].acsr.tipo == 2
            taux = state.march.input.celp[0].textern; // a temperatura em que esta estimativa sera feita
            // e dada pela temperatura da fonte, a pressao para o calculo das propriedades
            // e feita por pGSup
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            ///////////propriedades fisicas:
            double completionFraction = state.march.cells[0].acsr.injl.bet;
            completionFractionGuess = completionFraction;
            double visC = state.march.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
            double visP = state.march.cells[0].acsr.injl.FluidoPro.ViscOleo(pchute, taux);
            double visG = state.march.cells[0].acsr.injl.FluidoPro.ViscGas(pchute, taux);
            double visMis = (1 - completionFraction) * visP + completionFraction * visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
            double liquidDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspLiq(pchute, taux);
            double gasDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspGas(pchute, taux);
            rmis = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
            double rlpA = state.march.cells[0].acsr.injl.FluidoPro.MasEspLiq(1., 15.);
            double rlcA = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
            // estimativa grosseira da vazao massica do liquido complementar
            double massicC = rlcA * state.march.cells[0].acsr.injl.QLiq * state.march.cells[0].acsr.injl.bet / kSecondsPerDay;
            // estimativa grosseira da vazao massica de liquido e gas produzidos
            double massic;
            // massas especificas, condicao standard
            double Rhogs = state.march.cells[0].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
            double Rhols = (1000 * 141.5 / (131.5 + state.march.cells[0].acsr.injl.FluidoPro.API)) * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW) + 1000. * state.march.cells[0].acsr.injl.FluidoPro.Denag * state.march.cells[0].acsr.injl.FluidoPro.BSW;
            // multiplicador da vazao standard para se obter a Vazao massicagÃ¡s+liquido produzido
            double multiplicador = (Rhols + state.march.cells[0].acsr.injl.FluidoPro.RGO * Rhogs * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW));
            // massic *= multiplicador;//alteracao8
            massic = 1 * multiplicador * state.march.cells[0].acsr.injl.QLiq * (1. - state.march.cells[0].acsr.injl.bet) / kSecondsPerDay;
            // titulo do gas em relacao a misturaoleo+agua+gas
            double fracmasshidra = state.march.cells[0].acsr.injl.FluidoPro.FracMassHidra(pchute, taux);
            double massicP = (1. - fracmasshidra) * massic; // vazao massica de liquido produzido
            double massicG = fracmasshidra * massic;        // vazao massaica de gas
            /// estimativa grosseira da velocidade da mistura
            j = (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess) / state.march.cells[0].duto.area;
            // estimativa da fracao de vazio sem escorregamento
            double alfmis = (massicG / gasDensityAtGuess) / (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess);
            // propriedades de mistura
            rmis = (1 - alfmis) * rmis + alfmis * gasDensityAtGuess;
            visMis = (1 - alfmis) * visMis + alfmis * visG;
            // estimativa de numero de Reynolds
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            // estimativa de um fator de friccao
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);
        }
        double tmparea = state.march.cells[0].duto.area;
        double tmpperi = state.march.cells[0].duto.peri;
        // estimativa da perda por friccao no sistema, se o acessorio no inicio da tubulacao
        // nao for uma fonte de liquido, esta estimativa sera zero
        perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * tmpperi / tmparea;
        // se o acessorio for uma IPR, a estimativa serÃ¡ uma pressao prÃ³xima a uma vazao baixa
        // no fundo do poco
        // se nao for IPR, a estimativa continua, considerando a perda de carga e a hidrostatica,
        // que serÃ¡ avaliada neste laco
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchute, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchute, taux);
            // se o sistema nao for um anel principal de GL
            // para a estimativa de pressao, se utilizara so a hidrostatica de liquido, o que dara um
            // chute inicial de pressao muito alta
            double alfa = 0.;
            if (completionFractionGuess < 0.5) {
                double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
                alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
            }
            if (state.march.annulusDrift == 0)
                alfa = 1.; // se for o anel, so se tem gas
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            // avanco da pressao por meio da hidrostatica e da perda por friccao estimada
            pchute += ((rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2);
            // incremento de pressao devido a algum ganho de pressao constante em alguma celula
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchute += state.march.cells[i - 1].acsr.delp;
            // caso exista alguma IPR no meio do duto, verifica se a pressao estimada esta maior que
            // a pressao estatica, nÃ£o se trabalha com vazoes negativas neste solver
            // para evitar isto, se corrige a pressao para um valor proximo da
            // pressao estatica da IPR da celula
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchute) > -0.01 * state.march.cells[i - 1].acsr.ipr.Pres) && i == 1) {
                pchute = (10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute > 1.1 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 1.1 * state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute < 1.01 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            // tambÃ©m se evita um chute com valores maiores que o limite maximo de uma tabela
            // no caso de se estar usando PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
        }

        if (pchute > 1000)
            pchute = 1000; // pressao maxima de chute
    } else
        pchute = chute; // caso chute nao seja negativo utiliza a estimativa enviada na
}

}  // namespace

double searchReverseProductionBottomHolePressure(const SteadyStateSearchState &state, double chute) {
    state.reverseSteady = 1;
    state.march.convergenceMonitor = 1000.;
    // busca de dois chutes iniciais com valores com sinais opostos
    // para marchaProdPerm1 e assim iniciar o prpocesso de calculo de erro de funcao.
    double pchute = state.march.gasSurfacePressure; // inicializando o valor de pchute com o valor da pressao a jusante
    // do choke, pchute sera o valor de chute de fato no processo de busca
    // se o valor de chute>0, pchute=chute, senao, ele e estimado
    double taux; // valor de temperatura auxiliar para o eventual calculo
    // de pchute
    double j = 0.;
    double rmis = 0.;
    double frictionFactor = 0.;
    double perdafric = 0.;
    double completionFractionGuess = 0.;

    estimateInitialReverseBottomHolePressure(state, completionFractionGuess, perdafric, frictionFactor, rmis, j, taux, pchute, chute);
    // lista de parametro do metodo
    double pchute2;        // segundo chute de pressao da busca
    double pchuteAux = 0.; // auxiliar na busca de dos chutes de pressao
    // necessarios para se dar partida no metodo de calculo de raiz
    double marchResidual; // marchResidual recebe sempre o valor retornado pelo procedimento de marcha
    // neste caso, Ã© a pressao a montante do choke-a pressao da ultima celula, calculada pela
    // marcha
    marchResidual = marchReverseProductionSteady(state.march, pchute);
    state.march.searchOrigin = 1;
    // a marcha pode dar problemas, ou a pressao em algum ponto deu acima da pressao estÃ¡tica de
    // de alguma IPR no meio do caminho, retorna 1e10,
    // ou antes de atingir a ultima celula, a pressao ficou proximo de zero, ou a pressao
    // retorna -1e10, ou a pressao, para o caso PVTSim, ficou acima da pressao maxima da tabela
    // retorna 1e10
    classifyReverseMarchSentinel(state, marchResidual, pchuteAux, perdafric, taux);


    // para o caso de ter dado eerado a primeira marcha, com pchuteAux faz-se uma nova tentativa
    int kontaiter = 0; // contador para o laco em que se tentara uma nova estimativa
    // em que ao menos os valores 1e10 ou 1e-10 nÃ£o seja retornados
    double pchuteAux0 = pchuteAux;
    double abortValue;
    if (retryUntilReverseMarchCompletes(state, pchuteAux0, kontaiter, marchResidual, pchuteAux, pchute, chute, abortValue))
        return abortValue;
    if (kontaiter > 100) { // chegou ao limite da iteracao
        if ((*state.march.globals).chaverede == 0) {
            // neste caso, se finaliza a simulacao, nao tem um transiente
            // a ser feito a seguir e nem se estÃ¡ em uma rede
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
            else {
                // apresenta apenas um aviso
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                // se for em uma iteracao de rede, apos a primeira iteracao
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                // se logo apos tem uma simulacao transiente ou se esta na primeira iteracao de rede
                else
                    return 1.1e10;
            }
        } else {
            // se for em uma iteracao de rede, apos a primeira iteracao
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            // se logo apÃ³s tem uma simulacao transiente ou se esta na primeira iteracao de rede
            else
                return 1.1e10;
        }
    }

    // caso o laÃ§o anterior tenha levado a estimativa de um valor 1e10 para um valor -1e10
    // ou tenha levado para um valor de -1e10 para 1e10
    // ou seja, de um chute de pressao alto demais para finalizar a marcha
    // para um chute de pressao baixo demais para terminar a marcha ou vice versa
    while ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter < 100) {
        pchute = 0.5 * (pchute + pchuteAux); // busca-se um valor medio
        // entre o alto demais e o baixo demais
        if ((fabs(pchute - pchuteAux) / pchute) < 0.001 && marchResidual < -0.9e10) {
            pchute = 1.05 * pchute;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[0].acsr.ipr.Pres - pchute) > 0.00 * state.march.cells[0].acsr.ipr.Pres) {
                pchute = (10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute > 1.01 * state.march.cells[0].acsr.ipr.Pres)
                    pchute = 1.01 * state.march.cells[0].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
            pchuteAux = pchute;
        }
        marchResidual = marchReverseProductionSteady(state.march, pchute);
        kontaiter++;
    }
    if (kontaiter >= 100) {
        if ((*state.march.globals).chaverede == 0) {
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
            else {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        } else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    // caso tenha sido possivel obter um chute de pressao em que foi possivel finalizar a marcha
    kontaiter = 0;
    pchute2 = pchute;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    double guessLowerBound = 0.;
    state.march.input.buscaFC = fabs(state.march.input.buscaFC);
    double reduz = 1. - state.march.input.buscaFC;
    double amplifica = 1. + state.march.input.buscaFC;
    double aumenta = amplifica;
    if (fabs(marchResidual) < 1e-3)
        return pchute;
    else {
        return bracketReverseRoot(state, aumenta, reduz, guessLowerBound, positiveResidualGuess, negativeResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, pchute, chute);
    }
}

namespace {

/// Raises the guess until the march stops failing on too large an increment.
///
/// Inner loop of bracketFromLowGuess: the outer loop moves the guess, this one
/// backs off when the move overshot and the march returned a sentinel.
bool raisePressureUntilMarchCompletes(const SteadyStateSearchState &state, double &guessLowerBound, double &positiveResidualGuess, int &kontaiter, double &marchResidual, double pchuteAux, double &pchute2, double chute, int kontaTenta, double &abortValue) {
    // de ir ate o final sem que a pressao ficasse maior do que a pressao estatica
    // em um eventual IPR no meio da marcha ou maior que o limitre maximo de pressao
    // da tabela PVTSim, quando for este o caso
    // deve-se diminuir a estimativa de pressao baixa
    guessLowerBound = pchute2;
    pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a maior pressao
    // que leva a marchResidual>0 e a pressao alta demais
    marchResidual = marchProductionSteady(state.march, pchute2);
    if (fabs(marchResidual) > 1.01e10) {
        if (kontaTenta < 0)
            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                       "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                       "", "");
        cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
        if ((*state.march.globals).iterRede > 0)
            {
                abortValue = -1.1e10;
                return true;
            }
        else
            {
                abortValue = 1.1e10;
                return true;
            }
    }
    if (marchResidual > 0 && marchResidual < 0.9e10)
        positiveResidualGuess = pchute2;
    kontaiter++;
    if (kontaiter > 50 * 0.1 / state.march.input.buscaFC) {
        if ((*state.march.globals).chaverede == 0) {
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            } else {
                if (kontaTenta < 0)
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;

                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        } else {
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
    }
    return false;
}

/// Brackets the root upward when the first march came back positive.
///
/// One arm of the sign split that ends searchProductionBottomHolePressure. Its
/// twin is bracketFromHighGuess; the two are not each other's mirror, which is
/// why they are two functions and not one with a sign parameter.
bool bracketFromLowGuess(const SteadyStateSearchState &state, double &amplifica, double &reduz, int &reversao, double &val0, double &guessLowerBound, double &positiveResidualGuess, double &negativeResidualGuess, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double chute, int kontaTenta, double &abortValue) {
    // calculada pela
    // marcha, isto implica em pressao de chute baixa , deve-se agora buscar
    // uma pressao de chute alta para que marchResidual seja negativo e assim iniciar o processo
    // de calculo de zero de funcao
    positiveResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor positivo
    // a cada nova busca em que marchResidual se aproxima de zero, mas ainda positivo
    // positiveResidualGuess e atualizado com a Ãºltima pressao de chute
    int kontaReverso = 0;
    while (marchResidual > 0) {
        pchuteAux = pchute2;
        pchute2 *= amplifica; // Aumentando a pressao na busca de marchResidual>0
        int limpres = 0;
        if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
            (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
            if (chute < 0) {
                pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                    pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
                if (fabs(pchute2 - pchute) < 1e-15) {
                    pchute2 = -(1 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                    if (pchute2 < 0.9999 * state.march.cells[0].acsr.ipr.Pres)
                        pchute2 = 0.9999 * state.march.cells[0].acsr.ipr.Pres;
                }
            } else {
                pchute2 = pchuteAux * 1.0001;
            }

            limpres = 1;
        }
        marchResidual = marchProductionSteady(state.march, pchute2);
        if (fabs(marchResidual) > 1.01e10) {
            if (kontaTenta < 0)
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (marchResidual > 0. && limpres == 1) {
            int iterpres = 0;
            if (iterpres >= 10)
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                } else {
                    if (kontaTenta < 0)
                        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                   "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                   "", "");
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
        }
        if (marchResidual > 0 && marchResidual < 0.9e10) {
            positiveResidualGuess = pchute2; // atualizando o positiveResidualGuess
            if (marchResidual > val0 && state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance && state.march.input.lingas == 0) {

                kontaReverso++;
                if (kontaReverso == 2) {
                    reversao = 1;
                    marchResidual = -1;
                } else {
                    reduz = 1. + state.march.input.buscaFC;
                    amplifica = 1. - state.march.input.buscaFC;
                    val0 = marchResidual;
                }

            } else
                val0 = marchResidual;
        }
        kontaiter++;
        if (kontaiter > 50 * 0.1 / state.march.input.buscaFC) { // limite de iteracoes, falha na busca do segundo chute
            // fim da simulacao ou aviso de falha
            if ((*state.march.globals).chaverede == 0) {
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                } else {
                    if (kontaTenta < 0)
                        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                   "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                   "", "");
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            } else {
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        while (marchResidual > 0.9e10) { // o incremento de pressao foi demais e a marcha nao foi capaz
            if (raisePressureUntilMarchCompletes(state, guessLowerBound, positiveResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, chute, kontaTenta, abortValue))
                return true;
        }
    }
    negativeResidualGuess = pchute2;
    return false;
}

/// Brackets the root downward when the first march came back negative.
///
/// One arm of the sign split that ends searchProductionBottomHolePressure.
bool bracketFromHighGuess(const SteadyStateSearchState &state, double &amplifica, double &reduz, int &reversao, double &val0, double &guessLowerBound, double &negativeResidualGuess, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double chute, int kontaTenta, double &abortValue) {
    // marcha, isto implica em pressao de chute alta , deve-se agora buscar
    // uma pressao de chute baixa para que marchResidual seja positivo e assim iniciar o processo
    // de calculo de zero de funcao
    negativeResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor negativo
    // a cada nova busca em que marchResidual se aproxima de zero, mas ainda negativo
    // negativeResidualGuess Ã© atualizado com a Ãºltima pressao de chute
    int kontaReverso = 0;
    while (marchResidual < 0) {
        if (fabs(pchute2 - pchute) / pchute < (1. - reduz) / 10. && kontaiter > 50 * 0.1 / state.march.input.buscaFC) {
            pchute2 *= 0.5;
        }
        pchuteAux = pchute2;
        pchute2 *= reduz; // Diminuindo a pressao na busca de marchResidual>0
        if (pchute2 <= guessLowerBound)
            pchute2 = 0.5 * (pchuteAux + guessLowerBound); // guessLowerBound inicialmente
        // e zero, mas pode acontecer de baixar demais pchute2 ao ponto de marchResidual=-1e10
        //(pressao abaixo de 0.5 no meio da marcha), neste caso, guessLowerBound se torna este valor de
        // pchute2, pois, com isto, ja se sabe que nao se pode ir abaixo de guessLowerBound
        marchResidual = marchProductionSteady(state.march, pchute2);
        if (fabs(marchResidual) > 1.01e10) {
            if (kontaTenta < 0)
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (marchResidual < 0 && marchResidual > -0.9e10) {
            negativeResidualGuess = pchute2; // atualizando o negativeResidualGuess
            if (marchResidual < val0 && state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance && state.march.input.lingas == 0) {

                kontaReverso++;
                if (kontaReverso == 2) {
                    reversao = 1;
                    marchResidual = 1;
                } else {
                    reduz = 1. + state.march.input.buscaFC;
                    amplifica = 1. - state.march.input.buscaFC;
                    val0 = marchResidual;
                }
            } else
                val0 = marchResidual;
        }
        kontaiter++;
        if (kontaiter > 50 * 0.1 / state.march.input.buscaFC) { // limite de iteracoes, falha na busca do segundo chute
            // fim da simulacao ou aviso de falha
            if ((*state.march.globals).chaverede == 0) {
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                } else {
                    if (kontaTenta < 0)
                        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                   "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                   "", "");
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            } else {
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        while (marchResidual < -0.9e10) { // a reduÃ§Ã£o de pressao foi demais e a marcha nÃ£o foi capaz
            // de ir ate o final sem que a pressao ficasse inferior a 0.5kgf/cm2
            // deve-se aumentar a estimativa de pressao baixa
            guessLowerBound = pchute2;
            pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a menor pressao
            // que leva a marchResidual<0 e a pressao baixa demais
            marchResidual = marchProductionSteady(state.march, pchute2);
            if (fabs(marchResidual) > 1.01e10) {
                if (kontaTenta < 0)
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (marchResidual < 0 && marchResidual > -0.9e10) {
                negativeResidualGuess = pchute2;
            }
            kontaiter++;
            if (kontaiter > 50 * 0.1 / state.march.input.buscaFC) {
                if ((*state.march.globals).chaverede == 0) {
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        if (kontaTenta < 0)
                            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                       "", "");
                        else {
                            cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                            if ((*state.march.globals).iterRede > 0)
                                {
                                    abortValue = -1.1e10;
                                    return true;
                                }
                            else
                                {
                                    abortValue = 1.1e10;
                                    return true;
                                }
                        }
                } else {
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            }
        }
    }
    return false;
}

/// Walks the guess until the march stops returning a sentinel.
///
/// The march reports failure by returning 1e10 or -1e10 rather than by any other
/// means, so the search has to read the magnitude to know what happened.
bool retryUntilMarchCompletes(const SteadyStateSearchState &state, double pchuteAux0, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute, double chute, int kontaTenta, double &abortValue) {
    if ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter < 50) {
        double valtemp;
        valtemp = marchProductionSteady(state.march, pchuteAux); // marcha com pchuteAux
        if (fabs(valtemp) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if (kontaTenta < 0)
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        while (valtemp > 0.9e10 && marchResidual > 0.9e10) { // estimativa de pressao de fundo ainda alta
            pchuteAux *= 0.99;                     // reduzindo a estimativa
            valtemp = marchProductionSteady(state.march, pchuteAux);
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if (kontaTenta < 0)
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            kontaiter++; // 50 iteracoes no maximo
        }
        while (valtemp < -0.9e10 && marchResidual < -0.9e10 && kontaiter <= 50) { // estimativa de pressao de fundo ainda baixa
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection || state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection || state.march.cells[0].acsr.tipo == kAccessoryMultipleSource)
                pchuteAux *= 1.1;
            else
                pchuteAux *= 1.01; // aumentando a estimativa
            // verificando se este aumento ultrapassa o limite de pressao de uma eventual IPR no
            // fundo
            int limpres = 0;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[0].acsr.ipr.Pres) {

                pchuteAux = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[0].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[0].acsr.ipr.Pres;

                limpres = 1; // indicadpor de que este avanÃ§o de pressao esta muito alto
                // outros avancos devem ser feitos em um passo menor
            }
            // verificando se este aumento de estimativa foca acima do valor maximo de pressao
            // de uma eventual tabela PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
            valtemp = marchProductionSteady(state.march, pchuteAux); // nova tentativa
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if (kontaTenta < 0)
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (valtemp < -0.9e10 && limpres == 1) { // continua alto e existe o indicador
                // de que deve-se usar um passo de aumento de pressao menor
                int iterpres = 0;
                while (valtemp < -0.9e10 && iterpres < 10) { // novo laco com um passo menor,
                    // maximo de 10 iteracoes
                    pchuteAux *= 1.001;
                    valtemp = marchProductionSteady(state.march, pchuteAux);
                    if (fabs(valtemp) > 1.01e10) {
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                        if (kontaTenta < 0)
                            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                       "", "");
                        if ((*state.march.globals).iterRede > 0)
                            {
                                abortValue = -1.1e10;
                                return true;
                            }
                        else
                            {
                                abortValue = 1.1e10;
                                return true;
                            }
                    }
                    iterpres++;
                }
                if (iterpres >= 10) {
                    // caso em que atingiu o maximo de passos de incremento de pressa em uma situacao
                    // de passo pequeno, retorna um aviso que deu problema ou finaliza a simulacao
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                        // neste caso, se finaliza a simulacao, nao tem um transiente
                        // a ser feito a seguir e nem se estÃ¡ em uma rede
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                   "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                   "", "");
                    } else {
                        // apresenta apenas um aviso
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                        if (kontaTenta < 0)
                            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                                       "", "");
                        // se for em uma iteracao de rede, apos a primeira iteracao
                        if ((*state.march.globals).iterRede > 0)
                            {
                                abortValue = -1.1e10;
                                return true;
                            }
                        // se logo apos tem uma simulacao transiente ou se esta na primeira iteracao de rede
                        else
                            {
                                abortValue = 1.1e10;
                                return true;
                            }
                    }
                }
            }
            if (kontaiter == 50) {
                pchuteAux = 0.5 * pchuteAux0;
            }
            kontaiter++;
        }
        if (valtemp > -0.9e10 && valtemp < 0.9e10) { // conseguiu fazer a marcha atÃƒÂ© o final
            marchResidual = valtemp;
            pchute = pchuteAux;
        }
    }
    return false;
}

/// Moves the guess according to which sentinel the march returned.
///
/// -1e10 means the pressure fell near zero before the last cell, 1e10 means it
/// rose above a static pressure or above the table's maximum.
void classifyMarchSentinel(const SteadyStateSearchState &state, double marchResidual, double &pchuteAux, double completionFractionGuess, double perdafric, double &taux, double pchute) {
    if (marchResidual < -0.9e10) { // pressao muito baixa na marcha, deve ser aumentado o valor de chute
        // faz-se uma nova estimativa, sÃ³ que agora admitindo uma hidrostÃ¡tica de Ã¡gua,
        // o que darÃ¡ uma pressao de fundpo mais alta
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = 1000 + 0 * state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            double alfa = 0.;
            if (completionFractionGuess < 0.5) {
                double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
                alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
                alfa *= 0.5;
            }
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchuteAux -= state.march.cells[i - 1].acsr.delp;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = -(10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1000) {
            pchuteAux = 1000; // limite de pchuteAux
        }
    } else if (marchResidual > 0.9e10) { // pressao em algum ponto ficou acima de alguma pressao estatica
        // deve-se diminuir o valor do chute
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            // neste caso, utiliza-se uma fracao de vazio alta para a hidrostÃ¡tica
            double alfa = 0.8;
            if (completionFractionGuess < 0.5) {
                double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
                alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
                if (alfa > 0.9999)
                    alfa = 0.9999;
            }
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.2 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchuteAux -= state.march.cells[i - 1].acsr.delp;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = -(10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1000) {
            pchuteAux = 1000;
        }
    }
}

/// Estimates the bottom-hole pressure the search starts from.
///
/// With a negative guess the pressure is built from the surface pressure, the head
/// accessory and a friction estimate; otherwise the caller's guess stands. j, rmis
/// and frictionFactor became locals here rather than parameters: the range writes them and
/// nothing reads them afterwards.
void estimateInitialBottomHolePressure(const SteadyStateSearchState &state, double &completionFractionGuess, double &perdafric, double &frictionFactor, double &rmis, double &j, double &taux, double &pchute, double chute) {
    if (chute < 0) {

        if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection && fabs(state.march.cells[0].acsr.injl.QLiq) > 0.) {
            completionFractionGuess = state.march.cells[0].acsr.injl.bet;
            // este espaco faz uma estimativa de quanto deve ser a perda de carga media
            // a partir do valor da vazao no inicio da tubulacao, isto e feito apenas
            // se existir uma fonte de liquido, celula[0].acsr.tipo == 2
            taux = state.march.input.celp[0].textern; // a temperatura em que esta estimativa sera feita
            // e dada pela temperatura da fonte, a pressao para o calculo das propriedades
            // e feita por pGSup
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            ///////////propriedades fisicas:
            double completionFraction = state.march.cells[0].acsr.injl.bet;
            double visC = state.march.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
            double visP = state.march.cells[0].acsr.injl.FluidoPro.ViscOleo(pchute, taux);
            double visG = state.march.cells[0].acsr.injl.FluidoPro.ViscGas(pchute, taux);
            double visMis = (1 - completionFraction) * visP + completionFraction * visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
            double liquidDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspLiq(pchute, taux);
            double gasDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspGas(pchute, taux);
            rmis = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
            double rlcA = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
            // estimativa grosseira da vazao massica do liquido complementar
            double massicC = rlcA * state.march.cells[0].acsr.injl.QLiq * state.march.cells[0].acsr.injl.bet / kSecondsPerDay;
            // estimativa grosseira da vazao massica de liquido e gas produzidos
            double massic;
            // massas especificas, condicao standard
            double Rhogs = state.march.cells[0].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
            double Rhols = (1000 * 141.5 / (131.5 + state.march.cells[0].acsr.injl.FluidoPro.API)) * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW) + 1000. * state.march.cells[0].acsr.injl.FluidoPro.Denag * state.march.cells[0].acsr.injl.FluidoPro.BSW;
            // multiplicador da vazao standard para se obter a Vazao massicagÃ¡s+liquido produzido
            double multiplicador = (Rhols + state.march.cells[0].acsr.injl.FluidoPro.RGO * Rhogs * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW));
            // massic *= multiplicador;//alteracao8
            massic = 1 * multiplicador * state.march.cells[0].acsr.injl.QLiq * (1. - state.march.cells[0].acsr.injl.bet) / kSecondsPerDay;
            // titulo do gas em relacao a misturaoleo+agua+gas
            double fracmasshidra = state.march.cells[0].acsr.injl.FluidoPro.FracMassHidra(pchute, taux);
            double massicP = (1. - fracmasshidra) * massic; // vazao massica de liquido produzido
            double massicG = fracmasshidra * massic;        // vazao massaica de gas
            /// estimativa grosseira da velocidade da mistura
            j = (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess) / state.march.cells[0].duto.area;
            // estimativa da fracao de vazio sem escorregamento
            double alfmis = (massicG / gasDensityAtGuess) / (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess);
            // propriedades de mistura
            rmis = (1 - alfmis) * rmis + alfmis * gasDensityAtGuess;
            visMis = (1 - alfmis) * visMis + alfmis * visG;
            // estimativa de numero de Reynolds
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            // estimativa de um fator de friccao
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);
        }
        double tmparea = state.march.cells[0].duto.area;
        double tmpperi = state.march.cells[0].duto.peri;
        // estimativa da perda por friccao no sistema, se o acessorio no inicio da tubulacao
        // nao for uma fonte de liquido, esta estimativa sera zero
        perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * tmpperi / tmparea;
        // se o acessorio for uma IPR, a estimativa serÃ¡ uma pressao prÃ³xima a uma vazao baixa
        // no fundo do poco
        // se nao for IPR, a estimativa continua, considerando a perda de carga e a hidrostatica,
        // que serÃ¡ avaliada neste laco
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchute, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchute, taux);
            // se o sistema nao for um anel principal de GL
            // para a estimativa de pressao, se utilizara so a hidrostatica de liquido, o que dara um
            // chute inicial de pressao muito alta
            double alfa = 0.;
            if (completionFractionGuess < 0.5) {
                double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
                alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
            }
            if (state.march.annulusDrift == 0)
                alfa = 1.; // se for o anel, so se tem gas
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            // avanco da pressao por meio da hidrostatica e da perda por friccao estimada
            pchute += ((rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2);
            if (pchute < 0.8)
                pchute = 0.8;
            // incremento de pressao devido a algum ganho de pressao constante em alguma celula
            if (state.march.cells[i - 1].acsr.tipo == 7)
                pchute -= state.march.cells[i - 1].acsr.delp;
            // caso exista alguma IPR no meio do duto, verifica se a pressao estimada esta maior que
            // a pressao estatica, nÃ£o se trabalha com vazoes negativas neste solver
            // para evitar isto, se corrige a pressao para um valor proximo da
            // pressao estatica da IPR da celula
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchute) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres) && i == 1) {
                double flowRateGuess = 0.15 * state.march.cells[i - 1].duto.area * kSecondsPerDay;
                pchute = -(flowRateGuess / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute > 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute < 0.9 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 0.9 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            // tambÃ©m se evita um chute com valores maiores que o limite maximo de uma tabela
            // no caso de se estar usando PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
        }

        if (pchute > 1000) {
            pchute = 1000; // pressao maxima de chute
        }
    } else
        pchute = chute; // caso chute nao seja negativo utiliza a estimativa enviada na
}

}  // namespace

double searchProductionBottomHolePressure(const SteadyStateSearchState &state, double chute, int kontaTenta) {
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    // busca de dois chutes iniciais com valores com sinais opostos
    // para marchaProdPerm1 e assim iniciar o prpocesso de calculo de erro de funcao.
    double pchute = state.march.gasSurfacePressure; // inicializando o valor de pchute com o valor da pressao a jusante
    // do choke, pchute sera o valor de chute de fato no processo de busca
    // se o valor de chute>0, pchute=chute, senao, ele e estimado
    double taux; // valor de temperatura auxiliar para o eventual calculo
    // de pchute
    double j = 0.;
    double rmis = 0.;
    double frictionFactor = 0.;
    double perdafric = 0.;
    double completionFractionGuess = 0.;

    estimateInitialBottomHolePressure(state, completionFractionGuess, perdafric, frictionFactor, rmis, j, taux, pchute, chute);
    // lista de parametro do metodo
    double pchute2;        // segundo chute de pressao da busca
    double pchuteAux = 0.; // auxiliar na busca de dos chutes de pressao
    // necessarios para se dar partida no metodo de calculo de raiz
    double marchResidual; // marchResidual recebe sempre o valor retornado pelo procedimento de marcha
    // neste caso, Ã© a pressao a montante do choke-a pressao da ultima celula, calculada pela
    // marcha
    marchResidual = marchProductionSteady(state.march, pchute);
    state.march.searchOrigin = 1;
    // a marcha pode dar problemas, ou a pressao em algum ponto deu acima da pressao estÃ¡tica de
    // de alguma IPR no meio do caminho, retorna 1e10,
    // ou antes de atingir a ultima celula, a pressao ficou proximo de zero, ou a pressao
    // retorna -1e10, ou a pressao, para o caso PVTSim, ficou acima da pressao maxima da tabela
    // retorna 1e10
    classifyMarchSentinel(state, marchResidual, pchuteAux, completionFractionGuess, perdafric, taux, pchute);


    // para o caso de ter dado eerado a primeira marcha, com pchuteAux faz-se uma nova tentativa
    int kontaiter = 0; // contador para o laco em que se tentara uma nova estimativa
    // em que ao menos os valores 1e10 ou 1e-10 nao seja retornados
    double pchuteAux0 = pchuteAux;
    double abortValue;
    if (retryUntilMarchCompletes(state, pchuteAux0, kontaiter, marchResidual, pchuteAux, pchute, chute, kontaTenta, abortValue))
        return abortValue;
    if (kontaiter > 50) { // chegou ao limite da iteracao
        if ((*state.march.globals).chaverede == 0) {
            // neste caso, se finaliza a simulacao, nao tem um transiente
            // a ser feito a seguir e nem se estÃ¡ em uma rede
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            } else {
                // apresenta apenas um aviso
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if (kontaTenta < 0)
                    logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                               "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                               "", "");
                // se for em uma iteracao de rede, apos a primeira iteracao
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                // se logo apos tem uma simulacao transiente ou se esta na primeira iteracao de rede
                else
                    return 1.1e10;
            }
        } else {
            // se for em uma iteracao de rede, apos a primeira iteracao
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            // se logo apÃ³s tem uma simulacao transiente ou se esta na primeira iteracao de rede
            else
                return 1.1e10;
        }
    }

    // caso o laÃ§o anterior tenha levado a estimativa de um valor 1e10 para um valor -1e10
    // ou tenha levado para um valor de -1e10 para 1e10
    // ou seja, de um chute de pressao alto demais para finalizar a marcha
    // para um chute de pressao baixo demais para terminar a marcha ou vice versa
    while ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter < 50) {
        pchute = 0.5 * (pchute + pchuteAux); // busca-se um valor medio
        // entre o alto demais e o baixo demais
        if ((fabs(pchute - pchuteAux) / pchute) < 0.001 && marchResidual < -0.9e10) {
            pchute = 1.05 * pchute;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[0].acsr.ipr.Pres - pchute) < 0.00 * state.march.cells[0].acsr.ipr.Pres) {
                pchute = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute < 0.99 * state.march.cells[0].acsr.ipr.Pres)
                    pchute = 0.99 * state.march.cells[0].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
            pchuteAux = pchute;
        }
        marchResidual = marchProductionSteady(state.march, pchute);
        if (fabs(marchResidual) > 1.01e10) {
            if (kontaTenta < 0)
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
        kontaiter++;
    }
    if (kontaiter >= 50) {
        if ((*state.march.globals).chaverede == 0) {
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0) {
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                           "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes",
                           "", "");
            } else {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        } else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    // caso tenha sido possivel obter um chute de pressao em que foi possivel finalizar a marcha
    kontaiter = 0;
    pchute2 = pchute;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    double guessLowerBound = 0.;
    double val0 = marchResidual;
    int reversao = 0;
    double reduz = 1. - state.march.input.buscaFC;
    double amplifica = 1. + state.march.input.buscaFC;
    if (fabs(marchResidual) < 1e-3)
        return pchute;
    else {
        if (marchResidual < 0.) { // caso em que pressao a montante do choke < pressao da ultima celula, calculada pela
            double abortValue;
            if (bracketFromHighGuess(state, amplifica, reduz, reversao, val0, guessLowerBound, negativeResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, pchute, chute, kontaTenta, abortValue))
                return abortValue;
            positiveResidualGuess = pchute2;
        } else if (marchResidual > 0.) { // caso em que pressao a montante do choke > pressao da ultima celula,
            double abortValue;
            if (bracketFromLowGuess(state, amplifica, reduz, reversao, val0, guessLowerBound, positiveResidualGuess, negativeResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, pchute, chute, kontaTenta, abortValue))
                return abortValue;
        }

        if (reversao == 0)
            return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 0); // com as duas estimativas de pressao
        // em posicoes de sinal contrario da curva, inicia-se o processo de calculo de zero
        // de funcao
        else
            return searchReverseProductionBottomHolePressure(state, pchute * 1.);
    }
}

namespace {

/// Brackets the root when the choke passes less than the column delivers.
bool bracketFromLowGuessSecondary(const SteadyStateSearchState &state, double amplifica, double &positiveResidualGuess, double &negativeResidualGuess, double &guessLowerBound, int &kontaiter, double mult2, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double chute, double &abortValue) {
    // isto implica em pressao de chute baixa , deve-se agora buscar
    // uma pressao de chute alta para que marchResidual seja negativo e assim iniciar o processo
    // de calculo de zero de funcao
    positiveResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor positivo
    // a cada nova busca em que marchResidual se aproxima de zero, mas ainda positivo
    // positiveResidualGuess e atualizado com a Ãºltima pressao de chute
    while (marchResidual > 0) {
        pchuteAux = pchute2;
        if ((*state.march.globals).chaverede == 0)
            pchute2 *= amplifica; // Aumentando a pressao na busca de marchResidual>0
        else
            pchute2 *= mult2; // Aumentando a pressao na busca de marchResidual>0
        if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
            pchute2 = 0.9 * state.march.input.tabent.pmax;
        int limpres = 0;
        if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
            (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
            if (chute < 0) {
                pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                    pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
                if (fabs(pchute2 - pchute) < 1e-15) {
                    pchute2 = -(1 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                    if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                        pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
                }
            } else {
                pchute2 = pchuteAux * 1.0001;
            }
            limpres = 1;
        }
        marchResidual = marchProductionSteadySecondary(state.march, pchute2);
        if (fabs(marchResidual) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (marchResidual > 0. && limpres == 1) {
            int iterpres = 0;
            while (marchResidual > 0 && iterpres < 10) {
                pchute2 *= 1.001;
                marchResidual = marchProductionSteadySecondary(state.march, pchute2);
                if (fabs(marchResidual) > 1.01e10) {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
                iterpres++;
            }
            if (iterpres >= 10) {
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm2 atingiu maximo de iteracoes");
                else {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            }
        }
        if (marchResidual > 0 && marchResidual < 0.9e10)
            positiveResidualGuess = pchute2; // atualizando o positiveResidualGuess
        kontaiter++;
        if (kontaiter > 50) { // limite de iteracoes, falha na busca do segundo chute
            // fim da simulacao ou aviso de falha
            if ((*state.march.globals).chaverede == 0) {
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm2 atingiu maximo de iteracoes");
                else {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            } else {
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        while (marchResidual > 0.9e10) { // o incremento de pressao foi demais e a marcha nao foi capaz
            // de ir ate o final sem que a pressao ficasse maior do que a pressao estatica
            // em um eventual IPR no meio da marcha ou maior que o limitre maximo de pressao
            // da tabela PVTSim, quando for este o caso
            // deve-se diminuir a estimativa de pressao baixa
            guessLowerBound = pchute2;
            pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a maior pressao
            // que leva a marchResidual>0 e a pressao alta demais
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
                pchute2 = 0.9 * state.march.input.tabent.pmax;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
                pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                    pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
            }
            marchResidual = marchProductionSteadySecondary(state.march, pchute2);
            if (fabs(marchResidual) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (marchResidual > 0 && marchResidual < 0.9e10)
                positiveResidualGuess = pchute2;
            kontaiter++;

            if (kontaiter > 50) {
                if ((*state.march.globals).chaverede == 0) {
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm2 atingiu maximo de iteracoes");
                    else {
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                        if ((*state.march.globals).iterRede > 0)
                            {
                                abortValue = -1.1e10;
                                return true;
                            }
                        else
                            {
                                abortValue = 1.1e10;
                                return true;
                            }
                    }
                } else {
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            }
        }
    }
    negativeResidualGuess = pchute2;
    return false;
}

/// Brackets the root when the pressure guess was too high.
bool bracketFromHighGuessSecondary(const SteadyStateSearchState &state, double reduz, double &negativeResidualGuess, double &guessLowerBound, int &kontaiter, double mult1, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double chute, double &abortValue) {
    // Vazao no choke>Vazao da mistura na tubulaÃ§Ã£o
    // calculada pela marcha, isto implica em pressao de chute alta , deve-se agora buscar
    // uma pressao de chute baixa para que marchResidual seja positivo e assim iniciar o processo
    // de calculo de zero de funcao
    negativeResidualGuess = pchute; // armazenando o valor de chute de pressao que da o valor negativo
    // a cada nova busca em que marchResidual se aproxima de zero, mas ainda negativo
    // negativeResidualGuess Ã© atualizado com a Ãºltima pressao de chute
    while (marchResidual < 0) {
        if (fabs(pchute2 - pchute) / pchute < (1. - reduz) / 10. && kontaiter > 50) {
            pchute2 *= 0.5;
        }
        pchuteAux = pchute2;
        if ((*state.march.globals).chaverede == 0)
            pchute2 *= reduz; // Diminuindo a pressao na busca de marchResidual>0
        else
            pchute2 *= mult1; // Diminuindo a pressao na busca de marchResidual>0
        if (pchute2 <= guessLowerBound)
            pchute2 = 0.5 * (pchuteAux + guessLowerBound); // guessLowerBound inicialmente
        // e zero, mas pode acontecer de baixar demais pchute2 ao ponto de marchResidual=-1e10
        //(pressao abaixo de 0.5 no meio da marcha), neste caso, guessLowerBound se torna este valor de
        // pchute2, pois, com isto, ja se sabe que nao se pode ir abaixo de guessLowerBound
        if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
            pchute2 = 0.9 * state.march.input.tabent.pmax;
        if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
            (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
            pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
            if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
        }
        marchResidual = marchProductionSteadySecondary(state.march, pchute2);
        if (fabs(marchResidual) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (marchResidual < 0 && marchResidual > -0.9e10)
            negativeResidualGuess = pchute2; // atualizando o negativeResidualGuess
        kontaiter++;
        if (kontaiter > 50) { // limite de iteracoes, falha na busca do segundo chute
            // fim da simulacao ou aviso de falha
            if ((*state.march.globals).chaverede == 0) {
                if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                    NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm2 atingiu maximo de iteracoees");
                else {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            } else {
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        while (marchResidual < -0.9e10) { // a reduÃ§Ã£o de pressao foi demais e a marcha nÃ£o foi capaz
            // de ir ate o final sem que a pressao ficasse inferior a 0.5kgf/cm2
            // deve-se aumentar a estimativa de pressao baixa
            guessLowerBound = pchute2;
            pchute2 = 0.5 * (pchute2 + pchuteAux); // valor intermediario entre a menor pressao
            // que leva a marchResidual<0 e a pressao baixa demais
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
                pchute2 = 0.9 * state.march.input.tabent.pmax;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
                pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                    pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
            }
            marchResidual = marchProductionSteadySecondary(state.march, pchute2);
            if (fabs(marchResidual) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (marchResidual < 0 && marchResidual > -0.9e10)
                negativeResidualGuess = pchute2;
            kontaiter++;
            if (kontaiter > 50) {
                if ((*state.march.globals).chaverede == 0) {
                    if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm2 atingiu maximo de iteracoes");
                    else {
                        cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                        if ((*state.march.globals).iterRede > 0)
                            {
                                abortValue = -1.1e10;
                                return true;
                            }
                        else
                            {
                                abortValue = 1.1e10;
                                return true;
                            }
                    }
                } else {
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            }
        }
    }
    return false;
}

/// Walks the guess until the march stops returning a sentinel.
bool retryUntilMarchCompletesSecondary(const SteadyStateSearchState &state, double pchuteAux0, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute, double &abortValue) {
    if (marchResidual < -0.9e10 || marchResidual > 0.9e10) {
        double valtemp;
        valtemp = marchProductionSteadySecondary(state.march, pchuteAux); // marcha com pchuteAux
        if (fabs(valtemp) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        while (valtemp > 0.9e10 && marchResidual > 0.9e10 && kontaiter < 50) { // estimativa de pressao de fundo ainda alta
            pchuteAux *= 0.99;                                       // reduzindo a estimativa
            valtemp = marchProductionSteadySecondary(state.march, pchuteAux);
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            kontaiter++; // 50 iteracoes no maximo
        }
        while (valtemp < -0.9e10 && marchResidual < -0.9e10 && kontaiter <= 50) { // estimativa de pressao de fundo ainda baixa
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection || state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection || state.march.cells[0].acsr.tipo == kAccessoryMultipleSource)
                pchuteAux *= 1.1;
            else
                pchuteAux *= 1.01; // aumentando a estimativa
            // verificando se este aumento ultrapassa o limite de pressao de uma eventual IPR no
            // fundo
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[0].acsr.ipr.Pres) {
                pchuteAux = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[0].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[0].acsr.ipr.Pres;
            }
            // verificando se este aumento de estimativa foca acima do valor maximo de pressao
            // de uma eventual tabela PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
            valtemp = marchProductionSteadySecondary(state.march, pchuteAux); // nova tentativa
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (kontaiter == 50) {
                pchuteAux = 0.5 * pchuteAux0;
            }
            kontaiter++;
        }
        valtemp = marchProductionSteadySecondary(state.march, pchuteAux);
        if (fabs(valtemp) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (valtemp > -0.9e10 && valtemp < 0.9e10) { // conseguiu fazer a marcha atÃƒÂ© o final
            marchResidual = valtemp;
            pchute = pchuteAux;
        }
    }
    return false;
}

/// Moves the guess according to which sentinel the march returned.
void classifyMarchSentinelSecondary(const SteadyStateSearchState &state, double marchResidual, double &pchuteAux, double perdafric, double &taux) {
    if (marchResidual < -0.9e10) { // pressao muito baixa na marcha, deve ser aumentado o valor de chute
        // faz-se uma nova estimativa, sÃ³ que agora admitindo uma hidrostÃ¡tica de Ã¡gua,
        // o que darÃ¡ uma pressao de fundpo mais alta
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = 1000 + 0 * state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            double alfa = 0.;
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = -(10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;  // limite de pchuteAux
    } else if (marchResidual > 0.9e10) { // pressao em algum ponto ficou acima de alguma pressao estatica
        // deve-se diminuir o valor do chute
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            // neste caso, utiliza-se uma fraÃ§Ã£o de vazio alta para a hidrostÃ¡tica
            double alfa = 0.8;
            if (state.march.annulusDrift == 0)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.2 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchuteAux) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                pchuteAux = -(10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchuteAux < 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchuteAux = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;
    }
}

/// Estimates the bottom-hole pressure the secondary search starts from.
void estimateInitialBottomHolePressureSecondary(const SteadyStateSearchState &state, double &completionFractionGuess, double &perdafric, double &frictionFactor, double &rmis, double &j, double &taux, double &pchute, double chute) {
    if (chute < 0) {
        if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection && fabs(state.march.cells[0].acsr.injl.QLiq) > 0.) {
            // este espaco faz uma estimativa de quanto deve ser a perda de carga media
            // a partir do valor da vazao no inicio da tubulacao, isto e feito apenas
            // se existir uma fonte de liquido, celula[0].acsr.tipo == 2
            taux = state.march.input.celp[0].textern;
            // a temperatura em que esta estimativa sera feita
            // e dada pela temperatura da fonte, a pressao para o calculo das propriedades
            // e feita por pGSup
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            ///////////propriedades fisicas:
            double completionFraction = state.march.cells[0].acsr.injl.bet;
            completionFractionGuess = completionFraction;
            double visC = state.march.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
            double visP = state.march.cells[0].acsr.injl.FluidoPro.ViscOleo(pchute, taux);
            double visG = state.march.cells[0].acsr.injl.FluidoPro.ViscGas(pchute, taux);
            double visMis = (1 - completionFraction) * visP + completionFraction * visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
            double liquidDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspLiq(pchute, taux);
            double gasDensityAtGuess = state.march.cells[0].acsr.injl.FluidoPro.MasEspGas(pchute, taux);
            rmis = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
            double rlpA = state.march.cells[0].acsr.injl.FluidoPro.MasEspLiq(1., 15.);
            double rlcA = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
            // estimativa grosseira da vazao massica do liquido complementar
            double massicC = rlcA * state.march.cells[0].acsr.injl.QLiq * state.march.cells[0].acsr.injl.bet / kSecondsPerDay;
            // estimativa grosseira da vazao massica de liquido e gas produzidos
            double massic = rlpA * state.march.cells[0].acsr.injl.QLiq * (1. - state.march.cells[0].acsr.injl.bet) / kSecondsPerDay;
            // massas especificas, condicao standard
            double Rhogs = state.march.cells[0].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
            double Rhols = (1000 * 141.5 / (131.5 + state.march.cells[0].acsr.injl.FluidoPro.API)) * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW) + 1000. * state.march.cells[0].acsr.injl.FluidoPro.Denag * state.march.cells[0].acsr.injl.FluidoPro.BSW;
            // multiplicador da vazao standard para se obter a Vazao massicagÃ¡s+liquido produzido
            double multiplicador = (Rhols + state.march.cells[0].acsr.injl.FluidoPro.RGO * Rhogs * (1 - state.march.cells[0].acsr.injl.FluidoPro.BSW));
            // massic *= multiplicador;//alteracao8
            if ((*state.march.globals).chaverede == 1)
                massic = 1 * multiplicador * state.march.cells[0].acsr.injl.QLiq * (1. - state.march.cells[0].acsr.injl.bet) / kSecondsPerDay;
            else
                massic = 1 * multiplicador * state.march.cells[0].acsr.injl.QLiq * (1. - state.march.cells[0].acsr.injl.bet) / kSecondsPerDay;
            // titulo do gas em relacao a misturaoleo+agua+gas
            double fracmasshidra = state.march.cells[0].acsr.injl.FluidoPro.FracMassHidra(pchute, taux);
            double massicP = (1. - fracmasshidra) * massic; // vazao massica de liquido produzido
            double massicG = fracmasshidra * massic;        // vazao massaica de gas
            /// estimativa grosseira da velocidade da mistura
            j = (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess) / state.march.cells[0].duto.area;
            // estimativa da fracao de vazio sem escorregamento
            double alfmis = (massicG / gasDensityAtGuess) / (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess);
            // propriedades de mistura
            rmis = (1 - alfmis) * rmis + alfmis * gasDensityAtGuess;
            visMis = (1 - alfmis) * visMis + alfmis * visG;
            // estimativa de numero de Reynolds
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            // estimativa de um fator de friccao
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);
        }
        // estimativa da perda por friccao no sistema, se o acessorio no inicio da tubulacao
        // nao for uma fonte de liquido, esta estimativa sera zero
        perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.march.cells[0].duto.peri / state.march.cells[0].duto.area;
        // se o acessorio for uma IPR, a estimativa serÃ¡ uma pressao prÃ³xima a uma vazao baixa
        // no fundo do poco
        // se nao for IPR, a estimativa continua, considerando a perda de carga e a hidrostatica,
        // que serÃ¡ avaliada neste laco
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            if (taux < state.march.input.tmin)
                taux = state.march.input.tmin;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchute, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchute, taux);
            // se o sistema nao for um anel principal de GL
            // para a estimativa de pressao, se utilizara so a hidrostatica de liquido, o que dara um
            // chute inicial de pressao muito alta
            double alfa = 0.;
            if (completionFractionGuess < 0.5) {
                double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
                alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
            }
            if (state.march.annulusDrift == 0)
                alfa = 1.; // se for o anel, so se tem gas
            else if ((*state.march.globals).chaverede == 0 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            // avanco da pressao por meio da hidrostatica e da perda por friccao estimada
            pchute += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            // caso exista alguma IPR no meio do duto, verifica se a pressao estimada esta maior que
            // a pressao estatica, nÃ£o se trabalha com vazoes negativas neste solver
            // para evitar isto, se corrige a pressao para um valor proximo da
            // pressao estatica da IPR da celula
            if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                ((state.march.cells[i - 1].acsr.ipr.Pres - pchute) < 0.01 * state.march.cells[i - 1].acsr.ipr.Pres && i == 1)) {
                double flowRateGuess = 0.15 * state.march.cells[i - 1].duto.area * kSecondsPerDay;
                pchute = -(flowRateGuess / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute > 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute < 0.9 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 0.9 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            // tambÃ©m se evita um chute com valores maiores que o limite maximo de uma tabela
            // no caso de se estar usando PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
        }

        if (pchute > 1000)
            pchute = 1000.; // pressao maxima de chute
    } else
        pchute = chute; // caso chute nao seja negativo utiliza a estimativa enviada na
}

}  // namespace

double searchProductionBottomHolePressureSecondary(const SteadyStateSearchState &state, double chute, int kontaTenta) {
    // busca de dois chutes iniciais com valores com sinais opostos
    // para marchaProdPerm1 e assim iniciar o prpocesso de calculo de erro de funcao.
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    double pchute = state.march.gasSurfacePressure; // inicializando o valor de pchute com o valor da pressao a jusante
    // do choke, pchute sera o valor de chute de fato no processo de busca
    // se o valor de chute>0, pchute=chute, senao, ele e estimado
    double taux; // valor de temperatura auxiliar para o eventual calculo
    // de pchute
    double j = 0.;
    double rmis = 0.;
    double frictionFactor = 0.;
    double perdafric = 0.;
    double completionFractionGuess = 0.;

    estimateInitialBottomHolePressureSecondary(state, completionFractionGuess, perdafric, frictionFactor, rmis, j, taux, pchute, chute);
    // lista de parametro do metodo
    double pchute2;        // segundo chute de pressao da busca
    double pchuteAux = 0.; // auxiliar na busca de dos chutes de pressao
    // necessarios para se dar partida no metodo de calculo de raiz
    double marchResidual; // marchResidual recebe sempre o valor retornado pelo procedimento de marcha
    // neste caso, Ã© a pressao a montante do choke-a pressao da ultima celula, calculada pela
    // marcha
    marchResidual = marchProductionSteadySecondary(state.march, pchute);
    state.march.searchOrigin = 1;
    // a marcha pode dar problemas, ou a pressao em algum ponto deu acima da pressao estatica de
    // de alguma IPR no meio do caminho, retorna 1e10,
    // ou antes de atingir a ultima celula, a pressao ficou proximo de zero, ou a pressao
    // retorna -1e10, ou a pressao, para o caso PVTSim, ficou acima da pressao maxima da tabela
    // retorna 1e10
    classifyMarchSentinelSecondary(state, marchResidual, pchuteAux, perdafric, taux);

    double mult1 = 0.9;
    double mult2 = 1.1;

    // para o caso de ter dado eerado a primeira marcha, com pchuteAux faz-se uma nova tentativa
    int kontaiter = 0; // contador para o laco em que se tentara uma nova estimativa
    // em que ao menos os valores 1e10 ou 1e-10 nÃ£o seja retornados
    double pchuteAux0 = pchuteAux;
    double abortValue;
    if (retryUntilMarchCompletesSecondary(state, pchuteAux0, kontaiter, marchResidual, pchuteAux, pchute, abortValue))
        return abortValue;
    if (kontaiter > 50) { // chegou ao limite da iteracao
        if ((*state.march.globals).chaverede == 0) {
            // neste caso, se finaliza a simulacao, nao tem um transiente
            // a ser feito a seguir e nem se estÃ¡ em uma rede
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
            else {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        } else {
            // se for em uma iteracao de rede, apos a primeira iteracao
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            // se logo apÃ³s tem uma simulacao transiente ou se esta na primeira iteracao de rede
            else
                return 1.1e10;
        }
    }
    // caso o laÃ§o anterior tenha levado a estimativa de um valor 1e10 para um valor -1e10
    // ou tenha levado para um valor de -1e10 para 1e10
    // ou seja, de um chute de pressao alto demais para finalizar a marcha
    // para um chute de pressao baixo demais para terminar a marcha ou vice versa
    while ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter <= 50) {
        pchute = 0.5 * (pchute + pchuteAux); // busca-se um valor medio
        // entre o alto demais e o baixo demais
        if ((fabs(pchute - pchuteAux) / pchute) < 0.001 && marchResidual < -0.9e10) {
            pchute = 1.05 * pchute;
            if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                (state.march.cells[0].acsr.ipr.Pres - pchute) < 0.01 * state.march.cells[0].acsr.ipr.Pres) {
                pchute = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                if (pchute < 0.99 * state.march.cells[0].acsr.ipr.Pres)
                    pchute = 0.99 * state.march.cells[0].acsr.ipr.Pres;
            }
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
            pchuteAux = pchute;
        }

        marchResidual = marchProductionSteadySecondary(state.march, pchute);
        if (fabs(marchResidual) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
        kontaiter++;
    }
    if (kontaiter > 50) {
        if ((*state.march.globals).chaverede == 0) {
            if (state.march.input.transiente == 0 && chute < 0 && state.march.input.AP == 0)
                NumError("Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes");
            else {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        } else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }

    // caso tenha sido possivel obter um chute de pressao em que foi possivel finalizar a marcha
    kontaiter = 0;
    pchute2 = pchute;
    double guessLowerBound = 0.;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    double reduz = 1. - state.march.input.buscaFC;
    double amplifica = 1. + state.march.input.buscaFC;
    if (fabs(marchResidual) < 1e-3)
        return pchute;
    else {
        if (marchResidual < 0.) { // caso em que o chute de pressao foi alto
            double abortValue;
            if (bracketFromHighGuessSecondary(state, reduz, negativeResidualGuess, guessLowerBound, kontaiter, mult1, marchResidual, pchuteAux, pchute2, pchute, chute, abortValue))
                return abortValue;
            positiveResidualGuess = pchute2;
        } else if (marchResidual > 0.) { // caso em que Vazao no choke<Vazao da mistura na tubulaÃ§Ã£o,
            double abortValue;
            if (bracketFromLowGuessSecondary(state, amplifica, positiveResidualGuess, negativeResidualGuess, guessLowerBound, kontaiter, mult2, marchResidual, pchuteAux, pchute2, pchute, chute, abortValue))
                return abortValue;
        }
        return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 1); // com as duas estimativas de pressao
        // em posicoes de sinal contrario da curva, inicia-se o processo de calculo de zero
        // de funcao
    }
}

namespace {

/// Walks the column cell by cell for the tertiary search.
///
/// The same shape as advanceProductionCells in the march module, and not the
/// same function.
bool advanceTertiaryCells(const SteadyStateSearchState &state, int &i, double &abortValue) {
    
                advanceUpstreamSteadyPressure(state.march, i, 0); // avanco da marcha para obter a pressao na fronteira esquerda
                // da celula i
                // teste para ver se ocorreu algum problema:
                if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - state.march.cells[i].presaux) < (*state.march.globals).localtiny)
                    {
                        abortValue = 1e10;
                        return true;
                    }
                if (state.march.cells[i].presaux <= 0.1 ||
                    (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmin - state.march.cells[i].presaux) > (*state.march.globals).localtiny)) {
                    {
                        abortValue = -1e10;
                        return true;
                    }
                }
                refreshUpstreamProductionPeriphery(state.march, i); // atualizacao da pressao da fronteira esquerda,
                // caso exista alguma BCS ou incremento de pressao
                if (state.march.input.flashCompleto != 2)
                    advanceSteadyMass(state.march, i); // verifica se existe alguma fonte na celula anterior, com isto, atualiza
                // as vazoes massica na fronteira a esquerda, alÃ©m das propriedades dos fluidos,
                // densidade do gas, RGO, API, BSW, beta
                else
                    advanceCompositionalSteadyMass(state.march, i);
                state.march.updaters.advanceSteadyTemperature(i, 0); // faz o avanco da temperatura, da celula i-1 para a celula i
                // verifica se teve algum problema nos limites de temperatura
                // caso se esteja trabalhando com tabela PVTSim
                if (isnan(state.march.cells[i].temp))
                    NumError("Temperatrura na linha de producao com valor NaN");
                if (state.march.input.usaTabela == 1 && (state.march.cells[i].temp - state.march.input.tabent.tmin) < (*state.march.globals).localtiny)
                    state.march.cells[i].temp = state.march.input.tabent.tmin;
                state.march.updaters.updateProductionTemperaturePeriphery(i); // mera atualizacao de atributos de temperatura a esquerda e a direita
                advanceDownstreamSteadyPressure(state.march, i, 0); // evolui a pressao  fronteira a esquerda da celula i para o
                // seu centro de celula
                // verifica se ocorreu algum problema nesta evolucao de de pressao no centro da
                // celula
                if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - state.march.cells[i].pres) < (*state.march.globals).localtiny) {
                    {
                        abortValue = 1e10;
                        return true;
                    }
                }
                if (state.march.cells[i].pres <= 0.1 ||
                    (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmin - state.march.cells[i].pres) > (*state.march.globals).localtiny)) {
                    {
                        abortValue = -1e10;
                        return true;
                    }
                }
                refreshDownstreamProductionPeriphery(state.march, i); // mera atualizacao de atributos que guardam valores de pressao
                // das celulas a esquerda e a direita
                for (int j = 0; j < state.march.input.nvalvgas; j++) { // reavaliacao da vazao da valvula de gas lift, quando
                    // a celula tem uma.
                    // P.S. parece uma acao desnecessÃ¡ria e talvez atÃ© um complicador
                    // densecessario, em vavliacao
                    if (state.march.productionValveCellIndices[j] == i) {
                        int k = state.march.gasValveCellIndices[j];
                        state.march.updaters.computeSteadyGasFlowRate(k);
                    }
                }
                if (state.march.input.tipoFluido == 0)
                    advanceSteadyMassTransfer(state.march, i - 1);
                else
                    advanceSteadyGasMassTransfer(state.march, i - 1); // caso seja uma tabela PVTSim, calcula-se a
                // taxa de transferÃªncia de massa entre as fases para o uso no calculo de
                // calor latente da equacao de energia
                if (state.march.input.ordperm > 1) { // correcao de segunda ordem
                    double D0presaux = state.march.cells[i].presaux - state.march.cells[i - 1].pres;
                    double D0pres = state.march.cells[i].pres - state.march.cells[i].presaux;
                    double D0temp = state.march.cells[i].temp - state.march.cells[i - 1].temp;
                    advanceUpstreamSteadyPressure(state.march, i, 1);
                    refreshUpstreamProductionPeriphery(state.march, i);
                    advanceSteadyMass(state.march, i);
                    state.march.updaters.advanceSteadyTemperature(i, 1);
                    advanceDownstreamSteadyPressure(state.march, i, 1);
                    state.march.cells[i].pres = 0.5 * (state.march.cells[i].presaux + D0pres + state.march.cells[i].pres);
                    state.march.cells[i].presaux = 0.5 * (state.march.cells[i - 1].pres + D0presaux + state.march.cells[i].presaux);
                    state.march.cells[i].temp = 0.5 * (state.march.cells[i - 1].temp + D0temp + state.march.cells[i].temp);
                    refreshDownstreamProductionPeriphery(state.march, i);
                    refreshUpstreamProductionPeriphery(state.march, i);
                    state.march.updaters.updateProductionTemperaturePeriphery(i);
                    advanceSteadyMass(state.march, i);
                    if (state.march.input.tipoFluido == 0)
                        advanceSteadyMassTransfer(state.march, i - 1);
                    else
                        advanceSteadyGasMassTransfer(state.march, i - 1);
                }
    
                // apÃ³s se atingir a pressao no centro da celula i, primeira iteracao de marcha
                // verifica-se se existe uma VGL em i e faz-se uma estimativa inicial da Vazao de
                // GL (caso exista linha de gas). Observar que isto sÃ³ Ã© feito para a iteracao zero.
                if (state.march.input.lingas > 0 && state.march.input.nvalvgas > 0 && state.march.steadyIteration == 0)
                    state.march.updaters.initializeSteadyValveGasFlowRate(i);
                i++;
    
                // teste para verificar se a pressao do centro de celula ficou acima
                // da pressao estatica de uma eventual IPR
                if (state.march.cells[i - 1].pres <= 0.1 ||
                    (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmin - state.march.cells[i - 1].pres) > (*state.march.globals).localtiny)) {
                    {
                        abortValue = -1e10;
                        return true;
                    }
                } else if ((state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance &&
                            ((state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres) < (*state.march.globals).localtiny) && i == 1)) {
                    {
                        abortValue = 1e10;
                        return true;
                    }
                } else if ((state.march.cells[i - 1].acsr.tipo == kAccessoryRadialPorous &&
                            ((state.march.cells[i - 1].acsr.radialPoro.pRes[0] - state.march.cells[i - 1].pres) < (*state.march.globals).localtiny) && i == 1)) {
                    {
                        abortValue = 1e10;
                        return true;
                    }
                } else if ((state.march.cells[i - 1].acsr.tipo == kAccessoryPorous2D &&
                            ((state.march.cells[i - 1].acsr.poroso2D.dados.pRes - state.march.cells[i - 1].pres) < (*state.march.globals).localtiny) && i == 1)) {
                    {
                        abortValue = 1e10;
                        return true;
                    }
                }
    return false;
}

/// Marches the column to convergence for the tertiary search.
///
/// Not a search: it calls no solver and no march, and inlines a march of its
/// own.
bool marchTertiaryCellsUntilConverged(const SteadyStateSearchState &state, int &guessNeedsCorrection, int &i, double betini, double alfini, double pentrada, double &abortValue) {
    while (guessNeedsCorrection == 1) { // opcao antiga, ja nao tem mais efeito
        // efetivamente, este while sempre so e feito uma vez, quando a marcha consegue ir ate
        // a ultima celula sem problemas, caso ocorra algum problema, a marcha e finalizada e
        // sai do metodo retornando ou 1e10 ou -1e10

        // inicializando as pressoes e fracoes volumetricas das celulas iniciais, centro de celula
        // e fronteira de celula
        state.march.cells[0].presauxL = pentrada;
        state.march.cells[0].presLini = pentrada;
        state.march.cells[0].presL = pentrada;
        state.march.cells[0].pres = pentrada;
        state.march.cells[1].presL = pentrada;
        state.march.cells[0].presini = pentrada;
        state.march.cells[1].presLini = pentrada;
        state.march.cells[0].presaux = pentrada;
        state.march.cells[1].presauxL = pentrada;

        state.march.cells[0].alf = alfini;
        state.march.cells[0].alfini = alfini;
        state.march.cells[0].bet = betini;
        state.march.cells[0].betini = betini;
        state.march.cells[1].alfL = state.march.cells[0].alf;
        state.march.cells[1].alfLini = state.march.cells[0].alf;
        state.march.cells[0].alfPigD = state.march.cells[0].alf;
        state.march.cells[0].alfPigDini = state.march.cells[0].alf;
        state.march.cells[0].alfPigE = state.march.cells[0].alf;
        state.march.cells[0].alfPigEini = state.march.cells[0].alf;
        state.march.cells[1].betL = state.march.cells[0].bet;
        state.march.cells[1].betLini = state.march.cells[0].bet;
        state.march.cells[0].betPigD = state.march.cells[0].bet;
        state.march.cells[0].betPigDini = state.march.cells[0].bet;
        state.march.cells[0].betPigE = state.march.cells[0].bet;
        state.march.cells[0].betPigEini = state.march.cells[0].bet;
        state.march.cells[0].betI = state.march.cells[0].bet;
        state.march.cells[1].betLI = state.march.cells[0].bet;

        // verifica se ja existe algum problema no inicio da marcha
        if (state.march.cells[0].pres <= 0.1 ||
            (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmin - state.march.cells[0].pres) > (*state.march.globals).localtiny))
            {
                abortValue = -1e10;
                return true;
            }
        else if ((state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                  (state.march.cells[0].acsr.ipr.Pres - state.march.cells[0].pres) < (*state.march.globals).localtiny))
            {
                abortValue = 1e10;
                return true;
            }
        else if ((state.march.cells[0].acsr.tipo == kAccessoryRadialPorous &&
                  (state.march.cells[0].acsr.radialPoro.pRes[0] - state.march.cells[0].pres) < (*state.march.globals).localtiny))
            {
                abortValue = 1e10;
                return true;
            }
        else if ((state.march.cells[0].acsr.tipo == kAccessoryPorous2D &&
                  (state.march.cells[0].acsr.poroso2D.dados.pRes - state.march.cells[0].pres) < (*state.march.globals).localtiny))
            {
                abortValue = 1e10;
                return true;
            }
        // IniciaVazValvGasPerm e um metodo que faz uma estimativa da vazao na valvula de GL
        // quando ainda nao foi feita a marcha na linha de gas. Neste caso, ele recebe o
        // indice da celula de producao e verifica se nesta celula existe uma VGL, se existir,
        // caso a condicao na linha de gas seja vazao injetada, divide a vazao injetada pelo numero de
        // valvulas e indica este valor para a VGL relacionada a celula de producao
        // caso a condicao seja pressao de injecao, faz-se uma estimativa da pressao na linha de gas
        // na posicao da VGL por hidrotatica e com isto se calcula a vazao de injecao da VGL
        if (state.march.input.lingas > 0 && state.march.input.nvalvgas > 0 && state.march.steadyIteration == 0)
            state.march.updaters.initializeSteadyValveGasFlowRate(0);
        i = 1;
        // inicio da marcha propriamente dita
        while (i <= state.march.lastCell && state.march.cells[i - 1].pres >= 0.1 && fabs(pentrada - state.march.cells[0].pres) < (*state.march.globals).localtiny) {
            if (advanceTertiaryCells(state, i, abortValue))
                return true;
        }
        if (i == state.march.lastCell + 1)
            guessNeedsCorrection = 0; // fim da marcha
    }
    return false;
}


/// The marches per guess with convergence acceleration on: one more under a
/// pressure condition on the gas line, where it couples less easily with the
/// production line.
int acceleratedMarchesPerGuess(InjectionPressureCondition) {
    return 3;
}

int acceleratedMarchesPerGuess(InjectionFlowRateCondition) {
    return 2;
}

/// The gas line at steady state under its inlet condition, as the tertiary
/// search solves it: marched when the injection choke is not throttling and
/// searched (tertiary) when it is, under a pressure condition; searched
/// (secondary) under a flow-rate condition.
void solveGasLineForSearch(const SteadyStateSearchState &state, InjectionPressureCondition) {
    // marcha para o caso, pressao de injecao
    if (state.march.input.chokes.abertura[0] >= 0.2) { // choke de injecao inativo
        for (int iter = 0; iter < 1; iter++) {
            marchGasSteady(state.march);
        }
    } else
        searchGasPressureSteadyTertiary(state); // choke de injecao ativo
}

void solveGasLineForSearch(const SteadyStateSearchState &state, InjectionFlowRateCondition) {
    searchGasPressureSteadySecondary(state); // marcha na linha de gas para o caso de vazao de injecao
}

}  // namespace

double searchProductionBottomHolePressureTertiary(const SteadyStateSearchState &state, double pentrada) {
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    double alfini;
    double betini;

    // estimativa da fracao de vazio na primeira celula do sistema
    // fonte de massa no inicio da tubulacao
    // da mesma maneira que no caso de fonte de liquido, fracao de vazio= sem escorregamento
    state.march.cells[0].temp = state.march.cells[0].acsr.injm.temp;
    if (state.march.input.flashCompleto == 2) {
        if (state.march.input.tabelaDinamica == 0)
            state.march.cells[0].flui.atualizaPropComp(pentrada, state.march.cells[0].temp);
        state.march.cells[0].acsr.injm.FluidoPro.atualizaPropComp(pentrada, state.march.cells[0].temp);
    }
    state.march.cells[0].acsr.injm.condTermo = 0;
    state.march.cells[0].pres = pentrada;
    state.march.updaters.updateSource(0);
    double qgas = state.march.cells[0].acsr.injm.MassG /
                  state.march.cells[0].acsr.injm.FluidoPro.MasEspGas(pentrada, state.march.cells[0].temp);
    double qliq = state.march.cells[0].acsr.injm.MassP /
                      state.march.cells[0].acsr.injm.FluidoPro.MasEspLiq(pentrada, state.march.cells[0].temp) +
                  state.march.cells[0].acsr.injm.MassC /
                      state.march.cells[0].acsr.injm.fluidocol.MasEspFlu(pentrada, state.march.cells[0].temp);
    state.march.cells[0].acsr.injm.FluidoPro.MasEspLiq(pentrada, state.march.cells[0].temp);
    double qcomp = state.march.cells[0].acsr.injm.MassC /
                   state.march.cells[0].acsr.injm.fluidocol.MasEspFlu(pentrada, state.march.cells[0].temp);
    alfini = qgas / (qliq + qgas);
    if ((fabs(qcomp) + fabs(qliq)) > 1e-15)
        betini = fabs(qcomp) / (fabs(qcomp) + fabs(qliq));
    else
        betini = 0.;

    // esta marcha e feita para quando se tem alguma fonte no inicio da tubulacao,
    // portanto, admite-se que o duto esta fechado e coloca-se uma fonte no centro da
    // primeira celula. As vazoes na fronteira esquerda da celula sÃ£o portanto = 0
    state.march.cells[0].tempL = state.march.cells[0].temp;
    state.march.cells[1].tempL = state.march.cells[0].temp;
    state.march.cells[0].tempini = state.march.cells[0].temp;

    state.march.cells[0].ML = 0.;
    state.march.cells[0].MC = 0.;
    state.march.cells[1].ML = 0.;
    state.march.cells[0].MliqiniL = 0.;
    state.march.cells[0].Mliqini = 0.;
    state.march.cells[1].MliqiniL = 0.;
    state.march.cells[0].MComp = 0.;
    state.march.cells[0].QLL = 0.;
    state.march.cells[0].QL = 0;
    state.march.cells[1].QLL = 0.;
    state.march.cells[0].QG = 0.;
    double masfim = 1.;
    double masfim0 = -10;
    double presteste = 1.;
    double presteste0 = -10;
    state.march.steadyIteration = 0;

    int limIter = 200; // limite de iteracoes quando a opcao de aceleracao da convergencia esta desligado. desaconselhavel
    // desligar esta opcao, os ganhos sao pouco e a convergencia se torna instavel, principalmente
    // quando o tramo faz perte de um sistema de redes
    if (state.march.input.AceleraConvergPerm == 1) { // opcao aceleracao de convergencia ligada
        limIter = 2;                   // em geral faz-se apenas duas iteracoes de marcha para um determinado chute
        if (state.march.input.lingas == 1)
            limIter = withGasInletCondition(
                state.march.input.gasinj.tipoCC,
                [](auto condition) { return acceleratedMarchesPerGuess(condition); }); // no caso de se ter
        // uma condicao de contorno na injecao de gas = pressao, observou-se que o acoplamento dinamico
        // entre a linha de gas e de producao e mais difoicil, para se conseguir um sistema
        // melhor acoplado, deve-se fazer uma marcha iterativa a mais
    }
    // arq.CriterioConvergPerm Ã© um criterio de convergencia da marcha, so faz sentido
    // quando a aceleracao de convergencia esta desligada. Observe que esta convergencia nÃ£o
    // e de fato a conevregencia do problema, e apenas um repeticao de marcha para um determinado
    // chute de pressao ou de vazao no inicio da tubulacao. O que a convergencia busca de fato e
    // determinar qual a pressao ou vazao de fundo que satisfaz as condicoes de contorno no fim
    // da tubulacao, esta busca e feita nos metodos de busca.
    while ((fabs(masfim - masfim0) / fabs(masfim) > state.march.input.CriterioConvergPerm ||
            fabs(presteste - presteste0) / fabs(presteste) > state.march.input.CriterioConvergPerm) &&
           state.march.steadyIteration < limIter) {

        masfim0 = masfim;
        presteste0 = presteste;
        if (state.march.steadyIteration == 0 && state.march.input.lingas > 0 && state.march.input.nvalvgas > 0 && state.march.networkCoupled == 1)
            state.march.updaters.initializeTubingConnectionSteady(); // antes de iniciar a primeira iteracao de marcha,
        // faz-se uma estimativa inicial de como se da o acopamento termico entre a coluna e o anular,
        // caso se tenha linha de gas
        else if (state.march.input.lingas > 0 && state.march.input.nvalvgas > 0 && state.march.networkCoupled == 1)
            state.march.updaters.connectTubingSteady(); // acoplamento termico feito com valores de pressao e temperatura
        // obtidas na primeira iteracao de marcha
        int i;
        int guessNeedsCorrection = 1;
        double abortValue;
        if (marchTertiaryCellsUntilConverged(state, guessNeedsCorrection, i, betini, alfini, pentrada, abortValue))
            return abortValue;
        // apÃ³s o fim da marcha da linha de produÃ§Ã£o, Ã© feita a marcha da linha de gas
        // caso exista
        if (pentrada > 0 && state.march.input.lingas > 0 && state.march.input.nvalvgas > 0) {
            withGasInletCondition(state.march.gasCells[0].tipoCC,
                                  [&](auto condition) { solveGasLineForSearch(state, condition); });
        }

        if (state.march.input.acopColAnulPermForte > 0 && state.march.input.lingas > 0) {
            refreshProperties(state.march);
            refreshSteadyThermalVelocities(state.march);
            computePseudoTransientTimeStep(state.march);
            state.march.updaters.connectTubing();
            for (int kontaPseudo = 0; kontaPseudo < state.march.input.acopColAnulPermForte; kontaPseudo++) {
                for (int i = 0; i <= state.march.gasCellCount; i++)
                    state.march.gasCells[i].tempini = state.march.gasCells[i].temp;
                for (int i = 1; i <= state.march.gasCellCount; i++) {
                    state.march.updaters.computeGasTemperature(i, state.march.gasCells[i - 1].tempini, 1);
                }
                state.march.cells[0].tempini = state.march.cells[0].temp;
                state.march.cells[1].tempLini = state.march.cells[1].tempL;
                state.march.cells[1].tempL = state.march.cells[0].temp;
                for (int i = 1; i <= state.march.lastCell; i++) {
                    state.march.cells[i].tempini = state.march.cells[i].temp;
                }
                for (int i = 1; i <= state.march.lastCell; i++) {
                    state.march.updaters.computeTemperature(i, state.march.cells[i].tempini, 1);
                }
                refreshProperties(state.march);
                computePseudoTransientTimeStep(state.march);
                state.march.updaters.connectTubing();
            }
        }

        masfim = state.march.cells[state.march.lastCell - 1].MC; // guarda valor de vazao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        presteste = state.march.cells[state.march.lastCell].pres; // guarda valor de pressao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        state.march.steadyIteration++; // atualiza a ieteracao da marcha
        if (state.march.steadyIteration > 200 && state.march.input.AP == 0)
            NumError("ConvergÃƒÂªncia em marchaProdPerm1 atingiu maximo de iteracoes");
        else if (state.march.steadyIteration > 200)
            return 1.1e10;
    }

    double corrigePresF = 0.;
    if (((*state.march.globals).chaverede == 0 || state.march.endNode == 1 || (*state.march.globals).chaveRedeParalela == 1) && state.march.input.corrigeContSep == 1)
        corrigePresF = steadyPressureAtLastCell(state.march);

    state.march.gasSurfacePressure = (state.march.cells[state.march.lastCell].pres + corrigePresF);

    return state.march.cells[state.march.lastCell].pres; // caso a marcha tenha conseguido ir atÃ© a Ãºltima celula,
    // retorna a pressao da ultima celula calculada pela marcha
}

double searchReverseProductionPressureToPressure(const SteadyStateSearchState &state, double chute, double maximumFlowRate, int kontaiter) {
    state.reverseSteady = 1;
    state.march.convergenceMonitor = 1000.;
    if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
        state.march.cells[0].acsr.injg.FluidoPro = state.reverseNetworkFluid;
    } else if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
        state.march.cells[0].acsr.injl.FluidoPro = state.reverseNetworkFluid;
    } else if (state.march.cells[0].acsr.tipo == kAccessoryMultipleSource) {
        state.march.cells[0].acsr.injm.FluidoPro = state.reverseNetworkFluid;
    }
    double mchute = chute;

    double marchResidual = marchReverseProductionPressureToPressure(state.march, mchute);
    state.march.searchOrigin = 1;
    if (marchResidual < -0.9e10) {
        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
            NumError(
                "Busca de valores iniciais para callculo de zero de funcao em buscaProdPresPresPermRev com problemas, pressao abaixo de 1 kgf/cm2 na marcha reversa");
        else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    double mchuteAux;
    double mchute2 = mchute;
    state.march.input.buscaFC = fabs(state.march.input.buscaFC);

    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    int testaEscoa = 1;
    double val0 = marchResidual;
    int reversao = 0;
    if (fabs(marchResidual) < 1e-3)
        return mchute;
    else {
        if (marchResidual < 0.) {
            negativeResidualGuess = mchute;
            while (marchResidual < 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. + state.march.input.buscaFC);
                marchResidual = marchReverseProductionPressureToPressure(state.march, mchute2);
                if (marchResidual > -0.9e10) {
                    if (marchResidual < val0) {
                        reversao = 1;
                        marchResidual = 1;
                    } else
                        val0 = marchResidual;
                }
                if (marchResidual < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = mchute2;
                kontaiter++;
                double velocityGuess = mchute2;
                double completionFraction = state.march.cells[1].bet;
                double voidFraction = state.march.cells[1].alf;
                double firstCellPressure = state.march.cells[0].pres;
                double firstCellTemperature = state.march.cells[0].temp;
                double completionDensityAtGuess = state.march.cells[0].fluicol.MasEspFlu(firstCellPressure, firstCellTemperature);
                double liquidDensityAtGuess = state.march.cells[0].flui.MasEspLiq(firstCellPressure, firstCellTemperature);
                double gasDensityAtGuess = state.march.cells[0].flui.MasEspGas(firstCellPressure, firstCellTemperature);
                double rmisL = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
                double rmis = (1 - voidFraction) * rmisL + voidFraction * gasDensityAtGuess;
                double rCst = state.march.cells[0].fluicol.MasEspFlu(1., 20.);
                double rPst = state.march.cells[0].flui.MasEspLiq(1., 20.);
                double rGst = state.march.cells[0].flui.MasEspGas(1., 20.);
                double rmisLst = (1 - completionFraction) * rPst + completionFraction * rCst;
                double multiplica;
                if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection)
                    multiplica = rmisLst;
                else
                    multiplica = rGst;
                velocityGuess = mchute2 * multiplica / (rmis * state.march.cells[0].duto.area * kSecondsPerDay);
                if (kontaiter > 200 && fabs(velocityGuess) > 0.01) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPermRev atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (fabs(velocityGuess) <= 0.01) {
                    testaEscoa = 0;
                    marchResidual = 1.;
                }
            }
            positiveResidualGuess = mchute2;
        } else if (marchResidual > 0.) {
            positiveResidualGuess = mchute;
            while (marchResidual > 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. - state.march.input.buscaFC);
                marchResidual = marchReverseProductionPressureToPressure(state.march, mchute2);
                if (marchResidual > -0.9e10) {
                    if (marchResidual > val0) {
                        reversao = 1;
                        marchResidual = -1;
                    } else
                        val0 = marchResidual;
                }
                if (marchResidual > 0 && marchResidual > -0.9e10)
                    positiveResidualGuess = mchute2;
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm  atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
            negativeResidualGuess = mchute2;
        }
        if (testaEscoa == 1 && reversao == 0)
            return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 2, 0);
        else if (reversao == 0) {
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
                state.march.cells[0].acsr.injg.QGas = 0.;
            } else if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
                state.march.cells[0].acsr.injl.QLiq = 0.;
            }
            return 0.;
        } else {
            return searchProductionPressureToPressure(state, -chute, -maximumFlowRate, kontaiter);
        }
    }
}

double searchProductionPressureToPressure(const SteadyStateSearchState &state, double chute, double maximumFlowRate, int kontaiter) {
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    double mchute = chute;
    double marchResidual = marchProductionPressureToPressure(state.march, mchute);
    state.march.searchOrigin = 1;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10 && kontaiter < 800) {
            mchute *= 0.9;
            marchResidual = marchProductionPressureToPressure(state.march, mchute);
            if (mchute < 1e-5) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
        marchResidual = marchProductionPressureToPressure(state.march, mchute);
        kontaiter++;
    }
    double mchuteAux;
    double mchute2 = mchute;
    double guessLowerBound = 0.;

    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    int testaEscoa = 1;
    double sentido = 1.;
    double val0 = marchResidual;
    int reversao = 0;
    int kontaReverso;
    if (fabs(marchResidual) < 1e-3)
        return mchute;
    else {
        if (marchResidual * sentido < 0.) {
            negativeResidualGuess = mchute;
            kontaReverso = 0;
            while (marchResidual * sentido < 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. - state.march.input.buscaFC);
                marchResidual = marchProductionPressureToPressure(state.march, mchute2);
                if (marchResidual > -0.9e10) {
                    if (marchResidual < val0) {
                        kontaReverso++;
                        if (kontaReverso == 2) {
                            reversao = 1;
                            marchResidual = 1;
                        } else {
                            state.march.input.buscaFC = -fabs(state.march.input.buscaFC);
                            val0 = marchResidual;
                        }
                    } else
                        val0 = marchResidual;
                }
                if (marchResidual * sentido < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = mchute2;
                kontaiter++;
                double velocityGuess = mchute2;
                double completionFraction = state.march.cells[1].bet;
                double voidFraction = state.march.cells[1].alf;
                double firstCellPressure = state.march.cells[0].pres;
                double firstCellTemperature = state.march.cells[0].temp;
                double completionDensityAtGuess = state.march.cells[0].fluicol.MasEspFlu(firstCellPressure, firstCellTemperature);
                double liquidDensityAtGuess = state.march.cells[0].flui.MasEspLiq(firstCellPressure, firstCellTemperature);
                double gasDensityAtGuess = state.march.cells[0].flui.MasEspGas(firstCellPressure, firstCellTemperature);
                double rmisL = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
                double rmis = (1 - voidFraction) * rmisL + voidFraction * gasDensityAtGuess;
                double rCst = state.march.cells[0].fluicol.MasEspFlu(1., 20.);
                double rPst = state.march.cells[0].flui.MasEspLiq(1., 20.);
                double rGst = state.march.cells[0].flui.MasEspGas(1., 20.);
                double rmisLst = (1 - completionFraction) * rPst + completionFraction * rCst;
                double multiplica;
                if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection)
                    multiplica = rmisLst;
                else
                    multiplica = rGst;
                velocityGuess = mchute2 * multiplica / (rmis * state.march.cells[0].duto.area * kSecondsPerDay);
                if (kontaiter > 200 && velocityGuess > 0.01) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (velocityGuess <= 0.01) {
                    testaEscoa = 0;
                    marchResidual = 1.;
                }
                while (marchResidual < -0.9e10) {
                    guessLowerBound = mchute2;
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchProductionPressureToPressure(state.march, mchute2);
                    if (marchResidual * sentido < 0 && marchResidual > -0.9e10)
                        negativeResidualGuess = mchute2;
                    kontaiter++;
                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            positiveResidualGuess = mchute2;
        } else if (marchResidual * sentido > 0.) {
            positiveResidualGuess = mchute;
            kontaReverso = 0;
            while (marchResidual * sentido > 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. + state.march.input.buscaFC);
                marchResidual = marchProductionPressureToPressure(state.march, mchute2);
                if (marchResidual > -0.9e10) {
                    if (marchResidual > val0) {
                        kontaReverso++;
                        if (kontaReverso == 2) {
                            reversao = 1;
                            marchResidual = -1;
                        } else {
                            state.march.input.buscaFC = -fabs(state.march.input.buscaFC);
                            val0 = marchResidual;
                        }
                    } else
                        val0 = marchResidual;
                }
                if ((*state.march.globals).chaverede != 0 && mchute2 > maximumFlowRate && maximumFlowRate > 0) {
                    testaEscoa = 0;
                    marchResidual = -1.;
                }
                if (marchResidual * sentido > 0 && marchResidual > -0.9e10)
                    positiveResidualGuess = mchute2;
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm  atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                while (marchResidual < -0.9e10) {
                    guessLowerBound = mchute2;
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchProductionPressureToPressure(state.march, mchute2);
                    if (marchResidual * sentido > 0 && marchResidual > -0.9e10)
                        positiveResidualGuess = mchute2;
                    kontaiter++;

                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            negativeResidualGuess = mchute2;
        }
        if (testaEscoa == 1 && reversao == 0)
            return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 2, 0);
        else if (reversao == 0) {
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
                state.march.cells[0].acsr.injg.QGas = -2121212121;
            } else if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
                state.march.cells[0].acsr.injl.QLiq = -2121212121;
            }
            return 0.;
        } else {
            return searchReverseProductionPressureToPressure(state, -chute, -maximumFlowRate, kontaiter);
        }
    }
}

double searchProductionPressureToPressureSecondary(const SteadyStateSearchState &state, double chute, double maximumFlowRate) {
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    double mchute = chute;
    double marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
    state.march.searchOrigin = 1;
    int kontaiter = 0;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10 && kontaiter < 400) {
            mchute *= 0.9;
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
            if (mchute < 1e-5) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
            kontaiter++;
        }
    }
    double mchuteAux;
    double mchute2 = mchute;
    double guessLowerBound = 0.;
    kontaiter = 0;

    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;

    int testaEscoa = 1;
    if (fabs(marchResidual) < 1e-3)
        return mchute;
    else {
        if (marchResidual < 0.) {
            negativeResidualGuess = mchute;
            while (marchResidual < 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. - state.march.input.buscaFC);
                marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                if (marchResidual < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = mchute2;
                kontaiter++;
                double velocityGuess = mchute2;
                double completionFraction = state.march.cells[1].bet;
                double voidFraction = state.march.cells[1].alf;
                double firstCellPressure = state.march.cells[0].pres;
                double firstCellTemperature = state.march.cells[0].temp;
                double completionDensityAtGuess = state.march.cells[0].fluicol.MasEspFlu(firstCellPressure, firstCellTemperature);
                double liquidDensityAtGuess = state.march.cells[0].flui.MasEspLiq(firstCellPressure, firstCellTemperature);
                double gasDensityAtGuess = state.march.cells[0].flui.MasEspGas(firstCellPressure, firstCellTemperature);
                double rmisL = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
                double rmis = (1 - voidFraction) * rmisL + voidFraction * gasDensityAtGuess;
                double rCst = state.march.cells[0].fluicol.MasEspFlu(1., 20.);
                double rPst = state.march.cells[0].flui.MasEspLiq(1., 20.);
                double rGst = state.march.cells[0].flui.MasEspGas(1., 20.);
                double rmisLst = (1 - completionFraction) * rPst + completionFraction * rCst;
                double multiplica;
                if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection)
                    multiplica = rmisLst;
                else
                    multiplica = rGst;
                velocityGuess = mchute2 * multiplica / (rmis * state.march.cells[0].duto.area * kSecondsPerDay);
                if (kontaiter > 200 && velocityGuess > 0.01) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (velocityGuess <= 0.01) {
                    testaEscoa = 0;
                    marchResidual = 1.;
                }
                while (marchResidual < -0.9e10) {
                    guessLowerBound = mchute2;
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                    if (marchResidual < 0 && marchResidual > -0.9e10)
                        negativeResidualGuess = mchute2;
                    kontaiter++;
                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            positiveResidualGuess = mchute2;
        } else if (marchResidual > 0.) {
            positiveResidualGuess = mchute;
            while (marchResidual > 0) {
                mchuteAux = mchute2;
                mchute2 *= (1. + state.march.input.buscaFC);
                marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                if ((*state.march.globals).chaverede != 0 && mchute2 > maximumFlowRate) {
                    testaEscoa = 0;
                    marchResidual = -1.;
                }
                if (marchResidual > 0 && marchResidual < 0.9e10)
                    positiveResidualGuess = mchute2;
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm  atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                while (marchResidual < -0.9e10) {
                    guessLowerBound = mchute2;
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                    if (marchResidual > 0 && marchResidual > -0.9e10)
                        positiveResidualGuess = mchute2;
                    kontaiter++;

                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            negativeResidualGuess = mchute2;
        }
        if (testaEscoa == 1)
            return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 2, 1);
        else {
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
                state.march.cells[0].acsr.injg.QGas = 0.;
            } else if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
                state.march.cells[0].acsr.injl.QLiq = 0.;
            }
            return 0.;
        }
    }
}

namespace {

/// Brackets and solves the mass-flow root for the tertiary pressure-to-pressure search.
double bracketTertiaryPressureToPressureRoot(const SteadyStateSearchState &state, int &testaEscoa, double &positiveResidualGuess, double &negativeResidualGuess, double &guessLowerBound, double &mchute2, double &mchuteAux, int &kontaiter, double &marchResidual, double mchute, double maximumFlowRate) {
    if (marchResidual < 0.) {
        negativeResidualGuess = mchute;
        while (marchResidual < 0) {
            mchuteAux = mchute2;
            mchute2 *= (1. - state.march.input.buscaFC);
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
            if (marchResidual < 0 && marchResidual > -0.9e10)
                negativeResidualGuess = mchute2;
            kontaiter++;
            double velocityGuess = mchute2;
            double completionFraction = state.march.cells[1].bet;
            double voidFraction = state.march.cells[1].alf;
            double firstCellPressure = state.march.cells[0].pres;
            double firstCellTemperature = state.march.cells[0].temp;
            double completionDensityAtGuess = state.march.cells[0].fluicol.MasEspFlu(firstCellPressure, firstCellTemperature);
            double liquidDensityAtGuess = state.march.cells[0].flui.MasEspLiq(firstCellPressure, firstCellTemperature);
            double gasDensityAtGuess = state.march.cells[0].flui.MasEspGas(firstCellPressure, firstCellTemperature);
            double rmisL = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
            double rmis = (1 - voidFraction) * rmisL + voidFraction * gasDensityAtGuess;
            double rCst = state.march.cells[0].fluicol.MasEspFlu(1., 20.);
            double rPst = state.march.cells[0].flui.MasEspLiq(1., 20.);
            double rGst = state.march.cells[0].flui.MasEspGas(1., 20.);
            double rmisLst = (1 - completionFraction) * rPst + completionFraction * rCst;
            double multiplica;
            if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection)
                multiplica = rmisLst;
            else
                multiplica = rGst;
            velocityGuess = mchute2 * multiplica / (rmis * state.march.cells[0].duto.area * kSecondsPerDay);
            if (kontaiter > 200 && velocityGuess > 0.01) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            if (velocityGuess <= 0.01) {
                testaEscoa = 0;
                marchResidual = 1.;
            }
            while (marchResidual < -0.9e10) {
                guessLowerBound = mchute2;
                mchute2 = 0.5 * (mchute2 + mchuteAux);
                marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                if (marchResidual < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = mchute2;
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        positiveResidualGuess = mchute2;
    } else if (marchResidual > 0.) {
        positiveResidualGuess = mchute;
        while (marchResidual > 0) {
            mchuteAux = mchute2;
            mchute2 *= (1. + state.march.input.buscaFC);
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
            if ((*state.march.globals).chaverede != 0 && mchute2 > maximumFlowRate && maximumFlowRate > 0) {
                testaEscoa = 0;
                marchResidual = -1.;
            }
            if (marchResidual > 0 && marchResidual < 0.9e10)
                positiveResidualGuess = mchute2;
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm  atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            while (marchResidual > 0.9e10) {
                guessLowerBound = mchute2;
                mchute2 = 0.5 * (mchute2 + mchuteAux);
                marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute2);
                if (marchResidual > 0 && marchResidual < 0.9e10)
                    positiveResidualGuess = mchute2;
                kontaiter++;

                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        negativeResidualGuess = mchute2;
    }
    if (testaEscoa == 1)
        return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 2, 1);
    else {
        if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
            state.march.cells[0].acsr.injg.QGas = 0.;
        } else if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
            state.march.cells[0].acsr.injl.QLiq = 0.;
        }
        return 0.;
    }
}

}  // namespace

double searchProductionPressureToPressureTertiary(const SteadyStateSearchState &state, double chute, double maximumFlowRate) {

    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;
    double taux;
    double pchute = state.march.input.CCPres.pres[0];
    double totalSourceFlowRate = 0.;
    if (chute < 0.) {
        for (int i = 0; i < state.march.lastCell; i++) {
            taux = state.march.input.celp[i].textern;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchute, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchute, taux);
            // se o sistema nao for um anel principal de GL
            // para a estimativa de pressao, se utilizara so a hidrostatica de liquido, o que dara um
            // chute inicial de pressao muito alta
            double alfa = 0.1;
            if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i + 1].dx);
            // avanco da pressao por meio da hidrostatica e da perda por friccao estimada
            pchute -= ((rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed) / kPascalPerKgfPerCm2);
            // incremento de pressao devido a algum ganho de pressao constante em alguma celula
            if (state.march.cells[i].acsr.tipo == 7)
                pchute -= state.march.cells[i].acsr.delp;
            // caso exista alguma IPR no meio do duto, verifica se a pressao estimada esta menor que
            // a pressao estatica, nao se trabalha com vazoes negativas neste solver
            // para evitar isto, se corrige a pressao para um valor proximo da
            // pressao estatica da IPR da celula
            if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance && ((state.march.cells[i].acsr.ipr.Pres - pchute) > 0.0) && i == state.march.lastCell - 1) {

                pchute = -(10 / state.march.cells[i - 1].acsr.ipr.ip) + state.march.cells[i - 1].acsr.ipr.Pres;
                if (pchute < 0.99 * state.march.cells[i - 1].acsr.ipr.Pres)
                    pchute = 0.99 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
            double accessoryFlowRate = 0.;
            if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance)
                accessoryFlowRate = -(state.march.cells[i].acsr.ipr.Pres - pchute) * state.march.cells[i].acsr.ipr.ip;
            else if (state.march.cells[i].acsr.tipo == kAccessoryGasInjection && state.march.input.tipoFluido == 1)
                accessoryFlowRate = -state.march.cells[i].acsr.injg.QGas;
            else if (state.march.cells[i].acsr.tipo == kAccessoryLiquidInjection && state.march.input.tipoFluido == 0)
                accessoryFlowRate = -state.march.cells[i].acsr.injl.QLiq;
            else if (state.march.cells[i].acsr.tipo == kAccessoryLeak && state.march.cells[0].acsr.tipo == kAccessoryGasInjection) {
                accessoryFlowRate = 10000.;
            } else if (state.march.cells[i].acsr.tipo == kAccessoryLeak && state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection) {
                accessoryFlowRate = 1000.;
            }

            totalSourceFlowRate += accessoryFlowRate;
            if (totalSourceFlowRate < 0)
                totalSourceFlowRate -= accessoryFlowRate;
            // tambÃ©m se evita um chute com valores maiores que o limite maximo de uma tabela
            // no caso de se estar usando PVTSim
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
        }
        chute = totalSourceFlowRate;
    }
    double mchute = chute;
    double marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
    state.march.searchOrigin = 1;
    int kontaiter = 0;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10 && kontaiter < 400) {
            mchute *= (1. - state.march.input.buscaFC);
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
            if (mchute < 1e-5) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaProdPresPresPerm atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            marchResidual = marchProductionPressureToPressureSecondary(state.march, mchute);
            kontaiter++;
        }
    }
    double mchuteAux;
    double mchute2 = mchute;
    double guessLowerBound = 0.;
    kontaiter = 0;

    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;

    int testaEscoa = 1;
    if (fabs(marchResidual) < 1e-3)
        return mchute;
    else {
        return bracketTertiaryPressureToPressureRoot(state, testaEscoa, positiveResidualGuess, negativeResidualGuess, guessLowerBound, mchute2, mchuteAux, kontaiter, marchResidual, mchute, maximumFlowRate);
    }
}

double searchInjectionBottomHolePressure1(const SteadyStateSearchState &state, double chute) {
    double pavanc = state.march.gasSurfacePressure;
    double mchute = 0;
    double taux = state.march.cells[0].calor.Textern1;
    if (chute < 0) {
        if (state.march.input.flashCompleto < 1) {
            for (int i = 1; i <= state.march.lastCell; i++) {
                double rhol = state.march.cells[i].fluicol.MasEspFlu(pavanc, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pavanc -= rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed / kPascalPerKgfPerCm2;
                if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i].acsr.ipr.Pres - pavanc) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Poco injetor provavelmente produzindo quando se utiliza este valor de pressao na superficie");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= state.march.cells[i].acsr.ipr.ij * (state.march.cells[i].acsr.ipr.Pres - pavanc);
                    else
                        mchute -= (state.march.cells[i].acsr.ipr.ij * rhol) * (state.march.cells[i].acsr.ipr.Pres - pavanc) / state.march.cells[i].fluicol.MasEspFlu(1.01, 15.);
                }
                taux = state.march.cells[i].calor.Textern1;
            }
            if (state.march.input.condpocinj.CC == 3) {
                if (state.march.cells[0].acsr.injl.QLiq > 1e-5)
                    mchute = state.march.cells[0].acsr.injl.QLiq;
                else if (state.march.input.condpocinj.tipoFlui < 2)
                    mchute += 100.;
                else
                    mchute += 100. * state.march.cells[state.march.lastCell].fluicol.MasEspFlu(pavanc, taux) / state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.);
            }
        } else {
            for (int i = 1; i <= state.march.lastCell; i++) {
                double rhog = state.march.cells[i].flui.MasEspGas(pavanc, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pavanc -= rhog * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed / kPascalPerKgfPerCm2;
                if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i].acsr.ipr.Pres - pavanc) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Poco injetor provavelmente produzindo quando se utiliza este valor de pressao na superficie");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= state.march.cells[i].acsr.ipr.ij * (state.march.cells[i].acsr.ipr.Pres - pavanc);
                    else
                        mchute -= (state.march.cells[i].acsr.ipr.ij * rhog) * (state.march.cells[i].acsr.ipr.Pres - pavanc) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions);
                }
                taux = state.march.cells[i].calor.Textern1;
            }
            if (state.march.input.condpocinj.CC == 3) {
                if (state.march.cells[0].acsr.injg.QGas > 1e-5)
                    mchute = state.march.cells[0].acsr.injg.QGas;
                else
                    mchute += 100. * state.march.cells[state.march.lastCell].flui.MasEspGas(pavanc, taux) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions);
            }
        }
    } else
        mchute = chute;
    double mchute2;
    double mchuteAux;
    double marchResidual;
    marchResidual = marchInjectionSteady(state.march, mchute);
    int kontaiter = 0;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10) {
            mchute *= 0.9;
            marchResidual = marchInjectionSteady(state.march, mchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "PoÃ§o injetor com problemas na iteracao na busca de um chute de vazao inicial, perda de carga muito alta");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    } else if (marchResidual > 0.9e10) {
        while (marchResidual > 0.9e10) {
            mchute *= 1.1;
            marchResidual = marchInjectionSteady(state.march, mchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "PoÃ§o injetor com problemas na iteracao na busca de um chute de vazao inicial, Vazao muito baixa");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    }
    kontaiter = 0;
    mchute2 = mchute;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    if (fabs(marchResidual) < 1e-15)
        return mchute;
    else {
        if (marchResidual < 0.) {
            negativeResidualGuess = mchute;
            while (marchResidual < 0) {
                mchuteAux = mchute2;
                mchute2 *= 1.1;
                marchResidual = marchInjectionSteady(state.march, mchute2);
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm1 atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                while (marchResidual < -0.9e10) {
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchInjectionSteady(state.march, mchute2);
                    kontaiter++;
                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm1 atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            positiveResidualGuess = mchute2;
        } else if (marchResidual > 0.) {
            positiveResidualGuess = mchute;
            while (marchResidual > 0) {
                mchuteAux = mchute2;
                mchute2 *= 0.9;
                marchResidual = marchInjectionSteady(state.march, mchute2);
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm1 atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                while (marchResidual > 0.9e10) {
                    mchute2 = 0.5 * (mchute2 + mchuteAux);
                    marchResidual = marchInjectionSteady(state.march, mchute2);
                    kontaiter++;
                    if (kontaiter > 200) {
                        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                            NumError(
                                "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm1 atingiu maximo de iteracoes");
                        else {
                            if ((*state.march.globals).iterRede > 0)
                                return -1.1e10;
                            else
                                return 1.1e10;
                        }
                    }
                }
            }
            negativeResidualGuess = mchute2;
        }
        return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 0);
    }
}

namespace {

/// What the downward walk of a surface-pressure injection search does with a
/// guess below 1 kgf/cm2, the code's atmospheric pressure.
///
/// The second and the fifth injection searches bracket their root the same way
/// except here: the second refuses such a guess -- it steps back to 99 % of the
/// previous one instead of 90 %, and gives up if even that is below
/// atmospheric -- and the fifth hands it to the march. Nothing in the code says
/// why.
enum class SubAtmosphericGuess {
    refused,
    marched,
};

/// Brackets and solves the root for the two injection searches whose unknown is
/// the surface pressure: the second (condContorno 0) and the fifth (5).
double bracketInjectionPressureRoot(const SteadyStateSearchState &state, double &positiveResidualGuess, double &negativeResidualGuess, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, SubAtmosphericGuess subAtmosphericGuess) {
    if (marchResidual < 0.) {
        negativeResidualGuess = pchute;
        while (marchResidual < 0) {
            pchuteAux = pchute2;
            pchute2 *= 0.9;
            if (subAtmosphericGuess == SubAtmosphericGuess::refused && pchute2 < 1.) {
                pchute2 = pchuteAux;
                pchute2 *= 0.99;
                if (pchute2 < 1. && state.march.input.AP == 0)
                    NumError("Pressao de injecao abaixo da pressao atmosferica");
                else if (pchute2 < 1.)
                    return 1.1e10;
            }
            marchResidual = marchInjectionSteady(state.march, pchute2);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm2 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            while (marchResidual < -0.9e10) {
                pchute2 = 0.5 * (pchute2 + pchuteAux);
                marchResidual = marchInjectionSteady(state.march, pchute2);
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm2 atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        positiveResidualGuess = pchute2;
    } else if (marchResidual > 0.) {
        positiveResidualGuess = pchute;
        while (marchResidual > 0) {
            pchuteAux = pchute2;
            pchute2 *= 1.1;
            marchResidual = marchInjectionSteady(state.march, pchute2);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm2 atingiu maximo de iteracoes");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
            while (marchResidual > 0.9e10) {
                pchute2 = 0.5 * (pchute2 + pchuteAux);
                marchResidual = marchInjectionSteady(state.march, pchute2);
                kontaiter++;
                if (kontaiter > 200) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Busca de valores iniciais para calculo de zero de funcao em buscaInjPfundoPerm2 atingiu maximo de iteracoes");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
            }
        }
        negativeResidualGuess = pchute2;
    }
    return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 0);
}

}  // namespace

double searchInjectionBottomHolePressure2(const SteadyStateSearchState &state, double chute) {
    double pchute;
    if (chute < 0) {
        if (state.march.input.flashCompleto < 1) {
            if (state.march.input.condpocinj.tipoFlui < 2)
                pchute = state.march.cells[state.march.lastCell].acsr.ipr.Pres + ((state.march.cells[0].acsr.injl.QLiq) / (state.march.cells[state.march.lastCell].acsr.ipr.ij));
            else
                pchute = state.march.cells[state.march.lastCell].acsr.ipr.Pres + ((state.march.cells[0].acsr.injl.QLiq) * state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.) / (state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].fluicol.MasEspFlu(state.march.cells[state.march.lastCell].acsr.ipr.Pres, state.march.cells[state.march.lastCell].calor.Textern1)));
            double taux = state.march.cells[0].calor.Textern1;
            double rmis;
            double j;
            double frictionFactor;

            double visC = state.march.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
            double visMis = visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
            rmis = completionDensityAtGuess;
            double rlcA = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
            double massicC = rlcA * state.march.cells[0].acsr.injl.QLiq / kSecondsPerDay;
            j = (massicC / completionDensityAtGuess) / state.march.cells[0].duto.area;
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);

            double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.march.cells[0].duto.peri / state.march.cells[0].duto.area;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].fluicol.MasEspFlu(pchute, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pchute += (rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - pchute) > -(*state.march.globals).localtiny)
                    pchute = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
                taux = state.march.cells[i].calor.Textern1;
            }
        } else {
            pchute = state.march.cells[state.march.lastCell].acsr.ipr.Pres + ((state.march.cells[0].acsr.injg.QGas) * (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions) / (state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].flui.MasEspGas(state.march.cells[state.march.lastCell].acsr.ipr.Pres, state.march.cells[state.march.lastCell].calor.Textern1)));
            double taux = state.march.cells[0].calor.Textern1;
            double rmis;
            double j;
            double frictionFactor;

            double visC = state.march.cells[0].acsr.injg.FluidoPro.ViscGas(pchute, taux);
            double visMis = visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injg.FluidoPro.MasEspGas(pchute, taux);
            rmis = completionDensityAtGuess;
            double rlcA = (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions);
            double massicC = rlcA * state.march.cells[0].acsr.injg.QGas / kSecondsPerDay;
            j = (massicC / completionDensityAtGuess) / state.march.cells[0].duto.area;
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);

            double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.march.cells[0].duto.peri / state.march.cells[0].duto.area;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].flui.MasEspGas(pchute, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pchute += (rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - pchute) > -(*state.march.globals).localtiny)
                    pchute = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
                taux = state.march.cells[i].calor.Textern1;
            }
        }
    } else
        pchute = chute;
    if (pchute <= 1)
        pchute = 1.01;
    double pchute2;
    double pchuteAux = 0.;
    double marchResidual;
    marchResidual = marchInjectionSteady(state.march, pchute);
    int kontaiter = 0;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10) {
            pchute *= 1.1;
            marchResidual = marchInjectionSteady(state.march, pchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "PoÃ§o injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    } else if (marchResidual > 0.9e10) {
        while (marchResidual > 0.9e10) {
            pchute *= 0.9;
            marchResidual = marchInjectionSteady(state.march, pchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "PoÃ§o injetor com problemas na iteracao na busca de um chute de vazao inicial, pressao muito alta");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    }
    kontaiter = 0;
    pchute2 = pchute;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    if (fabs(marchResidual) < 1e-15)
        return pchute;
    else {
        return bracketInjectionPressureRoot(state, positiveResidualGuess, negativeResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, pchute, SubAtmosphericGuess::refused);
    }
}

double searchInjectionBottomHolePressure3(const SteadyStateSearchState &state, double chute) {
    state.march.cells[state.march.lastCell].pres = state.march.input.condpocinj.presfundo;
    if (state.march.cells[state.march.lastCell].pres < state.march.cells[state.march.lastCell].acsr.ipr.Pres) {
        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
            NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
        else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    double mchute;
    if (chute < 0.) {
        if (state.march.input.flashCompleto < 1) {
            if (state.march.input.condpocinj.tipoFlui < 2)
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij) * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            else
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].fluicol.MasEspFlu(state.march.cells[state.march.lastCell].pres,
                                                                                     state.march.cells[state.march.lastCell].calor.Textern1)) *
                         (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres) / state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.);
            double taux = state.march.cells[0].calor.Textern1;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].fluicol.MasEspFlu(state.march.cells[i].pres, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                state.march.cells[i - 1].pres = state.march.cells[i].pres + rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed / kPascalPerKgfPerCm2;
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i + 1].pres) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= (state.march.cells[i - 1].acsr.ipr.ij) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                    mchute -= (state.march.cells[i - 1].acsr.ipr.ij * state.march.cells[i - 1].fluicol.MasEspFlu(state.march.cells[i - 1].pres, 0.) / state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.)) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                }
                taux = state.march.cells[i].calor.Textern1;
            }
        } else {
            if (state.march.input.condpocinj.tipoFlui < 2)
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij) * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            else
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].flui.MasEspGas(state.march.cells[state.march.lastCell].pres,
                                                                                  state.march.cells[state.march.lastCell].calor.Textern1)) *
                         (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions);
            double taux = state.march.cells[0].calor.Textern1;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].flui.MasEspGas(state.march.cells[i].pres, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                state.march.cells[i - 1].pres = state.march.cells[i].pres + rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed / kPascalPerKgfPerCm2;
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i + 1].pres) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= (state.march.cells[i - 1].acsr.ipr.ij) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                    mchute -= (state.march.cells[i - 1].acsr.ipr.ij * state.march.cells[i - 1].flui.MasEspGas(state.march.cells[i - 1].pres, 0.) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions)) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                }
                taux = state.march.cells[i].calor.Textern1;
            }
        }
    } else
        mchute = chute;

    state.march.gasSurfacePressure = state.march.cells[0].pres;
    double marchResidual;
    marchResidual = marchInjectionSteady(state.march, mchute);
    if (marchResidual < -0.9e10) {
        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
            NumError(
                "PoÃ§o injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
        else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    if (marchResidual > 0.9e10) {
        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
            NumError(
                "PoÃ§o injetor com problemas na iteracao na busca de um chute de vazao inicial, pressao muito alta");
        else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    int konta = 0;
    while (fabs(state.march.cells[state.march.lastCell].pres - state.march.input.condpocinj.presfundo) / state.march.input.condpocinj.presfundo > 0.01 / 100 && konta < 200) {
        state.march.cells[state.march.lastCell].pres = state.march.input.condpocinj.presfundo;
        if (state.march.input.flashCompleto < 1) {
            if (state.march.input.condpocinj.tipoFlui < 2)
                mchute = -state.march.cells[state.march.lastCell].acsr.ipr.ij * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            else
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].fluicol.MasEspFlu(state.march.cells[state.march.lastCell].pres, state.march.cells[state.march.lastCell].temp) / state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.)) * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            for (int i = state.march.lastCell; i > 0; i--) {
                state.march.cells[i - 1].pres = state.march.cells[i].pres + state.march.updaters.steadyInjectionPressureDrop(i);
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i + 1].pres) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= state.march.cells[i - 1].acsr.ipr.ij * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                    else
                        mchute -= (state.march.cells[i - 1].acsr.ipr.ij * state.march.cells[i - 1].fluicol.MasEspFlu(state.march.cells[i - 1].pres, state.march.cells[i - 1].temp) / state.march.cells[state.march.lastCell].fluicol.MasEspFlu(1.01, 15.)) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                }
            }
        } else {
            if (state.march.input.condpocinj.tipoFlui < 2)
                mchute = -state.march.cells[state.march.lastCell].acsr.ipr.ij * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            else
                mchute = -(state.march.cells[state.march.lastCell].acsr.ipr.ij * state.march.cells[state.march.lastCell].flui.MasEspGas(state.march.cells[state.march.lastCell].pres, state.march.cells[state.march.lastCell].temp) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions)) * (state.march.cells[state.march.lastCell].acsr.ipr.Pres - state.march.cells[state.march.lastCell].pres);
            for (int i = state.march.lastCell; i > 0; i--) {
                state.march.cells[i - 1].pres = state.march.cells[i].pres + state.march.updaters.steadyInjectionPressureDrop(i);
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i + 1].pres) > -(*state.march.globals).localtiny) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance) {
                    if (state.march.input.condpocinj.tipoFlui < 2)
                        mchute -= state.march.cells[i - 1].acsr.ipr.ij * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                    else
                        mchute -= (state.march.cells[i - 1].acsr.ipr.ij * state.march.cells[i - 1].flui.MasEspGas(state.march.cells[i - 1].pres, state.march.cells[i - 1].temp) / (state.march.cells[0].flui.Deng * kAirDensityAtStandardConditions)) * (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres);
                }
            }
        }
        state.march.gasSurfacePressure = state.march.cells[0].pres;
        marchResidual = marchInjectionSteady(state.march, mchute);
        konta++;
        if (marchResidual < -0.9e10) {
            if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                NumError(
                    "PoÃ§o injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
            else {
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        }
        if (marchResidual > 0.9e10) {
            if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                NumError(
                    "PoÃ§o injetor com problemas na iteracao na busca de um chute de pressao inicial, pressao muito alta");
            else {
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
        }
    }
    if (konta >= 200) {
        if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
            NumError("PoÃ§o injetor com problemas na iteracao");
        else {
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    }
    return (state.march.cells[state.march.lastCell].pres - state.march.input.condpocinj.presfundo);
}

double searchInjectionBottomHolePressure4(const SteadyStateSearchState &state) {

    int guessNeedsCorrection = 1;

    if (state.march.input.flashCompleto < 1)
        state.march.cells[0].temp = state.march.cells[0].acsr.injl.temp;
    else
        state.march.cells[0].temp = state.march.cells[0].acsr.injg.temp;
    double alfini = 0.;
    double betini = 1.;

    if (state.march.input.flashCompleto >= 1) {
        if (state.march.cells[0].acsr.injg.seco == 1) {
            alfini = 1.;
            betini = 0.;
        } else {
            if (state.march.input.flashCompleto == 2) {
                if (state.march.input.tabelaDinamica == 0)
                    state.march.cells[0].flui.atualizaPropComp(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp, -1, NULL, NULL, state.march.input.pocinjec);
                state.march.cells[0].acsr.injg.FluidoPro.atualizaPropComp(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp, -1, NULL, NULL, state.march.input.pocinjec);
            }
            double masgas = state.march.cells[0].acsr.injg.VMas(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp);
            double quality;
            if (state.march.input.flashCompleto != 2)
                quality = state.march.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
            else
                quality = state.march.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
            double masT = masgas / quality;
            quality = state.march.cells[0].acsr.injg.FluidoPro.FracMassHidra(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp);
            double qgas = masT * quality /
                          state.march.cells[0].acsr.injg.FluidoPro.MasEspGas(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp);
            double qliq = masT * (1. - quality) /
                          state.march.cells[0].acsr.injg.FluidoPro.MasEspLiq(state.march.gasSurfacePressure, state.march.cells[0].acsr.injg.temp);
            alfini = qgas / (qliq + qgas);
            betini = 0.;
        }
    }

    state.march.cells[0].tempL = state.march.cells[0].temp;
    state.march.cells[1].tempL = state.march.cells[0].temp;
    state.march.cells[0].tempini = state.march.cells[0].temp;
    state.march.cells[0].ML = 0.;
    state.march.cells[0].MC = 0.;
    state.march.cells[1].ML = 0.;
    state.march.cells[0].MliqiniL = 0.;
    state.march.cells[0].Mliqini = 0.;
    state.march.cells[1].MliqiniL = 0.;
    state.march.cells[0].QLL = 0.;
    state.march.cells[0].QL = 0;
    state.march.cells[1].QLL = 0.;
    state.march.cells[0].QG = 0.;
    if (state.march.input.condpocinj.CC == 0)
        state.march.gasSurfacePressure = state.march.input.condpocinj.presinj;
    state.march.steadyIteration = 0;

    double masfim = 0.;

    int i;
    guessNeedsCorrection = 1;
    while (guessNeedsCorrection == 1) {
        state.march.cells[0].presauxL = state.march.gasSurfacePressure;
        state.march.cells[0].presLini = state.march.gasSurfacePressure;
        state.march.cells[0].presL = state.march.gasSurfacePressure;
        state.march.cells[0].pres = state.march.gasSurfacePressure;
        state.march.cells[1].presL = state.march.gasSurfacePressure;
        state.march.cells[0].presini = state.march.gasSurfacePressure;
        state.march.cells[1].presLini = state.march.gasSurfacePressure;
        state.march.cells[0].presaux = state.march.gasSurfacePressure;
        state.march.cells[1].presauxL = state.march.gasSurfacePressure;
        state.march.cells[0].alf = alfini;
        state.march.cells[0].alfini = alfini;
        state.march.cells[0].bet = betini;
        state.march.cells[0].betini = betini;
        state.march.cells[1].alfL = state.march.cells[0].alf;
        state.march.cells[1].alfLini = state.march.cells[0].alf;
        state.march.cells[0].alfPigD = state.march.cells[0].alf;
        state.march.cells[0].alfPigDini = state.march.cells[0].alf;
        state.march.cells[0].alfPigE = state.march.cells[0].alf;
        state.march.cells[0].alfPigEini = state.march.cells[0].alf;
        state.march.cells[1].betL = state.march.cells[0].bet;
        state.march.cells[1].betLini = state.march.cells[0].bet;
        state.march.cells[0].betPigD = state.march.cells[0].bet;
        state.march.cells[0].betPigDini = state.march.cells[0].bet;
        state.march.cells[0].betPigE = state.march.cells[0].bet;
        state.march.cells[0].betPigEini = state.march.cells[0].bet;
        state.march.cells[0].betI = state.march.cells[0].bet;
        state.march.cells[1].betLI = state.march.cells[0].bet;
        i = 1;
        while (i <= state.march.lastCell && state.march.cells[i - 1].pres >= 1. && masfim >= -(*state.march.globals).localtiny) {

            advanceUpstreamSteadyPressure(state.march, i, 0);
            refreshUpstreamProductionPeriphery(state.march, i);
            if (state.march.input.flashCompleto != 2)
                advanceSteadyMass(state.march, i);
            else
                advanceCompositionalSteadyMass(state.march, i);
            state.march.updaters.advanceSteadyTemperature(i, 0);
            state.march.updaters.updateProductionTemperaturePeriphery(i);
            advanceDownstreamSteadyPressure(state.march, i, 0);
            refreshDownstreamProductionPeriphery(state.march, i);
            advanceSteadyMassTransfer(state.march, i - 1);
            masfim += (state.march.cells[i - 1].fontemassCR + state.march.cells[i - 1].fontemassLR + state.march.cells[i - 1].fontemassGR);
            i++;
            if (state.march.cells[i - 1].pres <= 1) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("pressao menor que 1");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            } else if ((state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - state.march.cells[i - 1].pres) > -(*state.march.globals).localtiny)) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Valor de pressao de fundo menor que a pressao de reservatÃ³rio");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            } else if (i < (state.march.lastCell + 1) && masfim < -(*state.march.globals).localtiny) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError("Vazao massica de injeÃ§Ã£o menor que 0");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
        if (i == state.march.lastCell + 1)
            guessNeedsCorrection = 0;
    }
    state.march.updaters.updateSource(state.march.lastCell);
    masfim += (state.march.cells[state.march.lastCell].fontemassCR + state.march.cells[state.march.lastCell].fontemassLR + state.march.cells[state.march.lastCell].fontemassGR);

    return masfim;
}

double searchInjectionBottomHolePressure5(const SteadyStateSearchState &state, double chute) {

    double pchute = state.march.input.condpocinj.presfundo;
    double taux;
    double j = 0.;
    double rmis = 0.;
    double frictionFactor = 0.;
    taux = state.march.input.celp[0].textern;
    if (chute < 0.) {
        if (state.march.input.flashCompleto < 1) {
            double visC = state.march.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
            double visMis = visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
            rmis = completionDensityAtGuess;
            double rlcA = state.march.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
            double massicC = rlcA * state.march.cells[0].acsr.injl.QLiq / kSecondsPerDay;
            j = (massicC / completionDensityAtGuess) / state.march.cells[0].duto.area;
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);
            double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.march.cells[0].duto.peri / state.march.cells[0].duto.area;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].fluicol.MasEspFlu(pchute, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pchute += (rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
                if (pchute < 0.5) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Poco injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - pchute) > -(*state.march.globals).localtiny)
                    pchute = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
        } else {
            double visC = state.march.cells[0].acsr.injg.FluidoPro.ViscGas(pchute, taux);
            double visMis = visC;
            double completionDensityAtGuess = state.march.cells[0].acsr.injg.FluidoPro.MasEspGas(pchute, taux);
            rmis = completionDensityAtGuess;
            double rlcA = (state.march.cells[0].acsr.injg.FluidoPro.Deng * kAirDensityAtStandardConditions);
            double massicC = rlcA * state.march.cells[0].acsr.injg.QGas / kSecondsPerDay;
            j = (massicC / completionDensityAtGuess) / state.march.cells[0].duto.area;
            double reynolds;
            if (state.march.cells[0].duto.revest == 0)
                reynolds = state.march.cells[0].Rey(state.march.cells[0].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.march.cells[0].duto.area / state.march.cells[0].duto.peri;
                reynolds = state.march.cells[0].Rey(dhid, j, rmis, visMis);
            }
            frictionFactor = state.march.cells[0].fric(reynolds, state.march.cells[0].duto.rug / state.march.cells[0].duto.a);
            double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.march.cells[0].duto.peri / state.march.cells[0].duto.area;
            for (int i = state.march.lastCell; i > 0; i--) {
                double rhol = state.march.cells[i].flui.MasEspGas(pchute, taux);
                double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
                pchute += (rhol * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
                if (pchute < 0.5) {
                    if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                        NumError(
                            "Poco injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
                    else {
                        if ((*state.march.globals).iterRede > 0)
                            return -1.1e10;
                        else
                            return 1.1e10;
                    }
                }
                if (state.march.cells[i - 1].acsr.tipo == kAccessoryInflowPerformance && (state.march.cells[i - 1].acsr.ipr.Pres - pchute) > -(*state.march.globals).localtiny)
                    pchute = 1.01 * state.march.cells[i - 1].acsr.ipr.Pres;
            }
        }
    } else
        pchute = chute;
    double pchute2;
    double pchuteAux = 0.;
    double marchResidual;
    marchResidual = marchInjectionSteady(state.march, pchute);
    int kontaiter = 0;
    if (marchResidual < -0.9e10) {
        while (marchResidual < -0.9e10) {
            pchute *= 1.1;
            marchResidual = marchInjectionSteady(state.march, pchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Poco injetor com problemas na iteracao na busca de um chute de pressao inicial, perda de carga muito alta");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    } else if (marchResidual > 0.9e10) {
        while (marchResidual > 0.9e10) {
            pchute *= 0.9;
            marchResidual = marchInjectionSteady(state.march, pchute);
            kontaiter++;
            if (kontaiter > 200) {
                if ((*state.march.globals).chaverede == 0 && state.march.input.AP == 0)
                    NumError(
                        "Poco injetor com problemas na iteracao na busca de um chute de vazao inicial, pressao muito alta");
                else {
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
    }
    kontaiter = 0;
    pchute2 = pchute;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    if (fabs(marchResidual) < 1e-15)
        return pchute;
    else {
        return bracketInjectionPressureRoot(state, positiveResidualGuess, negativeResidualGuess, kontaiter, marchResidual, pchuteAux, pchute2, pchute, SubAtmosphericGuess::marched);
    }
}

namespace {

/// Brackets the branch root upward when the first march came back positive.
///
/// One arm of the sign split in bracketSecondaryBranchRoot.
bool bracketSecondaryBranchFromLowGuess(const SteadyStateSearchState &state, double amplifica, double &positiveResidualGuess, double &negativeResidualGuess, double &guessLowerBound, int &kontaiter, double mult2, double &marchResidual, double &pchuteAux, double &pchute2, double pchute, double &abortValue) {
    // Choke flow is below tubing flow, indicating a low pressure guess.
    // Increase the pressure until marchResidual becomes negative and brackets the root.
    positiveResidualGuess = pchute; // Store the pressure guess yielding marchResidual > 0.

    // Update the positive bound as marchResidual approaches zero.
    while (marchResidual > 0) {
        pchuteAux = pchute2;
        if ((*state.march.globals).chaverede == 0)
            pchute2 *= amplifica; // Increasing the pressure in search of marchResidual>0
        else
            pchute2 *= mult2; // Increasing the pressure in search of marchResidual>0
        if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
            pchute2 = 0.9 * state.march.input.tabent.pmax;
        int limpres = 0;
        marchResidual = marchProductionSteadySecondary(state.march, pchute2);
        if (fabs(marchResidual) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (marchResidual > 0. && limpres == 1) {
            int iterpres = 0;
            while (marchResidual > 0 && iterpres < 10) {
                pchute2 *= 1.001;
                marchResidual = marchProductionSteadySecondary(state.march, pchute2);
                if (fabs(marchResidual) > 1.01e10) {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
                iterpres++;
            }
            if (iterpres >= 10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        if (marchResidual > 0 && marchResidual < 0.9e10)
            positiveResidualGuess = pchute2; // updating positiveResidualGuess
        kontaiter++;
        if (kontaiter > 50) {
            // Iteration limit reached while searching for the second guess.
            // End the simulation or report the failure.
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        while (marchResidual > 0.9e10) {
            // The pressure step was too large, causing the marching process to exceed
            // the static or PVTSim pressure limit. Reduce the lower pressure estimate.
            guessLowerBound = pchute2;
            pchute2 = 0.5 * (pchute2 + pchuteAux); // Midpoint between the highest pressure yielding marchResidual > 0 and the upper pressure bound.
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
                pchute2 = 0.9 * state.march.input.tabent.pmax;
            marchResidual = marchProductionSteadySecondary(state.march, pchute2);
            if (fabs(marchResidual) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (marchResidual > 0 && marchResidual < 0.9e10)
                positiveResidualGuess = pchute2;
            kontaiter++;

            if (kontaiter > 50) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
    }
    negativeResidualGuess = pchute2;
    return false;
}

/// Brackets and solves the flow-rate root for the secondary branch.
double bracketSecondaryBranchRoot(const SteadyStateSearchState &state, double amplifica, double reduz, double &positiveResidualGuess, double &negativeResidualGuess, double &guessLowerBound, int &kontaiter, double mult2, double mult1, double &marchResidual, double &pchuteAux, double &pchute2, double pchute) {
    if (marchResidual < 0.) {
        // Choke flow exceeds tubing flow, indicating a high pressure guess.
        // Decrease the pressure until marchResidual becomes positive and brackets the root.
        negativeResidualGuess = pchute; // Store the pressure guess yielding marchResidual < 0.

        // Update the negative bound as marchResidual approaches zero.
        while (marchResidual < 0) {
            if (fabs(pchute2 - pchute) / pchute < (1. - reduz) / 10. && kontaiter > 50) {
                pchute2 *= 0.5;
            }
            pchuteAux = pchute2;
            if ((*state.march.globals).chaverede == 0)
                pchute2 *= reduz; // Decreasing the pressure in the search for marchResidual>0
            else
                pchute2 *= mult1; // Decreasing the pressure in the search for marchResidual>0
            if (pchute2 <= guessLowerBound)
                pchute2 = 0.5 * (pchuteAux + guessLowerBound);
            // guessLowerBound starts at zero. If pchute2 drops too far and marchResidual reaches -1e10,
            // set guessLowerBound to pchute2, establishing the minimum allowed pressure.
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
                pchute2 = 0.9 * state.march.input.tabent.pmax;
            marchResidual = marchProductionSteadySecondary(state.march, pchute2);
            if (fabs(marchResidual) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
            if (marchResidual < 0 && marchResidual > -0.9e10)
                negativeResidualGuess = pchute2; // updating negativeResidualGuess
            kontaiter++;
            if (kontaiter > 50) {
                // Iteration limit reached while searching for the second guess.
                // Stop the simulation or report the failure.
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    return -1.1e10;
                else
                    return 1.1e10;
            }
            while (marchResidual < -0.9e10) {
                // The pressure reduction was too large, causing the marching process
                // to fall below 0.5 kgf/cm². Increase the lower pressure estimate.
                guessLowerBound = pchute2;
                pchute2 = 0.5 * (pchute2 + pchuteAux); // Midpoint between the lowest pressure yielding marchResidual < 0 and the lower pressure bound.
                if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute2) < (*state.march.globals).localtiny)
                    pchute2 = 0.9 * state.march.input.tabent.pmax;
                if (state.march.cells[0].acsr.tipo == kAccessoryInflowPerformance &&
                    (state.march.cells[0].acsr.ipr.Pres - pchute2) < 0.001 * state.march.cells[0].acsr.ipr.Pres) {
                    pchute2 = -(10 / state.march.cells[0].acsr.ipr.ip) + state.march.cells[0].acsr.ipr.Pres;
                    if (pchute2 < 0.999 * state.march.cells[0].acsr.ipr.Pres)
                        pchute2 = 0.999 * state.march.cells[0].acsr.ipr.Pres;
                }
                marchResidual = marchProductionSteadySecondary(state.march, pchute2);
                if (fabs(marchResidual) > 1.01e10) {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
                if (marchResidual < 0 && marchResidual > -0.9e10)
                    negativeResidualGuess = pchute2;
                kontaiter++;
                if (kontaiter > 50) {
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    if ((*state.march.globals).iterRede > 0)
                        return -1.1e10;
                    else
                        return 1.1e10;
                }
            }
        }
        positiveResidualGuess = pchute2;
    } else if (marchResidual > 0.) {
        double abortValue;
        if (bracketSecondaryBranchFromLowGuess(state, amplifica, positiveResidualGuess, negativeResidualGuess, guessLowerBound, kontaiter, mult2, marchResidual, pchuteAux, pchute2, pchute, abortValue))
            return abortValue;
    }
    return solveSteadyRoot(state, negativeResidualGuess, positiveResidualGuess, 1, 1); // Find the root using pressure bounds with opposite signs.
}

/// Walks the guess until the branch march stops returning a sentinel.
bool retryUntilBranchMarchCompletes(const SteadyStateSearchState &state, double pchuteAux0, int &kontaiter, double &marchResidual, double &pchuteAux, double &pchute, double &abortValue) {
    if (marchResidual < -0.9e10 || marchResidual > 0.9e10) {
        double valtemp;
        valtemp = marchProductionSteadySecondary(state.march, pchuteAux); // marcha with pchuteAux
        if (fabs(valtemp) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        while (valtemp > 0.9e10 && marchResidual > 0.9e10 && kontaiter < 50) { // Estimated background pressure remains high
            pchuteAux *= 0.99;                                       // reducing the estimate
            valtemp = marchProductionSteadySecondary(state.march, pchuteAux);
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            kontaiter++; // 50 iterations maximum
        }
        while (valtemp < -0.9e10 && marchResidual < -0.9e10 && kontaiter <= 50) { // Estimated background pressure remains low
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection || state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection || state.march.cells[0].acsr.tipo == kAccessoryMultipleSource)
                pchuteAux *= 1.1;
            else
                pchuteAux *= 1.01; // increasing the estimate
            // Ensure the increased estimate does not exceed the PVTSim pressure limit.
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
            valtemp = marchProductionSteadySecondary(state.march, pchuteAux); // new try
            if (fabs(valtemp) > 1.01e10) {
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
                if ((*state.march.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
            if (kontaiter == 50) {
                pchuteAux = 0.5 * pchuteAux0;
            }
            kontaiter++;
        }
        valtemp = marchProductionSteadySecondary(state.march, pchuteAux);
        if (fabs(valtemp) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                {
                    abortValue = -1.1e10;
                    return true;
                }
            else
                {
                    abortValue = 1.1e10;
                    return true;
                }
        }
        if (valtemp > -0.9e10 && valtemp < 0.9e10) { // The marching process completed successfully.
            marchResidual = valtemp;
            pchute = pchuteAux;
        }
    }
    return false;
}

/// Moves the guess according to which sentinel the branch march returned.
void classifyBranchMarchSentinel(const SteadyStateSearchState &state, double marchResidual, double &pchuteAux, double perdafric, double &taux) {
    if (marchResidual < -0.9e10) {
        // Pressure is too low. Increase the guess using water hydrostatics
        // to obtain a higher bottom-hole pressure.
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = 1000 + 0 * state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            double alfa = 0.;
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.5 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess * 2.;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;  // limit of pchuteAux
    } else if (marchResidual > 0.9e10) { // pressure at some point exceeded some static pressure
        // must decrease the chute value
        pchuteAux = state.march.gasSurfacePressure;
        for (int i = state.march.lastCell; i > 0; i--) {
            taux = state.march.input.celp[i].textern;
            double rhol = state.march.cells[i].flui.MasEspLiq(pchuteAux, taux);
            double rhog = state.march.cells[i].flui.MasEspGas(pchuteAux, taux);
            // in this case, use a high void fraction for the hydrostatic calculation.
            double alfa = 0.8;
            if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection)
                alfa = 1.;
            else if ((*state.march.globals).chaverede == 1 && state.holdupGuess < 0.2 && state.holdupGuess > -1e-15)
                alfa = 1. - state.holdupGuess;
            double rhomix = (1. - alfa) * rhol + alfa * rhog;
            double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
            pchuteAux += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchuteAux) < (*state.march.globals).localtiny)
                pchuteAux = 0.9 * state.march.input.tabent.pmax;
        }
        if (pchuteAux > 1100)
            pchuteAux = 1100;
    }
}

}  // namespace

double searchSecondaryBranchFlowRate(const SteadyStateSearchState &state, double pPartida, int indPartida) {
    // Find two initial guesses with opposite signs to bracket the root.
    state.reverseSteady = 0;
    state.march.convergenceMonitor = 1000.;

    // Initialize the pressure guess with the downstream choke pressure.
    // If non-positive, estimate a new value.
    double pchute = pPartida;
    double taux; // Auxiliary temperature used to estimate pchute.
    double perdafric = 0.;
    double completionFractionGuess = 0.;
    if (state.march.cells[0].acsr.tipo == kAccessoryLiquidInjection)
        completionFractionGuess = state.march.cells[0].acsr.injl.bet;

    for (int i = indPartida; i > 0; i--) {
        taux = state.march.input.celp[i].textern;
        if (taux < state.march.input.tmin)
            taux = state.march.input.tmin;
        double rhol = state.march.cells[i].flui.MasEspLiq(pchute, taux);
        double rhog = state.march.cells[i].flui.MasEspGas(pchute, taux);
        // For systems other than the main gas-lift ring, estimate the pressure
        // using liquid hydrostatics, resulting in a high initial pressure guess.
        double alfa = 0.;
        if (completionFractionGuess < 0.5) {
            double quality = state.march.cells[i].flui.FracMassHidra(pchute, taux);
            alfa = quality * rhol / (rhog - quality * rhog + quality * rhol);
        }
        if (state.march.cells[0].acsr.tipo == kAccessoryGasInjection)
            alfa = 1.;
        double rhomix = (1. - alfa) * rhol + alfa * rhog;
        double dxmed = 0.5 * (state.march.cells[i].dx + state.march.cells[i - 1].dx);
        // Advance the pressure using hydrostatic head and estimated friction loss.
        pchute += (rhomix * 9.81 * sin(state.march.cells[i].duto.teta) * dxmed + perdafric * dxmed) / kPascalPerKgfPerCm2;
        if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
            pchute = 0.9 * state.march.input.tabent.pmax;
    }
    if (pchute > 1000)
        pchute = 1000.; // maximum pressure of maximum chute

    // Method parameters.
    double pchute2;        // Second pressure guess.
    double pchuteAux = 0.; // Helper used to bracket the root.
    double marchResidual;            // marchResidual stores the value returned by the marching procedure.
    // Difference between the upstream choke pressure and the final cell pressure.
    marchResidual = marchProductionSteadySecondary(state.march, pchute);
    state.march.searchOrigin = 1;
    // The marching process may fail if pressure exceeds an IPR static pressure
    // or the PVTSim table limit, returning 1e10, or approaches zero before
    // reaching the final cell, returning -1e10.
    classifyBranchMarchSentinel(state, marchResidual, pchuteAux, perdafric, taux);

    double mult1 = 0.9;
    double mult2 = 1.1;

    // Retry with pchuteAux if the first marching attempt fails.
    int kontaiter = 0; // Count attempts until the result differs from 1e10 or 1e-10.
    double pchuteAux0 = pchuteAux;
    double abortValue;
    if (retryUntilBranchMarchCompletes(state, pchuteAux0, kontaiter, marchResidual, pchuteAux, pchute, abortValue))
        return abortValue;
    if (kontaiter > 50) { // Iteration limit reached.
        // During a network iteration after the first one.
        if ((*state.march.globals).iterRede > 0)
            return -1.1e10;
        // If followed by a transient simulation or during the first network iteration.
        else
            return 1.1e10;
    }
    // The previous loop jumped between pressure guesses that were too high
    // and too low to complete the marching process.
    while ((marchResidual < -0.9e10 || marchResidual > 0.9e10) && kontaiter <= 50) {
        pchute = 0.5 * (pchute + pchuteAux); // Use the midpoint between the upper and lower bounds.
        if ((fabs(pchute - pchuteAux) / pchute) < 0.001 && marchResidual < -0.9e10) {
            pchute = 1.05 * pchute;
            if (state.march.input.usaTabela == 1 && (state.march.input.tabent.pmax - pchute) < (*state.march.globals).localtiny)
                pchute = 0.9 * state.march.input.tabent.pmax;
            pchuteAux = pchute;
        }

        marchResidual = marchProductionSteadySecondary(state.march, pchute);
        if (fabs(marchResidual) > 1.01e10) {
            cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
            if ((*state.march.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
        kontaiter++;
    }
    if (kontaiter > 50) {
        cout << "#################PERMANENTE FALHOU EM SUA CONVERGÃŠNCIA##############################" << endl;
        if ((*state.march.globals).iterRede > 0)
            return -1.1e10;
        else
            return 1.1e10;
    }

    // A pressure guess was found that allowed the marching process to complete.
    kontaiter = 0;
    pchute2 = pchute;
    double guessLowerBound = 0.;
    double negativeResidualGuess = 0.;
    double positiveResidualGuess = 0.;
    double reduz = 1. - state.march.input.buscaFC;
    double amplifica = 1. + state.march.input.buscaFC;
    if (fabs(marchResidual) < 1e-3)
        return pchute;
    else {
        return bracketSecondaryBranchRoot(state, amplifica, reduz, positiveResidualGuess, negativeResidualGuess, guessLowerBound, kontaiter, mult2, mult1, marchResidual, pchuteAux, pchute2, pchute);
    }
}

}  // namespace sisprod::steady
