
#define _USE_MATH_DEFINES // para M_PI

#include "Acidentes2.h"
#include "Bcsm2.h"
#include "BombaVol.h"
#include "Elem2D.h"
#include "Elem2DPoroso.h"
#include "Elem3DPoisson.h"
#include "FA_Hidratos.h"
#include "FA_Hidratos_Servico.h"
#include "FerramentasNumericas.h"
#include "FonteMas.h"
#include "FonteMasVap.h"
#include "FonteMassCHK.h"
#include "Geometria.h"
#include "GeometriaPoro.h"
#include "GradientCorrelations.h"
#include "Leitura.h"
#include "LeituraVapor.h"
#include "LerAP.h"
#include "LerRede.h"
#include "Log.h"
#include "Malha2D.h"
#include "Malha2DPoroso.h"
#include "Malha3DPoisson.h"
#include "Matriz.h"
#include "PorosoRad-Simples.h"
#include "PorosoRad.h"
#include "PropFlu.h"
#include "PropFluCol.h"
#include "PropFluColVF.h"
#include "PropVapor.h"
#include "SisProd_old.h"
#include "TrocaCalor.h"
#include "Vetor.h"
#include "acessorios.h"
#include "celRad-Simples.h"
#include "celRad.h"
#include "celula3.h"
#include "celulaGas.h"
#include "celulaVapor.h"
#include "chokegas.h"
#include "criterioIntermiSevera.h"
#include "dados1Poisson.h"
#include "dados1Poroso.h"
#include "dados3DPoisson.h"
#include "estrat.h"
#include "estruturaTabDin.h"
#include "estruturaUNV.h"
#include "estruturas.h"
#include "estruturasPoisson3D.h"
#include "estruturasPoroso.h"
#include "mapa.h"
#include "multiBCS.h"
#include "solver.h"
#include "solver3DPoisson.h"
#include "solverPoisson.h"
#include "solverPoroso.h"
#include "variaveisGlobais1D.h"
#include "versao.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <signal.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

using namespace std;

#ifdef linux
// barra para linux
#define BARRA "/"
#elif LINUX
// barra para linux
#define BARRA "/"
#elif Linux
// barra para linux
#define BARRA "/"
#elif UNIX
// barra para linux
#define BARRA "/"
#elif defined WIN32
// barra para windows
#define BARRA "\\"
#elif defined Win32
// barra para windows
#define BARRA "\\"
#elif defined win32
// barra para windows
#define BARRA "\\"
#else
// barra para outros
#define BARRA "/\\"
#endif

#define ARQUIVO_RELATORIO_PERFIS "resultado.log"

const char *saidaTexto[15] = {"                          Post Coitum Omine Animal Triste Est                   ",
                              "           'Ouca-me. O fim quase nunca esta longe, em nenhum momento!'          ",
                              "      So nos curamos de um sofrimento depois de o haver suportado ate o fim.    ",
                              "                   Infeliz e o espirito ansioso pelo futuro.                    ",
                              "                                    Memento Mori                                ",
                              " Somente um progresso calmo e constante, livre de precipitacao, conduz ao objetivo.",
                              "             Paciencia, nove mulheres nao conseguem gerar uma crianca em um mes. ",
                              "                  A necessidade e a mae da inovacao, mas a paciencia e o pai    ", /*, dando tempo e banana, qualquer macaco faz um projeto"*/
                              "O sucesso nao e uma linha reta, e um jogo de resistencia, e cada tropeco e apenas um degrau a mais para a vitoria!",
                              "                        Quem vive de navegar, o vento e quem lhe comanda                ",
                              "    Uma vez me perguntaram o que achava da passagem do tempo, e eu disse: sou contra    ",
                              "                 Nao importa o quanto voce va devagar, desde que nao pare                ",
                              "Um simulador que resolve uma parada de producao, comeca avancando pequenos incrementos de tempo",
                              "                            Nada e permanente, exceto a mudanca                           ",
                              "                  Uma jornada de mil quilometros come�a com um unico passo                "};
const char *saidaSubTexto[15] = {
    "                         Galeno de Pergamo do Transiente Longo                          ",
    "                     J. California Cooper depois da simulacao divergir                  ",
    "                                Marcel Proust no CrossFit                               ",
    "                              Seneca do Mindfulness                                     ",
    "                                   Zuleica da Funeraria                                 ",
    "                                       China In Box                                     ",
    "                                      Tiao do Linkedin                                  ",
    "                                      Marcao da Oficina                                 ",
    "                                    Mario Pascal do Insta                               ",
    "                          Seu Pereira na feira de artesanatos numericos                 ",
    "                        Luis Fernando Verissimo das Simulacoes Permanentes              ",
    "                             Confucio vendo a simulacao emperrar                       ",
    "                               Confucio das simulacoes sem fim                         ",
    "          Heraclito de Efeso vendo tudo mudar a cada incremento de tempo                ",
    "    Lao-Tse tomando coragem para simular um caso de parafinacao em dutos de producao    "};

detTempo tempVF;
detProp prop;
detMapProp mapprop;
detCI CI;
detCC CC;

int nthrdMatriz = 1;
string versao = "Marlim - " VERSAO;

time_t nowGlobIni;
tm *ltmGlobIni;

int diaIni;
int horaIni;
int minutoIni;
int segundoIni;

time_t nowGlobFim;
tm *ltmGlobFim;

int logRede = 0;
string nomeRedePrincipal;

struct fonteposic {
    int nmani;
    int mani[10];
    double posicmani;
};

struct noRede {
    int naflu = 0;
    int ncole = 0;
    int aflu[40];
    int cole[40];
    double normaP1 = 1000.;
    double normaP0 = 1000.;
    double normaMass1 = 1000.;
    double normaMass0 = 1000.;
    int cadastrado = 0;
};

struct tramoAtivo {
    int ncoleta;
    int nafluente;
    int nbloqueio;
    int ativo;
    std::vector<int> coleta;
    std::vector<int> afluente;
    int bloqueio;
    int permanente;
    int presimposta;
    double presMon;
    double presJus;
    int reverso;
    int derivaPrincipal;
    string tramoJson;
};
int redeLeitura = 0;

vector<tramoAtivo> tramos;
vector<tramoAtivo *> redesAlinhadas;

struct tramoPart {
    int ntram = 0;
    int part[40];
};

struct convergeNoPerm {
    double tL;
    double tH;
    double qostdTot;
    double qgstdTot;
    double qgstdTotNeg;
    double dengmist;
    double yco2mist;
    double apimist;
    double qwmist;
    double denagmist;
    double mliqmist;
    double mliqmistNeg;
    double mliqCmist;
    double mliqCmistNeg;
    double mgasmist;
    double mgasmistNeg;
    double moleomist;
    double moleomistNeg;
    double cpmist;
    double tempmist;
    double LVisL;
    double LVisH;
    double Qlmist;
    double betmist;
    double qlmistStd;
    double qlmistStdNeg;
    double qcmistStd;
    double qcmistStdNeg;
    double RGOmist;
    double bswmist;
    double moltot;
    double TRmist;
    ProFlu flu;
};

// criar relatorio dos arquivos de dados de saida da simulacao
ofstream arqRelatorioPerfis;

// criar objeto de logger
Logger logger("");

// criar string do path do arquivo de entrada
string pathArqEntrada("");

// criar string do path dos arquivos de entrada PVTSIM e SnapShot
string pathArqExtEntrada("");

// criar string do path e prefixo dos arquivos de saida para POCO_INJETOR
string pathPrefixoArqSaida("");

// criar string do arquivo de saida para o Snapshot, incluindo o path
string arqSaidaSnapShot("");

// diretorio contendo os arquivos de saida
string diretorioSaida("");

// criar ponteiro para a classe de sistema de producao
SProd *ptrSistemaProducao;

// criar ponteiro para a classe de Hidratos //chris
// FA_Hidrato solverHidrato(*ptrSistemaProducao);  // criação fora do loop

// definir o tipo de simulacao
tipoSimulacao_t tipoSimulacao = tipoSimulacao_t::transiente;

int match(Vcr<int> cadastrados, int ncadastro, int i) {
    int match = 0;
    for (int j = 0; j < ncadastro; j++)
        if (i == cadastrados[j])
            match = 1;
    return match;
}

int verficaConex(vector<tramoAtivo> &redeTeste, vector<int> &indice, Vcr<int> &cadastrados, int &ncadastro, int i, Vcr<int> &transporte) {
    int vmatch;
    vmatch = match(cadastrados, ncadastro, i);
    if (vmatch == 0) {
        redeTeste.push_back(tramos[i]);
        indice.push_back(i);
        cadastrados[ncadastro] = i;
        ncadastro++;
        int ref = redeTeste.size() - 1;
        transporte[i] = ref;
        for (int j = 0; j < tramos[i].nafluente; j++) {
            int nconex = tramos[i].afluente[j];
            verficaConex(redeTeste, indice, cadastrados, ncadastro, nconex, transporte);
        }
        for (int j = 0; j < tramos[i].ncoleta; j++) {
            int nconex = tramos[i].coleta[j];
            verficaConex(redeTeste, indice, cadastrados, ncadastro, nconex, transporte);
        }
    }
    return vmatch;
}

void descarteTramo(Rede &arqRede, int narq) {
    vector<tramoAtivo> redeAux;
    for (int i = 0; i < narq; i++) {
        tramoAtivo temp;
        temp.permanente = arqRede.malha[i].perm;
        temp.presimposta = arqRede.malha[i].presimposta;
        temp.tramoJson = arqRede.impfiles[i];
        temp.nafluente = arqRede.malha[i].nafluente;
        temp.presMon = arqRede.malha[i].presMon;
        temp.presJus = arqRede.malha[i].presJus;
        temp.reverso = arqRede.malha[i].reverso;
        temp.ativo = arqRede.malha[i].ativo;
        temp.derivaPrincipal = arqRede.malha[i].derivaPrincipal;
        if (temp.nafluente > 0) {
            if (arqRede.malha[i].ativo == 1) {
                for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
                    if (arqRede.malha[arqRede.malha[i].afluente[j]].ativo == 0) {
                        temp.nafluente--;
                    }
                }
                for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
                    if (arqRede.malha[arqRede.malha[i].afluente[j]].ativo == 1) {
                        temp.afluente.push_back(arqRede.malha[i].afluente[j]);
                    }
                }
            } else {
                temp.nafluente = 0;
            }
        }
        temp.ncoleta = arqRede.malha[i].ncoleta;
        if (temp.ncoleta > 0) {
            if (arqRede.malha[i].ativo == 1) {
                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                    if (arqRede.malha[arqRede.malha[i].coleta[j]].ativo == 0) {
                        temp.ncoleta--;
                    }
                }
                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                    if (arqRede.malha[arqRede.malha[i].coleta[j]].ativo == 1) {
                        temp.coleta.push_back(arqRede.malha[i].coleta[j]);
                    }
                }
            } else {
                temp.ncoleta = 0;
            }
        }
        if (arqRede.malha[i].ativo == 1) {
            temp.nbloqueio = arqRede.malha[i].nbloqueio;
            if (temp.nbloqueio > 0) {
                if (arqRede.malha[arqRede.malha[i].bloqueio].ativo == 1) {
                    temp.bloqueio = arqRede.malha[i].bloqueio;
                } else {
                    temp.nbloqueio = 0;
                    temp.bloqueio = -1;
                }
            } else
                temp.bloqueio = -1;
        } else {
            temp.nbloqueio = 0;
            temp.bloqueio = -1;
        }
        redeAux.push_back(temp);
    }
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].ativo == 1) {
            if (arqRede.malha[i].nafluente > 0)
                delete[] arqRede.malha[i].afluente;
            arqRede.malha[i].nafluente = redeAux[i].nafluente;
            if (arqRede.malha[i].nafluente > 0)
                arqRede.malha[i].afluente = new int[redeAux[i].nafluente];
            for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
                arqRede.malha[i].afluente[j] = redeAux[i].afluente[j];
            }
            if (arqRede.malha[i].ncoleta > 0)
                delete[] arqRede.malha[i].coleta;
            arqRede.malha[i].ncoleta = redeAux[i].ncoleta;
            if (arqRede.malha[i].ncoleta > 0)
                arqRede.malha[i].coleta = new int[redeAux[i].ncoleta];
            for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                arqRede.malha[i].coleta[j] = redeAux[i].coleta[j];
            }
            arqRede.malha[i].nbloqueio = redeAux[i].nbloqueio;
            arqRede.malha[i].bloqueio = redeAux[i].bloqueio;
        } else {
            arqRede.malha[i].perm = 0;
            arqRede.malha[i].presimposta = 0;
            if (arqRede.malha[i].nafluente > 0)
                delete[] arqRede.malha[i].afluente;
            arqRede.malha[i].nafluente = 0;
            if (arqRede.malha[i].ncoleta > 0)
                delete[] arqRede.malha[i].coleta;
            arqRede.malha[i].ncoleta = 0;
            arqRede.malha[i].bloqueio = -1;
            arqRede.malha[i].nbloqueio = 0;
        }
    }
}

void gravaRedeInterna(Rede &arqRede, int narq, int nredeAtiva, int nrede, int sizeRede, ostringstream &saidaT) {
    string tmp = saidaT.str();
    ofstream escreveIni(tmp.c_str(), ios_base::out);
    escreveIni << "{" << endl;
    escreveIni << "\"versao\":\"1.0\"," << endl;
    escreveIni << "\"configuracaoInicial\": {" << endl;
    escreveIni << "     " << "\"Rede\": true," << endl;
    escreveIni << "     " << "\"ParametroInicial\":" << arqRede.chutHol << "," << endl;
    escreveIni << "     " << "\"Relaxacao\":" << arqRede.relax << "," << endl;
    escreveIni << "     " << "\"limiteConvergencia\":" << arqRede.limConverge << "," << endl;
    if (arqRede.chute == 0)
        escreveIni << "     " << "\"ChuteNos\": false," << endl;
    else
        escreveIni << "     " << "\"ChuteNos\": true," << endl;
    if (arqRede.chaveredeT == 0) {
        escreveIni << "     " << "\"Transiente\": false," << endl;
        if (redesAlinhadas[nrede - 1][0].ativo == 0)
            escreveIni << "     " << "\"Transiente Falso por tramo inativo\": false," << endl;
    } else {
        if (redesAlinhadas[nrede - 1][0].ativo == 1)
            escreveIni << "     " << "\"Transiente\": true," << endl;
        else {
            escreveIni << "     " << "\"Transiente\": false," << endl;
            escreveIni << "     " << "\"Transiente Falso por tramo inativo\": false," << endl;
        }
    }
    escreveIni << "     " << "\"fluidoRede\":" << arqRede.fluidoRede << "," << endl;
    if (arqRede.tabelaDinamica == 0)
        escreveIni << "     " << "\"modeloTabelaDinamica\": false," << endl;
    else
        escreveIni << "     " << "\"modeloTabelaDinamica\": true," << endl;
    escreveIni << "     " << "\"TempoSimulacao\":" << arqRede.TmaxR << endl;
    escreveIni << " }," << endl;
    escreveIni << "\"Arquivos\": [" << endl;
    escreveIni << "     ";
    for (int i = 0; i < sizeRede; i++) {
        if (i < sizeRede - 1)
            escreveIni << "\"" << redesAlinhadas[nrede - 1][i].tramoJson << "\"" << ",";
        else
            escreveIni << "\"" << redesAlinhadas[nrede - 1][i].tramoJson << "\"" << endl;
    }
    escreveIni << " ]," << endl;
    escreveIni << "\"Conexao\": [" << endl;
    for (int i = 0; i < sizeRede; i++) {
        escreveIni << "{" << endl;
        escreveIni << "     " << "\"sisprod\": " << i << "," << endl;
        if (redesAlinhadas[nrede - 1][i].presimposta == 0)
            escreveIni << "     " << "\"PressaoImposta\": " << "false" << "," << endl;
        else
            escreveIni << "     " << "\"PressaoImposta\": " << "true" << "," << endl;
        escreveIni << "     " << "\"permanente\": " << redesAlinhadas[nrede - 1][i].permanente << "," << endl;
        escreveIni << "     " << "\"PressaoMontante\": " << redesAlinhadas[nrede - 1][i].presMon << "," << endl;
        escreveIni << "     " << "\"PressaoJusante\": " << redesAlinhadas[nrede - 1][i].presJus << "," << endl;
        escreveIni << "     " << "\"reverso\": " << redesAlinhadas[nrede - 1][i].reverso << "," << endl;
        escreveIni << "     " << "\"derivaPrincipal\": " << redesAlinhadas[nrede - 1][i].derivaPrincipal << "," << endl;
        escreveIni << "     " << "\"ncoletor\": " << redesAlinhadas[nrede - 1][i].ncoleta << "," << endl;
        escreveIni << "     " << "\"nafluente\": " << redesAlinhadas[nrede - 1][i].nafluente << "," << endl;
        escreveIni << "     " << "\"coletores\": [";
        for (int j = 0; j < redesAlinhadas[nrede - 1][i].ncoleta; j++) {
            if (j < redesAlinhadas[nrede - 1][i].ncoleta - 1)
                escreveIni << redesAlinhadas[nrede - 1][i].coleta[j] << ",";
            else
                escreveIni << redesAlinhadas[nrede - 1][i].coleta[j];
        }
        escreveIni << "]," << endl;
        escreveIni << "     " << "\"afluentes\": [";
        for (int j = 0; j < redesAlinhadas[nrede - 1][i].nafluente; j++) {
            if (j < redesAlinhadas[nrede - 1][i].nafluente - 1)
                escreveIni << redesAlinhadas[nrede - 1][i].afluente[j] << ",";
            else
                escreveIni << redesAlinhadas[nrede - 1][i].afluente[j];
        }
        escreveIni << "]," << endl;
        escreveIni << "     " << "\"bloqueio\": ";
        escreveIni << redesAlinhadas[nrede - 1][i].bloqueio;
        escreveIni << " " << endl;
        if (i < sizeRede - 1)
            escreveIni << "    }," << endl;
        else
            escreveIni << "    }" << endl;
    }
    escreveIni << "  ]" << endl;
    escreveIni << " }" << endl;
    escreveIni.close();
}

void preProcRede(Rede &arqRede, int narq) {

    for (int i = 0; i < narq; i++) {
        tramoAtivo temp;
        temp.permanente = arqRede.malha[i].perm;
        temp.presimposta = arqRede.malha[i].presimposta;
        temp.tramoJson = arqRede.impfiles[i];
        temp.nafluente = arqRede.malha[i].nafluente;
        temp.presMon = arqRede.malha[i].presMon;
        temp.presJus = arqRede.malha[i].presJus;
        temp.reverso = arqRede.malha[i].reverso;
        temp.ativo = arqRede.malha[i].ativo;
        temp.derivaPrincipal = arqRede.malha[i].derivaPrincipal;
        if (temp.nafluente > 0) {
            for (int j = 0; j < temp.nafluente; j++)
                temp.afluente.push_back(arqRede.malha[i].afluente[j]);
        }
        temp.ncoleta = arqRede.malha[i].ncoleta;
        if (temp.ncoleta > 0) {
            for (int j = 0; j < temp.ncoleta; j++)
                temp.coleta.push_back(arqRede.malha[i].coleta[j]);
        }
        temp.nbloqueio = arqRede.malha[i].nbloqueio;
        if (temp.nbloqueio > 0) {
            temp.bloqueio = arqRede.malha[i].bloqueio;
        } else
            temp.bloqueio = -1;
        tramos.push_back(temp);
    }
    int nativ = tramos.size();
    int nrede = 0;
    int nredeAtiva = 0;
    int partida = 0;
    int ncadastro = 0;
    vector<int> dimRede;
    Vcr<int> transporte(nativ);
    if (nativ > 0) {
        Vcr<int> cadastrados(nativ, -1);
        while (ncadastro < nativ) {
            vector<tramoAtivo> redeTeste;
            vector<int> indice;
            verficaConex(redeTeste, indice, cadastrados, ncadastro, partida, transporte);
            int sizeRede = redeTeste.size();
            if (sizeRede > 0) {
                tramoAtivo *alinhamento;
                dimRede.push_back(sizeRede);
                alinhamento = new tramoAtivo[sizeRede];
                redesAlinhadas.push_back(alinhamento);
                for (int i = 0; i < sizeRede; i++) {
                    alinhamento[i].permanente = redeTeste[i].permanente;
                    alinhamento[i].presimposta = redeTeste[i].presimposta;
                    alinhamento[i].presMon = redeTeste[i].presMon;
                    alinhamento[i].presJus = redeTeste[i].presJus;
                    alinhamento[i].reverso = redeTeste[i].reverso;
                    alinhamento[i].tramoJson = redeTeste[i].tramoJson;
                    alinhamento[i].nafluente = redeTeste[i].nafluente;
                    alinhamento[i].ncoleta = redeTeste[i].ncoleta;
                    alinhamento[i].nbloqueio = redeTeste[i].nbloqueio;
                    alinhamento[i].ativo = redeTeste[i].ativo;
                    alinhamento[i].derivaPrincipal = redeTeste[i].derivaPrincipal;
                    for (int j = 0; j < alinhamento[i].nafluente; j++)
                        alinhamento[i].afluente.push_back(transporte[redeTeste[i].afluente[j]]);
                    for (int j = 0; j < alinhamento[i].ncoleta; j++)
                        alinhamento[i].coleta.push_back(transporte[redeTeste[i].coleta[j]]);
                    alinhamento[i].bloqueio = -1;
                    for (int j = 0; j < sizeRede; j++) {
                        if (indice[j] == redeTeste[i].bloqueio) {
                            alinhamento[i].bloqueio = j;
                            break;
                        }
                    }
                }
            }
            int vmatch = 1;
            int konta = 0;
            while (vmatch == 1 && konta < nativ) {
                for (int j = 0; j < ncadastro; j++) {
                    if (konta == cadastrados[j]) {
                        vmatch = 1;
                        break;
                    } else
                        vmatch = 0;
                }
                if (vmatch == 1)
                    konta++;
            }
            partida = konta;
            nrede++;
            if (redesAlinhadas[nrede - 1][0].ativo == 1) {
                nredeAtiva++;
                ostringstream saidaT;
                ostringstream saidaT2;
                saidaT << "RedeInterna-" << nredeAtiva - 1 << ".json";
                saidaT2 << pathPrefixoArqSaida << "RedeInterna-" << nredeAtiva - 1 << ".json";

                gravaRedeInterna(arqRede, narq, nredeAtiva, nrede, sizeRede, saidaT2);
            }
        }
        redeLeitura = nredeAtiva;
        for (int i = 0; i < nrede; i++) {
            delete[] redesAlinhadas[i];
        }
    } else
        NumError("rede sem elementos ativos");
}

void WriteSnapShot(SProd &sis, double porc = 10, int tramo = -1) {

    srand(time(NULL));
    int frase = rand() % 15;

    if (sis.config().HISEP == 0) {
        ostringstream saidaS;
        int numero = round(porc);
        if (tramo < 0)
            saidaS << pathPrefixoArqSaida << "Saida" << numero << "-" << arqSaidaSnapShot;
        else
            saidaS << pathPrefixoArqSaida << "Tramo" << tramo << "-" << "Saida" << numero << "-" << arqSaidaSnapShot;
        string tmp = saidaS.str();
        ofstream rst(tmp.c_str(), ios_base::out);
        rst.width(16);
        rst.precision(12);
        rst << "*******************************************************************************" << endl;
        rst << saidaTexto[frase] << endl;
        rst << saidaSubTexto[frase] << endl;
        rst << "*******************************************************************************" << endl;
        rst << endl;
        rst << endl;
        rst << endl;
        /*  0->tempR; 1-> temp; 2-> tempL; 3-> tempini; 4-> presL; 5-> presLini; 6-> pres; 7-> presini; 8-> presauxL;
           *  9-> presaux; 10-> presauxR; 11-> presR; 12-> ML; 13-> MC; 14-> MR; 15-> MRini; 16-> MliqiniR; 17-> Mliqini;
           *  18-> MliqiniL; 19-> alfL; 20-> alfLini; 21-> alfR; 22-> alfRini; 23-> alf; 24-> alfini; 25-> betL; 26-> betLini;
           *  27-> betR; 28-> betRini; 29-> bet; 30-> betini; 31-> betLI; 32-> betI; 33-> betRI; 34-> FW; 35-> FWini;
           *  36-> massfonteCH; 37-> term1L; 38-> term2L; 39-> term1; 40-> term2; 41-> term1R; 42-> term2R; 43-> c0;
           *  44-> ud; 45-> c0Spare; 46-> udSpare; 47-> transmassL; 48-> transmassR; 49-> DTransDt1; 50-> DTransDt0;
           *  51-> DTransDxR; 52-> DTransDxL; 53-> coefTransBet; 54-> CoefDTR; 55-> CoefDTL; 56-> fontedissolv; 57-> TMModel;
           *  58-> TMModelL; 59-> fontemassLL; 60-> fontemassGL; 61-> fontemassCL; 62-> fontemassLR; 63-> fontemassGR; 64-> fontemassCR;
           *  65-> dpB; 66-> potB; 67-> transic; 68-> arranjo; 69-> arranjoR; 70-> chutearranjo; 71-> QLL; 72-> QL; 73-> QLR;
           *  74-> QG; 75-> rpL; 76-> rpC; 77-> rpR; 78-> rcL; 79-> rcC; 80-> rcR; 81-> rpLi; 82-> rpCi; 83-> rpRi; 84-> rcLi;
           *  85-> rcCi; 86-> rcRi; 87-> perdaEstratL; 88-> perdaEstratG; 89-> cinematico; 90-> VTemper; 91-> estadoPig; 92-> dtPig;
           *  93-> alfPigE; 94-> alfPigEini; 95-> alfPigD; 96-> alfPigDini; 97-> betPigE; 98-> betPigEini; 99-> betPigD;
           *  100-> betPigDini; 101-> velPig; 102-> velPigini; 103-> razPig; 104-> razPigini; 105-> VolLeveST; 106-> VolPesaST;
           *  107-> VolAguaST; 108-> DmasschokeG; 109-> DmasschokeL; 110-> DmasschokeC;
           *  111-> MliqiniR0; 112-> presBuf; 113-> presLiniBuf; 114-> presLBuf; 115-> presRBuf;116->presauxBuf; 117-> MCBuf;
           *  118-> MRBuf; 119->MliqiniRBuf; 120-> MRiniBuf; 121-> MLBuf; 122-> MliqiniBuf; 123-> MliqiniLBuf; 124->QLRini;
           *  125->presRini ; 126->presauxini ; 127->MCini ; 128->MRini;
           *  129->MLini ; 130->Mliqini0; 131->MliqiniL0 ;
           *  132->tempLini ; 133->tempRini;
           *
                flui= API; RGO; Deng; BSW; Denag; TempL; LVisL; TempH; LVisH; tipoemul; bswCorte; yco2; corrC;
                calor.Vint; calor.Tint; calor.Tint2; calor.kint; calor.cpint; calor.rhoint; calor.viscint;
                calor.Tcamada[i][j]
                *fluiL;
                *fluiR;
                fluicol;

                npseudo ; *fracmol ;
                *oCalculatedLiqComposition;*oCalculatedVapComposition;dCalculatedBeta;dCalculatedBubbleP;
                iCalculatedThermodynamicCondition;iCalculatedStockTankThermodynamicCondition;dVaporMassFraction;
                dStockTankVaporMassFraction;dStockTankLiquidDensity;dStockTankVaporDensity;
                */
        rst << "variavel de Sistema de Escoamento= " << " -> ";
        rst << sis.presfim << " ; ";
        rst << sis.pGSup << " ; ";
        rst << sis.masChkSup << " ; ";
        rst << sis.temperatura << " ; ";
        rst << sis.ncelGas << " ; ";
        rst << sis.presiniG << " ; ";
        rst << sis.tempiniG << " ; ";
        rst << sis.massfonte << " ; ";
        rst << sis.mult << " ; ";
        rst << sis.presMedMov << " ; ";
        rst << sis.jMedMov << " ; ";
        rst << sis.alfMedMov << " ; ";
        rst << sis.tMedMov << " ; ";
        rst << sis.aberto << " ; ";
        rst << sis.tempoaberto << " ; ";
        rst << sis.ktMedMov << " ; ";
        rst << sis.pTotal << " ; ";
        rst << sis.jTotal << " ; ";
        rst << sis.alfTotal << " ; ";
        rst << sis.presVet.size() << " ; ";
        for (int i = 0; i < sis.presVet.size(); i++) {
            rst << sis.presVet[i] << " ; ";
            rst << sis.jVet[i] << " ; ";
            rst << sis.alfVet[i] << " ; ";
            rst << sis.tVet[i] << " ; ";
        }
        rst << sis.celInter << " ; ";
        rst << sis.dtInter << " ; ";
        rst << sis.velInter << " ; ";
        rst << sis.celInterIni << " ; ";
        rst << sis.dtInterIni << " ; ";
        rst << sis.velInterIni << " ; ";
        rst << sis.fontemassPRBuf << " ; ";
        rst << sis.fontemassCRBuf << " ; ";
        rst << sis.fontemassGRBuf << " ; ";
        rst << endl;
        rst << endl;
        rst << endl;
        rst << endl;

        int ncel = sis.ncel + 1;
        for (int i = 0; i < ncel; i++) {
            rst << "Celula de producao= " << i << " -> ";
            rst << sis.cell(i).tempR << " ; " << sis.cell(i).temp << " ; " << sis.cell(i).tempL << " ; " << sis.cell(i).tempini << " ; " << sis.cell(i).presL << " ; " << sis.cell(i).presLini << " ; " << sis.cell(i).pres << " ; " << sis.cell(i).presini << " ; " << sis.cell(i).presauxL << " ; " << sis.cell(i).presaux << " ; " << sis.cell(i).presauxR << " ; " << sis.cell(i).presR << " ; " << sis.cell(i).ML << " ; " << sis.cell(i).MC << " ; " << sis.cell(i).MR << " ; " << sis.cell(i).MRini << " ; " << sis.cell(i).MliqiniR << " ; " << sis.cell(i).Mliqini << " ; " << sis.cell(i).MliqiniL << " ; " << sis.cell(i).alfL << " ; " << sis.cell(i).alfLini << " ; " << sis.cell(i).alfR << " ; " << sis.cell(i).alfRini << " ; " << sis.cell(i).alf << " ; " << sis.cell(i).alfini << " ; " << sis.cell(i).betL << " ; " << sis.cell(i).betLini << " ; " << sis.cell(i).betR << " ; " << sis.cell(i).betRini << " ; " << sis.cell(i).bet << " ; " << sis.cell(i).betini << " ; " << sis.cell(i).betLI << " ; " << sis.cell(i).betI << " ; " << sis.cell(i).betRI << " ; " << sis.cell(i).FW << " ; " << sis.cell(i).FWini << " ; " << sis.cell(i).massfonteCH << " ; " << sis.cell(i).term1L << " ; " << sis.cell(i).term2L << " ; " << sis.cell(i).term1 << " ; " << sis.cell(i).term2 << " ; " << sis.cell(i).term1R << " ; " << sis.cell(i).term2R << " ; " << sis.cell(i).c0 << " ; " << sis.cell(i).ud << " ; " << sis.cell(i).c0Spare << " ; " << sis.cell(i).udSpare << " ; " << sis.cell(i).transmassL << " ; " << sis.cell(i).transmassR << " ; " << sis.cell(i).DTransDt1 << " ; " << sis.cell(i).DTransDt0 << " ; " << sis.cell(i).DTransDxR << " ; " << sis.cell(i).DTransDxL << " ; " << sis.cell(i).coefTransBet << " ; " << sis.cell(i).CoefDTR << " ; " << sis.cell(i).CoefDTL << " ; " << sis.cell(i).fontedissolv << " ; " << sis.cell(i).TMModel << " ; " << sis.cell(i).TMModelL << " ; " << sis.cell(i).fontemassLL << " ; " << sis.cell(i).fontemassGL << " ; " << sis.cell(i).fontemassCL << " ; " << sis.cell(i).fontemassLR << " ; " << sis.cell(i).fontemassGR << " ; " << sis.cell(i).fontemassCR << " ; " << sis.cell(i).dpB << " ; " << sis.cell(i).potB << " ; " << sis.cell(i).transic << " ; " << sis.cell(i).arranjo << " ; " << sis.cell(i).arranjoR << " ; " << sis.cell(i).chutearranjo << " ; " << sis.cell(i).QLL << " ; " << sis.cell(i).QL << " ; " << sis.cell(i).QLR << " ; " << sis.cell(i).QG << " ; " << sis.cell(i).rpL << " ; " << sis.cell(i).rpC << " ; " << sis.cell(i).rpR << " ; " << sis.cell(i).rcL << " ; " << sis.cell(i).rcC << " ; " << sis.cell(i).rcR << " ; " << sis.cell(i).rpLi << " ; " << sis.cell(i).rpCi << " ; " << sis.cell(i).rpRi << " ; " << sis.cell(i).rcLi << " ; " << sis.cell(i).rcCi << " ; " << sis.cell(i).rcRi << " ; " << sis.cell(i).perdaEstratL << " ; " << sis.cell(i).perdaEstratG << " ; " << sis.cell(i).cinematico << " ; " << sis.cell(i).VTemper << " ; " << sis.cell(i).estadoPig << " ; " << sis.cell(i).dtPig << " ; " << sis.cell(i).alfPigE << " ; " << sis.cell(i).alfPigEini << " ; " << sis.cell(i).alfPigD << " ; " << sis.cell(i).alfPigDini << " ; " << sis.cell(i).betPigE << " ; " << sis.cell(i).betPigEini << " ; " << sis.cell(i).betPigD << " ; " << sis.cell(i).betPigDini << " ; " << sis.cell(i).velPig << " ; " << sis.cell(i).velPigini << " ; " << sis.cell(i).razPig << " ; " << sis.cell(i).razPigini << " ; " << sis.cell(i).VolLeveST << " ; " << sis.cell(i).VolPesaST << " ; " << sis.cell(i).VolAguaST << " ; " << sis.cell(i).DmasschokeG << " ; " << sis.cell(i).DmasschokeL << " ; " << sis.cell(i).DmasschokeC << " ; " << sis.cell(i).MliqiniR0 << " ; " << sis.cell(i).presBuf << " ; " << sis.cell(i).presLiniBuf << " ; " << sis.cell(i).presLBuf << " ; " << sis.cell(i).presRBuf << " ; " << sis.cell(i).presauxBuf << " ; " << sis.cell(i).MCBuf << " ; " << sis.cell(i).MRBuf << " ; " << sis.cell(i).MliqiniRBuf << " ; " << sis.cell(i).MRiniBuf << " ; " << sis.cell(i).MLBuf << " ; " << sis.cell(i).MliqiniBuf << " ; " << sis.cell(i).MliqiniLBuf << " ; " << sis.cell(i).QLRini << " ; " << sis.cell(i).presRini << " ; " << sis.cell(i).presauxini << " ; " << sis.cell(i).MCini << " ; " << sis.cell(i).MRini << " ; " << sis.cell(i).MLini << " ; " << sis.cell(i).Mliqini0 << " ; " << sis.cell(i).MliqiniL0 << " ; " << sis.cell(i).tempLini << " ; " << sis.cell(i).tempRini << " ; " << sis.cell(i).rgC << " ; " << sis.cell(i).rgL << " ; " << sis.cell(i).rgR << " ; " << sis.cell(i).rgCi << " ; " << sis.cell(i).rgLi << " ; " << sis.cell(i).rgRi << " ; " << sis.cell(i).flui.API << " ; " << sis.cell(i).flui.RGO << " ; " << sis.cell(i).flui.Deng << " ; " << sis.cell(i).flui.BSW << " ; " << sis.cell(i).flui.TempL << " ; " << sis.cell(i).flui.LVisL << " ; " << sis.cell(i).flui.TempH << " ; " << sis.cell(i).flui.LVisH << " ; " << sis.cell(i).flui.tipoemul << " ; " << sis.cell(i).flui.bswCorte << " ; " << sis.cell(i).flui.yco2 << " ; " << sis.cell(i).flui.corrC << " ; ";
            if (sis.config().modoParafina == 1) {
                rst << sis.cell(i).parafinado << " ; " << sis.cell(i).duto.espessuR[0] << " ; " << sis.cell(i).duto.cond[0] << " ; " << sis.cell(i).duto.cp[0] << " ; " << sis.cell(i).duto.rhoC[0] << " ; " << sis.cell(i).deltaPar << " ; ";
            }
            rst << sis.cell(i).calor.Vint << " ; " << sis.cell(i).calor.Tint << " ; " << sis.cell(i).calor.Tint2 << " ; " << sis.cell(i).calor.kint << " ; " << sis.cell(i).calor.cpint << " ; " << sis.cell(i).calor.rhoint << " ; " << sis.cell(i).calor.viscint << " ; " << sis.cell(i).calor.Vconf << " ; " << sis.cell(i).calor.fluxIni << " ; " << sis.cell(i).calor.fluxFim << " ; " << sis.cell(i).calor.ccon << " ; " << sis.cell(i).calor.ncon << " ; " << sis.cell(i).calor.mcon << " ; " << sis.cell(i).calor.npet << " ; " << sis.cell(i).calor.betint << " ; " << sis.cell(i).calor.betext << " ; " << sis.cell(i).calor.reyi << " ; " << sis.cell(i).calor.reye << " ; " << sis.cell(i).calor.grashi << " ; " << sis.cell(i).calor.grashe << " ; " << sis.cell(i).calor.nusi << " ; " << sis.cell(i).calor.nuse << " ; " << sis.cell(i).calor.pri << " ; " << sis.cell(i).calor.pre << " ; " << sis.cell(i).calor.hi << " ; " << sis.cell(i).calor.he;

            for (int j = 0; j < sis.cell(i).calor.geom.ncamadas; j++) {
                for (int k = 0; k <= sis.cell(i).calor.ncamada[j]; k++)
                    rst << " ; " << sis.cell(i).calor.Tcamada[j][k];
            }

            if (sis.config().flashCompleto == 2) {
                rst << " ; " << sis.cell(i).flui.npseudo;
                for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                    rst << " ; " << sis.cell(i).flui.fracMol[j];
                }
                for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                    rst << " ; " << sis.cell(i).flui.oCalculatedLiqComposition[j];
                }
                for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                    rst << " ; " << sis.cell(i).flui.oCalculatedVapComposition[j];
                }
                rst << " ; " << sis.cell(i).flui.dCalculatedBeta << " ; " << sis.cell(i).flui.dCalculatedBubbleP << " ; " << sis.cell(i).flui.iCalculatedThermodynamicCondition << " ; " << sis.cell(i).flui.iCalculatedStockTankThermodynamicCondition << " ; " << sis.cell(i).flui.dVaporMassFraction << " ; " << sis.cell(i).flui.dStockTankVaporMassFraction << " ; " << sis.cell(i).flui.dStockTankLiquidDensity << " ; " << sis.cell(i).flui.dStockTankVaporDensity;
            }

            rst << endl;
        }

        if (sis.config().lingas > 0) {
            rst << endl;
            rst << endl;
            rst << endl;
            rst << endl;
            ncel = sis.ncelGas + 1;
            for (int i = 0; i < ncel; i++) {
                rst << "Celula de servico= " << i << " -> ";
                rst << sis.gasCell(i).tempL << " ; " << sis.gasCell(i).temp << " ; " << sis.gasCell(i).tempR << " ; " << sis.gasCell(i).presL << " ; " << sis.gasCell(i).pres << " ; " << sis.gasCell(i).presini << " ; " << sis.gasCell(i).presR << " ; " << sis.gasCell(i).VGasL << " ; " << sis.gasCell(i).VGasR << " ; " << sis.gasCell(i).VGasRR << " ; " << sis.gasCell(i).u1LL << " ; " << sis.gasCell(i).u1L << " ; " << sis.gasCell(i).u1R << " ; " << sis.gasCell(i).massfonteCH << " ; " << sis.gasCell(i).fluxcal << " ; " << sis.gasCell(i).labelchk << " ; " << sis.gasCell(i).fechamon << " ; " << sis.gasCell(i).rpchk << " ; " << sis.gasCell(i).fonteM2 << " ; " << sis.gasCell(i).salinidade << " ; " << sis.gasCell(i).razInter << " ; " << sis.gasCell(i).razInterIni << " ; " << sis.gasCell(i).tempLini << " ; " << sis.gasCell(i).tempini << " ; " << sis.gasCell(i).tempRini << " ; " << sis.gasCell(i).presLini << " ; " << sis.gasCell(i).presini << " ; " << sis.gasCell(i).presRini << " ; " << sis.gasCell(i).VGasLini << " ; " << sis.gasCell(i).VGasRini << " ; " << sis.gasCell(i).VGasRRini << " ; " << sis.gasCell(i).u1LLini << " ; " << sis.gasCell(i).u1Lini << " ; " << sis.gasCell(i).u1Rini << " ; " << sis.gasCell(i).massfonteCHini << " ; " << sis.gasCell(i).fonteM2ini << " ; " << sis.gasCell(i).fechamonini << " ; " << sis.gasCell(i).posicini << " ; " << sis.gasCell(i).rpchkini << " ; " << sis.gasCell(i).calor.Vint << " ; " << sis.gasCell(i).calor.Tint << " ; " << sis.gasCell(i).calor.Tint2 << " ; " << sis.gasCell(i).calor.kint << " ; " << sis.gasCell(i).calor.cpint << " ; " << sis.gasCell(i).calor.rhoint << " ; " << sis.gasCell(i).calor.viscint << " ; " << sis.gasCell(i).calor.Vconf << " ; " << sis.gasCell(i).calor.fluxIni << " ; " << sis.gasCell(i).calor.fluxFim << " ; " << sis.gasCell(i).calor.ccon << " ; " << sis.gasCell(i).calor.ncon << " ; " << sis.gasCell(i).calor.mcon << " ; " << sis.gasCell(i).calor.npet << " ; " << sis.gasCell(i).calor.betext << " ; " << sis.gasCell(i).calor.betint << " ; " << sis.gasCell(i).calor.reyi << " ; " << sis.gasCell(i).calor.reye << " ; " << sis.gasCell(i).calor.grashi << " ; " << sis.gasCell(i).calor.grashe << " ; " << sis.gasCell(i).calor.nusi << " ; " << sis.gasCell(i).calor.nuse << " ; " << sis.gasCell(i).calor.pri << " ; " << sis.gasCell(i).calor.pre << " ; " << sis.gasCell(i).calor.hi << " ; " << sis.gasCell(i).calor.he;
                for (int j = 0; j < sis.gasCell(i).calor.geom.ncamadas; j++) {
                    for (int k = 0; k <= sis.gasCell(i).calor.ncamada[j]; k++)
                        rst << " ; " << sis.gasCell(i).calor.Tcamada[j][k];
                }
                rst << endl;
            }
        }

        rst << endl;
        rst << " Tempo =  " << (*sis.vg1dSP).lixo5;
        rst.close();

    } else {

        ostringstream saidaS;
        int numero = round(porc);
        if (tramo < 0)
            saidaS << pathPrefixoArqSaida << "Saida" << numero << "-" << arqSaidaSnapShot;
        else
            saidaS << pathPrefixoArqSaida << "Tramo" << tramo << "-" << "Saida" << numero << "-" << arqSaidaSnapShot;
        string tmp = saidaS.str();
        ofstream rst(tmp.c_str(), ios_base::out);
        rst.width(16);
        rst.precision(12);

        // t [s] \t H [] \t slipRatio [] \t Usg [m/s] \t Usl [m/s]

        int ncel = sis.ncel + 1;                                                                                                                                               // slipRatio=vg/vl
        rst << " t [s] " << ";" << " HL [] " << ";" << " slipRatio [] " << ";" << " Usg [m/s] " << ";" << " Usl [m/s] " << ";" << " Flow Pattern " << ";" << " cell " << "\n"; //<< i
        for (int i = 1; i < ncel; i++) {                                                                                                                                       // i=0 não possui velocidades (condição de parede fechada)
            // rst << " t [s] " << " HL [] " << " slipRatio [] " << " Usg [m/s] " << " Usl [m/s]" << " cell " << "\n"; //<< i
            rst << (*sis.vg1dSP).lixo5 << " ; " << 1 - sis.cell(i).alf << " ; " << (sis.cell(i).QG / sis.cell(i).duto.area * (sis.cell(i).alf + 0.001)) / (sis.cell(i).QL / (sis.cell(i).duto.area * ((1 - sis.cell(i).alf) + 0.001))) << " ; " << sis.cell(i).QG / sis.cell(i).duto.area << " ; " << sis.cell(i).QL / sis.cell(i).duto.area << " ; " << sis.cell(i).arranjo << " ; " << i << " ; ";

            rst << "\n";
        }

        rst.close();
    }
}

void ReadSnapShot(SProd &sis) {
    string dadosMR = sis.config().snapshotArqIn;
    ifstream rst(dadosMR.c_str(), ios_base::in);
    string chave;
    rst >> chave;
    while (chave != "->")
        rst >> chave;
    rst >> sis.presfim;
    rst >> chave;
    rst >> sis.pGSup;
    rst >> chave;
    rst >> sis.masChkSup;
    rst >> chave;
    rst >> sis.temperatura;
    rst >> chave;
    rst >> sis.ncelGas;
    rst >> chave;
    rst >> sis.presiniG;
    rst >> chave;
    rst >> sis.tempiniG;
    rst >> chave;
    rst >> sis.massfonte;
    rst >> chave;
    rst >> sis.mult;
    rst >> chave;
    rst >> sis.presMedMov;
    rst >> chave;
    rst >> sis.jMedMov;
    rst >> chave;
    rst >> sis.alfMedMov;
    rst >> chave;
    rst >> sis.tMedMov;
    rst >> chave;
    rst >> sis.aberto;
    rst >> chave;
    rst >> sis.tempoaberto;
    rst >> chave;
    rst >> sis.ktMedMov;
    rst >> chave;
    rst >> sis.pTotal;
    rst >> chave;
    rst >> sis.jTotal;
    rst >> chave;
    rst >> sis.alfTotal;
    rst >> chave;
    int nvet;
    rst >> nvet;
    rst >> chave;
    for (int i = 0; i < nvet; i++) {
        double temp;
        rst >> temp;
        rst >> chave;
        sis.presVet.push_back(temp);
        rst >> temp;
        rst >> chave;
        sis.jVet.push_back(temp);
        rst >> temp;
        rst >> chave;
        sis.alfVet.push_back(temp);
        rst >> temp;
        rst >> chave;
        sis.tVet.push_back(temp);
    }
    rst >> sis.celInter;
    rst >> chave;
    rst >> sis.dtInter;
    rst >> chave;
    rst >> sis.velInter;
    rst >> chave;
    rst >> sis.celInterIni;
    rst >> chave;
    rst >> sis.dtInterIni;
    rst >> chave;
    rst >> sis.velInterIni;
    rst >> chave;
    rst >> sis.fontemassPRBuf;
    rst >> chave;
    rst >> sis.fontemassCRBuf;
    rst >> chave;
    rst >> sis.fontemassGRBuf;
    rst >> chave;

    int ncel = sis.ncel + 1;
    for (int i = 0; i < ncel; i++) {
        rst >> chave;
        while (chave != "->")
            rst >> chave;
        rst >> sis.cell(i).tempR;
        rst >> chave;
        rst >> sis.cell(i).temp;
        rst >> chave;
        rst >> sis.cell(i).tempL;
        rst >> chave;
        rst >> sis.cell(i).tempini;
        rst >> chave;
        rst >> sis.cell(i).presL;
        rst >> chave;
        rst >> sis.cell(i).presLini;
        rst >> chave;
        rst >> sis.cell(i).pres;
        rst >> chave;
        rst >> sis.cell(i).presini;
        rst >> chave;
        rst >> sis.cell(i).presauxL;
        rst >> chave;
        rst >> sis.cell(i).presaux;
        rst >> chave;
        rst >> sis.cell(i).presauxR;
        rst >> chave;
        rst >> sis.cell(i).presR;
        rst >> chave;
        rst >> sis.cell(i).ML;
        rst >> chave;
        rst >> sis.cell(i).MC;
        rst >> chave;
        rst >> sis.cell(i).MR;
        rst >> chave;
        rst >> sis.cell(i).MRini;
        rst >> chave;
        rst >> sis.cell(i).MliqiniR;
        rst >> chave;
        rst >> sis.cell(i).Mliqini;
        rst >> chave;
        rst >> sis.cell(i).MliqiniL;
        rst >> chave;
        rst >> sis.cell(i).alfL;
        rst >> chave;
        rst >> sis.cell(i).alfLini;
        rst >> chave;
        rst >> sis.cell(i).alfR;
        rst >> chave;
        rst >> sis.cell(i).alfRini;
        rst >> chave;
        rst >> sis.cell(i).alf;
        rst >> chave;
        rst >> sis.cell(i).alfini;
        rst >> chave;
        rst >> sis.cell(i).betL;
        rst >> chave;
        rst >> sis.cell(i).betLini;
        rst >> chave;
        rst >> sis.cell(i).betR;
        rst >> chave;
        rst >> sis.cell(i).betRini;
        rst >> chave;
        rst >> sis.cell(i).bet;
        rst >> chave;
        rst >> sis.cell(i).betini;
        rst >> chave;
        rst >> sis.cell(i).betLI;
        rst >> chave;
        rst >> sis.cell(i).betI;
        rst >> chave;
        rst >> sis.cell(i).betRI;
        rst >> chave;
        rst >> sis.cell(i).FW;
        rst >> chave;
        rst >> sis.cell(i).FWini;
        rst >> chave;
        rst >> sis.cell(i).massfonteCH;
        rst >> chave;
        rst >> sis.cell(i).term1L;
        rst >> chave;
        rst >> sis.cell(i).term2L;
        rst >> chave;
        rst >> sis.cell(i).term1;
        rst >> chave;
        rst >> sis.cell(i).term2;
        rst >> chave;
        rst >> sis.cell(i).term1R;
        rst >> chave;
        rst >> sis.cell(i).term2R;
        rst >> chave;
        rst >> sis.cell(i).c0;
        rst >> chave;
        rst >> sis.cell(i).ud;
        rst >> chave;
        rst >> sis.cell(i).c0Spare;
        rst >> chave;
        rst >> sis.cell(i).udSpare;
        rst >> chave;
        rst >> sis.cell(i).transmassL;
        rst >> chave;
        rst >> sis.cell(i).transmassR;
        rst >> chave;
        rst >> sis.cell(i).DTransDt1;
        rst >> chave;
        rst >> sis.cell(i).DTransDt0;
        rst >> chave;
        rst >> sis.cell(i).DTransDxR;
        rst >> chave;
        rst >> sis.cell(i).DTransDxL;
        rst >> chave;
        rst >> sis.cell(i).coefTransBet;
        rst >> chave;
        rst >> sis.cell(i).CoefDTR;
        rst >> chave;
        rst >> sis.cell(i).CoefDTL;
        rst >> chave;
        rst >> sis.cell(i).fontedissolv;
        rst >> chave;
        rst >> sis.cell(i).TMModel;
        rst >> chave;
        rst >> sis.cell(i).TMModelL;
        rst >> chave;
        rst >> sis.cell(i).fontemassLL;
        rst >> chave;
        rst >> sis.cell(i).fontemassGL;
        rst >> chave;
        rst >> sis.cell(i).fontemassCL;
        rst >> chave;
        rst >> sis.cell(i).fontemassLR;
        rst >> chave;
        rst >> sis.cell(i).fontemassGR;
        rst >> chave;
        rst >> sis.cell(i).fontemassCR;
        rst >> chave;
        rst >> sis.cell(i).dpB;
        rst >> chave;
        rst >> sis.cell(i).potB;
        rst >> chave;
        rst >> sis.cell(i).transic;
        rst >> chave;
        rst >> sis.cell(i).arranjo;
        rst >> chave;
        rst >> sis.cell(i).arranjoR;
        rst >> chave;
        rst >> sis.cell(i).chutearranjo;
        rst >> chave;
        rst >> sis.cell(i).QLL;
        rst >> chave;
        rst >> sis.cell(i).QL;
        rst >> chave;
        rst >> sis.cell(i).QLR;
        rst >> chave;
        rst >> sis.cell(i).QG;
        rst >> chave;
        rst >> sis.cell(i).rpL;
        rst >> chave;
        rst >> sis.cell(i).rpC;
        rst >> chave;
        rst >> sis.cell(i).rpR;
        rst >> chave;
        rst >> sis.cell(i).rcL;
        rst >> chave;
        rst >> sis.cell(i).rcC;
        rst >> chave;
        rst >> sis.cell(i).rcR;
        rst >> chave;
        rst >> sis.cell(i).rpLi;
        rst >> chave;
        rst >> sis.cell(i).rpCi;
        rst >> chave;
        rst >> sis.cell(i).rpRi;
        rst >> chave;
        rst >> sis.cell(i).rcLi;
        rst >> chave;
        rst >> sis.cell(i).rcCi;
        rst >> chave;
        rst >> sis.cell(i).rcRi;
        rst >> chave;
        rst >> sis.cell(i).perdaEstratL;
        rst >> chave;
        rst >> sis.cell(i).perdaEstratG;
        rst >> chave;
        rst >> sis.cell(i).cinematico;
        rst >> chave;
        rst >> sis.cell(i).VTemper;
        rst >> chave;
        rst >> sis.cell(i).estadoPig;
        rst >> chave;
        rst >> sis.cell(i).dtPig;
        rst >> chave;
        rst >> sis.cell(i).alfPigE;
        rst >> chave;
        rst >> sis.cell(i).alfPigEini;
        rst >> chave;
        rst >> sis.cell(i).alfPigD;
        rst >> chave;
        rst >> sis.cell(i).alfPigDini;
        rst >> chave;
        rst >> sis.cell(i).betPigE;
        rst >> chave;
        rst >> sis.cell(i).betPigEini;
        rst >> chave;
        rst >> sis.cell(i).betPigD;
        rst >> chave;
        rst >> sis.cell(i).betPigDini;
        rst >> chave;
        rst >> sis.cell(i).velPig;
        rst >> chave;
        rst >> sis.cell(i).velPigini;
        rst >> chave;
        rst >> sis.cell(i).razPig;
        rst >> chave;
        rst >> sis.cell(i).razPigini;
        rst >> chave;
        rst >> sis.cell(i).VolLeveST;
        rst >> chave;
        rst >> sis.cell(i).VolPesaST;
        rst >> chave;
        rst >> sis.cell(i).VolAguaST;
        rst >> chave;
        rst >> sis.cell(i).DmasschokeG;
        rst >> chave;
        rst >> sis.cell(i).DmasschokeL;
        rst >> chave;
        rst >> sis.cell(i).DmasschokeC;
        rst >> chave;
        rst >> sis.cell(i).MliqiniR0;
        rst >> chave;
        rst >> sis.cell(i).presBuf;
        rst >> chave;
        rst >> sis.cell(i).presLiniBuf;
        rst >> chave;
        rst >> sis.cell(i).presLBuf;
        rst >> chave;
        rst >> sis.cell(i).presRBuf;
        rst >> chave;
        rst >> sis.cell(i).presauxBuf;
        rst >> chave;
        rst >> sis.cell(i).MCBuf;
        rst >> chave;
        rst >> sis.cell(i).MRBuf;
        rst >> chave;
        rst >> sis.cell(i).MliqiniRBuf;
        rst >> chave;
        rst >> sis.cell(i).MRiniBuf;
        rst >> chave;
        rst >> sis.cell(i).MLBuf;
        rst >> chave;
        rst >> sis.cell(i).MliqiniBuf;
        rst >> chave;
        rst >> sis.cell(i).MliqiniLBuf;
        rst >> chave;
        rst >> sis.cell(i).QLRini;
        rst >> chave;
        rst >> sis.cell(i).presRini;
        rst >> chave;
        rst >> sis.cell(i).presauxini;
        rst >> chave;
        rst >> sis.cell(i).MCini;
        rst >> chave;
        rst >> sis.cell(i).MRini;
        rst >> chave;
        rst >> sis.cell(i).MLini;
        rst >> chave;
        rst >> sis.cell(i).Mliqini0;
        rst >> chave;
        rst >> sis.cell(i).MliqiniL0;
        rst >> chave;
        rst >> sis.cell(i).tempLini;
        rst >> chave;
        rst >> sis.cell(i).tempRini;
        rst >> chave;
        rst >> sis.cell(i).rgC;
        rst >> chave;
        rst >> sis.cell(i).rgL;
        rst >> chave;
        rst >> sis.cell(i).rgR;
        rst >> chave;
        rst >> sis.cell(i).rgCi;
        rst >> chave;
        rst >> sis.cell(i).rgLi;
        rst >> chave;
        rst >> sis.cell(i).rgRi;
        rst >> chave;
        rst >> sis.cell(i).flui.API;
        rst >> chave;
        rst >> sis.cell(i).flui.RGO;
        rst >> chave;
        rst >> sis.cell(i).flui.Deng;
        rst >> chave;
        rst >> sis.cell(i).flui.BSW;
        rst >> chave;
        rst >> sis.cell(i).flui.TempL;
        rst >> chave;
        rst >> sis.cell(i).flui.LVisL;
        rst >> chave;
        rst >> sis.cell(i).flui.TempH;
        rst >> chave;
        rst >> sis.cell(i).flui.LVisH;
        rst >> chave;
        rst >> sis.cell(i).flui.tipoemul;
        rst >> chave;
        rst >> sis.cell(i).flui.bswCorte;
        rst >> chave;
        rst >> sis.cell(i).flui.yco2;
        rst >> chave;
        rst >> sis.cell(i).flui.corrC;
        rst >> chave;
        sis.cell(i).flui.RenovaFluido();
        if (i > 0)
            sis.cell(i - 1).fluiR = &(sis.cell(i).flui);
        if (i < ncel - 1)
            sis.cell(i + 1).fluiL = &(sis.cell(i).flui);
        if (sis.config().modoParafina == 1) {
            rst >> sis.cell(i).parafinado;
            rst >> chave;
            double espessura;
            double cpW;
            double kW;
            double rhoW;
            double rugosidade = 0.;
            if (sis.config().modoParafina == 1)
                rugosidade = sis.config().detalParafina.rug;
            rst >> espessura;
            rst >> chave;
            rst >> kW;
            rst >> chave;
            rst >> cpW;
            rst >> chave;
            rst >> rhoW;
            rst >> chave;
            rst >> sis.cell(i).deltaPar;
            if (sis.cell(i).parafinado == 1) {
                sis.cell(i).duto.atualizaCamada(espessura, rugosidade, cpW, kW, rhoW);
                sis.cell(i).calor.atualiza(sis.cell(i).duto, 1);
            }
        }
        rst >> sis.cell(i).calor.Vint;
        rst >> chave;
        rst >> sis.cell(i).calor.Tint;
        rst >> chave;
        rst >> sis.cell(i).calor.Tint2;
        rst >> chave;
        rst >> sis.cell(i).calor.kint;
        rst >> chave;
        rst >> sis.cell(i).calor.cpint;
        rst >> chave;
        rst >> sis.cell(i).calor.rhoint;
        rst >> chave;
        rst >> sis.cell(i).calor.viscint;
        rst >> chave;
        rst >> sis.cell(i).calor.Vconf;
        rst >> chave;
        rst >> sis.cell(i).calor.fluxIni;
        rst >> chave;
        rst >> sis.cell(i).calor.fluxFim;
        rst >> chave;
        rst >> sis.cell(i).calor.ccon;
        rst >> chave;
        rst >> sis.cell(i).calor.ncon;
        rst >> chave;
        rst >> sis.cell(i).calor.mcon;
        rst >> chave;
        rst >> sis.cell(i).calor.npet;
        rst >> chave;
        rst >> sis.cell(i).calor.betint;
        rst >> chave;
        rst >> sis.cell(i).calor.betext;
        rst >> chave;
        rst >> sis.cell(i).calor.reyi;
        rst >> chave;
        rst >> sis.cell(i).calor.reye;
        rst >> chave;
        rst >> sis.cell(i).calor.grashi;
        rst >> chave;
        rst >> sis.cell(i).calor.grashe;
        rst >> chave;
        rst >> sis.cell(i).calor.nusi;
        rst >> chave;
        rst >> sis.cell(i).calor.nuse;
        rst >> chave;
        rst >> sis.cell(i).calor.pri;
        rst >> chave;
        rst >> sis.cell(i).calor.pre;
        rst >> chave;
        rst >> sis.cell(i).calor.hi;
        rst >> chave;
        rst >> sis.cell(i).calor.he;
        rst >> chave;
        for (int j = 0; j < sis.cell(i).calor.geom.ncamadas; j++) {
            for (int k = 0; k <= sis.cell(i).calor.ncamada[j]; k++) {
                rst >> sis.cell(i).calor.Tcamada[j][k];
                rst >> chave;
            }
        }
        if (sis.config().flashCompleto == 2) {
            rst >> sis.cell(i).flui.npseudo;
            rst >> chave;
            for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                rst >> sis.cell(i).flui.fracMol[j];
                rst >> chave;
            }
            for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                rst >> sis.cell(i).flui.oCalculatedLiqComposition[j];
                rst >> chave;
            }
            for (int j = 0; j < sis.cell(i).flui.npseudo; j++) {
                rst >> sis.cell(i).flui.oCalculatedVapComposition[j];
                rst >> chave;
            }
            rst >> sis.cell(i).flui.dCalculatedBeta;
            rst >> chave;
            rst >> sis.cell(i).flui.dCalculatedBubbleP;
            rst >> chave;
            rst >> sis.cell(i).flui.iCalculatedThermodynamicCondition;
            rst >> chave;
            rst >> sis.cell(i).flui.iCalculatedStockTankThermodynamicCondition;
            rst >> chave;
            rst >> sis.cell(i).flui.dVaporMassFraction;
            rst >> chave;
            rst >> sis.cell(i).flui.dStockTankVaporMassFraction;
            rst >> chave;
            rst >> sis.cell(i).flui.dStockTankLiquidDensity;
            rst >> chave;
            rst >> sis.cell(i).flui.dStockTankVaporDensity;
            rst >> chave;

            sis.cell(i).flui.atualizaPropCompStandard();
            sis.cell(i).flui.atualizaPropComp(sis.cell(i).pres, sis.cell(i).temp, sis.cell(i).flui.dCalculatedBeta,
                                                sis.cell(i).flui.oCalculatedLiqComposition,
                                                sis.cell(i).flui.oCalculatedVapComposition, sis.config().pocinjec);
        }
    }

    if (sis.config().lingas > 0) {
        ncel = sis.ncelGas + 1;
        for (int i = 0; i < ncel; i++) {
            rst >> chave;
            while (chave != "->")
                rst >> chave;
            rst >> sis.gasCell(i).tempL;
            rst >> chave;
            rst >> sis.gasCell(i).temp;
            rst >> chave;
            rst >> sis.gasCell(i).tempR;
            rst >> chave;
            rst >> sis.gasCell(i).presL;
            rst >> chave;
            rst >> sis.gasCell(i).pres;
            rst >> chave;
            rst >> sis.gasCell(i).presini;
            rst >> chave;
            rst >> sis.gasCell(i).presR;
            rst >> chave;
            rst >> sis.gasCell(i).VGasL;
            rst >> chave;
            rst >> sis.gasCell(i).VGasR;
            rst >> chave;
            rst >> sis.gasCell(i).VGasRR;
            rst >> chave;
            rst >> sis.gasCell(i).u1LL;
            rst >> chave;
            rst >> sis.gasCell(i).u1L;
            rst >> chave;
            rst >> sis.gasCell(i).u1R;
            rst >> chave;
            rst >> sis.gasCell(i).massfonteCH;
            rst >> chave;
            rst >> sis.gasCell(i).fluxcal;
            rst >> chave;
            rst >> sis.gasCell(i).labelchk;
            rst >> chave;
            rst >> sis.gasCell(i).fechamon;
            rst >> chave;
            rst >> sis.gasCell(i).rpchk;
            rst >> chave;
            rst >> sis.gasCell(i).fonteM2;
            rst >> chave;
            rst >> sis.gasCell(i).salinidade;
            rst >> chave;
            rst >> sis.gasCell(i).razInter;
            rst >> chave;
            rst >> sis.gasCell(i).razInterIni;
            rst >> chave;
            sis.gasCell(i).celInter = &(sis.celInter);
            rst >> sis.gasCell(i).tempLini;
            rst >> chave;
            rst >> sis.gasCell(i).tempini;
            rst >> chave;
            rst >> sis.gasCell(i).tempRini;
            rst >> chave;
            rst >> sis.gasCell(i).presLini;
            rst >> chave;
            rst >> sis.gasCell(i).presini;
            rst >> chave;
            rst >> sis.gasCell(i).presRini;
            rst >> chave;
            rst >> sis.gasCell(i).VGasLini;
            rst >> chave;
            rst >> sis.gasCell(i).VGasRini;
            rst >> chave;
            rst >> sis.gasCell(i).VGasRRini;
            rst >> chave;
            rst >> sis.gasCell(i).u1LLini;
            rst >> chave;
            rst >> sis.gasCell(i).u1Lini;
            rst >> chave;
            rst >> sis.gasCell(i).u1Rini;
            rst >> chave;
            rst >> sis.gasCell(i).massfonteCHini;
            rst >> chave;
            rst >> sis.gasCell(i).fonteM2ini;
            rst >> chave;
            rst >> sis.gasCell(i).fechamonini;
            rst >> chave;
            sis.gasCell(i).celInterini = &(sis.celInterIni);
            rst >> sis.gasCell(i).posicini;
            rst >> chave;
            rst >> sis.gasCell(i).rpchkini;
            rst >> chave;
            rst >> sis.gasCell(i).calor.Vint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.Tint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.Tint2;
            rst >> chave;
            rst >> sis.gasCell(i).calor.kint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.cpint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.rhoint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.viscint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.Vconf;
            rst >> chave;
            rst >> sis.gasCell(i).calor.fluxIni;
            rst >> chave;
            rst >> sis.gasCell(i).calor.fluxFim;
            rst >> chave;
            rst >> sis.gasCell(i).calor.ccon;
            rst >> chave;
            rst >> sis.gasCell(i).calor.ncon;
            rst >> chave;
            rst >> sis.gasCell(i).calor.mcon;
            rst >> chave;
            rst >> sis.gasCell(i).calor.npet;
            rst >> chave;
            rst >> sis.gasCell(i).calor.betint;
            rst >> chave;
            rst >> sis.gasCell(i).calor.betext;
            rst >> chave;
            rst >> sis.gasCell(i).calor.reyi;
            rst >> chave;
            rst >> sis.gasCell(i).calor.reye;
            rst >> chave;
            rst >> sis.gasCell(i).calor.grashi;
            rst >> chave;
            rst >> sis.gasCell(i).calor.grashe;
            rst >> chave;
            rst >> sis.gasCell(i).calor.nusi;
            rst >> chave;
            rst >> sis.gasCell(i).calor.nuse;
            rst >> chave;
            rst >> sis.gasCell(i).calor.pri;
            rst >> chave;
            rst >> sis.gasCell(i).calor.pre;
            rst >> chave;
            rst >> sis.gasCell(i).calor.hi;
            rst >> chave;
            rst >> sis.gasCell(i).calor.he;
            rst >> chave;

            for (int j = 0; j < sis.gasCell(i).calor.geom.ncamadas; j++) {
                for (int k = 0; k <= sis.gasCell(i).calor.ncamada[j]; k++) {
                    rst >> sis.gasCell(i).calor.Tcamada[j][k];
                    rst >> chave;
                }
            }
        }
    }

    rst.close();
}

double FQpesado(Cel *celula, int i) {
    double bo;
    double fwi;
    double bswi;
    double beti;
    double qostd;
    if (celula[i].QL > 0.) {
        bo = celula[i - 1].flui.BOFunc(celula[i - 1].pres,
                                       celula[i - 1].temp);
        bswi = celula[i - 1].flui.BSW;
        fwi = bswi / (bo + bswi - bswi * bo);
        beti = celula[i - 1].bet;
    } else {
        bo = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp);
        bswi = celula[i].flui.BSW;
        fwi = bswi / (bo + bswi - bswi * bo);
        beti = celula[i].bet;
        //		rsi = celula[i].flui.RS(celula[i].pres, celula[i].temp)
    }
    return qostd = celula[i].QL * (1 - fwi) * (1. - beti) / bo;
}

double FQleve(Cel *celula, int i) {
    double bo;
    double fwi;
    double bswi;
    double beti;
    double rsi;
    double dgi;
    if (celula[i].QL > 0. && i > 0) {
        bo = celula[i - 1].flui.BOFunc(celula[i - 1].pres,
                                       celula[i - 1].temp);
        bswi = celula[i - 1].flui.BSW;
        fwi = bswi / (bo + bswi - bswi * bo);
        beti = celula[i - 1].bet;
        rsi = celula[i - 1].flui.RS(celula[i - 1].pres,
                                    celula[i - 1].temp) *
              (6.29 / 35.31467);
        dgi = celula[i - 1].flui.Deng;
    } else {
        bo = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp);
        bswi = celula[i].flui.BSW;
        fwi = bswi / (bo + bswi - bswi * bo);
        beti = celula[i].bet;
        rsi = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
        dgi = celula[i].flui.Deng;
    }
    double qostd = celula[i].QL * (1 - fwi) * (1. - beti) / bo;
    return qostd * rsi + (celula[i].MC - celula[i].Mliqini) / (dgi * 1.225);
}

int verificaFonteDuplaReversa(SProd *malha, int aux) {
    if (malha[aux].cell(0).acsr.tipo == 1) {
        if (malha[aux].cell(1).acsr.tipo == 1) {
            if (malha[aux].cell(1).acsr.injg.QGas * malha[aux].cell(0).acsr.injg.QGas < 0) {
                if (fabs(malha[aux].cell(0).acsr.injg.QGas) < fabs(malha[aux].cell(1).acsr.injg.QGas)) {
                    return -1;
                } else
                    return 1;
            } else
                return 1;
        } else
            return 1;
    } else if (malha[aux].cell(0).acsr.tipo == 2) {
        if (malha[aux].cell(1).acsr.tipo == 2) {
            if (malha[aux].cell(1).acsr.injl.QLiq * malha[aux].cell(0).acsr.injl.QLiq < 0) {
                if (fabs(malha[aux].cell(0).acsr.injl.QLiq) < fabs(malha[aux].cell(1).acsr.injl.QLiq)) {
                    return -1;
                } else
                    return 1;
            } else
                return 1;
        } else
            return 1;
    } else if (malha[aux].cell(1).acsr.tipo == 10) {
        double mass0 = malha[aux].cell(0).acsr.injm.MassC +
                       malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassP;
        double mass1 = malha[aux].cell(1).acsr.injm.MassC +
                       malha[aux].cell(1).acsr.injm.MassG + malha[aux].cell(1).acsr.injm.MassP;
        if (mass0 * mass1 < 0) {
            if (fabs(mass0) < fabs(mass1)) {
                return -1;
            } else
                return 1;
        } else
            return 1;
    } else
        return 1;
}

int trocaFonteColetor(SProd *malha, int aux) {
    if (malha[aux].cell(0).acsr.tipo == 1) {
        if (malha[aux].cell(1).acsr.tipo == 1) {
            if (malha[aux].cell(1).acsr.injg.QGas * malha[aux].cell(0).acsr.injg.QGas < 0) {
                if (fabs(malha[aux].cell(0).acsr.injg.QGas) < fabs(malha[aux].cell(1).acsr.injg.QGas)) {
                    InjGas reserva;
                    reserva = malha[aux].cell(0).acsr.injg;
                    malha[aux].cell(0).acsr.injg = malha[aux].cell(1).acsr.injg;
                    malha[aux].cell(1).acsr.injg = reserva;
                    return -1;
                } else
                    return 1;
            } else
                return 1;
        } else
            return 1;
    } else if (malha[aux].cell(0).acsr.tipo == 2) {
        if (malha[aux].cell(1).acsr.tipo == 2) {
            if (malha[aux].cell(1).acsr.injl.QLiq * malha[aux].cell(0).acsr.injl.QLiq < 0) {
                if (fabs(malha[aux].cell(0).acsr.injl.QLiq) < fabs(malha[aux].cell(1).acsr.injl.QLiq)) {
                    InjLiq reserva;
                    reserva = malha[aux].cell(0).acsr.injl;
                    malha[aux].cell(0).acsr.injl = malha[aux].cell(1).acsr.injl;
                    malha[aux].cell(1).acsr.injl = reserva;
                    return -1;
                } else
                    return 1;
            } else
                return 1;
        } else
            return 1;
    } else if (malha[aux].cell(1).acsr.tipo == 10) {
        double mass0 = malha[aux].cell(0).acsr.injm.MassC +
                       malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassP;
        double mass1 = malha[aux].cell(1).acsr.injm.MassC +
                       malha[aux].cell(1).acsr.injm.MassG + malha[aux].cell(1).acsr.injm.MassP;
        if (mass0 * mass1 < 0) {
            if (fabs(mass0) < fabs(mass1)) {
                InjMult reserva;
                reserva = malha[aux].cell(0).acsr.injm;
                malha[aux].cell(0).acsr.injm = malha[aux].cell(1).acsr.injm;
                malha[aux].cell(1).acsr.injm = reserva;
                return -1;
            } else
                return 1;
        } else
            return 1;
    } else
        return 1;
}

void retornaFonteColetor(SProd *malha, int aux) {
    if (malha[aux].cell(0).acsr.tipo == 1) {
        InjGas reserva;
        reserva = malha[aux].cell(0).acsr.injg;
        malha[aux].cell(0).acsr.injg = malha[aux].cell(1).acsr.injg;
        malha[aux].cell(1).acsr.injg = reserva;
    } else if (malha[aux].cell(0).acsr.tipo == 2) {
        InjLiq reserva;
        reserva = malha[aux].cell(0).acsr.injl;
        malha[aux].cell(0).acsr.injl = malha[aux].cell(1).acsr.injl;
        malha[aux].cell(1).acsr.injl = reserva;
    } else if (malha[aux].cell(1).acsr.tipo == 10) {
        InjMult reserva;
        reserva = malha[aux].cell(0).acsr.injm;
        malha[aux].cell(0).acsr.injm = malha[aux].cell(1).acsr.injm;
        malha[aux].cell(1).acsr.injm = reserva;
    }
}

void CicloRedeTrans(SProd *malha, Rede &arqRede,
                    Vcr<int> &inativo, Vcr<double> &alfrev,
                    Vcr<double> &betrev, Vcr<double> &titrev) {
    int narq = arqRede.nsisprod;
    Vcr<int> Resolv(narq, 0);
    Vcr<int> IndNorma(narq, 0);
    int ciclomalha = 0;
    ProFlu flumist;
    while (ciclomalha < narq - 1) {
        for (int i = 0; i < narq; i++) {
            if (arqRede.malha[i].nafluente == 0 && Resolv[i] == 0 && inativo[i] == 0) {
                Resolv[i] = 1;
                ciclomalha++;
            }
        }
        int resolvGlob = 0;
        for (int i = 0; i < narq; i++)
            resolvGlob += Resolv[i];

        while (resolvGlob < narq) {
            for (int i = 0; i < narq; i++) {
                int naflu = arqRede.malha[i].nafluente;
                int cicloaflu = 0;
                if (arqRede.malha[i].nafluente > 0 && Resolv[i] == 0) {
                    for (int j = 0; j < narq; j++) {
                        for (int k = 0; k < naflu; k++) {
                            if (arqRede.malha[i].afluente[k] == j && Resolv[j] == 1)
                                cicloaflu++;
                        }
                    }
                    if (cicloaflu == naflu && inativo[i] == 0) {
                        int afluAtivo = 0;
                        vector<int> indAflu;
                        int coleAtivo = 0;
                        vector<int> indCole;
                        for (int k = 0; k < naflu; k++) {
                            if (inativo[k] == 0) {
                                indAflu.push_back(arqRede.malha[i].afluente[k]);
                                afluAtivo++;
                            }
                        }
                        int ind = arqRede.malha[i].afluente[0];
                        int ncol = arqRede.malha[ind].ncoleta;
                        for (int k = 0; k < ncol; k++) {
                            if (inativo[k] == 0) {
                                indCole.push_back(arqRede.malha[ind].coleta[k]);
                                coleAtivo++;
                            }
                        }
                        int ntotal = coleAtivo + afluAtivo;

                        int nderiva = arqRede.malha[ind].ncoleta;
                        Vcr<int> ordCol(nderiva);
                        vector<double> dcol;
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            dcol.push_back(malha[aux].cell(0).duto.dia);
                            malha[aux].presE = -1.;
                        }
                        sort(dcol.begin(), dcol.end());
                        Vcr<int> carregado(nderiva, 0);
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            for (int icol2 = 0; icol2 < nderiva; icol2++) {
                                if (dcol[icol2] == malha[aux].cell(0).duto.dia && carregado[icol2] == 0) {
                                    ordCol[icol2] = aux;
                                    carregado[icol2] = 1;
                                    icol2 = nderiva;
                                }
                            }
                        }

                        int aux0 = ordCol[nderiva - 1];

                        Vcr<double> qostd(ntotal);
                        Vcr<double> RGO(ntotal);
                        Vcr<double> rhololeoST(ntotal);
                        Vcr<double> rhogST(ntotal);
                        Vcr<double> rholliqIS(ntotal);
                        Vcr<double> rhogIS(ntotal);
                        Vcr<double> API(ntotal);
                        Vcr<double> BSW(ntotal);
                        Vcr<double> Bet(ntotal);
                        Vcr<double> Rhocomp(ntotal);
                        Vcr<double> temp(ntotal);
                        Vcr<double> tempFim(ntotal);
                        Vcr<double> cpl(ntotal);
                        Vcr<double> cpg(ntotal);
                        Vcr<double> Mliq(ntotal);
                        Vcr<double> Qliq(ntotal);
                        Vcr<double> Mcomp(ntotal);
                        Vcr<double> Qcomp(ntotal);
                        Vcr<double> Mgas(ntotal);
                        Vcr<double> Qgas(ntotal);
                        Vcr<double> Deng(ntotal);
                        Vcr<double> yco2(ntotal);
                        Vcr<double> Denag(ntotal);

                        Vcr<double> MasstempP(ntotal);
                        Vcr<double> MasstempC(ntotal);
                        Vcr<double> MasstempG(ntotal);
                        Vcr<double> variaMassP(ntotal);
                        Vcr<double> variaMassC(ntotal);
                        Vcr<double> variaMassG(ntotal);

                        double tL = 0.;
                        if (malha[ind].config().flashCompleto == 1 && malha[i].config().tabent.tmin > tL)
                            tL = malha[ind].config().tabent.tmin;
                        double tH = 70.;
                        if (malha[ind].config().flashCompleto == 1 && malha[ind].config().tabent.tmax < tH)
                            tH = malha[ind].config().tabent.tmax;
                        double qostdTot = 0;
                        double qgstdTot = 0;
                        double dengmist = 0;
                        double yco2mist = 0;
                        double apimist = 0.;
                        double bswmist;
                        double RGOmist;
                        double qwmist = 0;
                        double denagmist = 0;
                        double mliqmistC = 0;
                        double mliqmist = 0;
                        double mgasmist = 0;
                        double cpmist = 0;
                        double tempmist = 0;
                        double LVisL = 0.;
                        double LVisH = 0.;
                        double Qlmist = 0.;
                        double QGmist = 0.;
                        double betmist = 0.;
                        double alfmist = 1.;
                        double titmist = 1.;
                        double qlmistStd = 0.;
                        double MassAfluP = 0.;
                        double MassAfluC = 0.;
                        double MassAfluG = 0.;

                        if ((*arqRede.vg1dSP).lixo5 >= 4999) {
                        }

                        int afluAlimenta = 0;
                        int coleAlimenta = 0;
                        for (int k = 0; k < afluAtivo; k++) {
                            int ind = indAflu[k];
                            int fim = malha[ind].ncel - 1;

                            if ((*arqRede.vg1dSP).iterRedeT == 0) {

                                MasstempC[k] = malha[ind].cell(fim + 1).fontemassCR;
                                MasstempP[k] = malha[ind].cell(fim + 1).fontemassLR;
                                MasstempG[k] = malha[ind].cell(fim + 1).fontemassGR;

                            } else {

                                MasstempC[k] = malha[ind].fontemassCRBuf;

                                MasstempP[k] = malha[ind].fontemassPRBuf;

                                MasstempG[k] = malha[ind].fontemassGRBuf;
                            }
                            double pres = malha[ind].pGSup;
                            malha[ind].calcTempFim();
                            tempFim[k] = malha[ind].tempSup;
                            temp[k] = malha[ind].tempSup;

                            if ((malha[ind].cell(fim).MliqiniR) > 0.) {
                                double bo = malha[ind].cell(fim).flui.BOFunc(pres, temp[k]);
                                RGO[k] = malha[ind].cell(fim).flui.RGO;
                                BSW[k] = malha[ind].cell(fim).flui.BSW;
                                API[k] = malha[ind].cell(fim).flui.API;
                                rhololeoST[k] = 1000. * 141.5 / (131.5 + API[k]);
                                Bet[k] = malha[ind].cell(fim).bet;
                                if (Bet[k] > 0) {
                                }
                                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;

                                double ba;
                                double rs;
                                if (malha[ind].cell(fim).flui.RGO < 1e6) {
                                    rs = malha[ind].cell(fim).flui.RS(malha[ind].cell(fim).pres,
                                                                        malha[ind].cell(fim).temp) *
                                         6.29 / 35.31467;
                                    bo = malha[ind].cell(fim).flui.BOFunc(malha[ind].cell(fim).pres,
                                                                            malha[ind].cell(fim).temp, rs);
                                    ba = malha[ind].cell(fim).flui.BAFunc(malha[ind].cell(fim).pres,
                                                                            malha[ind].cell(fim).temp);
                                } else {
                                    bo = 1;
                                    ba = 0.;
                                }
                                double FW = malha[ind].cell(fim).flui.BSW * ba /
                                            (bo + ba * malha[ind].cell(fim).flui.BSW -
                                             malha[ind].cell(fim).flui.BSW * bo);

                                qostd[k] = 0.;
                                if (malha[ind].cell(fim).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                                    qostd[k] = malha[ind].cell(fim).Mliqini *
                                               (1. - FW) * (1. - malha[ind].cell(fim).bet) /
                                               (bo * rholliqIS[k]);
                                Rhocomp[k] = rc;
                                cpl[k] = (1. - Bet[k]) * malha[ind].cell(fim).flui.CalorLiq(pres, temp[k]) +
                                         Bet[k] * malha[ind].cell(fim).fluicol.CalorLiq(pres, temp[k]);
                                cpg[k] = malha[ind].cell(fim).flui.CalorGas(pres, temp[k]);
                                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                                Qliq[k] = malha[ind].cell(fim + 1).QL;
                                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                                Qcomp[k] = Qliq[k] * Bet[k];
                                Denag[k] = malha[ind].cell(fim).flui.Denag;
                                Deng[k] = malha[ind].cell(fim).flui.Deng;
                                yco2[k] = malha[ind].cell(fim).flui.yco2;
                            }
                            if ((malha[ind].cell(fim).MR - malha[ind].cell(fim).MliqiniR) > 0.) {
                                rhogST[k] = malha[ind].cell(fim).flui.Deng * 1.225;
                                rhogIS[k] = malha[ind].cell(fim).flui.MasEspGas(pres, temp[k]);
                                cpg[k] = malha[ind].cell(fim).flui.CalorGas(pres, temp[k]);
                                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                                Qgas[k] = Mgas[k] / rhogIS[k];
                                Deng[k] = malha[ind].cell(fim).flui.Deng;
                                yco2[k] = malha[ind].cell(fim).flui.yco2;
                                if ((malha[ind].cell(fim).MliqiniR) < 0.)
                                    RGO[k] = 1400.;
                            }

                            qostdTot += qostd[k];
                            if ((malha[ind].cell(fim).MliqiniR) > 0.) {
                                double rst = malha[ind].cell(fim).flui.RS(pres, temp[k]);
                                qgstdTot += qostd[k] * rst;
                                dengmist += malha[ind].cell(fim).flui.rDgD * Deng[k] * qostd[k] * rst;
                                yco2mist += yco2[k] * qostd[k] * rst;
                            }
                            if ((malha[ind].cell(fim).MR - malha[ind].cell(fim).MliqiniR) > 0) {
                                double qgstdLocal = Mgas[k] / (Deng[k] * 1.225);
                                qgstdTot += qgstdLocal;
                                dengmist += malha[ind].cell(fim).flui.rDgL * Deng[k] * qgstdLocal;
                                yco2mist += yco2[k] * qgstdLocal;
                            }
                            apimist += (141.5 / (API[k] + 131.5)) * qostd[k];
                            double qw1;
                            if ((1. - BSW[k]) > 0)
                                qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                            else
                                qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                            qwmist += qw1;
                            denagmist += Denag[k] * qw1;
                            mliqmist += Mliq[k];
                            mgasmist += Mgas[k];
                            mliqmistC += Qcomp[k] * malha[i].cell(0).fluicol.MasEspFlu(pres, temp[k]);
                            cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
                            tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * tempFim[k];
                            LVisL += qostd[k] * malha[ind].cell(fim).flui.VisOM(tL);
                            LVisH += qostd[k] * malha[ind].cell(fim).flui.VisOM(tH);
                            Qlmist += Qliq[k];
                            QGmist += Qgas[k];
                            qlmistStd += (qw1 + qostd[k] +
                                          Qcomp[k] * malha[i].cell(0).fluicol.MasEspFlu(pres, temp[k]) /
                                              malha[i].cell(0).fluicol.MasEspFlu(1.001, 15.));
                            betmist += (Qcomp[k]);

                            afluAlimenta += 1;

                            MassAfluP -= (MasstempP[k] + 0 * variaMassP[k]);
                            MassAfluC -= (MasstempC[k] + 0 * variaMassC[k]);
                            MassAfluG -= (MasstempG[k] + 0 * variaMassG[k]);
                        }

                        for (int k = 0; k < coleAtivo; k++) {
                            int ind = indCole[k];
                            int ini = 0;
                            temp[k + afluAlimenta] = malha[ind].cell(ini).temp;
                            tempFim[k + afluAlimenta] = malha[ind].cell(ini).temp;
                            double pres;
                            if ((*arqRede.vg1dSP).iterRedeT == 0)
                                pres = malha[ind].cell(ini).pres;
                            else
                                pres = malha[ind].cell(ini).presBuf;

                            double rpini = malha[ind].cell(ini).flui.MasEspLiq(pres, temp[k + afluAlimenta]);
                            double rcini = malha[ind].cell(ini).fluicol.MasEspFlu(pres, temp[k + afluAlimenta]);
                            double betini = malha[ind].cell(ini).bet;
                            double rholini = (1. - betini) * rpini + betini * rcini;
                            double CoMassC;
                            double CoMassP;
                            double CoMassL;
                            double CoMassG;
                            if ((*arqRede.vg1dSP).iterRedeT == 0) {
                                CoMassC = (malha[ind].cell(ini).Mliqini / rholini) * betini * rcini;
                                CoMassP = (malha[ind].cell(ini).Mliqini / rholini) * (1. - betini) * rpini;
                                CoMassL = CoMassP + CoMassC;
                                CoMassG = malha[ind].cell(ini).MC - malha[ind].cell(ini).Mliqini;
                            } else {
                                CoMassC = (malha[ind].cell(ini).MliqiniBuf / rholini) * betini * rcini;
                                CoMassP = (malha[ind].cell(ini).MliqiniBuf / rholini) * (1. - betini) * rpini;
                                CoMassL = CoMassP + CoMassC;
                                CoMassG = malha[ind].cell(ini).MCBuf - malha[ind].cell(ini).MliqiniBuf;
                            }

                            if ((CoMassL) < 0.) {
                                double bo = malha[ind].cell(ini).flui.BOFunc(
                                    malha[ind].cell(ini).pres, malha[ind].cell(ini).temp);
                                double rp = malha[ind].cell(ini).flui.MasEspLiq(pres, temp[k + afluAlimenta]);
                                double rc = malha[ind].cell(ini).fluicol.MasEspFlu(pres, temp[k + afluAlimenta]);
                                RGO[k + afluAlimenta] = malha[ind].cell(ini).flui.RGO;
                                Bet[k + afluAlimenta] = (-CoMassC / rc) / (-(CoMassC / rc) - (CoMassP / rp));
                                BSW[k + afluAlimenta] = malha[ind].cell(ini).flui.BSW;
                                API[k + afluAlimenta] = malha[ind].cell(ini).flui.API;
                                rhololeoST[k + afluAlimenta] = 1000. * 141.5 / (131.5 + API[k + afluAlimenta]);

                                double ba;
                                double rs;
                                if (malha[ind].cell(ini).flui.RGO < 1e6) {
                                    rs = malha[ind].cell(ini).flui.RS(malha[ind].cell(ini).pres,
                                                                        malha[ind].cell(ini).temp) *
                                         6.29 / 35.31467;
                                    bo = malha[ind].cell(ini).flui.BOFunc(malha[ind].cell(ini).pres,
                                                                            malha[ind].cell(ini).temp, rs);
                                    ba = malha[ind].cell(ini).flui.BAFunc(malha[ind].cell(ini).pres,
                                                                            malha[ind].cell(ini).temp);
                                } else {
                                    bo = 1;
                                    ba = 0.;
                                }
                                double FW = malha[ind].cell(ini).flui.BSW * ba /
                                            (bo + ba * malha[ind].cell(ini).flui.BSW -
                                             malha[ind].cell(ini).flui.BSW * bo);

                                rholliqIS[k + afluAlimenta] = (1. - Bet[k + afluAlimenta]) * rp + Bet[k + afluAlimenta] * rc;
                                qostd[k + afluAlimenta] = 0.;
                                if (malha[ind].cell(ini).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                                    qostd[k + afluAlimenta] = CoMassL * (1. - FW) * (1. - Bet[k + afluAlimenta]) /
                                                              (bo * rholliqIS[k + afluAlimenta]);
                                Rhocomp[k + afluAlimenta] = rc;
                                cpl[k + afluAlimenta] = (1. - Bet[k + afluAlimenta]) * malha[ind].cell(ini).flui.CalorLiq(pres, temp[k + afluAlimenta]) +
                                                        Bet[k + afluAlimenta] * malha[ind].cell(ini).fluicol.CalorLiq(pres, temp[k + afluAlimenta]);
                                cpg[k + afluAlimenta] = malha[ind].cell(ini).flui.CalorGas(pres, temp[k + afluAlimenta]);
                                Mliq[k + afluAlimenta] = CoMassL;
                                Qliq[k + afluAlimenta] = CoMassL / rholliqIS[k + afluAlimenta];
                                Mcomp[k + afluAlimenta] = Mliq[k + afluAlimenta] * Bet[k + afluAlimenta] * Rhocomp[k + afluAlimenta] / rholliqIS[k + afluAlimenta];
                                Qcomp[k + afluAlimenta] = Qliq[k + afluAlimenta] * Bet[k + afluAlimenta];
                                Denag[k + afluAlimenta] = malha[ind].cell(ini).flui.Denag;
                                Deng[k + afluAlimenta] = malha[ind].cell(ini).flui.Deng;
                                yco2[k + afluAlimenta] = malha[ind].cell(ini).flui.yco2;
                            }
                            if ((CoMassG) < 0.) {
                                rhogIS[k + afluAlimenta] = malha[ind].cell(ini).flui.MasEspGas(pres, temp[k + afluAlimenta]);
                                cpg[k + afluAlimenta] = malha[ind].cell(ini).flui.CalorGas(pres, temp[k + afluAlimenta]);
                                Mgas[k + afluAlimenta] = CoMassG;
                                Qgas[k + afluAlimenta] = Mgas[k + afluAlimenta] / rhogIS[k + afluAlimenta];
                                Deng[k + afluAlimenta] = malha[ind].cell(ini).flui.Deng;
                                yco2[k + afluAlimenta] = malha[ind].cell(ini).flui.yco2;
                            }

                            qostdTot -= qostd[k + afluAlimenta];
                            if ((CoMassL) < 0.) {
                                double rst = malha[ind].cell(ini).flui.RS(pres, temp[k + afluAlimenta]);
                                qgstdTot -= qostd[k + afluAlimenta] * rst;
                                dengmist -= malha[ind].cell(ini).flui.rDgD * Deng[k + afluAlimenta] * qostd[k + afluAlimenta] * rst;
                                yco2mist -= yco2[k + afluAlimenta] * qostd[k + afluAlimenta] * rst;
                            }
                            if ((CoMassG) < 0.) {
                                double qgstdLocal = Mgas[k + afluAlimenta] / (Deng[k + afluAlimenta] * 1.225);
                                qgstdTot -= qgstdLocal;
                                dengmist -= malha[ind].cell(ini).flui.rDgL * Deng[k + afluAlimenta] * qgstdLocal;
                                yco2mist -= yco2[k + afluAlimenta] * qgstdLocal;
                            }
                            apimist -= (141.5 / (API[k + afluAlimenta] + 131.5)) * qostd[k + afluAlimenta];
                            double qw1;
                            if ((1. - BSW[k + afluAlimenta]) > 0)
                                qw1 = qostd[k + afluAlimenta] * BSW[k + afluAlimenta] / (1. - BSW[k + afluAlimenta]);
                            else
                                qw1 = Mliq[k + afluAlimenta] * (1. - Bet[k + afluAlimenta]) / (1000. * Denag[k + afluAlimenta]);
                            qwmist -= qw1;
                            denagmist -= Denag[k + afluAlimenta] * qw1;
                            mliqmist -= Mliq[k + afluAlimenta];
                            mgasmist -= Mgas[k + afluAlimenta];
                            mliqmistC -= Qcomp[k + afluAlimenta] *
                                         malha[i].cell(0).fluicol.MasEspFlu(pres, temp[k + afluAlimenta]);
                            cpmist -= (Mliq[k + afluAlimenta] * cpl[k + afluAlimenta] + Mgas[k + afluAlimenta] * cpg[k + afluAlimenta]);
                            tempmist -= (Mliq[k + afluAlimenta] * cpl[k + afluAlimenta] +
                                         Mgas[k + afluAlimenta] * cpg[k + afluAlimenta]) *
                                        tempFim[k + afluAlimenta];
                            LVisL -= qostd[k + afluAlimenta] * malha[ind].cell(ini).flui.VisOM(tL);
                            LVisH -= qostd[k + afluAlimenta] * malha[ind].cell(ini).flui.VisOM(tH);
                            Qlmist -= Qliq[k + afluAlimenta];
                            QGmist -= Qgas[k + afluAlimenta];
                            qlmistStd -= (qw1 + qostd[k + afluAlimenta] +
                                          Qcomp[k + afluAlimenta] * malha[i].cell(0).fluicol.MasEspFlu(pres, temp[k + afluAlimenta]) /
                                              malha[i].cell(0).fluicol.MasEspFlu(1.001, 15.));
                            betmist -= (Qcomp[k + afluAlimenta]);

                            coleAlimenta += 1;
                        }

                        if (mliqmist > (*arqRede.vg1dSP).localtiny * 1e-5 || mgasmist > (*arqRede.vg1dSP).localtiny * 1e-5) {
                            tempmist = tempmist / cpmist;

                            alfmist = malha[aux0].cell(0).alf;

                            if (qostdTot > (*arqRede.vg1dSP).localtiny * 1e-5) {
                                RGOmist = qgstdTot / qostdTot;
                                apimist = apimist / qostdTot;
                                apimist = 141.5 / apimist - 131.5;
                                LVisL = LVisL / qostdTot;
                                LVisH = LVisH / qostdTot;
                            } else {
                                RGOmist = malha[aux0].cell(0).flui.RGO;
                                apimist = malha[aux0].cell(0).flui.API;
                                LVisL = malha[aux0].cell(0).flui.VisOM(tL);
                                LVisH = malha[aux0].cell(0).flui.VisOM(tH);
                            }
                            if ((qgstdTot) > (*arqRede.vg1dSP).localtiny * 1e-5) {
                                dengmist = dengmist / qgstdTot;
                                yco2mist = yco2mist / qgstdTot;
                            } else {
                                dengmist = malha[aux0].cell(0).flui.Deng;
                                yco2mist = malha[aux0].cell(0).flui.yco2;
                            }
                            if ((qostdTot + qwmist) > (*arqRede.vg1dSP).localtiny * 1e-5)
                                bswmist = qwmist / (qostdTot + qwmist);
                            else
                                bswmist = malha[aux0].cell(0).flui.BSW;
                            if (qwmist > (*arqRede.vg1dSP).localtiny * 1e-5)
                                denagmist = denagmist / qwmist;
                            else
                                denagmist = malha[aux0].cell(0).flui.Denag;
                            if (Qlmist > (*arqRede.vg1dSP).localtiny * 1e-5)
                                betmist = betmist / Qlmist;
                            else
                                betmist = malha[aux0].cell(0).bet;

                            if (mgasmist > 0. && (mliqmist) > 0)
                                titmist = mgasmist / (mgasmist + mliqmist);
                            else {
                                double alfCol = malha[aux0].cell(0).alf;
                                double presCol = malha[aux0].cell(0).pres;
                                double tempCol = malha[aux0].cell(0).temp;
                                double betCol = malha[aux0].cell(0).bet;
                                double rhop = malha[aux0].cell(0).flui.MasEspLiq(presCol, tempCol);
                                double rhoc = malha[aux0].cell(0).fluicol.MasEspFlu(presCol, tempCol);
                                double rhog = malha[aux0].cell(0).flui.MasEspGas(presCol, tempCol);
                                double rhol = betCol * rhoc + (1. - betCol) * rhop;
                                double rhoCol = alfCol * rhog + (1. - alfCol) * rhol;
                                titmist = rhog * alfCol / ((1. - alfCol) * rhoCol + rhog * alfCol);
                            }

                        } else {
                            tempmist = malha[aux0].cell(0).temp;
                            RGOmist = malha[aux0].cell(0).flui.RGO;
                            apimist = malha[aux0].cell(0).flui.API;
                            LVisL = malha[aux0].cell(0).flui.VisOM(tL);
                            LVisH = malha[aux0].cell(0).flui.VisOM(tH);
                            dengmist = malha[aux0].cell(0).flui.Deng;
                            yco2mist = malha[aux0].cell(0).flui.yco2;
                            bswmist = malha[aux0].cell(0).flui.BSW;
                            denagmist = malha[aux0].cell(0).flui.Denag;
                            betmist = malha[aux0].cell(0).bet;
                            alfmist = malha[aux0].cell(0).alf;
                            double betCol = malha[aux0].cell(0).bet;
                            double alfCol = malha[aux0].cell(0).alf;
                            double presCol = malha[aux0].cell(0).pres;
                            double tempCol = malha[aux0].cell(0).temp;
                            double rhop = malha[aux0].cell(0).flui.MasEspLiq(presCol, tempCol);
                            double rhoc = malha[aux0].cell(0).fluicol.MasEspFlu(presCol, tempCol);
                            double rhog = malha[aux0].cell(0).flui.MasEspGas(presCol, tempCol);
                            double rhol = betCol * rhoc + (1. - betCol) * rhop;
                            double rhoCol = alfCol * rhog + (1. - alfCol) * rhol;
                            titmist = rhog * alfCol / ((1. - alfCol) * rhoCol + rhog * alfCol);
                        }

                        int naflu = arqRede.malha[i].nafluente;
                        if ((*arqRede.vg1dSP).iterRedeT == 0) {
                            for (int k = 0; k < naflu; k++) {
                                int iaflu = arqRede.malha[i].afluente[k];
                                titrev[iaflu] = titmist;
                                alfrev[iaflu] = alfmist;
                                betrev[iaflu] = betmist;
                            }
                        }

                        for (int icol = 0; icol < nderiva - 1; icol++) {
                            int aux = ordCol[icol];

                            double massP;
                            double massC;
                            double massG;

                            if ((*arqRede.vg1dSP).iterRedeT == 0) {

                                massC = malha[aux].cell(0).QL * malha[aux].cell(0).bet *
                                        malha[aux].cell(0).fluicol.MasEspFlu(malha[aux].cell(0).pres, malha[aux].cell(0).temp);

                                massP = malha[aux].cell(0).Mliqini - massC;

                                massG = malha[aux].cell(0).MC - malha[aux].cell(0).Mliqini;

                            } else {

                                double rhomistBuf = malha[aux].cell(0).bet *
                                                        malha[aux].cell(0).fluicol.MasEspFlu(malha[aux].cell(0).presBuf,
                                                                                               malha[aux].cell(0).temp) +
                                                    (1. - malha[aux].cell(0).bet) *
                                                        malha[aux].cell(0).flui.MasEspLiq(malha[aux].cell(0).presBuf,
                                                                                            malha[aux].cell(0).temp);
                                double QLbuf = malha[aux].cell(0).MliqiniBuf / rhomistBuf;
                                massC = QLbuf * malha[aux].cell(0).bet *
                                        malha[aux].cell(0).fluicol.MasEspFlu(malha[aux].cell(0).presBuf, malha[aux].cell(0).temp);

                                massP = malha[aux].cell(0).MliqiniBuf - massC;

                                massG = malha[aux].cell(0).MCBuf - malha[aux].cell(0).MliqiniBuf;
                            }

                            MassAfluP -= (massP);
                            MassAfluC -= (massC);
                            MassAfluG -= (massG);
                        }

                        for (int icol = nderiva - 1; icol >= 0; icol--) {
                            if (icol == nderiva - 1) {
                                int aux = ordCol[icol];
                                if ((*arqRede.vg1dSP).iterRedeT == 0) {
                                    malha[aux].cell(0).acsr.injm.temp = tempmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.RGO = RGOmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.Deng = dengmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.Denag = denagmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.yco2 = yco2mist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.API = apimist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.BSW = bswmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.TempL = tL;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.TempH = tH;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.LVisL = LVisL;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.LVisH = LVisH;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.RenovaFluido();
                                    malha[aux].cell(0).acsr.injm.temp = tempmist;
                                    if (malha[aux].config().corrDeng == 0) {
                                        malha[aux].cell(i).acsr.injm.FluidoPro.rDgD = 1.;
                                        malha[aux].cell(i).acsr.injm.FluidoPro.rDgL = 1.;
                                    } else {

                                        malha[aux].cell(i).acsr.injm.FluidoPro.razDegD(malha[aux].cell(i).pres,
                                                                                         malha[aux].cell(i).acsr.injm.temp);
                                        malha[aux].cell(i).acsr.injm.FluidoPro.rzDegL(malha[aux].cell(i).pres,
                                                                                        malha[aux].cell(i).acsr.injm.temp);
                                    }
                                    malha[aux].temperatura = tempmist;
                                }
                                malha[aux].cell(0).acsr.injm.MassP = MassAfluP;
                                malha[aux].cell(0).acsr.injm.MassC = MassAfluC;
                                malha[aux].cell(0).acsr.injm.MassG = MassAfluG;
                                if (MassAfluP < 0 && (malha[aux].cell(0).alf > 1 - (*arqRede.vg1dSP).localtiny)) {
                                    malha[aux].cell(0).acsr.injm.MassP = 0.;
                                }
                                if (MassAfluC < 0 && ((malha[aux].cell(0).alf > 1 - (*arqRede.vg1dSP).localtiny) || malha[aux].cell(0).bet < (*arqRede.vg1dSP).localtiny)) {
                                    malha[aux].cell(0).acsr.injm.MassC = 0.;
                                }
                                if (MassAfluG < 0 && (malha[aux].cell(0).alf < (*arqRede.vg1dSP).localtiny)) {
                                    malha[aux].cell(0).acsr.injm.MassG = 0.;
                                }

                                Resolv[aux] = 1;
                                resolvGlob += Resolv[aux];
                                ciclomalha++;
                                for (int iaflu = 0; iaflu < arqRede.malha[aux].nafluente; iaflu++) {
                                    int indaflu = arqRede.malha[aux].afluente[iaflu];
                                    malha[indaflu].titRev = alfmist;
                                    malha[indaflu].betaRev = betmist;

                                    if ((*arqRede.vg1dSP).iterRedeT == 0)
                                        malha[indaflu].tGSup = tempmist;
                                    if (malha[indaflu].cell(malha[indaflu].ncel).MC < -(*arqRede.vg1dSP).localtiny) {
                                        malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).flui;
                                    }
                                    if (arqRede.malha[indaflu].presimposta == 0) {
                                        if ((*arqRede.vg1dSP).iterRedeT == 0)
                                            malha[indaflu].pGSup = malha[aux].cell(0).pres;
                                        else
                                            malha[indaflu].pGSup = malha[aux].cell(0).presBuf;
                                        if (malha[indaflu].masChkSup == 0)
                                            malha[indaflu].presfim = malha[indaflu].pGSup;
                                    }
                                    malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                }
                            } else {
                                int aux = ordCol[icol];
                                malha[aux].config().ConContEntrada = 1;

                                malha[aux].cell(0).acsr.injm.MassP = 0.;
                                malha[aux].cell(0).acsr.injm.MassC = 0.;
                                malha[aux].cell(0).acsr.injm.MassG = 0.;

                                int auxfim = ordCol[nderiva - 1];
                                malha[aux].cell(0).acsr.injm.FluidoPro = malha[auxfim].cell(0).acsr.injm.FluidoPro;
                                malha[aux].cell(0).fluiL = &(malha[auxfim].cell(0).acsr.injm.FluidoPro);

                                malha[aux].tempE = tempmist;
                                malha[aux].titE = titmist;
                                malha[aux].betaE = betmist;

                                if ((*arqRede.vg1dSP).iterRedeT == 0) {
                                    malha[aux].presE = malha[auxfim].cell(0).pres;
                                } else {
                                    malha[aux].presE = malha[auxfim].cell(0).presBuf;
                                }
                                Resolv[aux] = 1;
                                resolvGlob += Resolv[aux];
                                ciclomalha++;
                            }
                        }
                    }
                }
            }
        }
    }
}

int buscaNoColetorMrestre(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int iaflu) {
    int ncol = arqRede.malha[iaflu].ncoleta;
    int icol = arqRede.malha[iaflu].coleta[0];

    double dmax = malha[icol].cell(0).duto.dia;
    int idmax = 0;
    for (int i = 1; i < ncol; i++) {
        icol = arqRede.malha[iaflu].coleta[i];
        if (malha[icol].cell(0).duto.dia > dmax && inativo[icol] == 0) {
            dmax = malha[icol].cell(0).duto.dia;
            idmax = i;
        }
    }
    return idmax;
}

void celAfluFinal(int naflu, int ncol, SProd *malha, Rede &arqRede, Vcr<int> &inativo) {
    int iaflu = malha[naflu].ncel;
    int icol = 0;
    int colet = arqRede.malha[naflu].ncoleta;
    int ncolM = buscaNoColetorMrestre(malha, arqRede, inativo, naflu);
    double Mt = 0;
    double Ml = 0;
    for (int j = 0; j < colet; j++) {
        int qtramo = arqRede.malha[naflu].coleta[j];
        if (j != ncolM) {
            Mt += malha[qtramo].cell(icol).MC;
            Ml += malha[qtramo].cell(icol).Mliqini;
        } else {
            Mt += malha[qtramo].cell(icol).MR;
            Ml += malha[qtramo].cell(icol).MliqiniR;
        }
    }

    malha[naflu].cell(iaflu).dutoR = malha[ncol].cell(icol).duto;
    malha[naflu].cell(iaflu).flui = malha[ncol].cell(icol).flui;
    malha[naflu].cell(iaflu).fluiR = malha[ncol].cell(icol).fluiR;
    malha[naflu].cell(iaflu).fluicol = malha[ncol].cell(icol).fluicol;

    malha[naflu].cell(iaflu).tempR = malha[ncol].cell(icol).tempR;
    if (malha[naflu].cell(iaflu).VTemper < 0.) {
        malha[naflu].cell(iaflu).temp = malha[ncol].cell(icol).acsr.injm.temp;
    }

    malha[naflu].pGSup = malha[ncol].cell(icol).pres;

    malha[naflu].cell(iaflu).presauxR = malha[ncol].cell(icol).presauxR;
    malha[naflu].cell(iaflu).presR = malha[ncol].cell(icol).presR;

    malha[naflu].cell(iaflu).MR = malha[naflu].cell(iaflu).MC;

    if (malha[naflu].masChkSup == 0 && malha[naflu].config().chkv == 0) {
        malha[naflu].cell(iaflu).alf = malha[naflu].titRev;
        malha[naflu].cell(iaflu).alfR = malha[ncol].cell(icol).alfR;
    }

    if (malha[naflu].masChkSup == 0 && malha[naflu].config().chkv == 0) {
        malha[naflu].cell(iaflu).bet = malha[naflu].betaRev;
        malha[naflu].cell(iaflu).betR = malha[ncol].cell(icol).betR;
    }

    malha[naflu].cell(iaflu).betI = malha[naflu].cell(iaflu).bet;
    malha[naflu].cell(iaflu).betRI = malha[naflu].cell(iaflu).bet;
    malha[naflu].cell(iaflu).FW = malha[ncol].cell(icol).FW;

    // malha[naflu].cell(iaflu).MliqiniR0 = Ml;// alteracao7
    malha[naflu].cell(iaflu).MliqiniR = malha[naflu].cell(iaflu).Mliqini;
    // malha[naflu].cell(iaflu).MliqiniR0 = malha[naflu].cell(iaflu).Mliqini;// alteracao7

    malha[naflu].cell(iaflu).QLR = malha[naflu].cell(iaflu).QL;
    // malha[naflu].cell(iaflu).QLRini=malha[naflu].cell(iaflu).QL;// alteracao7

    malha[naflu].cell(iaflu).rpR = malha[ncol].cell(icol).rpR;
    malha[naflu].cell(iaflu).rcR = malha[ncol].cell(icol).rcR;

    malha[naflu].cell(iaflu).rpRi = malha[ncol].cell(icol).rpRi;
    malha[naflu].cell(iaflu).rcRi = malha[ncol].cell(icol).rcRi;

    if (malha[naflu].masChkSup == 0 && malha[naflu].config().chkv == 0) {
        malha[naflu].cell(iaflu).alfPigE = malha[naflu].titRev;
        malha[naflu].cell(iaflu).alfPigER = malha[naflu].titRev;
        malha[naflu].cell(iaflu).alfPigD = malha[naflu].titRev;
        malha[naflu].cell(iaflu).betPigE = malha[naflu].betaRev;
        malha[naflu].cell(iaflu).betPigD = malha[naflu].betaRev;
    }
}

void corrigeVazNo(SProd *malha, int ind) {
    int ncel = malha[ind].ncel;

    if (((*malha[ind].vg1dSP).chaverede == 1 && malha[ind].noextremo != 1 && malha[ind].titRev >= 1. - (*malha[ind].vg1dSP).localtiny)) {
        if ((malha[ind].cell(ncel - 1).MliqiniR < 0)) {
            malha[ind].cell(ncel - 1).MR -= (malha[ind].cell(ncel - 1).MliqiniR);
            malha[ind].cell(ncel - 1).MliqiniR = 0;
            malha[ind].cell(ncel - 1).term1R = 0;
            malha[ind].cell(ncel - 1).term2R = 0;
            malha[ind].cell(ncel).Mliqini = 0;
            malha[ind].cell(ncel).MC = (malha[ind].cell(ncel - 1).MR);
            malha[ind].cell(ncel).Mliqini = 0;
            malha[ind].cell(ncel).term1 = 0;
            malha[ind].cell(ncel).term2 = 0;
        }
        if ((malha[ind].masChkSup == 0 || malha[ind].aberto == 1)) {
            malha[ind].cell(ncel - 1).QLR = 0;
            malha[ind].cell(ncel).QL = 0;
        }
    } else if (((*malha[ind].vg1dSP).chaverede == 1 && malha[ind].noextremo != 1 && malha[ind].titRev <= (*malha[ind].vg1dSP).localtiny)) {
        if (malha[ind].cell(ncel - 1).MliqiniR < 0) {
            malha[ind].cell(ncel - 1).MR = malha[ind].cell(ncel - 1).MliqiniR;
            malha[ind].cell(ncel - 1).term1R = 1;
            malha[ind].cell(ncel - 1).term2R = 0;
            malha[ind].cell(ncel).MC = malha[ind].cell(ncel - 1).MliqiniR;
            malha[ind].cell(ncel).term1 = 1;
            malha[ind].cell(ncel).term2 = 0;
        }
    }
}

void corrigeVazNoBuf(SProd *malha, int ind) {
    int ncel = malha[ind].ncel;

    if ((malha[ind].noextremo != 1 && malha[ind].titRev >= 1. - (*malha[ind].vg1dSP).localtiny)) {
        if ((malha[ind].fontemassPRBuf + malha[ind].fontemassCRBuf) < 0) {
            malha[ind].fontemassPRBuf = 0;
            malha[ind].fontemassCRBuf = 0;

            malha[ind].cell(ncel - 1).term1R = 0;
            malha[ind].cell(ncel - 1).term2R = 0;
            malha[ind].cell(ncel).term1 = 0;
            malha[ind].cell(ncel).term2 = 0;
        }
    } else if ((malha[ind].noextremo != 1 && malha[ind].titRev <= (*malha[ind].vg1dSP).localtiny)) {
        if (malha[ind].fontemassGRBuf < 0) {
            malha[ind].fontemassGRBuf = 0;
            malha[ind].cell(ncel - 1).term1R = 1;
            malha[ind].cell(ncel - 1).term2R = 0;
            malha[ind].cell(ncel).term1 = 1;
            malha[ind].cell(ncel).term2 = 0;
        }
    }
}

void SolveRedeTrans(SProd *malha, Rede &arqRede,
                    Vcr<int> &inativo, int indativo, int nrede) {

    int narq = arqRede.nsisprod;

    Vcr<int> colRev(narq);
    Vcr<double> alfRev(narq);
    Vcr<double> betRev(narq);
    Vcr<double> titRev(narq);
    // Vcr<double> razMast0(narq);//caso so Master
    // Vcr<double> razMast(narq);//caso so Master
    // Vcr<double> razMastCrit(narq);//caso so Master
    Vcr<int> celpos(narq);
    Vcr<int> kontasnp(narq, 1);

    for (int i = 0; i < narq; i++) {
        malha[i].modoPerm = 0;
        malha[i].config().tfinal = malha[i].tfinal = (*arqRede.vg1dSP).TmaxR;
        if (inativo[i] == 0) {
            malha[i].config().master1.razareaativ = 0.4;
            if (arqRede.malha[i].nafluente > 0) {
                if (malha[i].cell(0).acsr.tipo == 2) {
                    malha[i].cell(0).acsr.injm.FluidoPro = malha[i].cell(0).acsr.injl.FluidoPro;
                    malha[i].cell(0).acsr.injm.fluidocol = malha[i].cell(0).acsr.injl.fluidocol;
                } else if (malha[i].cell(0).acsr.tipo == 1) {
                    malha[i].cell(0).acsr.injm.FluidoPro = malha[i].cell(0).acsr.injg.FluidoPro;
                }
                int ind = arqRede.malha[i].afluente[0];

                int nderiva = arqRede.malha[ind].ncoleta;
                Vcr<int> ordCol(nderiva);
                vector<double> dcol;
                for (int icol = 0; icol < nderiva; icol++) {
                    int aux = arqRede.malha[ind].coleta[icol];
                    dcol.push_back(malha[aux].cell(0).duto.dia);
                }
                sort(dcol.begin(), dcol.end());
                Vcr<int> carregado(nderiva, 0);
                for (int icol = 0; icol < nderiva; icol++) {
                    Vcr<int> colRev(narq);
                    int aux = arqRede.malha[ind].coleta[icol];
                    for (int icol2 = 0; icol2 < nderiva; icol2++) {
                        if (dcol[icol2] == malha[aux].cell(0).duto.dia && carregado[icol2] == 0) {
                            ordCol[icol2] = aux;
                            carregado[icol2] = 1;
                            icol2 = nderiva;
                        }
                    }
                }
                if (i == ordCol[nderiva - 1]) {
                    if (malha[i].cell(0).acsr.tipo != 10) {
                        malha[i].cell(0).acsr.injm.MassC = malha[i].cell(0).QLR * malha[i].cell(0).bet *
                                                             malha[i].cell(0).fluicol.MasEspFlu(malha[i].cell(0).pres, malha[i].cell(0).temp);
                        malha[i].cell(0).acsr.injm.MassP = malha[i].cell(0).MliqiniR -
                                                             malha[i].cell(0).acsr.injm.MassC;
                        malha[i].cell(0).acsr.injm.MassG = malha[i].cell(0).MR - malha[i].cell(0).MliqiniR;
                        malha[i].cell(0).acsr.tipo = 10;
                    }
                    malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injm.FluidoPro;
                } else {
                    malha[i].config().ConContEntrada = 1;
                    malha[i].presE = malha[ordCol[nderiva - 1]].cell(0).pres +
                                     (malha[ordCol[nderiva - 1]].cell(0).pres -
                                      malha[ordCol[nderiva - 1]].cell(1).presaux);
                    malha[i].tempE = malha[i].cell(0).acsr.injl.temp;
                    malha[i].titE = malha[i].cell(0).fontemassGR / (malha[i].cell(0).fontemassLR +
                                                                      malha[i].cell(0).fontemassCR + malha[i].cell(0).fontemassGR);
                    malha[i].betaE = malha[ordCol[nderiva - 1]].cell(0).bet;
                    malha[i].cell(0).acsr.injl.QLiq = 0.;
                    malha[i].cell(0).acsr.injg.QGas = 0.;
                    malha[i].cell(0).acsr.injm.MassC = 0.;
                    malha[i].cell(0).acsr.injm.MassP = 0.;
                    malha[i].cell(0).acsr.injm.MassG = 0.;
                    malha[i].cell(0).fluiL = &malha[ordCol[nderiva - 1]].cell(0).acsr.injm.FluidoPro;
                    malha[i].cell(0).MC = malha[i].cell(0).fontemassLR +
                                            malha[i].cell(0).fontemassCR + malha[i].cell(0).fontemassGR;
                    malha[i].cell(0).MCini = malha[i].cell(0).MC;
                    malha[i].cell(0).Mliqini = malha[i].cell(0).fontemassLR + malha[i].cell(0).fontemassCR;
                    malha[i].cell(0).Mliqini0 = malha[i].cell(0).Mliqini;
                    double rholMix = malha[i].cell(0).bet * malha[i].cell(0).fluicol.MasEspFlu(
                                                                  malha[i].cell(0).pres, malha[i].cell(0).temp);
                    rholMix += (1. - malha[i].cell(0).bet) * malha[i].cell(0).flui.MasEspLiq(
                                                                   malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QL = malha[i].cell(0).Mliqini / rholMix;
                    malha[i].cell(0).QLini = malha[i].cell(0).QL;
                    double rhog = malha[i].cell(0).flui.MasEspGas(
                        malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QG = (malha[i].cell(0).MC - malha[i].cell(0).Mliqini) / rhog;
                    malha[i].cell(0).QGini = malha[i].cell(0).QG;
                }
            } else if (malha[i].config().ConContEntrada == 1) {

                if (malha[i].cell(0).acsr.tipo == 2) {
                    malha[i].cell(0).acsr.injl.QLiq = 0.;
                    malha[i].cell(0).MC = malha[i].cell(0).fontemassLR +
                                            malha[i].cell(0).fontemassCR + malha[i].cell(0).fontemassGR;
                    malha[i].cell(0).MCini = malha[i].cell(0).MC;
                    malha[i].cell(0).Mliqini = malha[i].cell(0).fontemassLR + malha[i].cell(0).fontemassCR;
                    malha[i].cell(0).Mliqini0 = malha[i].cell(0).Mliqini;
                    double rholMix = malha[i].cell(0).bet * malha[i].cell(0).fluicol.MasEspFlu(
                                                                  malha[i].cell(0).pres, malha[i].cell(0).temp);
                    rholMix += (1. - malha[i].cell(0).bet) * malha[i].cell(0).flui.MasEspLiq(
                                                                   malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QL = malha[i].cell(0).Mliqini / rholMix;
                    malha[i].cell(0).QLini = malha[i].cell(0).QL;
                    double rhog = malha[i].cell(0).flui.MasEspGas(
                        malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QG = (malha[i].cell(0).MC - malha[i].cell(0).Mliqini) / rhog;
                    malha[i].cell(0).QGini = malha[i].cell(0).QG;
                    malha[i].cell(0).fontemassLR = 0.;
                    malha[i].cell(0).fontemassCR = 0.;
                    malha[i].cell(0).fontemassGR = 0.;
                    malha[i].cell(0).fluiL = &malha[i].cell(0).flui;
                } else if (malha[i].cell(0).acsr.tipo == 1) {
                    malha[i].cell(0).acsr.injg.QGas = 0.;
                    malha[i].cell(0).MC = malha[i].cell(0).fontemassGR;
                    malha[i].cell(0).MCini = malha[i].cell(0).MC;
                    malha[i].cell(0).Mliqini = 0.;
                    malha[i].cell(0).Mliqini0 = malha[i].cell(0).Mliqini;
                    malha[i].cell(0).QL = 0;
                    malha[i].cell(0).QLini = malha[i].cell(0).QL;
                    double rhog = malha[i].cell(0).flui.MasEspGas(
                        malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QG = (malha[i].cell(0).MC - malha[i].cell(0).Mliqini) / rhog;
                    malha[i].cell(0).QGini = malha[i].cell(0).QG;
                    malha[i].cell(0).fontemassLR = 0.;
                    malha[i].cell(0).fontemassCR = 0.;
                    malha[i].cell(0).fontemassGR = 0.;
                    malha[i].cell(0).fluiL = &malha[i].cell(0).flui;
                }
            }
        }
    }

    Vcr<double> fonteG(narq);
    Vcr<double> fonteP(narq);
    Vcr<double> fonteC(narq);
    Vcr<double> injMG(narq);
    Vcr<double> injMP(narq);
    Vcr<double> injMC(narq);
    Vcr<double> rpR(narq);
    Vcr<double> rcR(narq);
    Vcr<double> rpRi(narq);
    Vcr<double> rcRi(narq);
    Vcr<double> temperatura(narq);
    Vcr<double> injTempIni(narq);
    Vcr<double> rDgD(narq);
    Vcr<double> rDgL(narq);
    Vcr<int> condcon(narq);
    ProFlu *fluiInjM;
    fluiInjM = new ProFlu[narq];
    ProFlu *fluiFim;
    fluiFim = new ProFlu[narq];
    ProFlu **fluiIni;
    fluiIni = new ProFlu *[narq];

    for (int i = 0; i < narq; i++) {
        if (inativo[i] == 0) {
            malha[i].titRev = malha[i].cell(malha[i].ncel).alf;
            malha[i].betaRev = malha[i].cell(malha[i].ncel).bet;
        }
    }

    Vcr<int> reverse(narq, 0);
    for (int i = 0; i < narq; i++) {
        if (reverse[i] == 0 && arqRede.malha[i].ncoleta > 0 && inativo[i] == 1) {
            int indCol = arqRede.malha[i].coleta[0];
            double qg = 0;
            double ql = 0;
            for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                int indCol2 = arqRede.malha[i].coleta[j];
                if (inativo[indCol2] == 1) {
                    if (malha[indCol2].cell(1).QG < 0)
                        qg -= malha[indCol2].cell(1).QG;
                    if (malha[indCol2].cell(1).QL < 0)
                        ql -= malha[indCol2].cell(1).QL;
                }
            }
            for (int j = 0; j < arqRede.malha[indCol].nafluente; j++) {
                int indAflu = arqRede.malha[i].afluente[j];
                int fim = malha[indAflu].ncel;
                if (inativo[indAflu] == 1) {
                    if (malha[indAflu].cell(fim).QG > 0)
                        qg += malha[indAflu].cell(fim).QG;
                    if (malha[indAflu].cell(fim).QL > 0)
                        ql += malha[indAflu].cell(fim).QL;
                }
            }
            for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                int indCol2 = arqRede.malha[i].coleta[j];
                if (fabs(qg) > 1e-15 || fabs(ql) > 1e-15)
                    alfRev[indCol2] = qg / (qg + ql);
                else
                    alfRev[indCol2] = 1.;
                reverse[indCol2] = 1;
            }
            for (int j = 0; j < arqRede.malha[indCol].nafluente; j++) {
                int indAflu = arqRede.malha[i].afluente[j];
                if (fabs(qg) > 1e-15 || fabs(ql) > 1e-15)
                    alfRev[indAflu] = qg / (qg + ql);
                else
                    alfRev[indAflu] = 1.;
                reverse[indAflu] = 1;
            }
        }
    }

    double dtteste;
    (*arqRede.vg1dSP).lixo5R = 0.;
    while ((*arqRede.vg1dSP).lixo5R < (*arqRede.vg1dSP).TmaxR) {

        for (int i = 0; i < narq; i++) {
            if (malha[i].config().calculaEnvelope == 1 && (*arqRede.vg1dSP).lixo5 <= (*arqRede.vg1dSP).TmaxR) { //*vg1dSP).lixo5>0 && //chris - Hidratos
                FA_Hidrato solverHidrato(malha[i]);
                solverHidrato.solverHidrato();
            }

            if (malha[i].config().lingas > 0 && malha[i].config().calculaEnvelope == 1 && (*arqRede.vg1dSP).lixo5 <= (*arqRede.vg1dSP).TmaxR) { //*vg1dSP).lixo5>0 && //chris - Hidratos
                FA_Hidrato_Servico solverHidratoG(malha[i]);
                solverHidratoG.solverHidratoG();
            }

            if (malha[i].config().flashCompleto == 2 && (*arqRede.vg1dSP).lixo5 < 1e-15) {
                malha[i].atualizaMiniTab();
            }

            if ((*arqRede.vg1dSP).lixo5 < 1e-15) {
                for (int j = 0; j < malha[i].config().ntendp; j++) {
                    malha[i].config().imprimeTrend(malha[i].celula, malha[i].MatTrendP[j], (*arqRede.vg1dSP).lixo5, i, malha[i].ntrend[j]);
                }
                malha[i].renovaTemp();
            }
        }

        dtteste = 1e6;
        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {
                malha[i].modeloCompleto = malha[i].config().correcaoMassaEspLiq;
                malha[i].determinaDT();
                malha[i].dtauxCFL = malha[i].dt;
                malha[i].dtauxFinal = malha[i].dt;
                if (malha[i].dt < dtteste)
                    dtteste = malha[i].dt;
            }
        }
        dtteste = dtteste / 1.;
        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {
                malha[i].dt = dtteste;
                for (int j = 0; j <= malha[i].ncel; j++) {
                    malha[i].cell(j).dt = malha[i].dt;
                    malha[i].cell(j).dt2 = malha[i].dt;
                    malha[i].cell(j).dtPig = malha[i].dt;
                }
            }
        }
        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {
                malha[i].reinicia = 0;
                celpos[i] = malha[i].config().master1.posic;
                // razMast0[i]=malha[i].cell(celpos[i]).acsr.chk.AreaGarg/malha[i].cell(celpos[i]).duto.area;//caso so Master
                malha[i].aberturaVal0(); // caso varias valvulas
                malha[i].solveLinGas();
                int presinterna = malha[i].noextremo;
                if (arqRede.malha[i].presimposta == 1)
                    presinterna = 1;
                malha[i].config().atualiza(malha[i].noinicial, presinterna, malha[i].derivaAnel,
                                      malha[i].chokeSup,
                                      malha[i].chokeInj, malha[i].celula,
                                      malha[i].celulaG, malha[i].pGSup, malha[i].temperatura,
                                      malha[i].presiniG, malha[i].tempiniG,
                                      malha[i].presE, malha[i].tempE, malha[i].titE, malha[i].betaE, (*arqRede.vg1dSP).lixo5R);
                malha[i].atualizaCC1();
                // razMastCrit[i]=0.5;//caso so Master
                for (int j = 0; j <= malha[i].config().nvalv; j++)
                    malha[i].vRazMastCrit[j] = 0.5; // caso varias valvulas
                // razMast[i]=malha[i].cell(celpos[i]).acsr.chk.AreaGarg/malha[i].cell(celpos[i]).duto.area;//caso so Master
                malha[i].aberturaVal1(); // caso varias valvulas
                if (arqRede.malha[i].ncoleta > 0 &&
                    (malha[i].chokeSup.AreaGarg / malha[i].cell(malha[i].ncel - 1).duto.area) < 0.6)
                    malha[i].config().chkv = 1;
                if ((*arqRede.vg1dSP).lixo5R >= 4000) {
                }
                for (int j = 0; j <= malha[i].config().nvalv; j++)
                    if (malha[i].vRazMast1[j] != malha[i].vRazMast0[j])
                        malha[i].modeloCompleto = 0; // caso varias valvulas
                if (malha[i].modeloCompleto == 1)
                    malha[i].avaliaVariaDpDt(); // caso varias valvulas
                if (malha[i].modeloCompleto == 0)
                    malha[i].config().cicloAcopTerm = 0;
                else
                    malha[i].config().cicloAcopTerm = 1;
            }
        }
        int modeloCompletoGlob = 1;
        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {
                if (malha[i].modeloCompleto == 0)
                    modeloCompletoGlob = 0;
            }
        }
        int ciclomax = modeloCompletoGlob;

        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {

                malha[i].pGSupIni = malha[i].pGSup;
                malha[i].tGSupIni = malha[i].tGSup;
                malha[i].presEini = malha[i].presE;
                malha[i].tempEini = malha[i].tempE;
                malha[i].titEini = malha[i].titE;
                malha[i].betaEini = malha[i].betaE;
                malha[i].alfEini = malha[i].alfE;
                malha[i].titRevini = malha[i].titRev;
                malha[i].betaRevini = malha[i].betaRev;
                malha[i].presfimini = malha[i].presfim;
                condcon[i] = malha[i].config().ConContEntrada;
                temperatura[i] = malha[i].temperatura;
                injTempIni[i] = malha[i].cell(0).acsr.injm.temp;
                rDgD[i] = malha[i].cell(0).acsr.injm.FluidoPro.rDgD;
                rDgL[i] = malha[i].cell(0).acsr.injm.FluidoPro.rDgL;
                if (malha[i].cell(0).acsr.tipo == 10) {
                    injMG[i] = malha[i].cell(0).acsr.injm.MassG;
                    injMP[i] = malha[i].cell(0).acsr.injm.MassP;
                    injMC[i] = malha[i].cell(0).acsr.injm.MassC;
                    fluiInjM[i] = malha[i].cell(0).acsr.injm.FluidoPro;
                }
                fluiFim[i] = malha[i].cell(malha[i].ncel).flui;
                fluiIni[i] = malha[i].cell(0).fluiL;
                rpR[i] = malha[i].cell(0).rpR;
                rcR[i] = malha[i].cell(0).rcR;
                rpRi[i] = malha[i].cell(0).rpRi;
                rcRi[i] = malha[i].cell(0).rcRi;
                malha[i].abertoini = malha[i].aberto;
                malha[i].tempoabertoini = malha[i].tempoaberto;
                malha[i].masChkSupini = malha[i].masChkSup;
                malha[i].mudaModoChkini = malha[i].mudaModoChk;
            }
        }

        for (int kontaAcop = 0; kontaAcop <= 1 * modeloCompletoGlob; kontaAcop++) {
            int reiniGlob = 0;

            if ((*arqRede.vg1dSP).lixo5R >= 1000) {
            }

            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0) {

                    if (malha[i].modeloCompleto == 0) {
                        for (int j = 0; j <= malha[i].ncel; j++)
                            malha[i].cell(j).m2d = 0.;
                    } else {

                        for (int j = 0; j <= malha[i].ncel; j++) {
                            double razDp = 0.1;
                            double razDT = 1;
                            if (j < celpos[i] && malha[i].cell(celpos[i]).acsr.chk.AreaGarg < 1e-15 * malha[i].cell(celpos[i]).acsr.chk.AreaTub) {
                                razDT = 1;
                            } else if (j == celpos[i] + 1 &&
                                       malha[i].cell(celpos[i]).acsr.chk.AreaGarg < 1e-15 * malha[i].cell(celpos[i]).acsr.chk.AreaTub) {
                                razDT = 1;
                            }
                            if ((fabs(malha[i].cell(j).dpdtIni) / malha[i].cell(j).pres < razDp) && fabs(malha[i].cell(j).dTdtIni) < razDT) {
                                if (malha[i].TransMassModel == 0)
                                    malha[i].cell(j).m2d = 1.;
                                else
                                    malha[i].cell(j).m2d = 0.;
                                malha[i].cell(j).mudaDT = 1.;
                            } else {
                                malha[i].cell(j).m2d = 0.;
                                malha[i].cell(j).mudaDT = 0.;
                            }
                        }
                    }
                    if (malha[i].config().estabCol == 1) {
                        for (int j = 0; j <= celpos[i]; j++) {
                            malha[i].cell(j).m2d = 0.;
                            malha[i].cell(j).mudaDT = 0.;
                            malha[i].cell(j).estabCol = 1;
                        }
                    }

                    malha[i].EvoluiFrac(alfRev[i], betRev[i], kontaAcop);
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).acsr.tipo == 15) {
                            malha[i].cell(j).acsr.radialPoro.avancoSW(malha[i].dt);
                            if (malha[i].cell(j).acsr.radialPoro.reinicia == -1) {
                                if (malha[i].reinicia > -1)
                                    malha[i].reinicia = -1;
                                // celula[i].acsr.radialPoro.reavaliaDT(Ndt)
                            }
                        } else if (malha[i].cell(j).acsr.tipo == 16) {
                            malha[i].cell(j).acsr.poroso2D.avancoSW(malha[i].dt);
                            if (malha[i].cell(j).acsr.poroso2D.reinicia == -1) {
                                if (malha[i].reinicia > -1)
                                    malha[i].reinicia = -1;
                                // celula[i].acsr.radialPoro.reavaliaDT(Ndt)
                            }
                        }
                    }
                    if (malha[i].config().correcaoMassaEspLiq == 1) {
                        for (int j = 0; j < malha[i].ncel; j++)
                            malha[i].cell(j + 1).mudaDTL = malha[i].cell(j).mudaDT;
                    }

                    if (kontaAcop == 0 && malha[i].config().controleDTvalv == 1) {
                        // caso so Master
                        // caso so Master
                        malha[i].restringeDTporValv(); // caso varias valvulas
                    }
                }
            }
            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0) {
                    if (malha[i].reinicia == -1)
                        reiniGlob = -1;
                }
            }
            if (reiniGlob < 0) {
                for (int i = 0; i < narq; i++) {
                    if (inativo[i] == 0) {
                        malha[i].ReiniEvolFrac0();
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).acsr.tipo == 15) {
                                malha[i].cell(j).acsr.radialPoro.reavaliaDT(malha[i].dt);
                            }
                            if (malha[i].cell(j).acsr.tipo == 16) {
                                malha[i].cell(j).acsr.poroso2D.reavaliaDT(malha[i].dt);
                            }
                        }
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).acsr.tipo == 15) {
                                malha[i].cell(j).acsr.radialPoro.reiniciaEvoluiSW(malha[i].dt);
                            }
                            if (malha[i].cell(j).acsr.tipo == 16) {
                                malha[i].cell(j).acsr.poroso2D.reiniciaEvoluiSW(malha[i].dt);
                            }
                        }
                        if (malha[i].dt < dtteste) {
                            dtteste = malha[i].dt;
                        }
                    }
                }
                for (int i = 0; i < narq; i++) {
                    if (inativo[i] == 0) {
                        malha[i].dt = dtteste;
                        malha[i].dtauxFinal = malha[i].dt;
                        malha[i].ReiniEvolFrac();
                        malha[i].EvoluiFrac(alfRev[i], betRev[i], kontaAcop);
                        malha[i].reinicia = 0;
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).acsr.tipo == 15) {
                                malha[i].cell(j).acsr.radialPoro.avancoSWcorrec();
                            }
                            if (malha[i].cell(j).acsr.tipo == 16) {
                                malha[i].cell(j).acsr.poroso2D.avancoSWcorrec();
                            }
                        }
                    }
                }
            }
            if ((*arqRede.vg1dSP).lixo5R <= 1e-15)
                (*arqRede.vg1dSP).dtrede = dtteste;

            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0 && arqRede.malha[i].ncoleta > 0) {
                    int colMaster = buscaNoColetorMrestre(malha, arqRede, inativo, i);
                    colRev[i] = arqRede.malha[i].coleta[colMaster];
                    betRev[i] = malha[colRev[i]].cell(0).bet;
                    double presRev = malha[colRev[i]].cell(0).pres;
                    double tempRev = malha[colRev[i]].cell(0).temp;
                    double rholp = malha[colRev[i]].cell(0).flui.MasEspLiq(presRev, tempRev);
                    double rholc = malha[colRev[i]].cell(0).fluicol.MasEspFlu(presRev, tempRev);
                    double rholmix = (1 - betRev[i]) * rholp + betRev[i] * rholc;
                    double rhog = malha[colRev[i]].cell(0).flui.MasEspGas(presRev, tempRev);
                    double romix = alfRev[i] * rhog + (1 - alfRev[i]) * rholmix;
                    titRev[i] = alfRev[i] * rhog / romix;

                    celAfluFinal(i, colRev[i], malha, arqRede, inativo);
                } else if (inativo[i] == 0) {
                    titRev[i] = 1.;
                    alfRev[i] = 1.;
                    betRev[i] = 0.;
                }
                if (inativo[i] == 0) {
                    malha[i].AtualizaPig();
                    if (kontaAcop == 0)
                        malha[i].dtCicMin = malha[i].dt;
                    if (kontaAcop == 1 * modeloCompletoGlob)
                        malha[i].atenuaDtMax();

                    if (modeloCompletoGlob == 1) {
                        fonteC[i] = malha[i].cell(malha[i].ncel).fontemassCR;
                        fonteP[i] = malha[i].cell(malha[i].ncel).fontemassLR;
                        fonteG[i] = malha[i].cell(malha[i].ncel).fontemassGR;
                    }

                    malha[i].calcCCpres(titRev[i], alfRev[i], betRev[i]);
                    int temaflu = 0;
                    if (arqRede.malha[i].ncoleta > 0) {
                        temaflu = 1;
                        if (malha[i].masChkSup == 0 && malha[i].config().chkv == 0) {
                            malha[i].cell(malha[i].ncel).alf = alfRev[i];
                            malha[i].cell(malha[i].ncel).bet = betRev[i];
                        }
                    }
                    if (malha[i].cell(malha[i].ncel).alf < 0.05 && malha[i].masChkSup == 1)
                        malha[i].cell(malha[i].ncel).alf = 0.05;
                    // malha[i].cell(celpos[i]).alf=0.05;//caso so Master
                    // caso varias valvulas
                    for (int j = 0; j <= malha[i].config().nvalv; j++) {
                        int celposAux;
                        if (j > 0)
                            celposAux = malha[i].config().valv[j - 1].posicP;
                        else
                            celposAux = celpos[i];
                        if (malha[i].cell(celposAux).alf < 0.05 && malha[i].vRazMast1[j] <= malha[i].config().master1.razareaativ)
                            malha[i].cell(celposAux).alf = 0.05;
                    }
                    // caso varias valvulas
                    malha[i].renovaterm(temaflu);
                }
            }

            (*arqRede.vg1dSP).iterRedeT = 0;
            Vcr<int> segundoAcop(narq, 0);
            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0 && (malha[i].noextremo == 1 ||
                                        (malha[i].noextremo == 0 &&
                                         ((malha[i].chokeSup.AreaGarg / malha[i].cell(malha[i].ncel - 1).duto.area) > 0.1 ||
                                          malha[i].masChkSup == 0))))
                    segundoAcop[i] = 1;
            }
            while ((*arqRede.vg1dSP).iterRedeT < 1) {
                CicloRedeTrans(malha, arqRede, inativo, alfRev, betRev, titRev);
                for (int i = 0; i < narq; i++) {
                    if (segundoAcop[i]) {
                        malha[i].SolveAcopPV();
                        malha[i].renovaBuffer();
                        if (arqRede.malha[i].ncoleta > 0)
                            malha[i].calcCCBuffer(titRev[i], alfRev[i], betRev[i]);

                    } else {
                        malha[i].renovaBufferCego();
                    }
                }
                (*arqRede.vg1dSP).iterRedeT++;
            }
            CicloRedeTrans(malha, arqRede, inativo, alfRev, betRev, titRev);
            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0) {
                    malha[i].calcCCpres(titRev[i], alfRev[i], betRev[i]);
                    if (arqRede.malha[i].ncoleta > 0) {
                        malha[i].renovatermAfluFim();
                    }
                    if (arqRede.malha[i].nafluente > 0 && malha[i].config().ConContEntrada > 0)
                        malha[i].renovatermColIni();
                    malha[i].SolveAcopPV();

                    if (kontaAcop < 1 * malha[i].modeloCompleto && malha[i].modeloCompleto == 1) {
                        for (int j = 0; j <= malha[i].ncel; j++)
                            malha[i].cell(j).dpdt =
                                1 * (malha[i].termolivreP[2 * j + 1] - malha[i].cell(j).pres) / malha[i].cell(j).dt;
                    }
                    malha[i].renova();
                    if (malha[i].config().cicloAcopTerm == 1) {
                        if (kontaAcop < 1 * malha[i].modeloCompleto)
                            for (int j = 0; j <= malha[i].ncel; j++)
                                malha[i].cell(j).dpdt = malha[i].cell(j).d2pdt2;
                        malha[i].marchaEnergTrans(kontaAcop, ciclomax);
                    }
                    if (arqRede.malha[i].nafluente > 0) {
                        malha[i].cell(0).dTdt = 0.;
                        malha[i].cell(0).dTdtL = 0.;
                        malha[i].cell(0).dpdt = 0.;
                        malha[i].cell(0).m2d = 0.;
                        malha[i].cell(0).mudaDT = 0.;
                        malha[i].cell(1).dTdt = 0.;
                        malha[i].cell(1).dTdtL = 0.;
                        malha[i].cell(1).dpdt = 0.;
                        malha[i].cell(1).m2d = 0.;
                        malha[i].cell(1).mudaDT = 0.;
                    }
                    double razChoke = (malha[i].chokeSup.AreaGarg / malha[i].cell(malha[i].ncel - 1).duto.area);
                    if ((razChoke > 1e-15 && malha[i].noextremo == 0)) {
                        malha[i].cell(malha[i].ncel).dTdt = 0.;
                        malha[i].cell(malha[i].ncel).dTdtL = 0.;
                        malha[i].cell(malha[i].ncel).dpdt = 0.;
                        malha[i].cell(malha[i].ncel).m2d = 0.;
                        malha[i].cell(malha[i].ncel).mudaDT = 0.;
                        malha[i].cell(malha[i].ncel - 1).dTdt = 0.;
                        malha[i].cell(malha[i].ncel - 1).dTdtL = 0.;
                        malha[i].cell(malha[i].ncel - 1).dpdt = 0.;
                        malha[i].cell(malha[i].ncel - 1).m2d = 0.;
                        malha[i].cell(malha[i].ncel - 1).mudaDT = 0.;
                    }
                    if (kontaAcop != 1 * modeloCompletoGlob) {
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            malha[i].cell(j).FeiticoDoTempo2();
                            if (malha[i].cell(j).acsr.tipo == 15) {
                                malha[i].cell(j).acsr.radialPoro.FeiticoDoTempoSW();
                            } else if (malha[i].cell(j).acsr.tipo == 16) {
                                malha[i].cell(j).acsr.poroso2D.FeiticoDoTempoSW();
                            }
                        }

                        if (malha[i].modeloCompleto == 0) {
                            for (int j = 0; j <= malha[i].ncel; j++) {
                                malha[i].cell(j).FeiticoDoTempo3();
                                if (malha[i].cell(j).acsr.tipo == 15) {
                                    malha[i].cell(j).acsr.radialPoro.FeiticoDoTempoPQ();
                                } else if (malha[i].cell(j).acsr.tipo == 16) {
                                    malha[i].cell(j).acsr.poroso2D.FeiticoDoTempoPQ();
                                }
                            }
                        }

                        malha[i].cell(malha[i].ncel).fontemassCR = fonteC[i];
                        malha[i].cell(malha[i].ncel).fontemassLR = fonteP[i];
                        malha[i].cell(malha[i].ncel).fontemassGR = fonteG[i];
                        malha[i].pGSup = malha[i].pGSupIni;
                        malha[i].tGSup = malha[i].tGSupIni;
                        malha[i].presfim = malha[i].presfimini;
                        malha[i].tempE = malha[i].tempEini;
                        malha[i].titE = malha[i].titEini;
                        malha[i].betaE = malha[i].betaEini;
                        malha[i].alfE = malha[i].alfEini;
                        malha[i].presE = malha[i].presEini;
                        malha[i].betaRev = malha[i].betaRevini;
                        malha[i].titRev = malha[i].titRevini;
                        malha[i].temperatura = temperatura[i];
                        malha[i].cell(malha[i].ncel).flui = fluiFim[i];
                        malha[i].cell(0).fluiL = fluiIni[i];
                        malha[i].cell(0).rpR = rpR[i];
                        malha[i].cell(0).rcR = rcR[i];
                        malha[i].cell(0).rpRi = rpRi[i];
                        malha[i].cell(0).rcRi = rcRi[i];
                        malha[i].config().ConContEntrada = condcon[i];
                        malha[i].aberto = malha[i].abertoini;
                        malha[i].tempoaberto = malha[i].tempoabertoini;
                        malha[i].masChkSup = malha[i].masChkSupini;
                        malha[i].mudaModoChk = malha[i].mudaModoChkini;
                        if (malha[i].cell(0).acsr.tipo == 10) {
                            malha[i].cell(0).acsr.injm.MassG = injMG[i];
                            malha[i].cell(0).acsr.injm.MassP = injMP[i];
                            malha[i].cell(0).acsr.injm.MassC = injMC[i];
                            malha[i].cell(0).acsr.injm.FluidoPro = fluiInjM[i];
                            malha[i].cell(0).acsr.injm.temp = injTempIni[i];
                            malha[i].cell(0).acsr.injm.FluidoPro.rDgD = rDgD[i];
                            malha[i].cell(0).acsr.injm.FluidoPro.rDgL = rDgL[i];
                        }
                    }
                }
            }
        }

        for (int i = 0; i < narq; i++) {
            if (reverse[i] == 0 && arqRede.malha[i].ncoleta > 0 && inativo[i] == 1) {
                int indCol = arqRede.malha[i].coleta[0];
                double qg = 0;
                double ql = 0;
                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                    int indCol2 = arqRede.malha[i].coleta[j];
                    if (inativo[indCol2] == 1) {
                        if (malha[indCol2].cell(1).QG < 0)
                            qg -= malha[indCol2].cell(1).QG;
                        if (malha[indCol2].cell(1).QL < 0)
                            ql -= malha[indCol2].cell(1).QL;
                    }
                }
                for (int j = 0; j < arqRede.malha[indCol].nafluente; j++) {
                    int indAflu = arqRede.malha[i].afluente[j];
                    int fim = malha[indAflu].ncel;
                    if (inativo[indAflu] == 1) {
                        if (malha[indAflu].cell(fim).QG > 0)
                            qg += malha[indAflu].cell(fim).QG;
                        if (malha[indAflu].cell(fim).QL > 0)
                            ql += malha[indAflu].cell(fim).QL;
                    }
                }
                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                    int indCol2 = arqRede.malha[i].coleta[j];
                    if (fabs(qg) > 1e-15 || fabs(ql) > 1e-15)
                        alfRev[indCol2] = qg / (qg + ql);
                    else
                        alfRev[indCol2] = malha[indCol].cell(0).alf;
                    reverse[indCol2] = 1;
                }
                for (int j = 0; j < arqRede.malha[indCol].nafluente; j++) {
                    int indAflu = arqRede.malha[i].afluente[j];
                    if (fabs(qg) > 1e-15 || fabs(ql) > 1e-15)
                        alfRev[indAflu] = qg / (qg + ql);
                    else
                        alfRev[indAflu] = malha[indCol].cell(0).alf;
                    reverse[indAflu] = 1;
                }
            }
        }
        for (int i = 0; i < narq; i++) {
            if (inativo[i] == 0) {
                if (malha[i].modeloCompleto == 0) {
                    malha[i].marchaEnergTrans(0, 0);
                }
                if (arqRede.malha[i].ncoleta > 0) {
                    int col = arqRede.malha[i].coleta[0];
                    malha[i].SolveTrans(titRev[i], alfRev[i], betRev[i], nrede, malha[col].cell(0).flui);
                } else
                    malha[i].SolveTrans(titRev[i], alfRev[i], betRev[i], nrede);
                for (int kontasnp = 0; kontasnp < malha[i].config().nsnp; kontasnp++) {
                    if (fabs((*arqRede.vg1dSP).lixo5 - malha[i].config().tempsnp[kontasnp]) < malha[i].dt) {
                        WriteSnapShot(malha[i], kontasnp, i); // registro do arquivo SNP
                    }
                }
            }
        }
        (*arqRede.vg1dSP).lixo5R += dtteste;
        (*arqRede.vg1dSP).dtrede = dtteste;
    }

    delete[] fluiInjM;
    delete[] fluiFim;
    delete[] fluiIni;

    for (int i = 0; i < narq; i++) {
        WriteSnapShot(malha[i], malha[i].config().nsnp, i);
    }
}

void aviso(int i) {
    cout << "*******************************************************************************" << endl;
    cout << "o tramo de indice " << i << " teve dificuldade de convergencia, por isto, sera desligado da rede" << endl;
    cout << "em certas situacoes, a dificuldade de convergencia de um tramo inserido em uma rede se deve ao" << endl;
    cout << "deste de fato nao conseguir escoar quando posto em competicao com outros tramos em uma rede de escoamento"
         << endl;
    cout << "talvez seja o caso de reavaliar as condicoes de escoamento deste tramo" << endl;
    cout << "*******************************************************************************" << endl;
    cout << endl;
}
void aviso2(int i) {
    cout << "*******************************************************************************" << endl;
    cout << "o tramo de indice " << i << " teve dificuldade de convergencia, por isto, sera� desligado da rede" << endl;
    cout << "Isto ocorreu logo na primeira iteracao, geralmente" << endl;
    cout << "em certas situacoes, a dificuldade de convergencia de um tramo inserido em uma rede se deve ao" << endl;
    cout << "deste de fato nao conseguir escoar quando posto em competicao com outros tramos em uma rede de escoamento"
         << endl;
    cout << "talvez seja o caso de reavaliar as condicoes de escoamento deste tramo" << endl;
    cout << "Porem, como o problema ocorreu na primeira iteracao, talvez seja um problema de chute inicial" << endl;
    cout << "Uma opcao para melhorar o chute e ir no arquivo de rede.json e diminuir a variavel ParametroInicial" << endl;
    cout << "*******************************************************************************" << endl;
    cout << endl;
}

void verificaTramoVazPres(int &ind, Rede &arqRede, SProd *malha) {
    int fimBusca = 0;
    while (fimBusca == 0) {
        if (arqRede.malha[ind].ncoleta > 0) {
            int aux = arqRede.malha[ind].coleta[0];
            if (arqRede.malha[ind].ncoleta != 1) {
                fimBusca = 1;
                ind = -1;
            } else if (malha[aux].config().ConContEntrada != 2) {
                ind = aux;
            } else {
                ind = aux;
                fimBusca = 1;
                string saida = "Erro na configuracao de redes. Sequencia de tramos interligados sem referencia de condicao de contorno a jusante,";
                saida += "um tramo a jusante sem divergencia de rede esta com condicao vazao-pressao, isto elimina a referencia de pressao a jusante, ";
                saida += "corrigir este tramo,  tramos  com condicao pressao-vazao : ";
                saida += arqRede.impfiles[aux];
                NumError(saida);
            }
        } else {
            fimBusca = 1;
        }
    }
}

void inativoColetor(int i, Rede &arqRede, Vcr<int> &inativo, int narq) {
    for (int j = 0; j < narq; j++) {
        if (arqRede.malha[j].nafluente > 0) {
            for (int k = 0; k < arqRede.malha[j].nafluente; k++) {
                if (arqRede.malha[j].afluente[k] == i) {
                    if (k == arqRede.malha[j].nafluente - 1)
                        arqRede.malha[j].nafluente -= 1;
                    else {
                        for (int k2 = k + 1; k2 < arqRede.malha[j].nafluente; k2++)
                            arqRede.malha[j].afluente[k2 - 1] = arqRede.malha[j].afluente[k2];
                        arqRede.malha[j].nafluente -= 1;
                    }
                }
            }
        }
    }
}

void inativoAfluente(int i, Rede &arqRede, Vcr<int> &inativo) {
    for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
        int k = arqRede.malha[i].afluente[j];
        if (arqRede.malha[k].ncoleta == 1) {
            inativo[k] = 1;
            if (arqRede.malha[k].nafluente > 0)
                inativoAfluente(k, arqRede, inativo);
        } else {
            for (int k1 = 0; k1 < arqRede.malha[j].ncoleta; k1++) {
                if (arqRede.malha[k].coleta[k1] == i) {
                    if (k1 == arqRede.malha[k].ncoleta - 1)
                        arqRede.malha[k].ncoleta -= 1;
                    else {
                        for (int k2 = k1 + 1; k2 < arqRede.malha[k].ncoleta; k2++)
                            arqRede.malha[k].coleta[k2 - 1] = arqRede.malha[k].coleta[k2];
                        arqRede.malha[k].coleta -= 1;
                    }
                }
            }
        }
    }
}
void avaliaBloq(int i, Rede &arqRede,
                Vcr<double> &nBloq, Vcr<double> &nBloq1, Vcr<double> &nBloq2,
                int &totBloq, int &totBloq1, int &totBloq2, int &col2, int naflu, vector<tramoPart> &bloq) {

    if (arqRede.malha[arqRede.malha[i].afluente[0]].ncoleta != 2) {
        for (int kaflu = 0; kaflu < naflu; kaflu++)
            nBloq[kaflu] = 1.;
        totBloq = naflu;
    } else {
        col2 = arqRede.malha[arqRede.malha[i].afluente[0]].coleta[0];
        if (col2 == i)
            col2 = arqRede.malha[arqRede.malha[i].afluente[0]].coleta[1];
        for (int kaflu = 0; kaflu < naflu; kaflu++) {
            if (bloq[i].part[kaflu] >= 0) {
                int match = 0;
                for (int kaflu2 = 0; kaflu2 < naflu; kaflu2++) {
                    if (bloq[col2].part[kaflu2] == bloq[i].part[kaflu]) {
                        nBloq[kaflu] = 1.;
                        totBloq++;
                        match = 1;
                        break;
                    }
                }
                if (match == 0) {
                    nBloq[kaflu] = 0.;
                }
            } else
                nBloq[kaflu] = 0.;
            if (bloq[i].part[kaflu] >= 0) {
                int match = 0;
                for (int kaflu2 = 0; kaflu2 < naflu; kaflu2++) {
                    if (bloq[col2].part[kaflu2] == bloq[i].part[kaflu]) {
                        match = 1;
                        break;
                    }
                }
                if (match == 0) {
                    nBloq1[kaflu] = 1.;
                    totBloq1++;
                } else
                    nBloq1[kaflu] = 0.;
            } else
                nBloq1[kaflu] = 0.;
            if (bloq[col2].part[kaflu] >= 0) {
                int match = 0;
                for (int kaflu2 = 0; kaflu2 < naflu; kaflu2++) {
                    if (bloq[col2].part[kaflu] == bloq[i].part[kaflu2]) {
                        match = 1;
                        break;
                    }
                }
                if (match == 0) {
                    nBloq2[kaflu] = 1.;
                    totBloq2++;
                } else
                    nBloq2[kaflu] = 0.;
            } else
                nBloq2[kaflu] = 0.;
        }
    }
}

int chutePresRede(int indprod, SProd *malha, Rede &arqRede, double chutehol,
                  Vcr<double> &razcolet, Vcr<double> &prescolet) {
    double vaz = (*arqRede.vg1dSP).somavaz * malha[indprod].cell(0).duto.area / (*arqRede.vg1dSP).somaarea;
    double vazG = (*arqRede.vg1dSP).somavazG * malha[indprod].cell(0).duto.area / (*arqRede.vg1dSP).somaarea;
    double presno = malha[indprod].hidroreverso(chutehol, vaz, vazG);
    int contarecur = 1;
    for (int i = 0; i < arqRede.malha[indprod].nafluente; i++) {
        int aflu = arqRede.malha[indprod].afluente[i];
        if (arqRede.malha[aflu].perm == 1) {
            if (arqRede.malha[aflu].presimposta == 0) {
                malha[aflu].pGSup = presno;
                razcolet[aflu] += 1.;
                prescolet[aflu] += presno;
            }
            if (arqRede.malha[aflu].nafluente > 0)
                contarecur += chutePresRede(aflu, malha, arqRede, chutehol, razcolet, prescolet);
        }
    }
    return contarecur;
}

void totalizaCicloRedeComp(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, int naflu, int i, int recb,
                           Vcr<double> nBloq, convergeNoPerm &noConv, int auxMaster, Vcr<int> coleta = Vcr<int>(), int ncoleta = 0) {
    Vcr<double> qostd(naflu + ncoleta);
    Vcr<double> rholliqIS(naflu + ncoleta);
    Vcr<double> rhogIS(naflu + ncoleta);
    Vcr<double> BSW(naflu + ncoleta);
    Vcr<double> Bet(naflu + ncoleta);
    Vcr<double> Rhocomp(naflu + ncoleta);
    Vcr<double> temp(naflu + ncoleta);
    Vcr<double> cpl(naflu + ncoleta);
    Vcr<double> cpg(naflu + ncoleta);
    Vcr<double> Mliq(naflu + ncoleta);
    Vcr<double> Qliq(naflu + ncoleta);
    Vcr<double> Mcomp(naflu + ncoleta);
    Vcr<double> Qcomp(naflu + ncoleta);
    Vcr<double> Mgas(naflu + ncoleta);
    Vcr<double> Denag(naflu + ncoleta);
    Vcr<double> titW(naflu + ncoleta);
    Vcr<double> vazMasLiqL(naflu + ncoleta);
    Vcr<double> vazMolV(naflu + ncoleta);
    for (int k = 0; k < naflu + ncoleta; k++) {
        vazMolV[k] = 0.;
    }
    noConv.flu = malha[i].cell(recb).flui;
    noConv.tL = 0.;
    noConv.tH = 70.;
    noConv.qostdTot = 0;
    noConv.moleomist = 0;
    noConv.moleomistNeg = 0;
    noConv.qwmist = 0;
    noConv.denagmist = 0;
    noConv.mliqmist = 0;
    noConv.mliqCmist = 0.;
    noConv.mgasmist = 0;
    noConv.cpmist = 0;
    noConv.tempmist = 0;
    noConv.LVisL = 0.;
    noConv.LVisH = 0.;
    noConv.Qlmist = 0.;
    noConv.betmist = 0.;
    noConv.qlmistStd = 0.;
    noConv.qgstdTot = 0.;
    noConv.moltot = 0.;
    noConv.mliqmistNeg = 0.;
    noConv.mliqCmistNeg = 0.;
    noConv.mgasmistNeg = 0.;
    noConv.qlmistStdNeg = 0.;
    noConv.qcmistStd = 0.;
    noConv.qcmistStdNeg = 0.;
    int kpos = 0;

    int aplicaFluiCol = 0;
    if ((*arqRede.vg1dSP).fluidoRede == 1)
        malha[i].cell(recb).acsr.injl.fluidocol = malha[i].cell(recb).fluicol;

    for (int k = 0; k < naflu; k++) {
        int ind = arqRede.malha[i].afluente[k];
        int fim = malha[ind].ncel - 1;
        if (nBloq[k] > 0. && malha[ind].cell(fim + 1).MC >= 0) {
            kpos++;
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                double ba = malha[ind].cell(fim).flui.BAFunc(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double fw = malha[ind].cell(fim).flui.BSW * ba /
                            (bo + ba * malha[ind].cell(fim).flui.BSW - malha[ind].cell(fim).flui.BSW * bo);
                double rhoO = malha[ind].cell(fim).flui.MasEspoleo(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double rhoW = malha[ind].cell(fim).flui.MasEspAgua(malha[ind].pGSup, malha[ind].cell(fim).temp);
                titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                    malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(fim).fluicol;
                    aplicaFluiCol = 1;
                }
                BSW[k] = malha[ind].cell(fim).flui.BSW;
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini * (1. - malha[ind].cell(fim).FW) *
                               (1. - malha[ind].cell(fim).bet) / (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                rhogIS[k] = malha[ind].cell(fim).flui.MasEspGas(pres, temp[k]);
                cpl[k] = (1. - Bet[k]) * malha[ind].cell(fim).flui.CalorLiq(pres, temp[k]) +
                         Bet[k] * malha[ind].cell(fim).fluicol.CalorLiq(pres, temp[k]);
                cpg[k] = malha[ind].cell(fim).flui.CalorGas(pres, temp[k]);
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;

                Denag[k] = malha[ind].cell(fim).flui.Denag;
                vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
                noConv.moleomist += vazMasLiqL[k];
                double pesoMolV = 0;
                for (int j = 0; j < malha[ind].cell(fim + 1).flui.npseudo; j++) {
                    pesoMolV += malha[ind].cell(fim + 1).flui.masMol[j] *
                                malha[ind].cell(fim + 1).flui.fracMol[j];
                }
                vazMolV[k] = (vazMasLiqL[k] + Mgas[k]) / pesoMolV;
                noConv.moltot += vazMolV[k];

                noConv.qostdTot += qostd[k];

                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.qwmist += qw1;
                noConv.denagmist += Denag[k] * qw1;
                noConv.mliqmist += Mliq[k];
                noConv.mliqCmist += Mcomp[k];
                noConv.TRmist += Mcomp[k] * malha[ind].cell(fim).fluicol.TR;
                noConv.mgasmist += Mgas[k];
                noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
                noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
                noConv.LVisL += qostd[k] * malha[ind].cell(fim).flui.VisOM(noConv.tL);
                noConv.LVisH += qostd[k] * malha[ind].cell(fim).flui.VisOM(noConv.tH);
                noConv.Qlmist += Qliq[k];
                noConv.qlmistStd += (qw1 + qostd[k] +
                                     Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                         malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStd += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                   malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        } else if (malha[ind].cell(fim + 1).MC < 0) {
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                BSW[k] = malha[ind].cell(fim).flui.BSW;
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini * (1. - malha[ind].cell(fim).FW) *
                               (1. - malha[ind].cell(fim).bet) / (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                Denag[k] = malha[ind].cell(fim).flui.Denag;
                double ba = malha[ind].cell(fim).flui.BAFunc(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double fw = malha[ind].cell(fim).flui.BSW * ba /
                            (bo + ba * malha[ind].cell(fim).flui.BSW - malha[ind].cell(fim).flui.BSW * bo);
                double rhoO = malha[ind].cell(fim).flui.MasEspoleo(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double rhoW = malha[ind].cell(fim).flui.MasEspAgua(malha[ind].pGSup, malha[ind].cell(fim).temp);
                titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
                vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
                noConv.moleomistNeg += vazMasLiqL[k];
                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.mliqmistNeg += Mliq[k];
                noConv.mliqCmistNeg += Mcomp[k];
                noConv.mgasmistNeg += Mgas[k];
                noConv.qlmistStdNeg += (qw1 + qostd[k] +
                                        Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                            malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStdNeg += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                        malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        }
    }
    for (int k = naflu; k < naflu + ncoleta; k++) {
        int ind = coleta[k - naflu];
        int iaflu = 0;
        int indAflu = arqRede.malha[i].afluente[iaflu];
        while (arqRede.malha[indAflu].perm == 0) {
            iaflu++;
            indAflu = arqRede.malha[i].afluente[iaflu];
        }
        if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1 && ind != auxMaster) {
            kpos++;
            int ini = 0;
            double bo = malha[ind].cell(ini).flui.BOFunc(
                malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double ba = malha[ind].cell(ini).flui.BAFunc(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double fw = malha[ind].cell(ini).flui.BSW * ba /
                        (bo + ba * malha[ind].cell(ini).flui.BSW - malha[ind].cell(ini).flui.BSW * bo);
            double rhoO = malha[ind].cell(ini).flui.MasEspoleo(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double rhoW = malha[ind].cell(ini).flui.MasEspAgua(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
            double pres = malha[indAflu].pGSup;
            malha[ind].calcTempFim();
            temp[k] = malha[ind].tempSup;
            Bet[k] = malha[ind].cell(ini).bet;

            if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(ini).fluicol;
                aplicaFluiCol = 1;
            }
            BSW[k] = malha[ind].cell(ini).flui.BSW;
            double rp = malha[ind].cell(ini).flui.MasEspLiq(pres, temp[k]);
            double rc = malha[ind].cell(ini).fluicol.MasEspFlu(pres, temp[k]);
            rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
            if (malha[ind].cell(ini + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                qostd[k] = -malha[ind].cell(ini + 1).Mliqini * (1. - malha[ind].cell(ini).FW) *
                           (1. - malha[ind].cell(ini).bet) / (bo * rholliqIS[k]);
            else
                qostd[k] = 0.;
            Rhocomp[k] = rc;
            rhogIS[k] = malha[ind].cell(ini).flui.MasEspGas(pres, temp[k]);
            cpl[k] = (1. - Bet[k]) * malha[ind].cell(ini).flui.CalorLiq(pres, temp[k]) +
                     Bet[k] * malha[ind].cell(ini).fluicol.CalorLiq(pres, temp[k]);
            cpg[k] = malha[ind].cell(ini).flui.CalorGas(pres, temp[k]);
            Mliq[k] = -malha[ind].cell(ini + 1).Mliqini;
            Qliq[k] = -malha[ind].cell(ini + 1).QL;
            Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
            Qcomp[k] = Qliq[k] * Bet[k];
            Mgas[k] = -malha[ind].cell(ini + 1).MC + malha[ind].cell(ini + 1).Mliqini;
            Denag[k] = malha[ind].cell(ini).flui.Denag;
            vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
            noConv.moleomist += vazMasLiqL[k];
            double pesoMolV = 0;
            for (int j = 0; j < malha[ind].cell(ini + 1).flui.npseudo; j++) {
                pesoMolV += malha[ind].cell(ini + 1).flui.masMol[j] *
                            malha[ind].cell(ini + 1).flui.fracMol[j];
            }
            vazMolV[k] = (vazMasLiqL[k] + Mgas[k]) / pesoMolV;
            noConv.moltot += vazMolV[k];

            noConv.qostdTot += qostd[k];

            double qw1;
            if ((1. - BSW[k]) > 0)
                qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
            else
                qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
            noConv.qwmist += qw1;
            noConv.denagmist += Denag[k] * qw1;
            noConv.mliqmist += Mliq[k];
            noConv.mliqCmist += Mcomp[k];
            noConv.TRmist += Mcomp[k] * malha[ind].cell(ini).fluicol.TR;
            noConv.mgasmist += Mgas[k];
            noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
            noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
            noConv.LVisL += qostd[k] * malha[ind].cell(ini).flui.VisOM(noConv.tL);
            noConv.LVisH += qostd[k] * malha[ind].cell(ini).flui.VisOM(noConv.tH);
            noConv.Qlmist += Qliq[k];
            noConv.qlmistStd += (qw1 + qostd[k] +
                                 Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.qcmistStd += (qw1 + qostd[k] +
                                 Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                               malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
        }
    }

    if (kpos > 0) {
        if (recb == 0 && noConv.moltot > 0) {
            for (int jgrup = 0; jgrup < malha[i].cell(recb).flui.npseudo; jgrup++) {
                malha[i].cell(recb).flui.fracMol[jgrup] = 0.;
                noConv.flu.fracMol[jgrup] = 0.;
            }
            for (int iaflu = 0; iaflu < naflu; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int fim = malha[ind].ncel - 1;
                for (int jgrup = 0; jgrup < malha[i].cell(recb).flui.npseudo; jgrup++) {
                    malha[i].cell(recb).flui.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            for (int iaflu = naflu; iaflu < naflu + ncoleta; iaflu++) {
                int ind = coleta[iaflu - naflu];
                int ini = 0;
                for (int jgrup = 0; jgrup < malha[i].cell(0).flui.npseudo; jgrup++) {
                    malha[i].cell(recb).flui.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            if (malha[i].config().tabelaDinamica == 0) {
                malha[i].cell(recb).flui.atualizaPropCompStandard();
                noConv.flu.atualizaPropCompStandard();
            }
        } else if (noConv.moltot > 0) {
            for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                noConv.flu.fracMol[jgrup] = 0.;
            }
            for (int iaflu = 0; iaflu < naflu; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int fim = malha[ind].ncel - 1;
                for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            for (int iaflu = naflu; iaflu < naflu + ncoleta; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int ini = 0;
                for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            if (malha[i].config().tabelaDinamica == 0)
                noConv.flu.atualizaPropCompStandard();
        } else {
            for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                noConv.flu.fracMol[jgrup] =
                    malha[auxMaster].fluiRevRede.fracMol[jgrup];
            }
            if (malha[i].config().tabelaDinamica == 0)
                noConv.flu.atualizaPropCompStandard();
        }

        if (noConv.qostdTot > (*arqRede.vg1dSP).localtiny) {
            noConv.LVisL = noConv.LVisL / noConv.qostdTot;
            noConv.LVisH = noConv.LVisH / noConv.qostdTot;
        } else {
            noConv.LVisL = malha[i].cell(recb).flui.VisOM(noConv.tL);
            noConv.LVisH = malha[i].cell(recb).flui.VisOM(noConv.tH);
        }
        if ((noConv.qostdTot + noConv.qwmist) > (*arqRede.vg1dSP).localtiny)
            noConv.bswmist = noConv.qwmist / (noConv.qostdTot + noConv.qwmist);
        else
            noConv.bswmist = malha[i].cell(recb).flui.BSW;
        if (noConv.qwmist > (*arqRede.vg1dSP).localtiny)
            noConv.denagmist = noConv.denagmist / noConv.qwmist;
        else
            noConv.denagmist = malha[i].cell(recb).flui.Denag;
        if (noConv.cpmist > (*arqRede.vg1dSP).localtiny)
            noConv.tempmist = noConv.tempmist / noConv.cpmist;
        else
            noConv.tempmist = 0.;
        if (noConv.qlmistStd > (*arqRede.vg1dSP).localtiny)
            noConv.betmist = noConv.betmist / noConv.qlmistStd;
        else
            noConv.betmist = 0.;
        if (noConv.mliqCmist > (*arqRede.vg1dSP).localtiny * 1e-15)
            noConv.TRmist /= noConv.mliqCmist;
        else
            noConv.TRmist = 0.;
    } else {
        for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
            noConv.flu.fracMol[jgrup] =
                malha[i].fluiRevRede.fracMol[jgrup];
        }
        if (malha[i].config().tabelaDinamica == 0)
            noConv.flu.atualizaPropCompStandard();
        noConv.LVisL = malha[i].fluiRevRede.LVisL;
        noConv.LVisH = malha[i].fluiRevRede.LVisH;
        noConv.bswmist = malha[i].fluiRevRede.BSW;
        noConv.denagmist = malha[i].fluiRevRede.Denag;
        noConv.betmist = 0.;
        int ind = arqRede.malha[i].afluente[0];
        int fim = malha[ind].ncel - 1;
        noConv.tempmist = malha[ind].cell(fim).temp;
        noConv.TRmist = malha[ind].cell(fim).fluicol.TR;
    }
    noConv.qlmistStd *= 86400;
    noConv.qlmistStdNeg *= 86400;
    noConv.qcmistStd *= 86400;
    noConv.qcmistStdNeg *= 86400;

    if (recb == 0) {
        malha[i].cell(recb).flui.BSW = noConv.bswmist;
        malha[i].cell(recb).flui.TempL = noConv.tL;
        malha[i].cell(recb).flui.TempH = noConv.tH;
        malha[i].cell(recb).flui.LVisL = noConv.LVisL;
        malha[i].cell(recb).flui.LVisH = noConv.LVisH;
    } else {
        noConv.flu.BSW = noConv.bswmist;
        noConv.flu.TempL = noConv.tL;
        noConv.flu.TempH = noConv.tH;
        noConv.flu.LVisL = noConv.LVisL;
        noConv.flu.LVisH = noConv.LVisH;
    }
}

double cicloRedeComp(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, vector<noRede> &normaEvol, vector<tramoPart> &bloq) {
    double norma = 0.;
    int nnos = 0;
    int narq = arqRede.nsisprod - 0 * indativo;
    (*arqRede.vg1dSP).relax = arqRede.relax;
    Vcr<int> Resolv(narq, 0);
    Vcr<int> IndNorma(narq, 0);
    Vcr<int> Retorna(narq, 0);
    if ((*arqRede.vg1dSP).iterRede >= 5) {
    }
    cout << "#################Iniciando ciclo iterativo: " << (*arqRede.vg1dSP).iterRede << endl;
    int ciclomalha = 0;
    double valor;
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].perm == 0) {
            Resolv[i] = 1;
            ciclomalha++;
        }
    }
    int testaNMalha = ciclomalha - narq + 1;
    while (ciclomalha < narq - 1) {
        int ntrdLocal = (*arqRede.vg1dSP).ntrd;
        if ((*arqRede.vg1dSP).ntrd > narq)
            ntrdLocal = narq;
#pragma omp parallel for num_threads(ntrdLocal)
        for (int i = 0; i < narq; i++) {
            int permAflu = 0;
            if (arqRede.malha[i].nafluente > 0) {
                for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                    int indaflu = arqRede.malha[i].afluente[iaflu];
                    permAflu += arqRede.malha[indaflu].perm;
                }
            }
            if ((arqRede.malha[i].nafluente == 0 || (permAflu == 0 && arqRede.malha[i].perm == 1)) &&
                Resolv[i] == 0 && inativo[i] == 0) {
                if (arqRede.malha[i].perm == 1) {
                    (*arqRede.vg1dSP).qualTramo = i;
                    if (malha[i].config().ConContEntrada == 0) {
                        int reverso = 0;
                        if (malha[i].cell(0).acsr.tipo == 1) {
                            if (malha[i].cell(0).acsr.injg.QGas < 0)
                                reverso = 1;
                        } else if (malha[i].cell(0).acsr.tipo == 2) {
                            if (malha[i].cell(0).acsr.injl.QLiq < 0)
                                reverso = 1;
                        } else if (malha[i].cell(0).acsr.tipo == 10) {
                            if ((malha[i].cell(0).acsr.injm.MassC +
                                 malha[i].cell(0).acsr.injm.MassG + malha[i].cell(0).acsr.injm.MassP) < 0)
                                reverso = 1;
                        }
                        malha[i].modoPerm = 1;
                        if (reverso == 0) {
                            if (malha[i].config().chokep.abertura[0] > 0.6) {
                                if ((*arqRede.vg1dSP).iterRede == 0) { // mudancaChute

                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 0 && malha[i].gasCell(0).tipoCC == 0)
                                        malha[i].config().gasinj.vazgas[0] = 150000 * malha[i].gasCell(0).duto.area / (*arqRede.vg1dSP).arearef;
                                    if (malha[i].config().lingas == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        double ciclo = 1.1e9;
                                        int konta = 0;
                                        double multVazGas;
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm();
                                        double testaPres1 = malha[i].gasCell(0).pres;
                                        malha[i].config().gasinj.vazgas[0] *= 1.05;
                                        malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        double testaPres2 = malha[i].gasCell(0).pres;
                                        if (testaPres1 < testaPres2) {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 1.05;
                                            else
                                                multVazGas = 0.95;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        } else {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 0.95;
                                            else
                                                multVazGas = 1.05;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        }
                                        while (ciclo > 0.9e9 && konta < 10) {
                                            malha[i].gasCell(0).tipoCC = 1;
                                            if (konta > 0)
                                                malha[i].buscaProdPfundoPerm();
                                            malha[i].gasCell(0).tipoCC = 0;
                                            ciclo = malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres, konta);
                                            if (ciclo > 0.9e9) {
                                                malha[i].config().gasinj.vazgas[0] *= multVazGas;
                                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                                            }
                                            konta++;
                                            if (konta > 10) {
                                                ciclo = 0.8e10;
                                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                                           "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################",
                                                           "", "");
                                            }
                                        }
                                    } // mudancaChute
                                    else
                                        valor = malha[i].buscaProdPfundoPerm();
                                } else { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        malha[i].gasCell(0).tipoCC = 0;
                                    }
                                    valor = malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                    //}//mudancaChute
                                }
                            } else {
                                if ((*arqRede.vg1dSP).iterRede == 0) { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 0 && malha[i].gasCell(0).tipoCC == 0)
                                        malha[i].config().gasinj.vazgas[0] = 150000 * malha[i].gasCell(0).duto.area / (*arqRede.vg1dSP).arearef;
                                    if (malha[i].config().lingas == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        double ciclo = 1.1e9;
                                        int konta = 0;
                                        double multVazGas;
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm2();
                                        double testaPres1 = malha[i].gasCell(0).pres;
                                        malha[i].config().gasinj.vazgas[0] *= 1.05;
                                        malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        double testaPres2 = malha[i].gasCell(0).pres;
                                        if (testaPres1 < testaPres2) {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 1.05;
                                            else
                                                multVazGas = 0.95;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        } else {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 0.95;
                                            else
                                                multVazGas = 1.05;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        }
                                        while (ciclo > 0.9e9 && konta < 10) {
                                            malha[i].gasCell(0).tipoCC = 1;
                                            if (konta > 0)
                                                malha[i].buscaProdPfundoPerm2();
                                            malha[i].gasCell(0).tipoCC = 0;
                                            ciclo = malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres, konta);
                                            if (ciclo > 0.9e9) {
                                                malha[i].config().gasinj.vazgas[0] *= multVazGas;
                                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                                            }
                                            konta++;
                                            if (konta > 10) {
                                                ciclo = 0.8e10;
                                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                                           "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################",
                                                           "", "");
                                            }
                                        }
                                    } // mudancaChute
                                    else
                                        valor = malha[i].buscaProdPfundoPerm2();
                                } else { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        malha[i].gasCell(0).tipoCC = 0;
                                    }
                                    valor = malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                    //}//mudancaChute
                                }
                            }
                        } else {
                            if ((*arqRede.vg1dSP).iterRede == 0)
                                valor = malha[i].buscaProdPfundoPermRev();
                            else
                                valor = malha[i].buscaProdPfundoPermRev(malha[i].cell(0).pres);
                        }
                    } else {
                        double chutemass = 0.;
                        if (malha[i].cell(0).acsr.tipo == 2)
                            chutemass = malha[i].cell(0).acsr.injl.QLiq;
                        else
                            chutemass = malha[i].cell(0).acsr.injg.QGas;
                        if (malha[i].config().chokep.abertura[0] > 0.6)
                            malha[i].buscaProdPresPresPerm(chutemass);
                        else
                            malha[i].buscaProdPresPresPerm2(chutemass);
                    }
                } else {
                    valor = 0;
                    int icol = arqRede.malha[i].coleta[0];
                    malha[i].pGSup = malha[icol].cell(0).pres;
                    malha[i].tGSup = malha[icol].cell(0).temp;
                    malha[i].tempSup = malha[icol].cell(0).temp;
                }
                if (valor < -1e9 || valor > 1e9) {
                    inativo[i] = 1;
                    inativoColetor(i, arqRede, inativo, narq);
                    if (valor < -1e9) {
                        aviso(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    } else {
                        aviso2(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    }
                    (*arqRede.vg1dSP).restartRede = 1;
                    Retorna[i] = -1;
                }
                Resolv[i] = 1;
#pragma omp atomic
                ciclomalha++;
                cout << "convergencia tramo: " << i << endl;
            }
        }
        int resolvGlob = 0;
        for (int i = 0; i < narq; i++) {
            resolvGlob += Resolv[i];
            if (Retorna[i] < 0)
                return 9000.;
        }

        while (resolvGlob < narq) {
            for (int i = 0; i < narq; i++) {
                int naflu = arqRede.malha[i].nafluente;
                int cicloaflu = 0;
                if (arqRede.malha[i].nafluente > 0 && Resolv[i] == 0) {
                    for (int j = 0; j < narq; j++) {
                        for (int k = 0; k < naflu; k++) {
                            if (arqRede.malha[i].afluente[k] == j && Resolv[j] == 1)
                                cicloaflu++;
                        }
                    }
                    if (cicloaflu == naflu && inativo[i] == 0) {

                        int ind = arqRede.malha[i].afluente[0];
                        int nderiva = arqRede.malha[ind].ncoleta;

                        int trocaDeriva = 0;
                        int novaDeriva;
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (arqRede.malha[aux].principal == 1) {
                                trocaDeriva = 1;
                                novaDeriva = aux;
                            }
                        }

                        Vcr<int> ordCol(nderiva);
                        vector<double> dcol;
                        Vcr<int> carregado(nderiva, 0);
                        if (trocaDeriva == 0) {
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                if (malha[aux].config().perm == 1) {
                                    double mult = 1.0;
                                    if (malha[aux].config().ConContEntrada == 2)
                                        mult = -1.;
                                    dcol.push_back(malha[aux].cell(0).duto.dia * mult);
                                }
                            }
                            sort(dcol.begin(), dcol.end());

                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                if (malha[aux].config().perm == 1) {
                                    double mult = 1.0;
                                    if (malha[aux].config().ConContEntrada == 2)
                                        mult = -1.;
                                    for (int icol2 = 0; icol2 < nderiva; icol2++) {
                                        if (fabs(dcol[icol2] - malha[aux].cell(0).duto.dia * mult) < 1.e-15 && carregado[icol2] == 0) {
                                            ordCol[icol2] = aux;
                                            carregado[icol2] = 1;
                                            icol2 = nderiva;
                                        }
                                    }
                                }
                            }
                        } else {
                            int icol2 = 0;
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                if (malha[aux].config().perm == 1) {
                                    if (aux != novaDeriva) {
                                        ordCol[icol2] = aux;
                                        carregado[icol2] = 1;
                                        icol2++;
                                    }
                                }
                            }
                            ordCol[nderiva - 1] = novaDeriva;
                            carregado[nderiva - 1] = 1;
                        }

                        convergeNoPerm noConv;
                        convergeNoPerm noConv1;
                        convergeNoPerm noConv2;
                        int totBloq = 0;
                        int totBloq1 = 0;
                        int totBloq2 = 0;
                        Vcr<double> nBloq(naflu);
                        Vcr<double> nBloq1(naflu);
                        Vcr<double> nBloq2(naflu);

                        int ncolneg = 0;
                        Vcr<int> colneg(nderiva);
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (malha[aux].config().perm == 1) {
                                if (malha[aux].cell(0).acsr.tipo == 1) {
                                    if (malha[aux].cell(0).acsr.injg.QGas < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                } else if (malha[aux].cell(0).acsr.tipo == 2) {
                                    if (malha[aux].cell(0).acsr.injl.QLiq < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                } else if (malha[aux].cell(0).acsr.tipo == 10) {
                                    if ((malha[aux].cell(0).acsr.injm.MassP + malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassC) < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                }
                            }
                        }

                        int auxMaster;
                        auxMaster = ordCol[nderiva - 1];

                        for (int iaflu = 0; iaflu < naflu; iaflu++) {
                            int indaflu = arqRede.malha[auxMaster].afluente[iaflu];
                            int celProp = 1 - verificaFonteDuplaReversa(malha, auxMaster);
                            malha[indaflu].fluiRevRede = malha[auxMaster].cell(celProp).flui;
                            malha[indaflu].tempRev = malha[auxMaster].cell(0).temp;
                            if ((*arqRede.vg1dSP).fluidoRede == 0)
                                malha[indaflu].config().razCompGasReves = malha[auxMaster].cell(celProp).acsr.injg.razCompGas;
                        }
                        for (int icol = 0; icol < nderiva - 1; icol++) {
                            int aux = ordCol[icol];
                            if (arqRede.malha[aux].ncoleta == 0)
                                malha[aux].fluiRevRede = malha[aux].config().flup[0];
                        }

                        int col2;
                        avaliaBloq(i, arqRede, nBloq, nBloq1, nBloq2, totBloq, totBloq1, totBloq2, col2, naflu, bloq);
                        if (totBloq > 0)
                            totalizaCicloRedeComp(malha, arqRede, inativo, indativo, naflu, i, 0, nBloq, noConv, auxMaster, colneg, ncolneg);
                        if (totBloq1 > 0) {
                            totalizaCicloRedeComp(malha, arqRede, inativo, indativo, naflu, i, 1, nBloq1, noConv1, auxMaster);
                            noConv1.qgstdTot = ((noConv1.mgasmist + noConv1.moleomist) /
                                                (noConv1.dengmist * 1.225)) *
                                               noConv1.flu.dStockTankVaporMassFraction;
                            noConv1.qgstdTotNeg = ((noConv1.mgasmistNeg + noConv1.moleomistNeg) /
                                                   (noConv1.dengmist * 1.225)) *
                                                  noConv1.flu.dStockTankVaporMassFraction;
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                malha[i].cell(1).acsr.tipo = 2;
                                if (malha[i].cell(1).acsr.injl.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injl.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injl.bet = noConv1.betmist;
                                malha[i].cell(1).acsr.injl.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injl.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injl.fluidocol.TR = noConv1.TRmist;

                                noConv1.qgstdTot *= 86400;

                                malha[i].cell(1).acsr.injl.QLiq = noConv1.qlmistStd + noConv1.qlmistStdNeg;
                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                malha[i].cell(1).acsr.tipo = 1;
                                if (malha[i].cell(1).acsr.injg.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injg.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injg.seco = 0;
                                malha[i].cell(1).acsr.injg.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injg.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injg.fluidocol.TR = noConv1.TRmist;
                                noConv1.qgstdTot *= 86400;
                                malha[i].cell(1).acsr.injg.QGas = noConv1.qgstdTot + noConv1.qgstdTotNeg * 86400.;
                                if (fabs(malha[i].cell(1).acsr.injg.QGas) > 1e-15)
                                    malha[i].cell(1).acsr.injg.razCompGas = (noConv1.qcmistStd + noConv1.qcmistStdNeg) /
                                                                              malha[i].cell(1).acsr.injg.QGas;
                                else
                                    malha[i].cell(1).acsr.injg.razCompGas = 0.;

                            } else {
                                malha[i].cell(1).acsr.tipo = 10;
                                if (malha[i].cell(1).acsr.injm.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injm.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injm.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injm.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injm.fluidocol.TR = noConv1.TRmist;
                                noConv1.qgstdTot *= 86400;

                                malha[i].cell(1).acsr.injm.MassP = noConv1.mliqmist - noConv1.mliqCmist +
                                                                     noConv1.mliqmistNeg - noConv1.mliqCmistNeg;
                                malha[i].cell(1).acsr.injm.MassC = noConv1.mliqCmist + noConv1.mliqCmistNeg;
                                malha[i].cell(1).acsr.injm.MassG = noConv1.mgasmist + noConv1.mgasmistNeg;
                            }
                        }
                        if (totBloq2 > 0) {
                            totalizaCicloRedeComp(malha, arqRede, inativo, indativo, naflu, col2, 1, nBloq2, noConv2, auxMaster);
                            noConv2.qgstdTot = ((noConv2.mgasmist + noConv2.moleomist) /
                                                (noConv2.dengmist * 1.225)) *
                                               noConv2.flu.dStockTankVaporMassFraction;
                            noConv2.qgstdTotNeg = ((noConv2.mgasmistNeg + noConv2.moleomistNeg) /
                                                   (noConv2.dengmist * 1.225)) *
                                                  noConv2.flu.dStockTankVaporMassFraction;
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                malha[col2].cell(1).acsr.tipo = 2;
                                if (malha[col2].cell(1).acsr.injl.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injl.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injl.bet = noConv2.betmist;
                                malha[col2].cell(1).acsr.injl.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injl.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injl.fluidocol.TR = noConv2.TRmist;

                                noConv2.qgstdTot *= 86400;

                                malha[col2].cell(1).acsr.injl.QLiq = noConv2.qlmistStd + noConv2.qlmistStdNeg;
                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                malha[col2].cell(1).acsr.tipo = 1;
                                if (malha[col2].cell(1).acsr.injg.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injg.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injg.seco = 0;
                                malha[col2].cell(1).acsr.injg.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injg.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injg.fluidocol.TR = noConv2.TRmist;
                                noConv2.qgstdTot *= 86400;
                                malha[col2].cell(1).acsr.injg.QGas = noConv2.qgstdTot + noConv2.qgstdTotNeg * 86400.;
                                if (fabs(malha[i].cell(1).acsr.injg.QGas) > 1e-15)
                                    malha[i].cell(1).acsr.injg.razCompGas = (noConv2.qcmistStd + noConv2.qcmistStdNeg) /
                                                                              malha[i].cell(1).acsr.injg.QGas;
                                else
                                    malha[i].cell(1).acsr.injg.razCompGas = 0.;
                            } else {
                                malha[col2].cell(1).acsr.tipo = 10;
                                if (malha[col2].cell(1).acsr.injm.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injm.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injm.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injm.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injm.fluidocol.TR = noConv2.TRmist;
                                noConv2.qgstdTot *= 86400;

                                malha[col2].cell(1).acsr.injm.MassP = noConv2.mliqmist - noConv2.mliqCmist +
                                                                        noConv2.mliqmistNeg - noConv2.mliqCmistNeg;
                                malha[col2].cell(1).acsr.injm.MassC = noConv2.mliqCmist + noConv2.mliqCmistNeg;
                                malha[col2].cell(1).acsr.injm.MassG = noConv2.mgasmist + noConv2.mgasmistNeg;
                            }
                        }

                        noConv.qgstdTot = ((noConv.mgasmist + noConv.moleomist) /
                                           (malha[ordCol[nderiva - 1]].cell(0).flui.Deng * 1.225)) *
                                          malha[ordCol[nderiva - 1]].cell(0).flui.dStockTankVaporMassFraction;
                        noConv.qgstdTotNeg = ((noConv.mgasmistNeg + noConv.moleomistNeg) /
                                              (malha[ordCol[nderiva - 1]].cell(0).flui.Deng * 1.225)) *
                                             malha[ordCol[nderiva - 1]].cell(0).flui.dStockTankVaporMassFraction;

                        for (int icol = 0; icol < nderiva; icol++) {
                            if (icol == nderiva - 1) {
                                int aux = ordCol[icol];
                                if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                    malha[aux].cell(0).acsr.tipo = 2;
                                    if (malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto == 0)
                                        malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto = 2;
                                    malha[aux].cell(0).acsr.injl.bet = noConv.betmist;
                                    malha[aux].cell(0).acsr.injl.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro = noConv.flu;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                    malha[aux].cell(0).acsr.injl.fluidocol.TR = noConv.TRmist;

                                    noConv.qgstdTot *= 86400;

                                    norma += pow((malha[aux].cell(0).acsr.injl.QLiq - (noConv.qlmistStd + noConv.qlmistStdNeg)) /
                                                     (noConv.qlmistStd + noConv.qlmistStdNeg),
                                                 2.);
                                    malha[aux].cell(0).acsr.injl.QLiq = noConv.qlmistStd + noConv.qlmistStdNeg;

                                } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                    malha[aux].cell(0).acsr.tipo = 1;
                                    if (malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto == 0)
                                        malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto = 2;
                                    malha[aux].cell(0).acsr.injg.seco = 0;
                                    malha[aux].cell(0).acsr.injg.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro = noConv.flu;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                    malha[aux].cell(0).acsr.injg.fluidocol.TR = noConv.TRmist;

                                    noConv.qgstdTot *= 86400;
                                    norma += pow((malha[aux].cell(0).acsr.injg.QGas - (noConv.qgstdTot + noConv.qgstdTotNeg * 86400.)) /
                                                     (noConv.qgstdTot + noConv.qgstdTotNeg * 86400.),
                                                 2.);
                                    malha[aux].cell(0).acsr.injg.QGas = noConv.qgstdTot + noConv.qgstdTotNeg * 86400.; //+qlmistStd;???????

                                    if (fabs(malha[aux].cell(0).acsr.injg.QGas) > 1e-15)
                                        malha[aux].cell(0).acsr.injg.razCompGas = (noConv.qcmistStd + noConv.qcmistStdNeg) /
                                                                                    malha[aux].cell(0).acsr.injg.QGas;
                                    else
                                        malha[aux].cell(0).acsr.injg.razCompGas = 0.;

                                } else {
                                    malha[aux].cell(0).acsr.tipo = 10;
                                    if (malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto == 0)
                                        malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto = 2;
                                    malha[aux].cell(0).acsr.injm.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro = noConv.flu;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                    malha[aux].cell(0).acsr.injm.fluidocol.TR = noConv.TRmist;

                                    noConv.qgstdTot *= 86400;

                                    double massT = malha[aux].cell(0).acsr.injm.MassP + malha[aux].cell(0).acsr.injm.MassC +
                                                   malha[aux].cell(0).acsr.injm.MassG;
                                    norma += pow((massT - (noConv.mliqmist + noConv.mgasmist + noConv.mliqmistNeg + noConv.mgasmistNeg)) /
                                                     (noConv.mliqmist + noConv.mgasmist + noConv.mliqmistNeg + noConv.mgasmistNeg),
                                                 2.);
                                    malha[aux].cell(0).acsr.injm.MassP = noConv.mliqmist - noConv.mliqCmist +
                                                                           noConv.mliqmistNeg - noConv.mliqCmistNeg;
                                    malha[aux].cell(0).acsr.injm.MassC = noConv.mliqCmist + noConv.mliqCmistNeg;
                                    malha[aux].cell(0).acsr.injm.MassG = noConv.mgasmist + noConv.mgasmistNeg;
                                }

                                int trocou = trocaFonteColetor(malha, aux);
                                (*arqRede.vg1dSP).qualTramo = aux;
                                double pini = malha[aux].cell(0).pres;
                                if (arqRede.malha[aux].perm == 1) {
                                    malha[aux].modoPerm = 1;
                                    int reverso = 0;
                                    if (malha[aux].cell(0).acsr.tipo == 1) {
                                        if (malha[aux].cell(0).acsr.injg.QGas < 0)
                                            reverso = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 2) {
                                        if (malha[aux].cell(0).acsr.injl.QLiq < 0)
                                            reverso = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 10) {
                                        if ((malha[aux].cell(0).acsr.injm.MassC +
                                             malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassP) < 0)
                                            reverso = 1;
                                    }
                                    if (reverso == 0) {
                                        if (malha[aux].config().chokep.abertura[0] > 0.6) {
                                            if ((*arqRede.vg1dSP).iterRede == 0)
                                                valor = malha[aux].buscaProdPfundoPerm();
                                            else
                                                valor = malha[aux].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        } else {
                                            if ((*arqRede.vg1dSP).iterRede == 0)
                                                valor = malha[aux].buscaProdPfundoPerm2();
                                            else
                                                valor = malha[aux].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        }
                                    } else {
                                        if ((*arqRede.vg1dSP).iterRede == 0)
                                            valor = malha[aux].buscaProdPfundoPermRev();
                                        else
                                            valor = malha[aux].buscaProdPfundoPermRev(malha[aux].cell(0).pres);
                                    }
                                } else
                                    valor = 0;
                                if (trocou == -1)
                                    retornaFonteColetor(malha, aux);
                                if (valor < -1e9 || valor > 1e9) {
                                    inativo[aux] = 1;
                                    inativoColetor(aux, arqRede, inativo, narq);
                                    inativoAfluente(aux, arqRede, inativo);
                                    if (valor < -1e9) {
                                        aviso(aux);
                                        (*arqRede.vg1dSP).iterRede = 200;
                                    } else {
                                        aviso2(aux);
                                        (*arqRede.vg1dSP).iterRede = 200;
                                    }
                                    (*arqRede.vg1dSP).restartRede = 1;
                                    return 9000.;
                                }
                                Resolv[aux] = 1;
                                resolvGlob += Resolv[aux];
                                ciclomalha++;
                                if (inativo[aux] == 0) {
                                    int qNo = -1;
                                    int kno = 0;
                                    while (qNo < 0) {
                                        noRede tempEvol;
                                        for (int kcol = 0; kcol < normaEvol[kno].ncole; kcol++) {
                                            if (normaEvol[kno].cole[kcol] == aux)
                                                qNo = kno;
                                        }
                                        kno++;
                                    }
                                    if (IndNorma[aux] == 0) {
                                        normaEvol[qNo].normaP0 = normaEvol[qNo].normaP1;
                                        normaEvol[qNo].normaP1 = pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                        norma += pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                        IndNorma[aux] = 1;
                                        nnos++;
                                    }
                                    for (int iaflu = 0; iaflu < arqRede.malha[aux].nafluente; iaflu++) {
                                        int indaflu = arqRede.malha[aux].afluente[iaflu];
                                        if ((*arqRede.vg1dSP).fluidoRede == 1)
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                        else if ((*arqRede.vg1dSP).fluidoRede == 0)
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                        else
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                        if (arqRede.malha[indaflu].perm == 1) {
                                            pini = malha[indaflu].pGSup;
                                            if (IndNorma[indaflu] == 0) {
                                                IndNorma[indaflu] = 1;
                                                if (arqRede.malha[indaflu].presimposta == 0) {
                                                    norma += pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                                    nnos++;
                                                }
                                            }
                                            if (arqRede.malha[indaflu].presimposta == 0)
                                                malha[indaflu].pGSup = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres + (1. - (*arqRede.vg1dSP).relax) * malha[indaflu].pGSup;
                                        } else {
                                            malha[indaflu].pGSup = malha[aux].cell(0).pres;
                                            malha[indaflu].tGSup = malha[aux].cell(0).temp;
                                            malha[indaflu].tempSup = malha[aux].cell(0).temp;
                                        }
                                    }
                                }
                                int iaflu = 0;
                                int indaflu = arqRede.malha[aux].afluente[iaflu];
                                while (arqRede.malha[indaflu].perm == 0) {
                                    iaflu++;
                                    indaflu = arqRede.malha[aux].afluente[iaflu];
                                }
                                for (int icol2 = 0; icol2 < nderiva - 1; icol2++) {
                                    norma += pow((malha[ordCol[icol2]].cell(0).pres - malha[aux].cell(0).pres) / malha[aux].cell(0).pres, 2.);
                                    nnos++;
                                    if ((*arqRede.vg1dSP).iterRede <= 1 && (*arqRede.vg1dSP).relax > 0.5) {
                                        malha[ordCol[icol2]].cell(0).pres = (0.5) * malha[aux].cell(0).pres + (1. - 0.5) * malha[ordCol[icol2]].cell(0).pres;
                                    } else {
                                        malha[ordCol[icol2]].cell(0).pres = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres +
                                                                              (1. - (*arqRede.vg1dSP).relax) * malha[ordCol[icol2]].cell(0).pres;
                                    }
                                    if ((*arqRede.vg1dSP).iterRede <= 1 && (*arqRede.vg1dSP).relax > 0.5) {
                                        malha[ordCol[icol2]].presE = (0.5) * malha[aux].cell(0).pres + (1. - 0.5) * malha[ordCol[icol2]].cell(0).pres;
                                    } else {
                                        malha[ordCol[icol2]].presE = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres +
                                                                     (1. - (*arqRede.vg1dSP).relax) * malha[ordCol[icol2]].cell(0).pres;
                                    }
                                    malha[ordCol[icol2]].tempE = noConv.tempmist;
                                    malha[ordCol[icol2]].betaE = noConv.betmist;
                                    malha[ordCol[icol2]].titE = noConv.mgasmist / (noConv.mgasmist + noConv.mliqmist);
                                }
                                cout << "convergencia tramo: " << aux << endl;
                            } else {
                                double chutemass = 0.;
                                double area = 0.;
                                int aux = ordCol[icol];
                                if ((*arqRede.vg1dSP).iterRede == 0) {
                                    int indAflu = arqRede.malha[aux].afluente[0];
                                    malha[aux].cell(0).pres = malha[indAflu].pGSup;
                                }
                                if (arqRede.malha[aux].perm == 1) {
                                    (*arqRede.vg1dSP).qualTramo = aux;
                                    if (malha[aux].config().ConContEntrada != 2) {
                                        if ((*arqRede.vg1dSP).fluidoRede != 0) {
                                            malha[aux].cell(0).acsr.injl.FluidoPro = noConv.flu;
                                            malha[aux].cell(0).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                            malha[aux].cell(0).acsr.tipo = 2;
                                            if (malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto == 0)
                                                malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto = 2;
                                            malha[aux].cell(0).acsr.injl.bet = noConv.betmist;
                                            malha[aux].cell(0).acsr.injl.temp = noConv.tempmist;
                                            malha[aux].cell(0).acsr.injl.fluidocol.TR = noConv.TRmist;

                                            malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                            if ((*arqRede.vg1dSP).iterRede > 0 && malha[aux].cell(0).acsr.injl.QLiq < 0) {
                                                malha[aux].cell(0).acsr.injl.FluidoPro = malha[aux].fluiRevRede;
                                            }

                                        } else {
                                            malha[aux].cell(0).acsr.injg.FluidoPro = noConv.flu;
                                            malha[aux].cell(0).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                            malha[aux].cell(0).acsr.tipo = 1;
                                            if (malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto == 0)
                                                malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto = 2;
                                            malha[aux].cell(0).acsr.injg.seco = 0;
                                            malha[aux].cell(0).acsr.injg.temp = noConv.tempmist;
                                            malha[aux].cell(0).acsr.injg.fluidocol.TR = noConv.TRmist;

                                            malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                            if ((*arqRede.vg1dSP).iterRede > 0 && malha[aux].cell(0).acsr.injg.QGas < 0) {
                                                malha[aux].cell(0).acsr.injg.FluidoPro = malha[aux].fluiRevRede;
                                            }
                                        }
                                    } else {
                                        malha[aux].cell(0).acsr.injm.FluidoPro = noConv.flu;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                        malha[aux].cell(0).acsr.tipo = 10;
                                        if (malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto == 0)
                                            malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto = 2;
                                        malha[aux].cell(0).acsr.injm.temp = noConv.tempmist;
                                        malha[aux].cell(0).acsr.injm.condTermo = 0;
                                        malha[aux].cell(0).acsr.injm.fluidocol.TR = noConv.TRmist;
                                        malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                    }

                                    if (malha[aux].config().ConContEntrada != 2) {
                                        for (int icol2 = icol; icol2 < nderiva; icol2++)
                                            area += malha[ordCol[icol2]].cell(0).duto.area;
                                        if ((*arqRede.vg1dSP).iterRede == 0) {
                                            if ((*arqRede.vg1dSP).fluidoRede != 0)
                                                chutemass = (noConv.qlmistStd + noConv.qlmistStdNeg) * malha[aux].cell(0).duto.area / area;
                                            else {
                                                chutemass = (noConv.qgstdTot * 86400 + noConv.qgstdTotNeg * 86400) *
                                                            malha[aux].cell(0).duto.area / area;
                                                if (fabs(chutemass) > 1e-15)
                                                    malha[aux].cell(0).acsr.injg.razCompGas = (noConv.qcmistStd + noConv.qcmistStdNeg) /
                                                                                                (noConv.qgstdTot * 86400. + noConv.qgstdTotNeg * 86400);
                                                else
                                                    malha[aux].cell(0).acsr.injg.razCompGas = 0.;
                                            }
                                        } else {
                                            if ((*arqRede.vg1dSP).fluidoRede == 1)
                                                chutemass = malha[aux].cell(0).acsr.injl.QLiq;
                                            else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                                chutemass = malha[aux].cell(0).acsr.injg.QGas;
                                                if (fabs(noConv.qgstdTot * 86400. + noConv.qgstdTotNeg * 86400) > 1e-15)
                                                    malha[aux].cell(0).acsr.injg.razCompGas = (noConv.qcmistStd + noConv.qcmistStdNeg) /
                                                                                                (noConv.qgstdTot * 86400. + noConv.qgstdTotNeg * 86400);
                                                else
                                                    malha[aux].cell(0).acsr.injg.razCompGas = 0.;
                                            }
                                        }
                                        if (arqRede.malha[aux].perm == 1) {
                                            malha[aux].modoPerm = 1;
                                            if (chutemass > 0.) {
                                                if (malha[aux].config().chokep.abertura[0] > 0.6)
                                                    valor = malha[aux].buscaProdPresPresPerm(chutemass, noConv.qgstdTot * 86400);
                                                else
                                                    valor = malha[aux].buscaProdPresPresPerm2(chutemass, noConv.qgstdTot * 86400);
                                            } else {
                                                malha[aux].cell(0).acsr.injg.razCompGas = malha[aux].config().razCompGasReves;
                                                valor = malha[aux].buscaProdPresPresPermRev(chutemass, -noConv.qgstdTot * 86400);
                                            }
                                            if ((*arqRede.vg1dSP).fluidoRede != 0) {
                                                if (malha[aux].cell(0).acsr.injl.QLiq > (noConv.qlmistStd + noConv.qlmistStdNeg)) {
                                                    malha[aux].cell(0).acsr.injl.QLiq = chutemass;
                                                    cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo: " << aux << endl;
                                                    arqRede.malha[aux].principal = 1;
                                                } else
                                                    arqRede.malha[aux].principal = 0;
                                            } else {
                                                if (malha[aux].cell(0).acsr.injg.QGas > (noConv.qgstdTot * 86400 + noConv.qgstdTotNeg * 86400.)) {
                                                    malha[aux].cell(0).acsr.injg.QGas = chutemass;
                                                    cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo: " << aux << endl;
                                                    arqRede.malha[aux].principal = 1;
                                                } else
                                                    arqRede.malha[aux].principal = 0;
                                            }
                                        } else
                                            valor = 0;

                                        if (valor < -1e9 || valor > 1e9) {
                                            inativo[aux] = 1;
                                            inativoColetor(aux, arqRede, inativo, narq);
                                            inativoAfluente(aux, arqRede, inativo);
                                            if (valor < -1e9) {
                                                aviso(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            } else {
                                                aviso2(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            }
                                            (*arqRede.vg1dSP).restartRede = 1;
                                            return 9000.;
                                        }
                                        Resolv[aux] = 1;
                                        resolvGlob += Resolv[aux];
                                        ciclomalha++;
                                        if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                            if (arqRede.malha[aux].perm == 1 && malha[aux].cell(0).acsr.injl.QLiq > 0.)
                                                noConv.qlmistStd -= malha[aux].cell(0).acsr.injl.QLiq;
                                            if (fabs(malha[aux].cell(0).acsr.injl.QLiq) < 1e-15)
                                                cout << "#################Sem vazao no tramo : " << aux << endl;
                                        } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                            if (arqRede.malha[aux].perm == 1 && malha[aux].cell(0).acsr.injg.QGas > 0.)
                                                noConv.qgstdTot -= malha[aux].cell(0).acsr.injg.QGas / 86400.;
                                            if (fabs(malha[aux].cell(0).acsr.injg.QGas) < 1e-15)
                                                cout << "#################Sem vazao no tramo : " << aux << endl;
                                        } else if ((*arqRede.vg1dSP).fluidoRede == 2) {
                                            if (arqRede.malha[aux].perm == 1 &&
                                                (malha[aux].cell(0).acsr.injl.QLiq) > 0) {
                                                double rlcA = malha[aux].cell(0).acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
                                                double vazmassc = rlcA * malha[aux].cell(0).acsr.injl.QLiq * malha[aux].cell(0).acsr.injl.bet / 86400;
                                                double massic = malha[aux].cell(0).acsr.injl.QLiq *
                                                                (1. - malha[aux].cell(0).acsr.injl.bet) / 86400;
                                                double Rhogs = malha[aux].cell(0).acsr.injl.FluidoPro.Deng * 1.225; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
                                                double Rhols = (1000 * 141.5 / (131.5 + malha[aux].cell(0).acsr.injl.FluidoPro.API)) * (1 - malha[aux].cell(0).acsr.injl.FluidoPro.BSW) + 1000. * malha[aux].cell(0).acsr.injl.FluidoPro.Denag * malha[aux].cell(0).acsr.injl.FluidoPro.BSW;
                                                double multiplicador = (Rhols + malha[aux].cell(0).acsr.injl.FluidoPro.RGO * Rhogs *
                                                                                    (1 - malha[aux].cell(0).acsr.injl.FluidoPro.BSW));
                                                massic *= multiplicador;
                                                double fracmasshidra = malha[aux].cell(0).acsr.injl.FluidoPro.FracMassHidra(malha[aux].cell(0).pres, malha[aux].cell(0).temp);
                                                noConv.mliqmist -= ((1. - fracmasshidra) * massic + vazmassc);
                                                noConv.mgasmist -= fracmasshidra * massic;
                                            }
                                            if (fabs(malha[aux].cell(0).acsr.injl.QLiq) < 1e-15)
                                                cout << "#################Sem vazao no tramo : " << aux << endl;
                                        }
                                    } else {
                                        if (arqRede.malha[aux].perm == 1) {
                                            int iaflu = 0;
                                            int indaflu = arqRede.malha[aux].afluente[iaflu];
                                            while (arqRede.malha[indaflu].perm == 0) {
                                                iaflu++;
                                                indaflu = arqRede.malha[aux].afluente[iaflu];
                                            }
                                            valor = malha[aux].buscaProdPfundoPerm3(malha[indaflu].pGSup);
                                        } else
                                            valor = 0.;
                                        if (valor < -1e9 || valor > 1e9) {
                                            inativo[aux] = 1;
                                            inativoColetor(aux, arqRede, inativo, narq);
                                            inativoAfluente(aux, arqRede, inativo);
                                            if (valor < -1e9) {
                                                aviso(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            } else {
                                                aviso2(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            }
                                            (*arqRede.vg1dSP).restartRede = 1;
                                            return 9000.;
                                        }
                                        Resolv[aux] = 1;
                                        resolvGlob += Resolv[aux];
                                        ciclomalha++;
                                        double rlcA = malha[aux].cell(0).acsr.injm.fluidocol.MasEspFlu(1.001, 15.);
                                        double Rhogs = malha[aux].cell(0).acsr.injm.FluidoPro.Deng * 1.225; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
                                        double Rhols = (1000 * 141.5 / (131.5 + malha[aux].cell(0).acsr.injm.FluidoPro.API)) * (1 - malha[aux].cell(0).acsr.injm.FluidoPro.BSW) + 1000. * malha[aux].cell(0).acsr.injm.FluidoPro.Denag *
                                                                                                                                                                                          malha[aux].cell(0).acsr.injm.FluidoPro.BSW;
                                        double multiplicador = (Rhols + malha[aux].cell(0).acsr.injm.FluidoPro.RGO * Rhogs *
                                                                            (1 - malha[aux].cell(0).acsr.injm.FluidoPro.BSW));
                                        double fracmasshidra = malha[aux].cell(0).acsr.injm.FluidoPro.FracMassHidra(malha[aux].cell(0).pres,
                                                                                                                      malha[aux].cell(0).temp);
                                        if ((*arqRede.vg1dSP).fluidoRede == 1 && malha[aux].cell(0).acsr.injm.MassP > 0.) {
                                            if (arqRede.malha[aux].perm == 1) {
                                                noConv.qlmistStd -= (malha[aux].cell(0).acsr.injm.MassP * 86400. /
                                                                         ((1. - fracmasshidra) * multiplicador) +
                                                                     malha[aux].cell(0).acsr.injm.MassC * 86400. / rlcA);
                                            }
                                        } else if ((*arqRede.vg1dSP).fluidoRede == 0 && malha[aux].cell(0).acsr.injm.MassG > 0.) {
                                            double tit = malha[aux].cell(0).acsr.injm.FluidoPro.dStockTankVaporMassFraction;
                                            noConv.qgstdTot -= malha[aux].cell(0).acsr.injm.MassG * tit /
                                                               (malha[aux].cell(0).acsr.injm.FluidoPro.Deng * 1.225 * fracmasshidra);
                                        } else if ((malha[aux].cell(0).acsr.injm.MassC + malha[aux].cell(0).acsr.injm.MassP +
                                                    malha[aux].cell(0).acsr.injm.MassG) > 0.) {
                                            noConv.mliqmist -= (malha[aux].cell(0).acsr.injm.MassC + malha[aux].cell(0).acsr.injm.MassP);
                                            noConv.mgasmist -= malha[aux].cell(0).acsr.injm.MassG;
                                        }
                                        if (fabs(malha[aux].cell(0).acsr.injm.MassC +
                                                 malha[aux].cell(0).acsr.injm.MassP +
                                                 malha[aux].cell(0).acsr.injm.MassG) < 1e-15)
                                            cout << "#################Sem vazao no tramo : " << aux << endl;
                                    }
                                    cout << "convergencia tramo: " << aux << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (fabs(testaNMalha) > 0) {
        (*arqRede.vg1dSP).erroRede = fabs(sqrt(norma) / nnos - (*arqRede.vg1dSP).norma0) / (*arqRede.vg1dSP).norma0;
        (*arqRede.vg1dSP).norma0 = sqrt(norma) / nnos;
        return sqrt(norma) / nnos;
    } else
        return 0;
}

void atualizaComp(SProd &sistem1, int i) {
    double titF = 0.;
    ProFlu fluF;

    double boF = 1.;
    double baF = 1.;
    double fwF = 1.;
    double rhoOF = 900.;
    double rhoWF = 1000.;

    if (sistem1.cell(i - 1).acsr.tipo == 1) {
        sistem1.cell(i - 1).acsr.injg.FluidoPro.tabelaDinamica = 0;
        if (i > 0)
            sistem1.cell(i - 1).acsr.injg.FluidoPro.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition, sistem1.cell(i - 1).flui.oCalculatedVapComposition,
                sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.injg.FluidoPro.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.injg.FluidoPro;
        fwF = 0.;
        titF = 1.;
    } else if (sistem1.cell(i - 1).acsr.tipo == 2) {
        sistem1.cell(i - 1).acsr.injl.FluidoPro.tabelaDinamica = 0;
        if (i > 0)
            sistem1.cell(i - 1).acsr.injl.FluidoPro.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.injl.FluidoPro.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.injl.FluidoPro;
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (sistem1.cell(i - 1).acsr.tipo == 3) {
        sistem1.cell(i - 1).acsr.ipr.FluidoPro.tabelaDinamica = 0;
        if (i > 0)
            sistem1.cell(i - 1).acsr.ipr.FluidoPro.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.ipr.FluidoPro.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.ipr.FluidoPro;
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (sistem1.cell(i - 1).acsr.tipo == 15) {
        sistem1.cell(i - 1).acsr.radialPoro.flup.tabelaDinamica = 0;
        int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
        double tRes = sistem1.cell(i).acsr.radialPoro.tRes;
        for (int k = 0; k < ncelRad; k++) {
            sistem1.cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
        }
        if (i > 0)
            sistem1.cell(i - 1).acsr.radialPoro.flup.atualizaPropComp(
                sistem1.cell(i - 1).pres, tRes, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.radialPoro.flup.atualizaPropComp(
                sistem1.cell(i).pres, tRes, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.radialPoro.flup;
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (sistem1.cell(i - 1).acsr.tipo == 16) {
        sistem1.cell(i - 1).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
        int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
        double tRes = sistem1.cell(i).acsr.poroso2D.dados.transfer.tRes;
        for (int k = 0; k < ncelRad; k++) {
            sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
        }
        int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
        tRes = sistem1.cell(i).acsr.poroso2D.dados.tRes;
        for (int k = 0; k < ncelEle; k++) {
            sistem1.cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
        }
        if (i > 0)
            sistem1.cell(i - 1).acsr.poroso2D.dados.flup.atualizaPropComp(
                sistem1.cell(i - 1).pres, tRes, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.poroso2D.dados.flup.atualizaPropComp(
                sistem1.cell(i).pres, tRes, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.poroso2D.dados.flup;
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (sistem1.cell(i - 1).acsr.tipo == 10) {
        sistem1.cell(i - 1).acsr.injm.FluidoPro.tabelaDinamica = 0;
        if (i > 0)
            sistem1.cell(i - 1).acsr.injm.FluidoPro.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.injm.FluidoPro.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        fluF = sistem1.cell(i - 1).acsr.injm.FluidoPro;
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (sistem1.cell(i - 1).acsr.tipo == 9 && sistem1.cell(i - 1).acsr.fontechk.abertura > 1e-6) {
        sistem1.cell(i - 1).acsr.fontechk.fluidoP.tabelaDinamica = 0;
        sistem1.cell(i - 1).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
        if (i > 0)
            sistem1.cell(i - 1).acsr.fontechk.fluidoP.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp,
                sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.fontechk.fluidoP.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        if (i > 0)
            sistem1.cell(i - 1).acsr.fontechk.fluidoPamb.atualizaPropComp(
                sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, sistem1.cell(i - 1).flui.dCalculatedBeta,
                sistem1.cell(i - 1).flui.oCalculatedLiqComposition,
                sistem1.cell(i - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        else
            sistem1.cell(i - 1).acsr.fontechk.fluidoPamb.atualizaPropComp(
                sistem1.cell(i).pres, sistem1.cell(i).temp, -1, NULL, NULL, sistem1.config().pocinjec);
        if (sistem1.cell(i - 1).acsr.fontechk.presT > sistem1.cell(i - 1).acsr.fontechk.pamb) {
            fluF = sistem1.cell(i - 1).acsr.fontechk.fluidoP;
        } else {
            fluF = sistem1.cell(i - 1).acsr.fontechk.fluidoPamb;
        }
        boF = fluF.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        baF = fluF.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rhoWF = fluF.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    }

    // sistem1.renovaFonte(i - 1);//metodo que verifica se existe uma fonte na celula e calcula o valor das
    // vazoes massicas de liquido produzido (oleo+agua), gas e liquido complementar
    // relacao entre fonte a esquerda e a direita de uma celula:
    sistem1.cell(i).fontemassLL = sistem1.cell(i - 1).fontemassLR;
    sistem1.cell(i).fontemassCL = sistem1.cell(i - 1).fontemassCR;
    sistem1.cell(i).fontemassGL = sistem1.cell(i - 1).fontemassGR;

    double vazMasLiqL = sistem1.cell(i).MliqiniL;
    double vazMasGasL = sistem1.cell(i).ML - sistem1.cell(i).MliqiniL;
    double fonteMasLiqL = sistem1.cell(i).fontemassLL;
    double fonteMasGasL = sistem1.cell(i).fontemassGL;

    // definicao de temperaturas maximas e minimas para uma eventual reavaliacao do metodo
    // ASTM quando ocorre mistura de fluidos e se deseja atualizar o modelo
    // de viscosidade de oleo morto se este for o ASTM que trabalha com um par de temperatura

    double bo;
    double ba;
    double rs;
    double titV = 0.;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (sistem1.cell(i - 1).flui.RGO < 1e6) {
        rs = sistem1.cell(i - 1).flui.RS(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        bo = sistem1.cell(i - 1).flui.BOFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp, rs);
        ba = sistem1.cell(i - 1).flui.BAFunc(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
        rs = rs * 6.29 / 35.31467;
    } else {
        bo = 1;
        rs = 0;
        ba = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    sistem1.cell(i - 1).FW = sistem1.cell(i - 1).flui.BSW * ba / (bo + ba * sistem1.cell(i - 1).flui.BSW - sistem1.cell(i - 1).flui.BSW * bo);
    sistem1.cell(i - 1).FWini = sistem1.cell(i - 1).FW;
    // temperatura na fronteira entre a celula i-1 e a celula i, da segunda iteracao em diante, pode-se
    // usar a temperatura da celula i, pois já a tem calculada, mas isto pode ser um complicador de
    // convergencia, mais seguro manter o criterio em todas as iteracoes, neste caso, a temperatura na
    // fronteira e admitida = a temperatura da celula i-1
    //	  tmed = (sistem1.cell(i).dx * sistem1.cell(i).temp + sistem1.cell(i).dxL * sistem1.cell(i).tempL)

    double fwV = sistem1.cell(i - 1).FW;
    double rhoOV = sistem1.cell(i - 1).flui.MasEspoleo(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
    double rhoWV = sistem1.cell(i - 1).flui.MasEspAgua(sistem1.cell(i - 1).pres, sistem1.cell(i - 1).temp);
    titV = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
    vazMasLiqL *= titV;
    fonteMasLiqL *= titF;
    double pesoMolV = 0;
    double pesoMolF = 0;
    for (int j = 0; j < sistem1.cell(i).flui.npseudo; j++) {
        pesoMolV += sistem1.cell(i).flui.masMol[j] * sistem1.cell(i).flui.fracMol[j];
        pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
    }
    double vazMolV = (vazMasLiqL + vazMasGasL) / pesoMolV;
    double vazMolF = (fonteMasLiqL + fonteMasGasL) / pesoMolF;
    double razMolV = vazMolV / (vazMolV + vazMolF);
    if (vazMolF > 0.) {
        for (int j = 0; j < sistem1.cell(i).flui.npseudo; j++) {
            sistem1.cell(i).flui.fracMol[j] = razMolV * sistem1.cell(i - 1).flui.fracMol[j] +
                                                (1. - razMolV) * fluF.fracMol[j];
        }
    } else {
        for (int j = 0; j < sistem1.cell(i).flui.npseudo; j++)
            sistem1.cell(i).flui.fracMol[j] = sistem1.cell(i - 1).flui.fracMol[j];
    }
}

void atualizaCel(SProd &sistem1, int i) {
    for (int j = 0; j < sistem1.cell(i).flui.npseudo; j++)
        sistem1.cell(i).flui.fracMol[j] = sistem1.cell(i - 1).flui.fracMol[j];
    sistem1.cell(i).flui.iCalculatedStockTankThermodynamicCondition =
        sistem1.cell(i - 1).flui.iCalculatedStockTankThermodynamicCondition;
    sistem1.cell(i).flui.dStockTankVaporMassFraction =
        sistem1.cell(i - 1).flui.dStockTankVaporMassFraction;
    sistem1.cell(i).flui.dStockTankLiquidDensity =
        sistem1.cell(i - 1).flui.dStockTankLiquidDensity;
    sistem1.cell(i).flui.dStockTankVaporDensity =
        sistem1.cell(i - 1).flui.dStockTankVaporDensity;

    if (sistem1.cell(i).flui.dStockTankLiquidDensity > 0.01) {
        sistem1.cell(i).flui.API = 141.5 / (sistem1.cell(i).flui.dStockTankLiquidDensity / 1000.) - 131.5;
    } else
        sistem1.cell(i).flui.API = 50;
    sistem1.cell(i).flui.Deng = sistem1.cell(i).flui.dStockTankVaporDensity / 1.225;
    sistem1.cell(i).flui.RGO = sistem1.cell(i - 1).flui.RGO;
    sistem1.cell(i).flui.IRGO = sistem1.cell(i - 1).flui.IRGO;

    double tmed = (sistem1.cell(i).dx * sistem1.cell(i).temp + sistem1.cell(i).dxL * sistem1.cell(i).tempL) / (sistem1.cell(i).dx + sistem1.cell(i).dxL);
    double pmed = sistem1.cell(i).presaux + sistem1.cell(i - 1).dpB / 98066.5;

    int veriI = i - 1;
    while ((sistem1.cell(veriI).flui.dCalculatedBeta > 1 - (0.0 + 1e-15) ||
            sistem1.cell(veriI).flui.dCalculatedBeta < (0.0 + 1e-15)) &&
           (veriI > 0))
        veriI--;
    if ((veriI == 0) && (sistem1.cell(veriI).flui.dCalculatedBeta > 1 - (0.0 + 1e-15) || sistem1.cell(veriI).flui.dCalculatedBeta < (0.0 + 1e-15)))
        veriI = i - 1;
    sistem1.cell(i).flui.atualizaPropComp(pmed, tmed, sistem1.cell(veriI).flui.dCalculatedBeta,
                                            sistem1.cell(veriI).flui.oCalculatedLiqComposition,
                                            sistem1.cell(veriI).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
}

int ranqueiaCol(Rede &arqRede, int i) {
    int rank = 0;
    if (arqRede.malha[i].ncoleta > 0) {
        for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
            int icol = arqRede.malha[i].coleta[j];
            if (arqRede.malha[icol].ncoleta == 0)
                rank = 1;
            else if (arqRede.malha[icol].ncoleta > 0)
                rank += (ranqueiaCol(arqRede, icol) + 1);
        }
    }
    return rank;
}

void preparaTabDin(SProd &sistem1) {
    for (int i = 0; i <= sistem1.ncel; i++) {
        sistem1.cell(i).flui.tabelaDinamica = 0;
        if (sistem1.cell(i).acsr.tipo == 1) {
            if (sistem1.cell(i).acsr.tipo == 1 && sistem1.cell(i).acsr.injg.QGas < 0) {
                sistem1.cell(i).acsr.injg.FluidoPro = sistem1.cell(0).flui;
            }
            sistem1.cell(i).acsr.injg.FluidoPro.tabelaDinamica = 0;
            sistem1.cell(i).acsr.injg.FluidoPro.atualizaPropCompStandard();
        } else if (sistem1.cell(i).acsr.tipo == 2) {
            if (sistem1.cell(i).acsr.tipo == 2 && sistem1.cell(i).acsr.injl.QLiq < 0) {
                sistem1.cell(i).acsr.injl.FluidoPro = sistem1.cell(0).flui;
            }
            sistem1.cell(i).acsr.injl.FluidoPro.tabelaDinamica = 0;
            sistem1.cell(i).acsr.injl.FluidoPro.atualizaPropCompStandard();
        } else if (sistem1.cell(i).acsr.tipo == 3) {
            sistem1.cell(i).acsr.ipr.FluidoPro.tabelaDinamica = 0;
            sistem1.cell(i).acsr.ipr.FluidoPro.atualizaPropCompStandard();
        } else if (sistem1.cell(i).acsr.tipo == 15) {
            sistem1.cell(i).acsr.radialPoro.flup.tabelaDinamica = 0;
            int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
            for (int k = 0; k < ncelRad; k++) {
                sistem1.cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
            }
        } else if (sistem1.cell(i).acsr.tipo == 16) {
            sistem1.cell(i).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
            int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
            for (int k = 0; k < ncelRad; k++) {
                sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
            }
            int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
            for (int k = 0; k < ncelEle; k++) {
                sistem1.cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
            }
        } else if (sistem1.cell(i).acsr.tipo == 10) {
            if (sistem1.cell(i).acsr.tipo == 10 &&
                sistem1.cell(i).acsr.injm.MassC + sistem1.cell(i).acsr.injm.MassG + sistem1.cell(i).acsr.injm.MassP < 0) {
                sistem1.cell(i).acsr.injm.FluidoPro = sistem1.cell(0).flui;
            }
            sistem1.cell(i).acsr.injm.FluidoPro.atualizaPropCompStandard();
        } else if (sistem1.cell(i).acsr.tipo == 9 && sistem1.cell(i).acsr.fontechk.abertura > 1e-6) {
            sistem1.cell(i).acsr.fontechk.fluidoP.tabelaDinamica = 0;
            sistem1.cell(i).acsr.fontechk.fluidoP.atualizaPropCompStandard();
            sistem1.cell(i).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
            sistem1.cell(i).acsr.fontechk.fluidoPamb.atualizaPropCompStandard();
        }
    }
    for (int i = 0; i < sistem1.ntabDin; i++) {
        int i0 = sistem1.tabDin[i].celIni;
        int i1 = sistem1.tabDin[i].celFim;
        sistem1.tabDin[i].pmax = sistem1.cell(i0).pres;
        sistem1.tabDin[i].pmin = sistem1.cell(i0).pres;
        sistem1.tabDin[i].tmax = sistem1.cell(i0).temp;
        sistem1.tabDin[i].tmin = sistem1.cell(i0).temp;
        for (int j = i0 + 1; j <= i1; j++) {
            if (sistem1.cell(j).pres < sistem1.tabDin[i].pmin)
                sistem1.tabDin[i].pmin = sistem1.cell(j).pres;
            else if (sistem1.cell(j).pres > sistem1.tabDin[i].pmax)
                sistem1.tabDin[i].pmax = sistem1.cell(j).pres;
            if (sistem1.cell(j).temp < sistem1.tabDin[i].tmin)
                sistem1.tabDin[i].tmin = sistem1.cell(j).temp;
            else if (sistem1.cell(j).temp > sistem1.tabDin[i].tmax)
                sistem1.tabDin[i].tmax = sistem1.cell(j).temp;
        }
        sistem1.tabDin[i].pmax *= 1.6;
        sistem1.tabDin[i].pmin *= 0.6;
        sistem1.tabDin[i].tmax = sistem1.tabDin[i].tmax + 20;
        sistem1.tabDin[i].tmin = sistem1.tabDin[i].tmin - 20;
        if (sistem1.tabDin[i].tmin < -30)
            sistem1.tabDin[i].tmin = -30;
        sistem1.tabDin[i].delP = 5.;
        sistem1.tabDin[i].delT = 5.;
        sistem1.tabDin[i].npontosT = floor((sistem1.tabDin[i].tmax - sistem1.tabDin[i].tmin) / (sistem1.tabDin[i].delT)) + 1;
        if (sistem1.tabDin[i].npontosT <= 3)
            sistem1.tabDin[i].npontosT = 4;
        sistem1.tabDin[i].npontosP = floor((sistem1.tabDin[i].pmax - sistem1.tabDin[i].pmin) / (sistem1.tabDin[i].delP)) + 1;
        if (sistem1.tabDin[i].npontosP <= 3)
            sistem1.tabDin[i].npontosP = 4;
        sistem1.tabDin[i].delP = (sistem1.tabDin[i].pmax - sistem1.tabDin[i].pmin) /
                                 (sistem1.tabDin[i].npontosP - 1);
        sistem1.tabDin[i].delT = (sistem1.tabDin[i].tmax - sistem1.tabDin[i].tmin) /
                                 (sistem1.tabDin[i].npontosT - 1);
    }
    int i0 = sistem1.tabDin[0].celIni;
    int i1 = sistem1.tabDin[0].celFim;
    sistem1.cell(i0).flui.atualizaPropCompStandard();
    double pmed = sistem1.cell(i0).presaux;
    double tmed = sistem1.cell(i0).temp;
    sistem1.cell(i0).flui.atualizaPropComp(pmed, tmed, -1, NULL, NULL, sistem1.config().pocinjec);
    for (int i = i0 + 1; i <= i1; i++) {
        for (int k = 0; k < sistem1.cell(i).flui.npseudo; k++)
            sistem1.cell(i).flui.fracMol[k] = sistem1.cell(i - 1).flui.fracMol[k];
        atualizaCel(sistem1, i);
    }

    for (int j = 1; j < sistem1.ntabDin; j++) {
        int i0 = sistem1.tabDin[j].celIni;
        int i1 = sistem1.tabDin[j].celFim;
        atualizaComp(sistem1, i0);
        sistem1.cell(i0).flui.atualizaPropCompStandard();
        tmed = (sistem1.cell(i0).dx * sistem1.cell(i0).temp + sistem1.cell(i0).dxL * sistem1.cell(i0).tempL) / (sistem1.cell(i0).dx + sistem1.cell(i0).dxL);
        pmed = sistem1.cell(i0).presaux + sistem1.cell(i0 - 1).dpB / 98066.5;
        sistem1.cell(i0).flui.atualizaPropComp(pmed, tmed, sistem1.cell(i0 - 1).flui.dCalculatedBeta,
                                                 sistem1.cell(i0 - 1).flui.oCalculatedLiqComposition,
                                                 sistem1.cell(i0 - 1).flui.oCalculatedVapComposition, sistem1.config().pocinjec);
        for (int i = i0 + 1; i <= i1; i++) {
            for (int k = 0; k < sistem1.cell(i).flui.npseudo; k++)
                sistem1.cell(i).flui.fracMol[k] = sistem1.cell(i - 1).flui.fracMol[k];
            atualizaCel(sistem1, i);
        }
    }
    for (int i = 0; i < sistem1.ntabDin; i++) {
        sistem1.tabDin[i].rholF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].rhogF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].DrholDpF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].DrhogDpF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].DrholDtF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].DrhogDtF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].valBO = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].tit = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].rs = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].cplF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].cpgF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].valZ = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].HlF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].HgF = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].valdZdT = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].valdZdP = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].viscO = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].viscG = new double *[sistem1.tabDin[i].npontosP + 1];
        sistem1.tabDin[i].PBF = new double[sistem1.tabDin[i].npontosT];
        sistem1.tabDin[i].TBF = new double[sistem1.tabDin[i].npontosT];
        for (int k = 0; k < sistem1.tabDin[i].npontosP + 1; k++) {
            sistem1.tabDin[i].rholF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].rhogF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].DrholDpF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].DrhogDpF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].DrholDtF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].DrhogDtF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].valBO[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].tit[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].rs[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].cplF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].cpgF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].valZ[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].HlF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].HgF[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].valdZdT[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].valdZdP[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].viscO[k] = new double[sistem1.tabDin[i].npontosT + 1];
            sistem1.tabDin[i].viscG[k] = new double[sistem1.tabDin[i].npontosT + 1];
        }
        for (int k = 0; k < sistem1.tabDin[i].npontosP + 1; k++) {
            for (int j = 0; j < sistem1.tabDin[i].npontosT + 1; j++) {
                sistem1.tabDin[i].rholF[k][j] = 0;
                sistem1.tabDin[i].rhogF[k][j] = 0;
                sistem1.tabDin[i].DrholDpF[k][j] = 0;
                sistem1.tabDin[i].DrhogDpF[k][j] = 0;
                sistem1.tabDin[i].DrholDtF[k][j] = 0;
                sistem1.tabDin[i].DrhogDtF[k][j] = 0;
                sistem1.tabDin[i].valBO[k][j] = 0;
                sistem1.tabDin[i].tit[k][j] = 0;
                sistem1.tabDin[i].rs[k][j] = 0;
                sistem1.tabDin[i].cplF[k][j] = 0;
                sistem1.tabDin[i].cpgF[k][j] = 0;
                sistem1.tabDin[i].valZ[k][j] = 0;
                sistem1.tabDin[i].HlF[k][j] = 0;
                sistem1.tabDin[i].HgF[k][j] = 0;
                sistem1.tabDin[i].valdZdT[k][j] = 0;
                sistem1.tabDin[i].valdZdP[k][j] = 0;
                sistem1.tabDin[i].viscO[k][j] = 0;
                sistem1.tabDin[i].viscG[k][j] = 0;
            }
        }
        for (int k = 0; k < sistem1.tabDin[i].npontosT; k++) {
            sistem1.tabDin[i].PBF[k] = 0;
            sistem1.tabDin[i].TBF[k] = 0;
        }
        sistem1.tabDin[i].TBF[0] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].rholF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].rhogF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].DrholDpF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].DrhogDpF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].DrholDtF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].DrhogDtF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].valBO[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].tit[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].rs[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].cplF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].cpgF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].valZ[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].HlF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].HgF[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].valdZdT[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].valdZdP[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].viscO[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].viscG[0][1] = sistem1.tabDin[i].tmin;
        sistem1.tabDin[i].rholF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].rhogF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].DrholDpF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].DrhogDpF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].DrholDtF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].DrhogDtF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].valBO[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].tit[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].rs[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].cplF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].cpgF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].valZ[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].HlF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].HgF[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].valdZdT[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].valdZdP[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].viscO[1][0] = sistem1.tabDin[i].pmin;
        sistem1.tabDin[i].viscG[1][0] = sistem1.tabDin[i].pmin;
        for (int k = 1; k < sistem1.tabDin[i].npontosT; k++) {
            sistem1.tabDin[i].TBF[k] = sistem1.tabDin[i].TBF[k - 1] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].rholF[0][k + 1] = sistem1.tabDin[i].rholF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].rhogF[0][k + 1] = sistem1.tabDin[i].rhogF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].DrholDpF[0][k + 1] = sistem1.tabDin[i].DrholDpF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].DrhogDpF[0][k + 1] = sistem1.tabDin[i].DrhogDpF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].DrholDtF[0][k + 1] = sistem1.tabDin[i].DrholDtF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].DrhogDtF[0][k + 1] = sistem1.tabDin[i].DrhogDtF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].valBO[0][k + 1] = sistem1.tabDin[i].valBO[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].tit[0][k + 1] = sistem1.tabDin[i].tit[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].rs[0][k + 1] = sistem1.tabDin[i].rs[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].cplF[0][k + 1] = sistem1.tabDin[i].cplF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].cpgF[0][k + 1] = sistem1.tabDin[i].cpgF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].valZ[0][k + 1] = sistem1.tabDin[i].valZ[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].HlF[0][k + 1] = sistem1.tabDin[i].HlF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].HgF[0][k + 1] = sistem1.tabDin[i].HgF[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].valdZdT[0][k + 1] = sistem1.tabDin[i].valdZdT[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].valdZdP[0][k + 1] = sistem1.tabDin[i].valdZdP[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].viscO[0][k + 1] = sistem1.tabDin[i].viscO[0][k] + sistem1.tabDin[i].delT;
            sistem1.tabDin[i].viscG[0][k + 1] = sistem1.tabDin[i].viscG[0][k] + sistem1.tabDin[i].delT;
        }
        for (int k = 1; k < sistem1.tabDin[i].npontosP; k++) {
            sistem1.tabDin[i].rholF[k + 1][0] = sistem1.tabDin[i].rholF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].rhogF[k + 1][0] = sistem1.tabDin[i].rhogF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].DrholDpF[k + 1][0] = sistem1.tabDin[i].DrholDpF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].DrhogDpF[k + 1][0] = sistem1.tabDin[i].DrhogDpF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].DrholDtF[k + 1][0] = sistem1.tabDin[i].DrholDtF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].DrhogDtF[k + 1][0] = sistem1.tabDin[i].DrhogDtF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].valBO[k + 1][0] = sistem1.tabDin[i].valBO[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].tit[k + 1][0] = sistem1.tabDin[i].tit[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].rs[k + 1][0] = sistem1.tabDin[i].rs[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].cplF[k + 1][0] = sistem1.tabDin[i].cplF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].cpgF[k + 1][0] = sistem1.tabDin[i].cpgF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].valZ[k + 1][0] = sistem1.tabDin[i].valZ[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].HlF[k + 1][0] = sistem1.tabDin[i].HlF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].HgF[k + 1][0] = sistem1.tabDin[i].HgF[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].valdZdT[k + 1][0] = sistem1.tabDin[i].valdZdT[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].valdZdP[k + 1][0] = sistem1.tabDin[i].valdZdP[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].viscO[k + 1][0] = sistem1.tabDin[i].viscO[k][0] + sistem1.tabDin[i].delP;
            sistem1.tabDin[i].viscG[k + 1][0] = sistem1.tabDin[i].viscG[k][0] + sistem1.tabDin[i].delP;
        }
#pragma omp parallel for num_threads((*sistem1.vg1dSP).ntrd)
        for (int k = 1; k < sistem1.tabDin[i].npontosP + 1; k++) {
            ProFlu fluF;
            for (int j = 1; j < sistem1.tabDin[i].npontosT + 1; j++) {
                double pres = sistem1.tabDin[i].rholF[k][0];
                double temp = sistem1.tabDin[i].rholF[0][j];
                int indIni = sistem1.tabDin[i].celIni;
                fluF = sistem1.cell(indIni).flui;
                if (j == 1)
                    fluF.atualizaPropComp(pres, temp, -1, NULL, NULL, sistem1.config().pocinjec);
                else
                    fluF.atualizaPropComp(pres, temp, fluF.dCalculatedBeta, fluF.oCalculatedLiqComposition,
                                          fluF.oCalculatedVapComposition, sistem1.config().pocinjec);
                sistem1.tabDin[i].rholF[k][j] = fluF.MasEspoleo(pres, temp);
                sistem1.tabDin[i].rhogF[k][j] = fluF.MasEspGas(pres, temp);
                sistem1.tabDin[i].DrhogDpF[k][j] = fluF.drhodp(pres, temp);
                sistem1.tabDin[i].DrholDtF[k][j] = fluF.DrholDT(pres, temp);
                sistem1.tabDin[i].DrhogDtF[k][j] = fluF.drhodt(pres, temp);
                sistem1.tabDin[i].valBO[k][j] = fluF.BOFunc(pres, temp);
                sistem1.tabDin[i].valZ[k][j] = fluF.Zdran(pres, temp);
                sistem1.tabDin[i].valdZdT[k][j] = fluF.DZDT(pres, temp);
                sistem1.tabDin[i].valdZdP[k][j] = fluF.DZDP(pres, temp);
                sistem1.tabDin[i].tit[k][j] = fluF.FracMass(pres, temp);
                sistem1.tabDin[i].rs[k][j] = fluF.RS(pres, temp);
                sistem1.tabDin[i].cplF[k][j] = fluF.CalorLiq(pres, temp);
                sistem1.tabDin[i].cpgF[k][j] = fluF.CalorGas(pres, temp);
                sistem1.tabDin[i].HlF[k][j] = fluF.EntalpLiq(pres, temp);
                sistem1.tabDin[i].HgF[k][j] = fluF.EntalpGas(pres, temp);
                sistem1.tabDin[i].viscO[k][j] = fluF.ViscOleo(pres, temp);
                if (k == sistem1.tabDin[i].npontosP)
                    sistem1.tabDin[i].PBF[j - 1] = fluF.PB(pres, temp);
            }
        }
    }

    for (int i = 0; i < sistem1.ntabDin; i++) {
        int i0 = sistem1.tabDin[i].celIni;
        int i1 = sistem1.tabDin[i].celFim;
        for (int j = i0; j <= i1; j++) {
            sistem1.cell(j).flui.tabDin.TBF = sistem1.tabDin[i].TBF;
            sistem1.cell(j).flui.tabDin.PBF = sistem1.tabDin[i].PBF;
            sistem1.cell(j).flui.tabDin.rholF = sistem1.tabDin[i].rholF;
            sistem1.cell(j).flui.tabDin.rhogF = sistem1.tabDin[i].rhogF;
            sistem1.cell(j).flui.tabDin.DrholDpF = sistem1.tabDin[i].DrholDpF;
            sistem1.cell(j).flui.tabDin.DrhogDpF = sistem1.tabDin[i].DrhogDpF;
            sistem1.cell(j).flui.tabDin.DrholDtF = sistem1.tabDin[i].DrholDtF;
            sistem1.cell(j).flui.tabDin.DrhogDtF = sistem1.tabDin[i].DrhogDtF;
            sistem1.cell(j).flui.tabDin.valBO = sistem1.tabDin[i].valBO;
            sistem1.cell(j).flui.tabDin.tit = sistem1.tabDin[i].tit;
            sistem1.cell(j).flui.tabDin.rs = sistem1.tabDin[i].rs;
            sistem1.cell(j).flui.tabDin.cplF = sistem1.tabDin[i].cplF;
            sistem1.cell(j).flui.tabDin.cpgF = sistem1.tabDin[i].cpgF;
            sistem1.cell(j).flui.tabDin.valZ = sistem1.tabDin[i].valZ;
            sistem1.cell(j).flui.tabDin.HlF = sistem1.tabDin[i].HlF;
            sistem1.cell(j).flui.tabDin.HgF = sistem1.tabDin[i].HgF;
            sistem1.cell(j).flui.tabDin.valdZdT = sistem1.tabDin[i].valdZdT;
            sistem1.cell(j).flui.tabDin.valdZdP = sistem1.tabDin[i].valdZdP;
            sistem1.cell(j).flui.tabDin.viscO = sistem1.tabDin[i].viscO;
            sistem1.cell(j).flui.tabDin.viscG = sistem1.tabDin[i].viscG;
            sistem1.cell(j).flui.tabDin.delP = sistem1.tabDin[i].delP;
            sistem1.cell(j).flui.tabDin.delT = sistem1.tabDin[i].delT;
            sistem1.cell(j).flui.tabDin.npontosP = sistem1.tabDin[i].npontosP;
            sistem1.cell(j).flui.tabDin.npontosT = sistem1.tabDin[i].npontosT;
            sistem1.cell(j).flui.tabDin.pmax = sistem1.tabDin[i].pmax;
            sistem1.cell(j).flui.tabDin.pmin = sistem1.tabDin[i].pmin;
            sistem1.cell(j).flui.tabDin.tmax = sistem1.tabDin[i].tmax;
            sistem1.cell(j).flui.tabDin.tmin = sistem1.tabDin[i].tmin;
        }
    }

    for (int i = 0; i <= sistem1.ncel; i++) {
        sistem1.cell(i).flui.tabelaDinamica = 1;
        if (sistem1.cell(i).acsr.tipo == 15) {
            sistem1.cell(i).acsr.radialPoro.flup.tabelaDinamica = 1;
            sistem1.cell(i).acsr.radialPoro.preparaTabDin();
        }
        if (sistem1.cell(i).acsr.tipo == 16) {
            sistem1.cell(i).acsr.poroso2D.dados.flup.tabelaDinamica = 1;
            sistem1.cell(i).acsr.poroso2D.preparaTabDin();
        }
    }
}

void totalizaCicloRedeCompCego(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, int naflu, int i, int recb,
                               Vcr<double> nBloq, convergeNoPerm &noConv, int auxMaster, Vcr<int> coleta = Vcr<int>(), int ncoleta = 0) {
    Vcr<double> qostd(naflu + ncoleta);
    Vcr<double> rholliqIS(naflu + ncoleta);
    Vcr<double> rhogIS(naflu + ncoleta);
    Vcr<double> BSW(naflu + ncoleta);
    Vcr<double> Bet(naflu + ncoleta);
    Vcr<double> Rhocomp(naflu + ncoleta);
    Vcr<double> temp(naflu + ncoleta);
    Vcr<double> cpl(naflu + ncoleta);
    Vcr<double> cpg(naflu + ncoleta);
    Vcr<double> Mliq(naflu + ncoleta);
    Vcr<double> Qliq(naflu + ncoleta);
    Vcr<double> Mcomp(naflu + ncoleta);
    Vcr<double> Qcomp(naflu + ncoleta);
    Vcr<double> Mgas(naflu + ncoleta);
    Vcr<double> Denag(naflu + ncoleta);
    Vcr<double> titW(naflu + ncoleta);
    Vcr<double> vazMasLiqL(naflu + ncoleta);
    Vcr<double> vazMolV(naflu + ncoleta);
    for (int k = 0; k < naflu + ncoleta; k++) {
        vazMolV[k] = 0.;
    }
    noConv.flu = malha[i].cell(recb).flui;
    noConv.qostdTot = 0;
    noConv.moleomist = 0;
    noConv.moleomistNeg = 0;
    noConv.qwmist = 0;
    noConv.denagmist = 0;
    noConv.mliqmist = 0;
    noConv.mliqCmist = 0.;
    noConv.mgasmist = 0;
    noConv.cpmist = 0;
    noConv.tempmist = 0;
    noConv.Qlmist = 0.;
    noConv.betmist = 0.;
    noConv.qlmistStd = 0.;
    noConv.qgstdTot = 0.;
    noConv.moltot = 0.;
    noConv.mliqmistNeg = 0.;
    noConv.mliqCmistNeg = 0.;
    noConv.mgasmistNeg = 0.;
    noConv.qlmistStdNeg = 0.;
    noConv.qcmistStd = 0.;
    noConv.qcmistStdNeg = 0.;
    int kpos = 0;

    int aplicaFluiCol = 0;
    if ((*arqRede.vg1dSP).fluidoRede == 1)
        malha[i].cell(recb).acsr.injl.fluidocol = malha[i].cell(recb).fluicol;

    for (int k = 0; k < naflu; k++) {
        int ind = arqRede.malha[i].afluente[k];
        int fim = malha[ind].ncel - 1;
        if (nBloq[k] > 0. && malha[ind].cell(fim + 1).MC >= 0) {
            kpos++;
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                double ba = malha[ind].cell(fim).flui.BAFunc(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double fw = malha[ind].cell(fim).flui.BSW * ba /
                            (bo + ba * malha[ind].cell(fim).flui.BSW - malha[ind].cell(fim).flui.BSW * bo);
                double rhoO = malha[ind].cell(fim).flui.MasEspoleo(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double rhoW = malha[ind].cell(fim).flui.MasEspAgua(malha[ind].pGSup, malha[ind].cell(fim).temp);
                titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                    malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(fim).fluicol;
                    aplicaFluiCol = 1;
                }
                BSW[k] = malha[ind].cell(fim).flui.BSW;
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini * (1. - malha[ind].cell(fim).FW) *
                               (1. - malha[ind].cell(fim).bet) / (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                rhogIS[k] = malha[ind].cell(fim).flui.MasEspGas(pres, temp[k]);
                cpl[k] = (1. - Bet[k]) * malha[ind].cell(fim).flui.CalorLiq(pres, temp[k]) +
                         Bet[k] * malha[ind].cell(fim).fluicol.CalorLiq(pres, temp[k]);
                cpg[k] = malha[ind].cell(fim).flui.CalorGas(pres, temp[k]);
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                Denag[k] = malha[ind].cell(fim).flui.Denag;
                vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
                noConv.moleomist += vazMasLiqL[k];
                double pesoMolV = 0;
                for (int j = 0; j < malha[ind].cell(fim + 1).flui.npseudo; j++) {
                    pesoMolV += malha[ind].cell(fim + 1).flui.masMol[j] *
                                malha[ind].cell(fim + 1).flui.fracMol[j];
                }
                vazMolV[k] = (vazMasLiqL[k] + Mgas[k]) / pesoMolV;
                noConv.moltot += vazMolV[k];

                noConv.qostdTot += qostd[k];

                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.qwmist += qw1;
                noConv.denagmist += Denag[k] * qw1;
                noConv.mliqmist += Mliq[k];
                noConv.mliqCmist += Mcomp[k];
                noConv.TRmist += Mcomp[k] * malha[ind].cell(fim).fluicol.TR;
                noConv.mgasmist += Mgas[k];
                noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
                noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
                noConv.Qlmist += Qliq[k];
                noConv.qlmistStd += (qw1 + qostd[k] +
                                     Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                         malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStd += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                   malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        } else if (malha[ind].cell(fim + 1).MC < 0) {
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                BSW[k] = malha[ind].cell(fim).flui.BSW;
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini * (1. - malha[ind].cell(fim).FW) *
                               (1. - malha[ind].cell(fim).bet) / (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                Denag[k] = malha[ind].cell(fim).flui.Denag;
                double ba = malha[ind].cell(fim).flui.BAFunc(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double fw = malha[ind].cell(fim).flui.BSW * ba /
                            (bo + ba * malha[ind].cell(fim).flui.BSW - malha[ind].cell(fim).flui.BSW * bo);
                double rhoO = malha[ind].cell(fim).flui.MasEspoleo(malha[ind].pGSup, malha[ind].cell(fim).temp);
                double rhoW = malha[ind].cell(fim).flui.MasEspAgua(malha[ind].pGSup, malha[ind].cell(fim).temp);
                titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
                vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
                noConv.moleomistNeg += vazMasLiqL[k];
                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.mliqmistNeg += Mliq[k];
                noConv.mliqCmistNeg += Mcomp[k];
                noConv.mgasmistNeg += Mgas[k];
                noConv.qlmistStdNeg += (qw1 + qostd[k] +
                                        Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                            malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStdNeg += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                        malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        }
    }
    for (int k = naflu; k < naflu + ncoleta; k++) {
        int ind = coleta[k - naflu];
        if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1 && ind != auxMaster) {
            kpos++;
            int ini = 0;
            int indAflu = arqRede.malha[i].afluente[0];
            double bo = malha[ind].cell(ini).flui.BOFunc(
                malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double ba = malha[ind].cell(ini).flui.BAFunc(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double fw = malha[ind].cell(ini).flui.BSW * ba /
                        (bo + ba * malha[ind].cell(ini).flui.BSW - malha[indAflu].cell(ini).flui.BSW * bo);
            double rhoO = malha[ind].cell(ini).flui.MasEspoleo(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            double rhoW = malha[ind].cell(ini).flui.MasEspAgua(malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            titW[k] = (1 - fw) * rhoO / ((1 - fw) * rhoO + fw * rhoW);
            double pres = malha[indAflu].pGSup;
            malha[ind].calcTempFim();
            temp[k] = malha[ind].tempSup;
            Bet[k] = malha[ind].cell(ini).bet;

            if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(ini).fluicol;
                aplicaFluiCol = 1;
            }
            BSW[k] = malha[ind].cell(ini).flui.BSW;
            double rp = malha[ind].cell(ini).flui.MasEspLiq(pres, temp[k]);
            double rc = malha[ind].cell(ini).fluicol.MasEspFlu(pres, temp[k]);
            rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
            if (malha[ind].cell(ini + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                qostd[k] = -malha[ind].cell(ini + 1).Mliqini * (1. - malha[ind].cell(ini).FW) *
                           (1. - malha[ind].cell(ini).bet) / (bo * rholliqIS[k]);
            else
                qostd[k] = 0.;
            Rhocomp[k] = rc;
            rhogIS[k] = malha[ind].cell(ini).flui.MasEspGas(pres, temp[k]);
            cpl[k] = (1. - Bet[k]) * malha[ind].cell(ini).flui.CalorLiq(pres, temp[k]) +
                     Bet[k] * malha[ind].cell(ini).fluicol.CalorLiq(pres, temp[k]);
            cpg[k] = malha[ind].cell(ini).flui.CalorGas(pres, temp[k]);
            Mliq[k] = -malha[ind].cell(ini + 1).Mliqini;
            Qliq[k] = -malha[ind].cell(ini + 1).QL;
            Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
            Qcomp[k] = Qliq[k] * Bet[k];
            Mgas[k] = -malha[ind].cell(ini + 1).MC + malha[ind].cell(ini + 1).Mliqini;
            Denag[k] = malha[ind].cell(ini).flui.Denag;
            vazMasLiqL[k] = titW[k] * (Mliq[k] - Mcomp[k]);
            noConv.moleomist += vazMasLiqL[k];
            double pesoMolV = 0;
            for (int j = 0; j < malha[ind].cell(ini + 1).flui.npseudo; j++) {
                pesoMolV += malha[ind].cell(ini + 1).flui.masMol[j] *
                            malha[ind].cell(ini + 1).flui.fracMol[j];
            }
            vazMolV[k] = (vazMasLiqL[k] + Mgas[k]) / pesoMolV;
            noConv.moltot += vazMolV[k];

            noConv.qostdTot += qostd[k];

            double qw1;
            if ((1. - BSW[k]) > 0)
                qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
            else
                qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
            noConv.qwmist += qw1;
            noConv.denagmist += Denag[k] * qw1;
            noConv.mliqmist += Mliq[k];
            noConv.mliqCmist += Mcomp[k];
            noConv.TRmist += Mcomp[k] * malha[ind].cell(ini).fluicol.TR;
            noConv.mgasmist += Mgas[k];
            noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
            noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
            noConv.Qlmist += Qliq[k];
            noConv.qlmistStd += (qw1 + qostd[k] +
                                 Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.qcmistStd += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                 malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                               malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
        }
    }

    if (kpos > 0) {
        if (recb == 0 && noConv.moltot > 0) {
            for (int jgrup = 0; jgrup < malha[i].cell(recb).flui.npseudo; jgrup++) {
                malha[i].cell(recb).flui.fracMol[jgrup] = 0.;
                noConv.flu.fracMol[jgrup] = 0.;
            }
            for (int iaflu = 0; iaflu < naflu; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int fim = malha[ind].ncel - 1;
                for (int jgrup = 0; jgrup < malha[i].cell(0).flui.npseudo; jgrup++) {
                    malha[i].cell(recb).flui.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            for (int iaflu = naflu; iaflu < naflu + ncoleta; iaflu++) {
                int ind = coleta[iaflu - naflu];
                int ini = 0;
                for (int jgrup = 0; jgrup < malha[i].cell(0).flui.npseudo; jgrup++) {
                    malha[i].cell(recb).flui.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            malha[i].cell(recb).flui.atualizaPropCompStandard();
            noConv.flu.atualizaPropCompStandard();
        } else if (noConv.moltot > 0) {
            for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++)
                noConv.flu.fracMol[jgrup] = 0.;
            for (int iaflu = 0; iaflu < naflu; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int fim = malha[ind].ncel - 1;
                for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(fim + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            for (int iaflu = naflu; iaflu < naflu + ncoleta; iaflu++) {
                int ind = arqRede.malha[i].afluente[iaflu];
                int ini = 0;
                for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                    noConv.flu.fracMol[jgrup] +=
                        vazMolV[iaflu] * malha[ind].cell(ini + 1).flui.fracMol[jgrup] / noConv.moltot;
                }
            }
            noConv.flu.atualizaPropCompStandard();
        } else {
            for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
                noConv.flu.fracMol[jgrup] =
                    malha[i].fluiRevRede.fracMol[jgrup];
            }
            noConv.flu.atualizaPropCompStandard();
        }
        if (noConv.mliqCmist > (*arqRede.vg1dSP).localtiny * 1e-15)
            noConv.TRmist /= noConv.mliqCmist;
        else
            noConv.TRmist = 0.;
    } else {
        for (int jgrup = 0; jgrup < noConv.flu.npseudo; jgrup++) {
            noConv.flu.fracMol[jgrup] =
                malha[i].fluiRevRede.fracMol[jgrup];
        }
        noConv.flu.atualizaPropCompStandard();
        noConv.LVisL = malha[i].fluiRevRede.LVisL;
        noConv.LVisH = malha[i].fluiRevRede.LVisH;
        noConv.bswmist = malha[i].fluiRevRede.BSW;
        noConv.denagmist = malha[i].fluiRevRede.Denag;
        noConv.betmist = 0.;
        int ind = arqRede.malha[i].afluente[0];
        int fim = malha[ind].ncel - 1;
        noConv.tempmist = malha[ind].cell(fim).temp;
        noConv.TRmist = malha[ind].cell(fim).fluicol.TR;
    }
    noConv.qlmistStd *= 86400.;
    noConv.qlmistStdNeg *= 86400.;
    noConv.qcmistStd *= 86400.;
    noConv.qcmistStdNeg *= 86400.;
    kpos = 0;
}

void cicloRedeCompCego(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, vector<tramoPart> &bloq) {
    int narq = arqRede.nsisprod - 0 * indativo;
    (*arqRede.vg1dSP).relax = arqRede.relax;
    Vcr<int> Resolv(narq, 0);
    Vcr<int> IndNorma(narq, 0);
    Vcr<int> prepTab(narq, 0);
    int ciclomalha = 0;
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].perm == 0) {
            Resolv[i] = 1;
            ciclomalha++;
        }
    }
    while (ciclomalha < narq - 1) {
        int ntrdLocal = (*arqRede.vg1dSP).ntrd;
        if ((*arqRede.vg1dSP).ntrd > narq)
            ntrdLocal = narq;
#pragma omp parallel for num_threads(ntrdLocal)
        for (int i = 0; i < narq; i++) {
            int permAflu = 0;
            if (arqRede.malha[i].nafluente > 0) {
                for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                    int indaflu = arqRede.malha[i].afluente[iaflu];
                    permAflu += arqRede.malha[indaflu].perm;
                }
            }
            if ((arqRede.malha[i].nafluente == 0 || (permAflu == 0 && arqRede.malha[i].perm == 1)) &&
                Resolv[i] == 0 && inativo[i] == 0) {
                if (malha[i].cell(0).acsr.tipo == 1 && malha[i].cell(0).acsr.injg.QGas >= 0.) {
                    preparaTabDin(malha[i]);
                    prepTab[i] = 1;
                } else if (malha[i].cell(0).acsr.tipo == 2 && malha[i].cell(0).acsr.injl.QLiq >= 0.) {
                    preparaTabDin(malha[i]);
                    prepTab[i] = 1;
                } else if (malha[i].cell(0).acsr.tipo == 10 && malha[i].cell(0).acsr.injm.MassP >= 0.) {
                    preparaTabDin(malha[i]);
                    prepTab[i] = 1;
                }
                Resolv[i] = 1;
                ciclomalha++;
            }
        }
        int resolvGlob = 0;
        for (int i = 0; i < narq; i++)
            resolvGlob += Resolv[i];

        while (resolvGlob < narq) {
            for (int i = 0; i < narq; i++) {

                int naflu = arqRede.malha[i].nafluente;
                int cicloaflu = 0;
                if (arqRede.malha[i].nafluente > 0 && Resolv[i] == 0) {
                    for (int j = 0; j < narq; j++) {
                        for (int k = 0; k < naflu; k++) {
                            if (arqRede.malha[i].afluente[k] == j && Resolv[j] == 1)
                                cicloaflu++;
                        }
                    }
                    if (cicloaflu == naflu && inativo[i] == 0) {

                        int ind = arqRede.malha[i].afluente[0];
                        int nderiva = arqRede.malha[ind].ncoleta;

                        int trocaDeriva = 0;
                        int novaDeriva;
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (arqRede.malha[aux].principal == 1) {
                                trocaDeriva = 1;
                                novaDeriva = aux;
                            }
                        }

                        Vcr<int> ordCol(nderiva);
                        vector<double> dcol;
                        Vcr<int> carregado(nderiva, 0);
                        if (trocaDeriva == 0) {
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                double mult = 1.0;
                                if (malha[aux].config().ConContEntrada == 2)
                                    mult = -1.;
                                dcol.push_back(malha[aux].cell(0).duto.dia * mult);
                            }
                            sort(dcol.begin(), dcol.end());

                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                double mult = 1.0;
                                if (malha[aux].config().ConContEntrada == 2)
                                    mult = -1.;
                                for (int icol2 = 0; icol2 < nderiva; icol2++) {
                                    if (fabs(dcol[icol2] - malha[aux].cell(0).duto.dia * mult) < 1.e-15 && carregado[icol2] == 0) {
                                        ordCol[icol2] = aux;
                                        carregado[icol2] = 1;
                                        icol2 = nderiva;
                                    }
                                }
                            }
                        } else {
                            int icol2 = 0;
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                if (aux != novaDeriva) {
                                    ordCol[icol2] = aux;
                                    carregado[icol2] = 1;
                                    icol2++;
                                }
                            }
                            ordCol[nderiva - 1] = novaDeriva;
                            carregado[nderiva - 1] = 1;
                        }

                        convergeNoPerm noConv;
                        convergeNoPerm noConv1;
                        convergeNoPerm noConv2;
                        int totBloq = 0;
                        int totBloq1 = 0;
                        int totBloq2 = 0;
                        Vcr<double> nBloq(naflu);
                        Vcr<double> nBloq1(naflu);
                        Vcr<double> nBloq2(naflu);

                        int ncolneg = 0;
                        Vcr<int> colneg(nderiva);
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (malha[aux].cell(0).acsr.tipo == 1) {
                                if (malha[aux].cell(0).acsr.injg.QGas < 0) {
                                    colneg[ncolneg] = aux;
                                    ncolneg++;
                                }
                            } else if (malha[aux].cell(0).acsr.tipo == 2) {
                                if (malha[aux].cell(0).acsr.injl.QLiq < 0) {
                                    colneg[ncolneg] = aux;
                                    ncolneg++;
                                }
                            } else if (malha[aux].cell(0).acsr.tipo == 10) {
                                if ((malha[aux].cell(0).acsr.injm.MassP + malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassC) < 0) {
                                    colneg[ncolneg] = aux;
                                    ncolneg++;
                                }
                            }
                        }

                        int auxMaster;
                        auxMaster = ordCol[nderiva - 1];

                        for (int iaflu = 0; iaflu < naflu; iaflu++) {
                            int indaflu = arqRede.malha[auxMaster].afluente[iaflu];
                            malha[indaflu].fluiRevRede = malha[auxMaster].cell(0).flui;
                            malha[indaflu].tempRev = malha[auxMaster].cell(0).temp;
                            if ((*arqRede.vg1dSP).fluidoRede == 0)
                                malha[indaflu].config().razCompGasReves = malha[auxMaster].cell(0).acsr.injg.razCompGas;
                        }
                        for (int icol = 0; icol < nderiva - 1; icol++) {
                            int aux = ordCol[icol];
                            if (arqRede.malha[aux].ncoleta == 0)
                                malha[aux].fluiRevRede = malha[aux].config().flup[0];
                        }

                        int col2;
                        avaliaBloq(i, arqRede, nBloq, nBloq1, nBloq2, totBloq, totBloq1, totBloq2, col2, naflu, bloq);
                        if (totBloq > 0)
                            totalizaCicloRedeCompCego(malha, arqRede, inativo, indativo, naflu, i, 0, nBloq, noConv, auxMaster, colneg, ncolneg);
                        if (totBloq1 > 0) {
                            totalizaCicloRedeCompCego(malha, arqRede, inativo, indativo, naflu, i, 1, nBloq1, noConv1, auxMaster);
                            noConv1.qgstdTot = ((noConv1.mgasmist + noConv1.moleomist) /
                                                (noConv1.dengmist * 1.225)) *
                                               noConv1.flu.dStockTankVaporMassFraction;
                            noConv1.qgstdTotNeg = ((noConv1.mgasmistNeg + noConv1.moleomistNeg) /
                                                   (noConv1.dengmist * 1.225)) *
                                                  noConv1.flu.dStockTankVaporMassFraction;
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                malha[i].cell(1).acsr.tipo = 2;
                                if (malha[i].cell(1).acsr.injl.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injl.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injl.bet = noConv1.betmist;
                                malha[i].cell(1).acsr.injl.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injl.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injl.fluidocol.TR = noConv1.TRmist;

                                noConv1.qgstdTot *= 86400;

                                malha[i].cell(1).acsr.injl.QLiq = noConv1.qlmistStd + noConv1.qlmistStdNeg;
                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                malha[i].cell(1).acsr.tipo = 1;
                                if (malha[i].cell(1).acsr.injg.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injg.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injg.seco = 0;
                                malha[i].cell(1).acsr.injg.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injg.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injg.fluidocol.TR = noConv1.TRmist;
                                noConv1.qgstdTot *= 86400;
                                malha[i].cell(1).acsr.injg.QGas = noConv1.qgstdTot + noConv1.qgstdTotNeg * 86400.;

                                if (fabs(malha[i].cell(1).acsr.injg.QGas) > 1e-15)
                                    malha[i].cell(1).acsr.injg.razCompGas = (noConv1.qcmistStd + noConv1.qcmistStdNeg) /
                                                                              malha[i].cell(1).acsr.injg.QGas;
                                else
                                    malha[i].cell(1).acsr.injg.razCompGas = 0.;

                            } else {
                                malha[i].cell(1).acsr.tipo = 10;
                                if (malha[i].cell(1).acsr.injm.FluidoPro.flashCompleto == 0)
                                    malha[i].cell(1).acsr.injm.FluidoPro.flashCompleto = 2;
                                malha[i].cell(1).acsr.injm.temp = noConv1.tempmist;
                                malha[i].cell(1).acsr.injm.FluidoPro = noConv1.flu;
                                malha[i].cell(1).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                malha[i].cell(1).acsr.injm.fluidocol.TR = noConv1.TRmist;
                                noConv1.qgstdTot *= 86400;

                                malha[i].cell(1).acsr.injm.MassP = noConv1.mliqmist - noConv1.mliqCmist +
                                                                     noConv1.mliqmistNeg - noConv1.mliqCmistNeg;
                                malha[i].cell(1).acsr.injm.MassC = noConv1.mliqCmist + noConv1.mliqCmistNeg;
                                malha[i].cell(1).acsr.injm.MassG = noConv1.mgasmist + noConv1.mgasmistNeg;
                            }
                        }
                        if (totBloq2 > 0) {
                            totalizaCicloRedeCompCego(malha, arqRede, inativo, indativo, naflu, col2, 1, nBloq2, noConv2, auxMaster);
                            noConv2.qgstdTot = ((noConv2.mgasmist + noConv2.moleomist) /
                                                (noConv2.dengmist * 1.225)) *
                                               noConv2.flu.dStockTankVaporMassFraction;
                            noConv2.qgstdTotNeg = ((noConv2.mgasmistNeg + noConv2.moleomistNeg) /
                                                   (noConv2.dengmist * 1.225)) *
                                                  noConv2.flu.dStockTankVaporMassFraction;
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                malha[col2].cell(1).acsr.tipo = 2;
                                if (malha[col2].cell(1).acsr.injl.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injl.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injl.bet = noConv2.betmist;
                                malha[col2].cell(1).acsr.injl.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injl.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injl.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injl.fluidocol.TR = noConv2.TRmist;

                                noConv2.qgstdTot *= 86400;

                                malha[col2].cell(1).acsr.injl.QLiq = noConv2.qlmistStd + noConv2.qlmistStdNeg;
                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                malha[col2].cell(1).acsr.tipo = 1;
                                if (malha[col2].cell(1).acsr.injg.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injg.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injg.seco = 0;
                                malha[col2].cell(1).acsr.injg.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injg.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injg.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injg.fluidocol.TR = noConv2.TRmist;
                                noConv2.qgstdTot *= 86400;
                                malha[col2].cell(1).acsr.injg.QGas = noConv2.qgstdTot + noConv2.qgstdTotNeg * 86400.;
                                if (fabs(malha[col2].cell(1).acsr.injg.QGas) > 1e-15)
                                    malha[col2].cell(1).acsr.injg.razCompGas = (noConv2.qcmistStd + noConv2.qcmistStdNeg) /
                                                                                 malha[col2].cell(1).acsr.injg.QGas;
                                else
                                    malha[col2].cell(1).acsr.injg.razCompGas = 0.;
                            } else {
                                malha[col2].cell(1).acsr.tipo = 10;
                                if (malha[col2].cell(1).acsr.injm.FluidoPro.flashCompleto == 0)
                                    malha[col2].cell(1).acsr.injm.FluidoPro.flashCompleto = 2;
                                malha[col2].cell(1).acsr.injm.temp = noConv2.tempmist;
                                malha[col2].cell(1).acsr.injm.FluidoPro = noConv2.flu;
                                malha[col2].cell(1).acsr.injm.FluidoPro.tabelaDinamica = 0;
                                malha[col2].cell(1).acsr.injm.fluidocol.TR = noConv2.TRmist;
                                noConv2.qgstdTot *= 86400;

                                malha[col2].cell(1).acsr.injm.MassP = noConv2.mliqmist - noConv2.mliqCmist +
                                                                        noConv2.mliqmistNeg - noConv2.mliqCmistNeg;
                                malha[col2].cell(1).acsr.injm.MassC = noConv2.mliqCmist + noConv2.mliqCmistNeg;
                                malha[col2].cell(1).acsr.injm.MassG = noConv2.mgasmist + noConv2.mgasmistNeg;
                            }
                        }

                        for (int icol = 0; icol < nderiva; icol++) {
                            if (icol == nderiva - 1) {
                                int aux = ordCol[icol];
                                if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                    malha[aux].cell(0).acsr.injl.FluidoPro = noConv.flu;
                                } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                    malha[aux].cell(0).acsr.injg.FluidoPro = noConv.flu;
                                } else {
                                    malha[aux].cell(0).acsr.injm.FluidoPro = noConv.flu;
                                }
                                malha[aux].cell(0).flui = noConv.flu;
                            } else {
                                int aux = ordCol[icol];
                                if (malha[aux].config().ConContEntrada != 2) {
                                    if ((*arqRede.vg1dSP).fluidoRede != 0) {
                                        malha[aux].cell(0).acsr.injl.FluidoPro = noConv.flu;
                                        malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                        if ((malha[aux].cell(0).acsr.injl.QLiq) < 0) {
                                            malha[aux].cell(0).acsr.injl.FluidoPro = malha[aux].fluiRevRede;
                                        }

                                    } else {
                                        malha[aux].cell(0).acsr.injg.FluidoPro = noConv.flu;
                                        malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                        if ((malha[aux].cell(0).acsr.injg.QGas) < 0) {
                                            malha[aux].cell(0).acsr.injg.FluidoPro = malha[aux].fluiRevRede;
                                        }
                                    }
                                } else {
                                    malha[aux].cell(0).acsr.injm.FluidoPro = noConv.flu;
                                    malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                }
                                malha[aux].cell(0).flui = noConv.flu;
                            }
                        }

                        for (int icol = 0; icol < nderiva; icol++) {
                            if (icol == nderiva - 1) {
                                int aux = ordCol[icol];

                                if (nderiva > 1) {
                                    if (arqRede.malha[aux].perm == 1) {
                                        preparaTabDin(malha[aux]);
                                        prepTab[aux] = 1;
                                    }
                                } else {
                                    if (malha[aux].cell(0).acsr.tipo == 1 &&
                                        (malha[aux].cell(0).acsr.injg.QGas >= 0. || arqRede.malha[aux].ncoleta == 0)) {
                                        if (malha[aux].cell(0).acsr.injg.QGas >= 0.)
                                            preparaTabDin(malha[aux]);
                                        else {
                                            for (int kontaCel = 0; kontaCel <= malha[aux].ncel; kontaCel++) {
                                                malha[aux].cell(kontaCel).flui = malha[aux].fluiRevRede;
                                            }
                                            preparaTabDin(malha[aux]);
                                            malha[aux].cell(0).acsr.injg.FluidoPro = malha[aux].fluiRevRede;
                                        }
                                        prepTab[aux] = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 2 &&
                                               (malha[aux].cell(0).acsr.injl.QLiq >= 0. || arqRede.malha[aux].ncoleta == 0)) {
                                        if (malha[aux].cell(0).acsr.injl.QLiq >= 0.)
                                            preparaTabDin(malha[aux]);
                                        else {
                                            for (int kontaCel = 0; kontaCel <= malha[aux].ncel; kontaCel++) {
                                                malha[aux].cell(kontaCel).flui = malha[aux].fluiRevRede;
                                            }
                                            malha[aux].cell(0).acsr.injl.FluidoPro = malha[aux].fluiRevRede;
                                            preparaTabDin(malha[aux]);
                                        }
                                        prepTab[aux] = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 10 &&
                                               (malha[aux].cell(0).acsr.injm.MassP >= 0. || arqRede.malha[aux].ncoleta == 0)) {
                                        if (malha[aux].cell(0).acsr.injm.MassP >= 0.)
                                            preparaTabDin(malha[aux]);
                                        else {
                                            for (int kontaCel = 0; kontaCel <= malha[aux].ncel; kontaCel++) {
                                                malha[aux].cell(kontaCel).flui = malha[aux].fluiRevRede;
                                            }
                                            malha[aux].cell(0).acsr.injm.FluidoPro = malha[aux].fluiRevRede;
                                            preparaTabDin(malha[aux]);
                                        }
                                        prepTab[aux] = 1;
                                    }
                                }

                                Resolv[aux] = 1;
                                resolvGlob += Resolv[aux];
                                ciclomalha++;
                            } else {
                                int aux = ordCol[icol];

                                if (malha[aux].config().ConContEntrada != 2) {
                                    if (arqRede.malha[aux].perm == 1) {
                                        preparaTabDin(malha[aux]);
                                        prepTab[aux] = 1;
                                    }
                                    Resolv[aux] = 1;
                                    resolvGlob += Resolv[aux];
                                    ciclomalha++;
                                } else {
                                    if (arqRede.malha[aux].perm == 1) {
                                        preparaTabDin(malha[aux]);
                                        prepTab[aux] = 1;
                                    }
                                    Resolv[aux] = 1;
                                    resolvGlob += Resolv[aux];
                                    ciclomalha++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int semprepara = 1;
    while (semprepara > 0) {
        semprepara = 0;
        for (int i = 0; i < narq; i++) {
            if (prepTab[i] == 0 && arqRede.malha[i].perm == 1) {
                if (arqRede.malha[i].ncoleta > 0) {
                    int fimPrepara = 0;
                    for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                        int indCol = arqRede.malha[i].coleta[j];
                        if ((malha[indCol].cell(0).acsr.tipo == 1 && (malha[indCol].cell(0).acsr.injg.QGas +
                                                                        malha[indCol].cell(1).acsr.injg.QGas) >= 0) ||
                            (malha[indCol].cell(0).acsr.tipo == 2 && (malha[indCol].cell(0).acsr.injl.QLiq +
                                                                        malha[indCol].cell(1).acsr.injl.QLiq) >= 0) ||
                            (malha[indCol].cell(0).acsr.tipo == 10 && (malha[indCol].cell(0).acsr.injm.MassP +
                                                                         malha[indCol].cell(1).acsr.injm.MassP) >= 0) ||
                            prepTab[indCol] == 1) {
                            for (int k = 0; k < malha[i].ncel; k++) {
                                malha[i].cell(k).flui = malha[indCol].cell(0).flui;
                                if (malha[i].cell(k).acsr.tipo == 1 && malha[i].cell(k).acsr.injg.QGas < 0) {
                                    malha[i].cell(k).acsr.injg.FluidoPro = malha[indCol].cell(0).flui;
                                } else if (malha[i].cell(k).acsr.tipo == 2 && malha[i].cell(k).acsr.injl.QLiq < 0) {
                                    malha[i].cell(k).acsr.injl.FluidoPro = malha[indCol].cell(0).flui;
                                } else if (malha[i].cell(k).acsr.tipo == 10 && malha[i].cell(k).acsr.injm.MassP < 0) {
                                    malha[i].cell(k).acsr.injm.FluidoPro = malha[indCol].cell(0).flui;
                                }
                            }
                            preparaTabDin(malha[i]);
                            prepTab[i] = 1;
                            fimPrepara = 1;
                        }
                        if (fimPrepara == 1)
                            break;
                    }
                    if (fimPrepara == 0)
                        semprepara++;
                }
            }
        }
    }
}

void totalizaCicloRede(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, int naflu, int i, int recb,
                       Vcr<double> nBloq, convergeNoPerm &noConv, int auxMaster, Vcr<int> coleta = Vcr<int>(), int ncoleta = 0) {
    Vcr<double> qostd(naflu + ncoleta);
    Vcr<double> RGO(naflu + ncoleta);
    Vcr<double> rhololeoST(naflu + ncoleta);
    Vcr<double> rhogST(naflu + ncoleta);
    Vcr<double> rholliqIS(naflu + ncoleta);
    Vcr<double> rhogIS(naflu + ncoleta);
    Vcr<double> API(naflu + ncoleta);
    Vcr<double> BSW(naflu + ncoleta);
    Vcr<double> Bet(naflu + ncoleta);
    Vcr<double> Rhocomp(naflu + ncoleta);
    Vcr<double> temp(naflu + ncoleta);
    Vcr<double> cpl(naflu + ncoleta);
    Vcr<double> cpg(naflu + ncoleta);
    Vcr<double> Mliq(naflu + ncoleta);
    Vcr<double> Qliq(naflu + ncoleta);
    Vcr<double> Mcomp(naflu + ncoleta);
    Vcr<double> Qcomp(naflu + ncoleta);
    Vcr<double> Mgas(naflu + ncoleta);
    Vcr<double> Deng(naflu + ncoleta);
    Vcr<double> yco2(naflu + ncoleta);
    Vcr<double> Denag(naflu + ncoleta);
    noConv.flu = malha[i].cell(recb).flui;
    noConv.tL = 0.;
    if (malha[i].config().flashCompleto == 1 && malha[i].config().tabent.tmin > noConv.tL)
        noConv.tL = malha[i].config().tabent.tmin + 1.;
    noConv.tH = 70.;
    if (malha[i].config().flashCompleto == 1 && malha[i].config().tabent.tmax < noConv.tH)
        noConv.tH = malha[i].config().tabent.tmax - 1.;
    noConv.qostdTot = 0;
    noConv.qgstdTot = 0;
    noConv.qgstdTotNeg = 0;
    noConv.dengmist = 0.;
    noConv.yco2mist = 0;
    noConv.apimist = 0.;
    noConv.qwmist = 0;
    noConv.denagmist = 0;
    noConv.mliqmist = 0;
    noConv.mliqmistNeg = 0;
    noConv.mliqCmist = 0.;
    noConv.mliqCmistNeg = 0.;
    noConv.mgasmist = 0;
    noConv.mgasmistNeg = 0;
    noConv.cpmist = 0;
    noConv.tempmist = 0;
    noConv.LVisL = 0.;
    noConv.LVisH = 0.;
    noConv.Qlmist = 0.;
    noConv.betmist = 0.;
    noConv.qlmistStd = 0.;
    noConv.qlmistStdNeg = 0.;
    noConv.qcmistStd = 0.;
    noConv.qcmistStdNeg = 0.;
    noConv.TRmist = 0.;
    int aplicaFluiCol = 0;
    if ((*arqRede.vg1dSP).fluidoRede == 1)
        malha[i].cell(recb).acsr.injl.fluidocol = malha[i].cell(recb).fluicol;
    int kpos = 0;
    for (int k = 0; k < naflu; k++) {
        int ind = arqRede.malha[i].afluente[k];
        int fim = malha[ind].ncel - 1;
        if (nBloq[k] > 0. && malha[ind].cell(fim + 1).MC >= 0) {
            kpos++;
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                RGO[k] = malha[ind].cell(fim).flui.RGO;
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                    malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(fim).fluicol;
                    aplicaFluiCol = 1;
                }
                BSW[k] = malha[ind].cell(fim).flui.BSW;
                rhogST[k] = malha[ind].cell(fim).flui.Deng * 1.225;
                API[k] = malha[ind].cell(fim).flui.API;
                rhololeoST[k] = 1000. * 141.5 / (131.5 + API[k]);
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini *
                               (1. - malha[ind].cell(fim).FW) * (1. - malha[ind].cell(fim).bet) /
                               (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                rhogIS[k] = malha[ind].cell(fim).flui.MasEspGas(pres, temp[k]);
                cpl[k] = (1. - Bet[k]) * malha[ind].cell(fim).flui.CalorLiq(pres, temp[k]) +
                         Bet[k] * malha[ind].cell(fim).fluicol.CalorLiq(pres, temp[k]);
                cpg[k] = malha[ind].cell(fim).flui.CalorGas(pres, temp[k]);
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                Deng[k] = malha[ind].cell(fim).flui.Deng;
                yco2[k] = malha[ind].cell(fim).flui.yco2;
                Denag[k] = malha[ind].cell(fim).flui.Denag;

                noConv.qostdTot += qostd[k];
                double qgtemp;
                if ((*arqRede.vg1dSP).tipoFluidoRedeGlob == 0)
                    qgtemp = qostd[k] * RGO[k];
                else
                    qgtemp = Mgas[k] / rhogST[k] + qostd[k] * malha[ind].cell(fim).flui.RS(pres, temp[k]) * 6.29 / 35.31467;
                if (qgtemp <= 1e-15)
                    qgtemp = malha[ind].cell(fim + 1).MC / rhogST[k];
                noConv.qgstdTot += qgtemp;
                noConv.dengmist += Deng[k] * qgtemp;
                noConv.yco2mist += yco2[k] * qgtemp;
                noConv.apimist += (141.5 / (API[k] + 131.5)) * qostd[k];
                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.qwmist += qw1;
                noConv.denagmist += Denag[k] * qw1;
                noConv.mliqmist += Mliq[k];
                noConv.mliqCmist += Mcomp[k];
                noConv.TRmist += Mcomp[k] * malha[ind].cell(fim).fluicol.TR;
                noConv.mgasmist += Mgas[k];
                noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
                noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
                noConv.LVisL += qostd[k] * malha[ind].cell(fim).flui.VisOM(noConv.tL);
                noConv.LVisH += qostd[k] * malha[ind].cell(fim).flui.VisOM(noConv.tH);
                noConv.Qlmist += Qliq[k];
                noConv.qlmistStd += (qw1 + qostd[k] +
                                     Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                         malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStd += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                   malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        } else if (malha[ind].cell(fim + 1).MC < 0) {
            if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1) {
                double bo = malha[ind].cell(fim).flui.BOFunc(
                    malha[ind].pGSup, malha[ind].cell(fim).temp);
                RGO[k] = malha[ind].cell(fim).flui.RGO;
                double pres = malha[ind].pGSup;
                malha[ind].calcTempFim();
                temp[k] = malha[ind].tempSup;
                Bet[k] = malha[ind].cell(fim).bet;

                BSW[k] = malha[ind].cell(fim).flui.BSW;
                rhogST[k] = malha[ind].cell(fim).flui.Deng * 1.225;
                double rp = malha[ind].cell(fim).flui.MasEspLiq(pres, temp[k]);
                double rc = malha[ind].cell(fim).fluicol.MasEspFlu(pres, temp[k]);
                rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
                if (malha[ind].cell(fim + 1).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                    qostd[k] = malha[ind].cell(fim + 1).Mliqini *
                               (1. - malha[ind].cell(fim).FW) * (1. - malha[ind].cell(fim).bet) /
                               (bo * rholliqIS[k]);
                else
                    qostd[k] = 0.;
                Rhocomp[k] = rc;
                Mliq[k] = malha[ind].cell(fim + 1).Mliqini;
                Qliq[k] = malha[ind].cell(fim + 1).QL;
                Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
                Qcomp[k] = Qliq[k] * Bet[k];
                Mgas[k] = malha[ind].cell(fim + 1).MC - malha[ind].cell(fim + 1).Mliqini;
                Denag[k] = malha[ind].cell(fim).flui.Denag;

                double qgtemp;
                if ((*arqRede.vg1dSP).tipoFluidoRedeGlob == 0)
                    qgtemp = qostd[k] * RGO[k];
                else
                    qgtemp = Mgas[k] / rhogST[k] + qostd[k] * malha[ind].cell(fim).flui.RS(pres, temp[k]) * 6.29 / 35.31467;
                if (qgtemp <= 1e-15)
                    qgtemp = malha[ind].cell(fim + 1).MC / rhogST[k];
                noConv.qgstdTotNeg += qgtemp;
                double qw1;
                if ((1. - BSW[k]) > 0)
                    qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
                else
                    qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
                noConv.mliqmistNeg += Mliq[k];
                noConv.mliqCmistNeg += Mcomp[k];
                noConv.mgasmistNeg += Mgas[k];
                noConv.qlmistStdNeg += (qw1 + qostd[k] +
                                        Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                            malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
                noConv.qcmistStdNeg += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                        malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            }
        }
    }
    for (int k = naflu; k < naflu + ncoleta; k++) {
        int ind = coleta[k - naflu];
        int iaflu = 0;
        int indAflu = arqRede.malha[i].afluente[iaflu];
        while (arqRede.malha[indAflu].perm == 0) {
            iaflu++;
            indAflu = arqRede.malha[i].afluente[iaflu];
        }
        if (inativo[ind] == 0 && arqRede.malha[ind].perm == 1 && ind != auxMaster) {
            kpos++;
            int ini = 0;
            double bo = malha[ind].cell(ini).flui.BOFunc(
                malha[indAflu].pGSup, malha[ind].cell(ini).temp);
            RGO[k] = malha[ind].cell(ini).flui.RGO;
            double pres = malha[indAflu].pGSup;
            malha[ind].calcTempFim();
            temp[k] = malha[ind].tempSup;
            Bet[k] = malha[ind].cell(ini).bet;

            if (Bet[k] > (*arqRede.vg1dSP).localtiny && aplicaFluiCol == 0 && (*arqRede.vg1dSP).fluidoRede == 1) {
                malha[i].cell(recb).acsr.injl.fluidocol = malha[ind].cell(ini).fluicol;
                aplicaFluiCol = 1;
            }
            BSW[k] = malha[ind].cell(ini).flui.BSW;
            rhogST[k] = malha[ind].cell(ini).flui.Deng * 1.225;
            API[k] = malha[ind].cell(ini).flui.API;
            rhololeoST[k] = 1000. * 141.5 / (131.5 + API[k]);
            double rp = malha[ind].cell(ini).flui.MasEspLiq(pres, temp[k]);
            double rc = malha[ind].cell(ini).fluicol.MasEspFlu(pres, temp[k]);
            rholliqIS[k] = (1. - Bet[k]) * rp + Bet[k] * rc;
            if (malha[ind].cell(ini).flui.dStockTankVaporMassFraction < 1. - 1e-15)
                qostd[k] = -malha[ind].cell(ini + 1).Mliqini *
                           (1. - malha[ind].cell(ini).FW) * (1. - malha[ind].cell(ini).bet) /
                           (bo * rholliqIS[k]);
            else
                qostd[k] = 0.;
            Rhocomp[k] = rc;
            rhogIS[k] = malha[ind].cell(ini).flui.MasEspGas(pres, temp[k]);
            cpl[k] = (1. - Bet[k]) * malha[ind].cell(ini).flui.CalorLiq(pres, temp[k]) +
                     Bet[k] * malha[ind].cell(ini).fluicol.CalorLiq(pres, temp[k]);
            cpg[k] = malha[ind].cell(ini).flui.CalorGas(pres, temp[k]);
            Mliq[k] = -malha[ind].cell(ini + 1).Mliqini;
            Qliq[k] = -malha[ind].cell(ini + 1).QL;
            Mcomp[k] = Mliq[k] * Bet[k] * Rhocomp[k] / rholliqIS[k];
            Qcomp[k] = Qliq[k] * Bet[k];
            Mgas[k] = -malha[ind].cell(ini + 1).MC + malha[ind].cell(ini + 1).Mliqini;
            Deng[k] = malha[ind].cell(ini).flui.Deng;
            yco2[k] = malha[ind].cell(ini).flui.yco2;
            Denag[k] = malha[ind].cell(ini).flui.Denag;

            noConv.qostdTot += qostd[k];
            double qgtemp;
            if ((*arqRede.vg1dSP).tipoFluidoRedeGlob == 0)
                qgtemp = qostd[k] * RGO[k];
            else
                qgtemp = Mgas[k] / rhogST[k] + qostd[k] * malha[ind].cell(ini).flui.RS(pres, temp[k]) * 6.29 / 35.31467;
            if (qgtemp <= 1e-15)
                qgtemp = -malha[ind].cell(ini + 1).MC / rhogST[k];
            noConv.qgstdTot += qgtemp;
            noConv.dengmist += Deng[k] * qgtemp;
            noConv.yco2mist += yco2[k] * qgtemp;
            noConv.apimist += (141.5 / (API[k] + 131.5)) * qostd[k];
            double qw1;
            if ((1. - BSW[k]) > 0)
                qw1 = qostd[k] * BSW[k] / (1. - BSW[k]);
            else
                qw1 = Mliq[k] * (1. - Bet[k]) / (1000. * Denag[k]);
            noConv.qwmist += qw1;
            noConv.denagmist += Denag[k] * qw1;
            noConv.mliqmist += Mliq[k];
            noConv.mliqCmist += Mcomp[k];
            noConv.TRmist += Mcomp[k] * malha[ind].cell(ini).fluicol.TR;
            noConv.mgasmist += Mgas[k];
            noConv.cpmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]);
            noConv.tempmist += (Mliq[k] * cpl[k] + Mgas[k] * cpg[k]) * temp[k];
            noConv.LVisL += qostd[k] * malha[ind].cell(ini).flui.VisOM(noConv.tL);
            noConv.LVisH += qostd[k] * malha[ind].cell(ini).flui.VisOM(noConv.tH);
            noConv.Qlmist += Qliq[k];
            noConv.qlmistStd += (qw1 + qostd[k] +
                                 Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                     malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.qcmistStd += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                                 malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
            noConv.betmist += (Qcomp[k] * malha[i].cell(recb).fluicol.MasEspFlu(pres, temp[k]) /
                               malha[i].cell(recb).fluicol.MasEspFlu(1.001, 15.));
        }
    }
    if (kpos > 0) {
        if (noConv.qostdTot > (*arqRede.vg1dSP).localtiny) {
            noConv.RGOmist = noConv.qgstdTot / noConv.qostdTot;
            noConv.apimist = noConv.apimist / noConv.qostdTot;
            noConv.apimist = 141.5 / noConv.apimist - 131.5;
            noConv.LVisL = noConv.LVisL / noConv.qostdTot;
            noConv.LVisH = noConv.LVisH / noConv.qostdTot;
        } else {
            noConv.RGOmist = 140000.;
            noConv.apimist = malha[i].cell(0).flui.API;
            noConv.LVisL = malha[i].cell(recb).flui.VisOM(noConv.tL);
            noConv.LVisH = malha[i].cell(recb).flui.VisOM(noConv.tH);
        }
        if (noConv.qgstdTot > (*arqRede.vg1dSP).localtiny) {
            noConv.dengmist = noConv.dengmist / noConv.qgstdTot;
            noConv.yco2mist = noConv.yco2mist / noConv.qgstdTot;
        } else {
            noConv.dengmist = malha[i].cell(recb).flui.Deng;
            noConv.yco2mist = malha[i].cell(recb).flui.yco2;
        }
        if ((noConv.qostdTot + noConv.qwmist) > (*arqRede.vg1dSP).localtiny)
            noConv.bswmist = noConv.qwmist / (noConv.qostdTot + noConv.qwmist);
        else
            noConv.bswmist = malha[i].cell(recb).flui.BSW;
        if (noConv.qwmist > (*arqRede.vg1dSP).localtiny)
            noConv.denagmist = noConv.denagmist / noConv.qwmist;
        else
            noConv.denagmist = malha[i].cell(recb).flui.Denag;
        if (noConv.cpmist > (*arqRede.vg1dSP).localtiny)
            noConv.tempmist = noConv.tempmist / noConv.cpmist;
        else
            noConv.tempmist = 0.;
        if (noConv.qlmistStd > (*arqRede.vg1dSP).localtiny)
            noConv.betmist = noConv.betmist / noConv.qlmistStd;
        else
            noConv.betmist = 0.;
        if (noConv.mliqCmist > (*arqRede.vg1dSP).localtiny * 1e-15)
            noConv.TRmist /= noConv.mliqCmist;
        else
            noConv.TRmist = 0.;
    } else {
        noConv.RGOmist = malha[i].fluiRevRede.RGO;
        noConv.apimist = malha[i].fluiRevRede.API;
        noConv.LVisL = malha[i].fluiRevRede.LVisL;
        noConv.LVisH = malha[i].fluiRevRede.LVisH;
        noConv.dengmist = malha[i].fluiRevRede.Deng;
        noConv.yco2mist = malha[i].fluiRevRede.yco2;
        noConv.bswmist = malha[i].fluiRevRede.BSW;
        noConv.denagmist = malha[i].fluiRevRede.Denag;
        noConv.betmist = 0.;
        int ind = arqRede.malha[i].afluente[0];
        int fim = malha[ind].ncel - 1;
        noConv.tempmist = malha[ind].cell(fim).temp;
        noConv.TRmist = malha[ind].cell(fim).fluicol.TR;
    }
    noConv.qlmistStd *= 86400;
    noConv.qlmistStdNeg *= 86400;
    noConv.qcmistStd *= 86400;
    noConv.qcmistStdNeg *= 86400;

    noConv.flu.RGO = noConv.RGOmist;
    noConv.flu.Deng = noConv.dengmist;
    noConv.flu.Denag = noConv.denagmist;
    noConv.flu.yco2 = noConv.yco2mist;
    noConv.flu.API = noConv.apimist;
    noConv.flu.BSW = noConv.bswmist;
    noConv.flu.TempL = noConv.tL;
    noConv.flu.TempH = noConv.tH;
    noConv.flu.LVisL = noConv.LVisL;
    noConv.flu.LVisH = noConv.LVisH;
    noConv.flu.RenovaFluido();
}

double cicloRede(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo, vector<noRede> &normaEvol, vector<tramoPart> &bloq) {
    double norma = 0.;
    int nnos = 0;
    int narq = arqRede.nsisprod - 0 * indativo;
    (*arqRede.vg1dSP).relax = arqRede.relax;
    Vcr<int> Resolv(narq, 0);
    Vcr<int> IndNorma(narq, 0);
    Vcr<int> Retorna(narq, 0);
    if ((*arqRede.vg1dSP).iterRede >= 5) {
    }
    cout << "#################Iniciando ciclo iterativo: " << (*arqRede.vg1dSP).iterRede << endl;
    int ciclomalha = 0;
    double valor;
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].perm == 0) {
            Resolv[i] = 1;
            ciclomalha++;
        }
    }
    int testaNMalha = ciclomalha - narq + 1;
    while (ciclomalha < narq - 1) {
        int ntrdLocal = (*arqRede.vg1dSP).ntrd;
        if ((*arqRede.vg1dSP).ntrd > narq)
            ntrdLocal = narq;
#pragma omp parallel for num_threads(ntrdLocal)
        for (int i = 0; i < narq; i++) {
            int permAflu = 0;
            if (arqRede.malha[i].nafluente > 0) {
                for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                    int indaflu = arqRede.malha[i].afluente[iaflu];
                    permAflu += arqRede.malha[indaflu].perm;
                }
            }
            if ((arqRede.malha[i].nafluente == 0 || (permAflu == 0 && arqRede.malha[i].perm == 1)) &&
                Resolv[i] == 0 && inativo[i] == 0) {
                if (arqRede.malha[i].perm == 1) {
                    (*arqRede.vg1dSP).qualTramo = i;
                    if (malha[i].config().ConContEntrada == 0) {
                        int reverso = 0;
                        if (malha[i].cell(0).acsr.tipo == 1) {
                            if (malha[i].cell(0).acsr.injg.QGas < 0)
                                reverso = 1;
                        } else if (malha[i].cell(0).acsr.tipo == 2) {
                            if (malha[i].cell(0).acsr.injl.QLiq < 0)
                                reverso = 1;
                        } else if (malha[i].cell(0).acsr.tipo == 10) {
                            if ((malha[i].cell(0).acsr.injm.MassC +
                                 malha[i].cell(0).acsr.injm.MassG + malha[i].cell(0).acsr.injm.MassP) < 0)
                                reverso = 1;
                        }
                        malha[i].modoPerm = 1;
                        if (reverso == 0) {
                            if (malha[i].config().chokep.abertura[0] > 0.6) {
                                if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0) { // mudancaChute

                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 0 && malha[i].gasCell(0).tipoCC == 0)
                                        malha[i].config().gasinj.vazgas[0] = 150000 * malha[i].gasCell(0).duto.area / (*arqRede.vg1dSP).arearef;
                                    if (malha[i].config().lingas == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        double ciclo = 1.1e9;
                                        int konta = 0;
                                        double multVazGas;
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm();
                                        double testaPres1 = malha[i].gasCell(0).pres;
                                        malha[i].config().gasinj.vazgas[0] *= 1.05;
                                        malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        double testaPres2 = malha[i].gasCell(0).pres;
                                        if (testaPres1 < testaPres2) {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 1.05;
                                            else
                                                multVazGas = 0.95;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        } else {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 0.95;
                                            else
                                                multVazGas = 1.05;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        }
                                        while (ciclo > 0.9e9 && konta < 10) {
                                            malha[i].gasCell(0).tipoCC = 1;
                                            if (konta > 0)
                                                malha[i].buscaProdPfundoPerm();
                                            malha[i].gasCell(0).tipoCC = 0;
                                            ciclo = malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres, konta);
                                            if (ciclo > 0.9e9) {
                                                malha[i].config().gasinj.vazgas[0] *= multVazGas;
                                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                                            }
                                            konta++;
                                            if (konta > 10) {
                                                ciclo = 0.8e10;
                                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                                if (malha[i].config().transiente == 0)
                                                    NumError(
                                                        "Falha na busca de solucao permanente para sistema com GL e condicao de pressao");
                                            }
                                        }
                                    } // mudancaChute
                                    else {
                                        valor = malha[i].buscaProdPfundoPerm();
                                    }
                                } else { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        malha[i].gasCell(0).tipoCC = 0;
                                    }
                                    valor = malha[i].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                    //}//mudancaChute
                                }
                            } else {
                                if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0) { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 0 && malha[i].gasCell(0).tipoCC == 0)
                                        malha[i].config().gasinj.vazgas[0] = 150000 * malha[i].gasCell(0).duto.area / (*arqRede.vg1dSP).arearef;
                                    if (malha[i].config().lingas == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        double ciclo = 1.1e9;
                                        int konta = 0;
                                        double multVazGas;
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm2();
                                        double testaPres1 = malha[i].gasCell(0).pres;
                                        malha[i].config().gasinj.vazgas[0] *= 1.05;
                                        malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        double testaPres2 = malha[i].gasCell(0).pres;
                                        if (testaPres1 < testaPres2) {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 1.05;
                                            else
                                                multVazGas = 0.95;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        } else {
                                            if (malha[i].gasCell(0).pres > testaPres1)
                                                multVazGas = 0.95;
                                            else
                                                multVazGas = 1.05;
                                            malha[i].config().gasinj.vazgas[0] /= 1.05;
                                        }
                                        while (ciclo > 0.9e9 && konta < 10) {
                                            malha[i].gasCell(0).tipoCC = 1;
                                            if (konta > 0)
                                                malha[i].buscaProdPfundoPerm2();
                                            malha[i].gasCell(0).tipoCC = 0;
                                            ciclo = malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres, konta);
                                            if (ciclo > 0.9e9) {
                                                malha[i].config().gasinj.vazgas[0] *= multVazGas;
                                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                                            }
                                            konta++;
                                            if (konta > 10) {
                                                ciclo = 0.8e10;
                                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                                if (malha[i].config().transiente == 0)
                                                    NumError(
                                                        "Falha na busca de solucao permanente para sistema com GL e condicao de pressao");
                                            }
                                        }
                                    } // mudancaChute
                                    else
                                        valor = malha[i].buscaProdPfundoPerm2();
                                } else { // mudancaChute
                                    if (malha[i].config().lingas == 1 && malha[i].config().gasinj.chuteVaz == 1 && malha[i].gasCell(0).tipoCC == 0) {
                                        malha[i].gasCell(0).tipoCC = 1;
                                        malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        malha[i].gasCell(0).tipoCC = 0;
                                    }
                                    valor = malha[i].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                    //}//mudancaChute
                                }
                            }
                        } else {
                            if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0)
                                valor = malha[i].buscaProdPfundoPermRev();
                            else
                                valor = malha[i].buscaProdPfundoPermRev(malha[i].cell(0).pres);
                        }
                    } // alteracao4XXXXXXXXXXXX
                    else {
                        double chutemass = 0.;
                        if (malha[i].cell(0).acsr.tipo == 2) {
                            chutemass = malha[i].cell(0).acsr.injl.QLiq;
                            if (fabs(chutemass) < 1e-5)
                                chutemass = 100.;
                        } else {
                            chutemass = malha[i].cell(0).acsr.injg.QGas;
                            if (fabs(chutemass) < 1e-5)
                                chutemass = 1000.;
                        }
                        if (malha[i].config().chokep.abertura[0] > 0.6)
                            malha[i].buscaProdPresPresPerm(chutemass);
                        else
                            malha[i].buscaProdPresPresPerm2(chutemass);
                    }
                } else {
                    valor = 0;
                    int icol = arqRede.malha[i].coleta[0];
                    malha[i].pGSup = malha[icol].cell(0).pres;
                    malha[i].tGSup = malha[icol].cell(0).temp;
                    malha[i].tempSup = malha[icol].cell(0).temp;
                }
                if (valor < -1e9 || valor > 1e9) {
                    inativo[i] = 1;
                    inativoColetor(i, arqRede, inativo, narq);
                    if (valor < -1e9) {
                        aviso(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    } else {
                        aviso2(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    }
                    (*arqRede.vg1dSP).restartRede = 1;
                    Retorna[i] = -1;
                }
                Resolv[i] = 1;
#pragma omp atomic
                ciclomalha++;
                cout << "convergencia tramo: " << i << endl;
            }
        }
        int resolvGlob = 0;
        for (int i = 0; i < narq; i++) {
            resolvGlob += Resolv[i];
            if (Retorna[i] < 0)
                return 9000.;
        }

        while (resolvGlob < narq) {
            for (int i = 0; i < narq; i++) {
                int naflu = arqRede.malha[i].nafluente;
                int cicloaflu = 0;
                if (arqRede.malha[i].nafluente > 0 && Resolv[i] == 0) {
                    for (int j = 0; j < narq; j++) {
                        for (int k = 0; k < naflu; k++) {
                            if (arqRede.malha[i].afluente[k] == j && Resolv[j] == 1)
                                cicloaflu++;
                        }
                    }
                    if (cicloaflu == naflu && inativo[i] == 0) {

                        int ind = arqRede.malha[i].afluente[0];
                        int nderiva = arqRede.malha[ind].ncoleta;

                        int trocaDeriva = 0;
                        int novaDeriva;
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (arqRede.malha[aux].principal == 1) {
                                trocaDeriva = 1;
                                novaDeriva = aux;
                            }
                        }

                        Vcr<int> ordCol(nderiva);
                        vector<double> dcol;
                        Vcr<int> carregado(nderiva, 0);
                        if (trocaDeriva == 0) {
                            int derivaP = 0;
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                derivaP += arqRede.malha[aux].derivaPrincipal;
                            }
                            if (derivaP == 0) {
                                for (int icol = 0; icol < nderiva; icol++) {
                                    int aux = arqRede.malha[ind].coleta[icol];
                                    if (malha[aux].config().perm == 1) {
                                        double mult = 1.0;
                                        if (malha[aux].config().ConContEntrada == 2)
                                            mult = -1.;
                                        dcol.push_back(malha[aux].cell(0).duto.dia * mult);
                                    }
                                }
                                sort(dcol.begin(), dcol.end());

                                for (int icol = 0; icol < nderiva; icol++) {
                                    int aux = arqRede.malha[ind].coleta[icol];
                                    if (malha[aux].config().perm == 1) {
                                        double mult = 1.0;
                                        if (malha[aux].config().ConContEntrada == 2)
                                            mult = -1.;
                                        for (int icol2 = 0; icol2 < nderiva; icol2++) {
                                            if (fabs(dcol[icol2] - malha[aux].cell(0).duto.dia * mult) < 1.e-15 && carregado[icol2] == 0) {
                                                ordCol[icol2] = aux;
                                                carregado[icol2] = 1;
                                                icol2 = nderiva;
                                            }
                                        }
                                    }
                                }
                            } else {
                                int icol;
                                int aux;
                                for (icol = 0; icol < nderiva; icol++) {
                                    aux = arqRede.malha[ind].coleta[icol];
                                    if (arqRede.malha[aux].derivaPrincipal == 1)
                                        break;
                                }
                                int icol2 = 0;
                                while (icol2 < nderiva - 1) {
                                    int aux2 = arqRede.malha[ind].coleta[icol2];
                                    if (malha[aux2].config().perm == 1) {
                                        if (aux2 != aux) {
                                            ordCol[icol2] = aux2;
                                            carregado[icol2] = 1;
                                            icol2++;
                                        }
                                    }
                                }
                                ordCol[nderiva - 1] = aux;
                                carregado[nderiva - 1] = 1;
                            }
                        } else {
                            int icol2 = 0;
                            for (int icol = 0; icol < nderiva; icol++) {
                                int aux = arqRede.malha[ind].coleta[icol];
                                if (malha[aux].config().perm == 1) {
                                    if (aux != novaDeriva) {
                                        ordCol[icol2] = aux;
                                        carregado[icol2] = 1;
                                        icol2++;
                                    }
                                }
                            }
                            ordCol[nderiva - 1] = novaDeriva;
                            carregado[nderiva - 1] = 1;
                        }

                        convergeNoPerm noConv;
                        convergeNoPerm noConv1;
                        convergeNoPerm noConv2;
                        int totBloq = 0;
                        int totBloq1 = 0;
                        int totBloq2 = 0;
                        Vcr<double> nBloq(naflu);
                        Vcr<double> nBloq1(naflu);
                        Vcr<double> nBloq2(naflu);
                        int ncolneg = 0;
                        Vcr<int> colneg(nderiva);
                        for (int icol = 0; icol < nderiva; icol++) {
                            int aux = arqRede.malha[ind].coleta[icol];
                            if (malha[aux].config().perm == 1) {
                                if (malha[aux].cell(0).acsr.tipo == 1) {
                                    if (malha[aux].cell(0).acsr.injg.QGas < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                } else if (malha[aux].cell(0).acsr.tipo == 2) {
                                    if (malha[aux].cell(0).acsr.injl.QLiq < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                } else if (malha[aux].cell(0).acsr.tipo == 10) {
                                    if ((malha[aux].cell(0).acsr.injm.MassP + malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassC) < 0) {
                                        colneg[ncolneg] = aux;
                                        ncolneg++;
                                    }
                                }
                            }
                        }

                        int auxMaster;
                        auxMaster = ordCol[nderiva - 1];

                        for (int iaflu = 0; iaflu < naflu; iaflu++) {
                            int celProp = 1 - verificaFonteDuplaReversa(malha, auxMaster);
                            int indaflu = arqRede.malha[auxMaster].afluente[iaflu];
                            malha[indaflu].fluiRevRede = malha[auxMaster].cell(celProp).flui;
                            malha[indaflu].tempRev = malha[auxMaster].cell(0).temp;
                            if ((*arqRede.vg1dSP).fluidoRede == 0)
                                malha[indaflu].config().razCompGasReves = malha[auxMaster].cell(celProp).acsr.injg.razCompGas;
                        }
                        for (int icol = 0; icol < nderiva - 1; icol++) {
                            int aux = ordCol[icol];
                            if (arqRede.malha[aux].ncoleta == 0) {
                                malha[aux].fluiRevRede = malha[aux].config().flup[0];
                            }
                        }

                        int col2;
                        avaliaBloq(i, arqRede, nBloq, nBloq1, nBloq2, totBloq, totBloq1, totBloq2, col2, naflu, bloq);
                        if (totBloq > 0)
                            totalizaCicloRede(malha, arqRede, inativo, indativo, naflu, i, 0, nBloq, noConv, auxMaster, colneg, ncolneg);
                        int recb = 1;
                        if (totBloq == 0)
                            recb = 0;
                        if (totBloq1 > 0) {
                            totalizaCicloRede(malha, arqRede, inativo, indativo, naflu, i, recb, nBloq1, noConv1, auxMaster);
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                if (malha[i].cell(recb).acsr.injl.FluidoPro.flashCompleto == 2)
                                    malha[i].cell(recb).acsr.injl.FluidoPro.flashCompleto = 0;
                                malha[i].cell(recb).acsr.tipo = 2;
                                malha[i].cell(recb).acsr.injl.bet = noConv1.betmist;
                                malha[i].cell(recb).acsr.injl.temp = noConv1.tempmist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.RGO = noConv1.RGOmist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.Deng = noConv1.dengmist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.Denag = noConv1.denagmist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.yco2 = noConv1.yco2mist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.API = noConv1.apimist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.BSW = noConv1.bswmist;
                                malha[i].cell(recb).acsr.injl.FluidoPro.TempL = noConv1.tL;
                                malha[i].cell(recb).acsr.injl.FluidoPro.TempH = noConv1.tH;
                                malha[i].cell(recb).acsr.injl.FluidoPro.LVisL = noConv1.LVisL;
                                malha[i].cell(recb).acsr.injl.FluidoPro.LVisH = noConv1.LVisH;
                                malha[i].cell(recb).acsr.injl.FluidoPro.RenovaFluido();
                                malha[i].cell(recb).acsr.injl.fluidocol.TR = noConv1.TRmist;

                                malha[i].cell(recb).acsr.injl.QLiq = noConv1.qlmistStd + noConv1.qlmistStdNeg;

                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                if (malha[i].cell(recb).acsr.injg.FluidoPro.flashCompleto == 2)
                                    malha[i].cell(recb).acsr.injg.FluidoPro.flashCompleto = 0;
                                malha[i].cell(recb).acsr.tipo = 1;
                                malha[i].cell(recb).acsr.injg.seco = 0;
                                malha[i].cell(recb).acsr.injg.temp = noConv1.tempmist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.RGO = noConv1.RGOmist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.API = noConv1.apimist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.BSW = noConv1.bswmist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.TempL = noConv1.tL;
                                malha[i].cell(recb).acsr.injg.FluidoPro.TempH = noConv1.tH;
                                malha[i].cell(recb).acsr.injg.FluidoPro.LVisL = noConv1.LVisL;
                                malha[i].cell(recb).acsr.injg.FluidoPro.LVisH = noConv1.LVisH;
                                malha[i].cell(recb).acsr.injg.FluidoPro.Deng = noConv1.dengmist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.yco2 = noConv1.yco2mist;
                                malha[i].cell(recb).acsr.injg.FluidoPro.RenovaFluido();
                                malha[i].cell(recb).acsr.injg.fluidocol.TR = noConv1.TRmist;

                                noConv1.qgstdTot *= 86400.;
                                malha[i].cell(recb).acsr.injg.QGas = noConv1.qgstdTot + noConv1.qgstdTotNeg * 86400.;
                                if (fabs(malha[i].cell(recb).acsr.injg.QGas) > 1e-15)
                                    malha[i].cell(recb).acsr.injg.razCompGas = (noConv1.qcmistStd + noConv1.qcmistStdNeg) /
                                                                                 malha[i].cell(recb).acsr.injg.QGas;
                                else
                                    malha[i].cell(recb).acsr.injg.razCompGas = 0.;

                            } else {
                                if (malha[i].cell(recb).acsr.injm.FluidoPro.flashCompleto == 2)
                                    malha[i].cell(recb).acsr.injm.FluidoPro.flashCompleto = 0;
                                malha[i].cell(recb).acsr.tipo = 10;
                                malha[i].cell(recb).acsr.injm.temp = noConv1.tempmist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.RGO = noConv1.RGOmist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.Deng = noConv1.dengmist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.Denag = noConv1.denagmist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.yco2 = noConv1.yco2mist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.API = noConv1.apimist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.BSW = noConv1.bswmist;
                                malha[i].cell(recb).acsr.injm.FluidoPro.TempL = noConv1.tL;
                                malha[i].cell(recb).acsr.injm.FluidoPro.TempH = noConv1.tH;
                                malha[i].cell(recb).acsr.injm.FluidoPro.LVisL = noConv1.LVisL;
                                malha[i].cell(recb).acsr.injm.FluidoPro.LVisH = noConv1.LVisH;
                                malha[i].cell(recb).acsr.injm.FluidoPro.RenovaFluido();
                                malha[i].cell(recb).acsr.injm.fluidocol.TR = noConv1.TRmist;

                                malha[i].cell(recb).acsr.injm.MassP = noConv1.mliqmist - noConv1.mliqCmist +
                                                                        noConv1.mliqmistNeg - noConv1.mliqCmistNeg;
                                malha[i].cell(recb).acsr.injm.MassC = noConv1.mliqCmist + noConv1.mliqCmistNeg;
                                malha[i].cell(recb).acsr.injm.MassG = noConv1.mgasmist + noConv1.mgasmistNeg;
                            }
                        }
                        if (totBloq2 > 0) {
                            totalizaCicloRede(malha, arqRede, inativo, indativo, naflu, col2, recb, nBloq2, noConv2, auxMaster);
                            if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                if (malha[col2].cell(recb).acsr.injl.FluidoPro.flashCompleto == 2)
                                    malha[col2].cell(recb).acsr.injl.FluidoPro.flashCompleto = 0;
                                malha[col2].cell(recb).acsr.tipo = 2;
                                malha[col2].cell(recb).acsr.injl.bet = noConv2.betmist;
                                malha[col2].cell(recb).acsr.injl.temp = noConv2.tempmist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.RGO = noConv2.RGOmist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.Deng = noConv2.dengmist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.Denag = noConv2.denagmist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.yco2 = noConv2.yco2mist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.API = noConv2.apimist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.BSW = noConv2.bswmist;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.TempL = noConv2.tL;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.TempH = noConv2.tH;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.LVisL = noConv2.LVisL;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.LVisH = noConv2.LVisH;
                                malha[col2].cell(recb).acsr.injl.FluidoPro.RenovaFluido();
                                malha[col2].cell(recb).acsr.injl.fluidocol.TR = noConv2.TRmist;

                                malha[col2].cell(recb).acsr.injl.QLiq = noConv2.qlmistStd + noConv2.qlmistStdNeg;

                            } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                if (malha[col2].cell(recb).acsr.injg.FluidoPro.flashCompleto == 2)
                                    malha[col2].cell(recb).acsr.injg.FluidoPro.flashCompleto = 0;
                                malha[col2].cell(recb).acsr.tipo = 1;
                                malha[col2].cell(recb).acsr.injg.seco = 0;
                                malha[col2].cell(recb).acsr.injg.temp = noConv2.tempmist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.RGO = noConv2.RGOmist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.API = noConv2.apimist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.BSW = noConv2.bswmist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.TempL = noConv2.tL;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.TempH = noConv2.tH;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.LVisL = noConv2.LVisL;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.LVisH = noConv2.LVisH;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.Deng = noConv2.dengmist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.yco2 = noConv2.yco2mist;
                                malha[col2].cell(recb).acsr.injg.FluidoPro.RenovaFluido();
                                malha[col2].cell(recb).acsr.injg.fluidocol.TR = noConv2.TRmist;

                                noConv2.qgstdTot *= 86400.;
                                malha[col2].cell(recb).acsr.injg.QGas = noConv2.qgstdTot + noConv2.qgstdTotNeg * 86400.;
                                if (fabs(malha[col2].cell(recb).acsr.injg.QGas) > 1e-15)
                                    malha[col2].cell(recb).acsr.injg.razCompGas = (noConv2.qcmistStd + noConv2.qcmistStdNeg) /
                                                                                    malha[col2].cell(recb).acsr.injg.QGas;
                                else
                                    malha[col2].cell(recb).acsr.injg.razCompGas = 0.;

                            } else {
                                if (malha[col2].cell(recb).acsr.injm.FluidoPro.flashCompleto == 2)
                                    malha[col2].cell(recb).acsr.injm.FluidoPro.flashCompleto = 0;
                                malha[col2].cell(recb).acsr.tipo = 10;
                                malha[col2].cell(recb).acsr.injm.temp = noConv2.tempmist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.RGO = noConv2.RGOmist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.Deng = noConv2.dengmist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.Denag = noConv2.denagmist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.yco2 = noConv2.yco2mist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.API = noConv2.apimist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.BSW = noConv2.bswmist;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.TempL = noConv2.tL;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.TempH = noConv2.tH;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.LVisL = noConv2.LVisL;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.LVisH = noConv2.LVisH;
                                malha[col2].cell(recb).acsr.injm.FluidoPro.RenovaFluido();
                                malha[col2].cell(recb).acsr.injm.fluidocol.TR = noConv2.TRmist;

                                malha[col2].cell(recb).acsr.injm.MassP = noConv2.mliqmist - noConv2.mliqCmist +
                                                                           noConv2.mliqmistNeg - noConv2.mliqCmistNeg;
                                malha[col2].cell(recb).acsr.injm.MassC = noConv2.mliqCmist + noConv2.mliqCmistNeg;
                                malha[col2].cell(recb).acsr.injm.MassG = noConv2.mgasmist + noConv2.mgasmistNeg;
                            }
                        }

                        for (int icol = 0; icol < nderiva; icol++) {
                            if (icol == nderiva - 1) {
                                int aux = ordCol[icol];
                                (*arqRede.vg1dSP).qualTramo = aux;
                                if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                    if (malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto == 2)
                                        malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto = 0;
                                    malha[aux].cell(0).acsr.tipo = 2;
                                    malha[aux].cell(0).acsr.injl.bet = noConv.betmist;
                                    malha[aux].cell(0).acsr.injl.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.RGO = noConv.RGOmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.Deng = noConv.dengmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.Denag = noConv.denagmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.yco2 = noConv.yco2mist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.API = noConv.apimist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.BSW = noConv.bswmist;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.TempL = noConv.tL;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.TempH = noConv.tH;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.LVisL = noConv.LVisL;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.LVisH = noConv.LVisH;
                                    malha[aux].cell(0).acsr.injl.FluidoPro.RenovaFluido();
                                    malha[aux].cell(0).acsr.injl.fluidocol.TR = noConv.TRmist;

                                    norma += pow((malha[aux].cell(0).acsr.injl.QLiq - (noConv.qlmistStd + noConv.qlmistStdNeg)) /
                                                     (noConv.qlmistStd + noConv.qlmistStdNeg),
                                                 2.);
                                    malha[aux].cell(0).acsr.injl.QLiq = noConv.qlmistStd + noConv.qlmistStdNeg;

                                    malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                    if (malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto == 2)
                                        malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto = 0;
                                    malha[aux].cell(0).acsr.tipo = 1;
                                    malha[aux].cell(0).acsr.injg.seco = 0;
                                    malha[aux].cell(0).acsr.injg.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.RGO = noConv.RGOmist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.API = noConv.apimist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.BSW = noConv.bswmist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.TempL = noConv.tL;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.TempH = noConv.tH;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.LVisL = noConv.LVisL;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.LVisH = noConv.LVisH;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.Deng = noConv.dengmist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.yco2 = noConv.yco2mist;
                                    malha[aux].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                                    malha[aux].cell(0).acsr.injg.fluidocol.TR = noConv.TRmist;

                                    noConv.qgstdTot *= 86400.;
                                    norma += pow((malha[aux].cell(0).acsr.injg.QGas - (noConv.qgstdTot +
                                                                                         noConv.qgstdTotNeg * 86400.)) /
                                                     (noConv.qgstdTot + noConv.qgstdTotNeg * 86400.),
                                                 2.);
                                    malha[aux].cell(0).acsr.injg.QGas = noConv.qgstdTot +
                                                                          noConv.qgstdTotNeg * 86400.;

                                    if (fabs(malha[aux].cell(0).acsr.injg.QGas) > 1e-15)
                                        malha[aux].cell(0).acsr.injg.razCompGas = (noConv.qcmistStd + noConv.qcmistStdNeg) /
                                                                                    malha[aux].cell(0).acsr.injg.QGas;
                                    else
                                        malha[aux].cell(0).acsr.injg.razCompGas = 0.;

                                    malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                } else {
                                    if (malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto == 2)
                                        malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto = 0;
                                    malha[aux].cell(0).acsr.tipo = 10;
                                    malha[aux].cell(0).acsr.injm.temp = noConv.tempmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.RGO = noConv.RGOmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.Deng = noConv.dengmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.Denag = noConv.denagmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.yco2 = noConv.yco2mist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.API = noConv.apimist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.BSW = noConv.bswmist;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.TempL = noConv.tL;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.TempH = noConv.tH;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.LVisL = noConv.LVisL;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.LVisH = noConv.LVisH;
                                    malha[aux].cell(0).acsr.injm.FluidoPro.RenovaFluido();
                                    malha[aux].cell(0).acsr.injm.fluidocol.TR = noConv.TRmist;

                                    double massT = malha[aux].cell(0).acsr.injm.MassP + malha[aux].cell(0).acsr.injm.MassC +
                                                   malha[aux].cell(0).acsr.injm.MassG;
                                    norma += pow((massT - (noConv.mliqmist + noConv.mgasmist + noConv.mliqmistNeg + noConv.mgasmistNeg)) /
                                                     (noConv.mliqmist + noConv.mgasmist + noConv.mliqmistNeg + noConv.mgasmistNeg),
                                                 2.);
                                    malha[aux].cell(0).acsr.injm.MassP = noConv.mliqmist - noConv.mliqCmist +
                                                                           noConv.mliqmistNeg - noConv.mliqCmistNeg;
                                    malha[aux].cell(0).acsr.injm.MassC = noConv.mliqCmist + noConv.mliqCmistNeg;
                                    malha[aux].cell(0).acsr.injm.MassG = noConv.mgasmist + noConv.mgasmistNeg;

                                    malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                }

                                int trocou = trocaFonteColetor(malha, aux);

                                double pini = malha[aux].cell(0).pres;
                                if (arqRede.malha[aux].perm == 1) {
                                    malha[aux].modoPerm = 1;
                                    int reverso = 0;
                                    if (malha[aux].cell(0).acsr.tipo == 1) {
                                        if (malha[aux].cell(0).acsr.injg.QGas < 0)
                                            reverso = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 2) {
                                        if (malha[aux].cell(0).acsr.injl.QLiq < 0)
                                            reverso = 1;
                                    } else if (malha[aux].cell(0).acsr.tipo == 10) {
                                        if ((malha[aux].cell(0).acsr.injm.MassC +
                                             malha[aux].cell(0).acsr.injm.MassG + malha[aux].cell(0).acsr.injm.MassP) < 0)
                                            reverso = 1;
                                    }
                                    if (reverso == 0) {
                                        if (malha[aux].config().chokep.abertura[0] > 0.6) {
                                            if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0)
                                                valor = malha[aux].buscaProdPfundoPerm();
                                            else
                                                valor = malha[aux].buscaProdPfundoPerm(malha[i].cell(0).pres);
                                        } else {
                                            if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0)
                                                valor = malha[aux].buscaProdPfundoPerm2();
                                            else
                                                valor = malha[aux].buscaProdPfundoPerm2(malha[i].cell(0).pres);
                                        }
                                    } else {
                                        if ((*arqRede.vg1dSP).iterRede == 0 && arqRede.chute == 0)
                                            valor = malha[aux].buscaProdPfundoPermRev();
                                        else
                                            valor = malha[aux].buscaProdPfundoPermRev(malha[aux].cell(0).pres);
                                    }
                                } else
                                    valor = 0;
                                if (trocou == -1)
                                    retornaFonteColetor(malha, aux);
                                if (valor < -1e9 || valor > 1e9) {
                                    inativo[aux] = 1;
                                    inativoColetor(aux, arqRede, inativo, narq);
                                    inativoAfluente(aux, arqRede, inativo);
                                    if (valor < -1e9) {
                                        aviso(aux);
                                        (*arqRede.vg1dSP).iterRede = 200;
                                    } else {
                                        aviso2(aux);
                                        (*arqRede.vg1dSP).iterRede = 200;
                                    }
                                    (*arqRede.vg1dSP).restartRede = 1;
                                    return 9000.;
                                }
                                Resolv[aux] = 1;
                                resolvGlob += Resolv[aux];
                                ciclomalha++;
                                if (inativo[aux] == 0) {
                                    int qNo = -1;
                                    int kno = 0;
                                    while (qNo < 0) {
                                        for (int kcol = 0; kcol < normaEvol[kno].ncole; kcol++) {
                                            if (normaEvol[kno].cole[kcol] == aux)
                                                qNo = kno;
                                        }
                                        kno++;
                                    }
                                    if (IndNorma[aux] == 0) {
                                        normaEvol[qNo].normaP0 = normaEvol[qNo].normaP1;
                                        normaEvol[qNo].normaP1 = pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                        norma += pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                        IndNorma[aux] = 1;
                                        nnos++;
                                    }
                                    for (int iaflu = 0; iaflu < arqRede.malha[aux].nafluente; iaflu++) {
                                        int indaflu = arqRede.malha[aux].afluente[iaflu];
                                        if ((*arqRede.vg1dSP).fluidoRede == 1)
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                        else if ((*arqRede.vg1dSP).fluidoRede == 0)
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                        else
                                            malha[indaflu].cell(malha[indaflu].ncel).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                        if (arqRede.malha[indaflu].perm == 1) {
                                            pini = malha[indaflu].pGSup;
                                            if (IndNorma[indaflu] == 0) {
                                                IndNorma[indaflu] = 1;
                                                if (arqRede.malha[indaflu].presimposta == 0) {
                                                    norma += pow((malha[aux].cell(0).pres - pini) / pini, 2.);
                                                    nnos++;
                                                }
                                            }
                                            if (arqRede.malha[indaflu].presimposta == 0) {
                                                if ((*arqRede.vg1dSP).iterRede <= 1 && (*arqRede.vg1dSP).relax > 0.5)
                                                    malha[indaflu].pGSup = (0.5) * malha[aux].cell(0).pres + (1. - 0.5) * malha[indaflu].pGSup;
                                                else
                                                    malha[indaflu].pGSup = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres +
                                                                           (1. - (*arqRede.vg1dSP).relax) * malha[indaflu].pGSup;
                                            }
                                        } else {
                                            malha[indaflu].pGSup = malha[aux].cell(0).pres;
                                            malha[indaflu].tGSup = malha[aux].cell(0).temp;
                                            malha[indaflu].tempSup = malha[aux].cell(0).temp;
                                        }
                                    }
                                }
                                int iaflu = 0;
                                int indaflu = arqRede.malha[aux].afluente[iaflu];
                                while (arqRede.malha[indaflu].perm == 0) {
                                    iaflu++;
                                    indaflu = arqRede.malha[aux].afluente[iaflu];
                                }
                                for (int icol2 = 0; icol2 < nderiva - 1; icol2++) {
                                    norma += pow((malha[ordCol[icol2]].cell(0).pres - malha[aux].cell(0).pres) / malha[aux].cell(0).pres, 2.);
                                    nnos++;
                                    if ((*arqRede.vg1dSP).iterRede <= 1 && (*arqRede.vg1dSP).relax > 0.5) {
                                        malha[ordCol[icol2]].cell(0).pres = (0.5) * malha[aux].cell(0).pres + (1. - 0.5) * malha[ordCol[icol2]].cell(0).pres;
                                    } else {
                                        malha[ordCol[icol2]].cell(0).pres = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres +
                                                                              (1. - (*arqRede.vg1dSP).relax) * malha[ordCol[icol2]].cell(0).pres;
                                    }
                                    if ((*arqRede.vg1dSP).iterRede <= 1 && (*arqRede.vg1dSP).relax > 0.5) {
                                        malha[ordCol[icol2]].presE = (0.5) * malha[aux].cell(0).pres + (1. - 0.5) * malha[ordCol[icol2]].cell(0).pres;
                                    } else {
                                        malha[ordCol[icol2]].presE = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres +
                                                                     (1. - (*arqRede.vg1dSP).relax) * malha[ordCol[icol2]].cell(0).pres;
                                    }
                                    malha[ordCol[icol2]].tempE = noConv.tempmist;
                                    malha[ordCol[icol2]].betaE = noConv.betmist;
                                    malha[ordCol[icol2]].titE = noConv.mgasmist / (noConv.mgasmist + noConv.mliqmist);
                                }
                                cout << "convergencia tramo: " << aux << endl;
                            } else {
                                double chutemass = 0.;
                                double area = 0.;
                                int aux = ordCol[icol];
                                if ((*arqRede.vg1dSP).iterRede == 0) {
                                    int indAflu = arqRede.malha[aux].afluente[0];
                                    malha[aux].cell(0).pres = malha[indAflu].pGSup;
                                }
                                (*arqRede.vg1dSP).qualTramo = aux;
                                if (malha[aux].config().ConContEntrada != 2) {
                                    if (malha[aux].config().perm == 1) {
                                        if ((*arqRede.vg1dSP).fluidoRede != 0) {
                                            if (malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto == 2)
                                                malha[aux].cell(0).acsr.injl.FluidoPro.flashCompleto = 0;
                                            malha[aux].cell(0).acsr.tipo = 2;
                                            malha[aux].cell(0).acsr.injl.bet = noConv.betmist;
                                            malha[aux].cell(0).acsr.injl.temp = noConv.tempmist;
                                            if ((*arqRede.vg1dSP).iterRede == 0 || malha[aux].cell(0).acsr.injl.QLiq > 0) {
                                                malha[aux].cell(0).acsr.injl.FluidoPro.RGO = noConv.RGOmist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.Deng = noConv.dengmist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.Denag = noConv.denagmist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.yco2 = noConv.yco2mist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.API = noConv.apimist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.BSW = noConv.bswmist;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.TempL = noConv.tL;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.TempH = noConv.tH;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.LVisL = noConv.LVisL;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.LVisH = noConv.LVisH;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.RenovaFluido();
                                                malha[aux].cell(0).acsr.injl.fluidocol.TR = noConv.TRmist;
                                            } else {
                                                malha[aux].cell(0).acsr.injl.FluidoPro = malha[aux].fluiRevRede;
                                                malha[aux].cell(0).acsr.injl.FluidoPro.RenovaFluido();
                                            }

                                            malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injl.FluidoPro;
                                        } else {
                                            if (malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto == 2)
                                                malha[aux].cell(0).acsr.injg.FluidoPro.flashCompleto = 0;
                                            malha[aux].cell(0).acsr.tipo = 1;
                                            malha[aux].cell(0).acsr.injg.seco = 0;
                                            malha[aux].cell(0).acsr.injg.temp = noConv.tempmist;
                                            if ((*arqRede.vg1dSP).iterRede == 0 || malha[aux].cell(0).acsr.injg.QGas > 0) {
                                                malha[aux].cell(0).acsr.injg.FluidoPro.RGO = noConv.RGOmist;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.API = noConv.apimist;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.BSW = noConv.bswmist;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.TempL = noConv.tL;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.TempH = noConv.tH;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.LVisL = noConv.LVisL;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.LVisH = noConv.LVisH;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.Deng = noConv.dengmist;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.yco2 = noConv.yco2mist;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                                                malha[aux].cell(0).acsr.injg.fluidocol.TR = noConv.TRmist;
                                            } else {
                                                malha[aux].cell(0).acsr.injg.FluidoPro = malha[aux].fluiRevRede;
                                                malha[aux].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                                            }

                                            malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injg.FluidoPro;
                                        }
                                    }
                                } else {
                                    if (malha[aux].config().perm == 1) {
                                        if (malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto == 2)
                                            malha[aux].cell(0).acsr.injm.FluidoPro.flashCompleto = 0;
                                        malha[aux].cell(0).acsr.tipo = 10;

                                        malha[aux].cell(0).acsr.injm.temp = noConv.tempmist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.RGO = noConv.RGOmist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.Deng = noConv.dengmist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.Denag = noConv.denagmist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.yco2 = noConv.yco2mist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.API = noConv.apimist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.BSW = noConv.bswmist;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.TempL = noConv.tL;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.TempH = noConv.tH;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.LVisL = noConv.LVisL;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.LVisH = noConv.LVisH;
                                        malha[aux].cell(0).acsr.injm.FluidoPro.RenovaFluido();
                                        malha[aux].cell(0).acsr.injm.fluidocol.TR = noConv.TRmist;

                                        malha[aux].cell(0).flui = malha[aux].cell(0).acsr.injm.FluidoPro;
                                    }
                                }

                                if (malha[aux].config().perm == 1) {
                                    if (malha[aux].config().ConContEntrada != 2) {
                                        for (int icol2 = icol; icol2 < nderiva; icol2++)
                                            area += malha[ordCol[icol2]].cell(0).duto.area;
                                        int indicadorVP = aux;
                                        indicadorVP = -1;
                                        if (indicadorVP == -1) {
                                            if ((*arqRede.vg1dSP).iterRede == 0) {
                                                if ((*arqRede.vg1dSP).fluidoRede != 0)
                                                    chutemass = (noConv.qlmistStd + noConv.qlmistStdNeg) * malha[aux].cell(0).duto.area / area;
                                                else {
                                                    chutemass = (noConv.qgstdTot * 86400. +
                                                                 noConv.qgstdTotNeg * 86400) *
                                                                malha[aux].cell(0).duto.area / area;
                                                    if (fabs(chutemass) > 1e-15)
                                                        malha[aux].cell(0).acsr.injg.razCompGas = (noConv.qcmistStd + noConv.qcmistStdNeg) /
                                                                                                    (noConv.qgstdTot * 86400. + noConv.qgstdTotNeg * 86400);
                                                    else
                                                        malha[aux].cell(0).acsr.injg.razCompGas = 0.;
                                                }
                                            } else {
                                                if ((*arqRede.vg1dSP).fluidoRede == 1)
                                                    chutemass = malha[aux].cell(0).acsr.injl.QLiq;
                                                else if ((*arqRede.vg1dSP).fluidoRede == 0)
                                                    chutemass = malha[aux].cell(0).acsr.injg.QGas;
                                            }
                                        }
                                        if (arqRede.malha[aux].perm == 1) {
                                            malha[aux].modoPerm = 1;
                                            if (malha[aux].config().ConContEntrada != 2) {
                                                if (chutemass > 0.) {
                                                    if (aux == 32) {
                                                    }
                                                    if (malha[aux].config().chokep.abertura[0] > 0.6)
                                                        valor = malha[aux].buscaProdPresPresPerm(chutemass, 1 * noConv.qgstdTot * 86400);
                                                    else
                                                        valor = malha[aux].buscaProdPresPresPerm2(chutemass, 1 * noConv.qgstdTot * 86400);
                                                } else {
                                                    malha[aux].cell(0).acsr.injg.razCompGas = malha[aux].config().razCompGasReves;
                                                    valor = malha[aux].buscaProdPresPresPermRev(chutemass, -1 * noConv.qgstdTot * 86400);
                                                }
                                                if ((*arqRede.vg1dSP).fluidoRede != 0) {
                                                    if (malha[aux].cell(0).acsr.injl.QLiq > (noConv.qlmistStd + noConv.qlmistStdNeg)) {
                                                        malha[aux].cell(0).acsr.injl.QLiq = chutemass;
                                                        cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo: " << aux << endl;
                                                        arqRede.malha[aux].principal = 1;
                                                    } else
                                                        arqRede.malha[aux].principal = 0;
                                                } else {
                                                    if (malha[aux].cell(0).acsr.injg.QGas > (noConv.qgstdTot * 86400 +
                                                                                               noConv.qgstdTotNeg * 86400)) {
                                                        malha[aux].cell(0).acsr.injg.QGas = chutemass;
                                                        cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo: " << aux << endl;
                                                        arqRede.malha[aux].principal = 1;
                                                    } else
                                                        arqRede.malha[aux].principal = 0;
                                                }
                                            } else
                                                valor = 0;
                                        }

                                        if (valor < -1e9 || valor > 1e9) {
                                            inativo[aux] = 1;
                                            inativoColetor(aux, arqRede, inativo, narq);
                                            inativoAfluente(aux, arqRede, inativo);
                                            if (valor < -1e9) {
                                                aviso(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            } else {
                                                aviso2(aux);
                                                (*arqRede.vg1dSP).iterRede = 200;
                                            }
                                            (*arqRede.vg1dSP).restartRede = 1;
                                            return 9000.;
                                        }
                                        Resolv[aux] = 1;
                                        resolvGlob += Resolv[aux];
                                        ciclomalha++;
                                        if ((*arqRede.vg1dSP).fluidoRede == 1) {
                                            if (arqRede.malha[aux].perm == 1 && (malha[aux].cell(0).acsr.injl.QLiq) > 0)
                                                noConv.qlmistStd -= malha[aux].cell(0).acsr.injl.QLiq;
                                            if (fabs(malha[aux].cell(0).acsr.injl.QLiq + 2121212121) < 1e-15)
                                                cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo  : " << aux << endl;
                                        } else if ((*arqRede.vg1dSP).fluidoRede == 0) {
                                            if (arqRede.malha[aux].perm == 1 && (malha[aux].cell(0).acsr.injg.QGas) > 0)
                                                noConv.qgstdTot -= malha[aux].cell(0).acsr.injg.QGas / 86400;
                                            if (fabs(malha[aux].cell(0).acsr.injg.QGas + 2121212121) < 1e-15)
                                                cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo : " << aux << endl;
                                        } else if ((*arqRede.vg1dSP).fluidoRede == 2) {
                                            if (arqRede.malha[aux].perm == 1 &&
                                                (malha[aux].cell(0).acsr.injl.QLiq) > 0) {
                                                double rlcA = malha[aux].cell(0).acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
                                                double vazmassc = rlcA * malha[aux].cell(0).acsr.injl.QLiq * malha[aux].cell(0).acsr.injl.bet / 86400;
                                                double massic = malha[aux].cell(0).acsr.injl.QLiq *
                                                                (1. - malha[aux].cell(0).acsr.injl.bet) / 86400;
                                                double Rhogs = malha[aux].cell(0).acsr.injl.FluidoPro.Deng * 1.225; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
                                                double Rhols = (1000 * 141.5 / (131.5 + malha[aux].cell(0).acsr.injl.FluidoPro.API)) * (1 - malha[aux].cell(0).acsr.injl.FluidoPro.BSW) + 1000. * malha[aux].cell(0).acsr.injl.FluidoPro.Denag * malha[aux].cell(0).acsr.injl.FluidoPro.BSW;
                                                double multiplicador = (Rhols + malha[aux].cell(0).acsr.injl.FluidoPro.RGO * Rhogs *
                                                                                    (1 - malha[aux].cell(0).acsr.injl.FluidoPro.BSW));
                                                massic *= multiplicador;
                                                double fracmasshidra = malha[aux].cell(0).acsr.injl.FluidoPro.FracMassHidra(malha[aux].cell(0).pres, malha[aux].cell(0).temp);
                                                noConv.mliqmist -= ((1. - fracmasshidra) * massic + vazmassc);
                                                noConv.mgasmist -= fracmasshidra * massic;
                                            }
                                            if (fabs(malha[aux].cell(0).acsr.injl.QLiq + 2121212121) < 1e-15)
                                                cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo : " << aux << endl;
                                        }
                                    } else {
                                        if (arqRede.malha[aux].perm == 1) {
                                            int iaflu = 0;
                                            int indaflu = arqRede.malha[aux].afluente[iaflu];
                                            while (arqRede.malha[indaflu].perm == 0) {
                                                iaflu++;
                                                indaflu = arqRede.malha[aux].afluente[iaflu];
                                            }
                                            valor = malha[aux].buscaProdPfundoPerm3(malha[indaflu].pGSup);
                                        } else
                                            valor = 0.;
                                        if (arqRede.malha[aux].perm == 1) {
                                            if (valor < -1e9 || valor > 1e9) {
                                                inativo[aux] = 1;
                                                inativoColetor(aux, arqRede, inativo, narq);
                                                inativoAfluente(aux, arqRede, inativo);
                                                if (valor < -1e9) {
                                                    aviso(aux);
                                                    (*arqRede.vg1dSP).iterRede = 200;
                                                } else {
                                                    aviso2(aux);
                                                    (*arqRede.vg1dSP).iterRede = 200;
                                                }
                                                (*arqRede.vg1dSP).restartRede = 1;
                                                return 9000.;
                                            }
                                            Resolv[aux] = 1;
                                            resolvGlob += Resolv[aux];
                                            ciclomalha++;
                                            double rlcA = malha[aux].cell(0).acsr.injm.fluidocol.MasEspFlu(1.001, 15.);
                                            double Rhogs = malha[aux].cell(0).acsr.injm.FluidoPro.Deng * 1.225; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
                                            double Rhols = (1000 * 141.5 / (131.5 + malha[aux].cell(0).acsr.injm.FluidoPro.API)) * (1 - malha[aux].cell(0).acsr.injm.FluidoPro.BSW) + 1000. * malha[aux].cell(0).acsr.injm.FluidoPro.Denag *
                                                                                                                                                                                              malha[aux].cell(0).acsr.injm.FluidoPro.BSW;
                                            double multiplicador = (Rhols + malha[aux].cell(0).acsr.injm.FluidoPro.RGO * Rhogs *
                                                                                (1 - malha[aux].cell(0).acsr.injm.FluidoPro.BSW));
                                            double fracmasshidra = malha[aux].cell(0).acsr.injm.FluidoPro.FracMassHidra(malha[aux].cell(0).pres,
                                                                                                                          malha[aux].cell(0).temp);
                                            if ((*arqRede.vg1dSP).fluidoRede == 1 && malha[aux].cell(0).acsr.injm.MassP > 0.) {
                                                if (arqRede.malha[aux].perm == 1) {
                                                    noConv.qlmistStd -= (malha[aux].cell(0).acsr.injm.MassP * 86400. /
                                                                             ((1. - fracmasshidra) * multiplicador) +
                                                                         malha[aux].cell(0).acsr.injm.MassC * 86400. / rlcA);
                                                }
                                            } else if ((*arqRede.vg1dSP).fluidoRede == 0 && malha[aux].cell(0).acsr.injm.MassG > 0.) {
                                                double tit = malha[aux].cell(0).acsr.injm.FluidoPro.FracMassHidra(1., 20.);
                                                noConv.qgstdTot -= malha[aux].cell(0).acsr.injm.MassG * tit /
                                                                   (malha[aux].cell(0).acsr.injm.FluidoPro.Deng * 1.225 * fracmasshidra);
                                            } else if ((malha[aux].cell(0).acsr.injm.MassC + malha[aux].cell(0).acsr.injm.MassP +
                                                        malha[aux].cell(0).acsr.injm.MassG) > 0.) {
                                                noConv.mliqmist -= (malha[aux].cell(0).acsr.injm.MassC + malha[aux].cell(0).acsr.injm.MassP);
                                                noConv.mgasmist -= malha[aux].cell(0).acsr.injm.MassG;
                                            }
                                            if (fabs(malha[aux].cell(0).acsr.injm.MassC +
                                                     malha[aux].cell(0).acsr.injm.MassP +
                                                     malha[aux].cell(0).acsr.injm.MassG + 2121212121) < 1e-15)
                                                cout << "Tramo de bifurcacao com problemas-aceitando mais vazao do que a soma das correntes-revisar condiçao de contorno no tramo : " << aux << endl;
                                        }
                                    }
                                    if (arqRede.malha[aux].perm == 1)
                                        cout << "convergencia tramo: " << aux << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (fabs(testaNMalha) > 0) {
        (*arqRede.vg1dSP).erroRede = fabs(sqrt(norma) / nnos - (*arqRede.vg1dSP).norma0) / (*arqRede.vg1dSP).norma0;
        (*arqRede.vg1dSP).norma0 = sqrt(norma) / nnos;
        return sqrt(norma) / nnos;
    } else {
        return 0;
    }
}

void alteraModoFluidoCompBlack(SProd *malha, int narq, Vcr<int> &calclat0, Vcr<int> &tipoFluido0) {
    for (int j = 0; j < narq; j++) {
        malha[j].config().flashCompleto = 0;
        calclat0[j] = malha[j].CalcLat;
        malha[j].CalcLat = 0;
        tipoFluido0[j] = malha[j].config().tipoFluido;
        (*malha[j].vg1dSP).blackOilTemp = 1;
        malha[j].config().tipoFluido = 0;
        malha[j].fluiRevRede.flashCompleto = 0;
        malha[j].fluiRevRede.modoBlackTemp = 0;
        for (int i = 0; i <= malha[j].ncel; i++) {
            malha[j].cell(i).flui.flashCompleto = 0;
            malha[j].cell(i).flui.modoBlackTemp = 0;
            if (malha[j].cell(i).acsr.tipo == 1) {
                malha[j].cell(i).acsr.injg.FluidoPro.flashCompleto = 0;
                malha[j].cell(i).acsr.injg.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 2) {
                malha[j].cell(i).acsr.injl.FluidoPro.flashCompleto = 0;
                malha[j].cell(i).acsr.injl.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 3) {
                malha[j].cell(i).acsr.ipr.FluidoPro.flashCompleto = 0;
                malha[j].cell(i).acsr.ipr.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 15) {
                malha[j].cell(i).acsr.radialPoro.flup.flashCompleto = 0;
                malha[j].cell(i).acsr.radialPoro.flup.modoBlackTemp = 0;
                int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                for (int k = 0; k < ncelRad; k++) {
                    malha[j].cell(i).acsr.radialPoro.celula[k].flup.flashCompleto = 0;
                    malha[j].cell(i).acsr.radialPoro.celula[k].flup.modoBlackTemp = 0;
                }
            } else if (malha[j].cell(i).acsr.tipo == 16) {
                malha[j].cell(i).acsr.poroso2D.dados.flup.flashCompleto = 0;
                malha[j].cell(i).acsr.poroso2D.dados.flup.modoBlackTemp = 0;
                int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                for (int k = 0; k < ncelRad; k++) {
                    malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.flashCompleto = 0;
                    malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.modoBlackTemp = 0;
                }
                int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                for (int k = 0; k < ncelEle; k++) {
                    malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.flashCompleto = 0;
                    malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.modoBlackTemp = 0;
                }
            } else if (malha[j].cell(i).acsr.tipo == 10) {
                malha[j].cell(i).acsr.injm.FluidoPro.flashCompleto = 0;
                malha[j].cell(i).acsr.injm.FluidoPro.modoBlackTemp = 0;
            }
        }
    }
}

void alteraModoFluidoBlackComp(SProd *malha, int narq, Vcr<int> calclat0, Vcr<int> tipoFluido0) {
    for (int j = 0; j < narq; j++) {
        malha[j].config().flashCompleto = 2;
        malha[j].CalcLat = calclat0[j];
        malha[j].config().tipoFluido = tipoFluido0[j];
        (*malha[j].vg1dSP).blackOilTemp = 0;
        malha[j].fluiRevRede.flashCompleto = 2;
        malha[j].fluiRevRede.modoBlackTemp = 0;
        for (int i = 0; i <= malha[j].ncel; i++) {
            malha[j].cell(i).flui.flashCompleto = 2;
            malha[j].cell(i).flui.modoBlackTemp = 0;
            if (malha[j].cell(i).acsr.tipo == 1) {
                malha[j].cell(i).acsr.injg.FluidoPro.flashCompleto = 2;
                malha[j].cell(i).acsr.injg.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 2) {
                malha[j].cell(i).acsr.injl.FluidoPro.flashCompleto = 2;
                malha[j].cell(i).acsr.injl.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 3) {
                malha[j].cell(i).acsr.ipr.FluidoPro.flashCompleto = 2;
                malha[j].cell(i).acsr.ipr.FluidoPro.modoBlackTemp = 0;
            } else if (malha[j].cell(i).acsr.tipo == 15) {
                malha[j].cell(i).acsr.radialPoro.flup.flashCompleto = 2;
                malha[j].cell(i).acsr.radialPoro.flup.modoBlackTemp = 0;
                int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                for (int k = 0; k < ncelRad; k++) {
                    malha[j].cell(i).acsr.radialPoro.celula[k].flup.flashCompleto = 2;
                    malha[j].cell(i).acsr.radialPoro.celula[k].flup.modoBlackTemp = 0;
                }
            } else if (malha[j].cell(i).acsr.tipo == 16) {
                malha[j].cell(i).acsr.poroso2D.dados.flup.flashCompleto = 2;
                malha[j].cell(i).acsr.poroso2D.dados.flup.modoBlackTemp = 0;
                int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                for (int k = 0; k < ncelRad; k++) {
                    malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.flashCompleto = 2;
                    malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.modoBlackTemp = 0;
                }
                int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                for (int k = 0; k < ncelEle; k++) {
                    malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.flashCompleto = 2;
                    malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.modoBlackTemp = 0;
                }
            } else if (malha[j].cell(i).acsr.tipo == 10) {
                malha[j].cell(i).acsr.injm.FluidoPro.flashCompleto = 2;
                malha[j].cell(i).acsr.injm.FluidoPro.modoBlackTemp = 0;
            }
        }
        malha[j].cell(0).flui.npseudo = malha[j].cell(1).flui.npseudo;
        malha[j].cell(malha[j].ncel).flui.npseudo = malha[j].cell(malha[j].ncel - 1).flui.npseudo;
        malha[j].buscaIni = 0;
    }
}

int convergeRede(SProd *malha, Rede &arqRede, int narq, Vcr<int> &inativo, int &indativo, vector<noRede> &normaEvol, vector<tramoPart> &bloq) {
    double norma = 1000.;
    double limConvergeLocal = arqRede.limConverge;
    int convergeaux = 0;
    if (malha[0].config().flashCompleto == 0 && arqRede.fluidoRede == 0) {
        limConvergeLocal *= 2.1;
        convergeaux = 1;
    }
    while (norma > limConvergeLocal && indativo < narq && (*arqRede.vg1dSP).restartRede == 0 && (*arqRede.vg1dSP).iterRede < 200) {
        if ((*arqRede.vg1dSP).erroRede < 0.01) {
            cout << "convergencia da Rede se tornando lenta ou perdida, "
                    "um motivo provável para isto sao os limites de valores de tabela, caso a opçao PVTSim puro esteja sendo usada"
                 << endl;
            (*arqRede.vg1dSP).relax *= 0.5;
            if ((*arqRede.vg1dSP).relax < 0.1)
                (*arqRede.vg1dSP).relax = 0.1;
        }
        indativo = 0;
        for (int i = 0; i < narq; i++)
            indativo += inativo[i];
        if (malha[0].config().flashCompleto != 2) {
            (*arqRede.vg1dSP).verificaTramo = 0;
            norma = cicloRede(malha, arqRede, inativo, indativo, normaEvol, bloq);
            if (convergeaux == 1 && (*arqRede.vg1dSP).iterRede > 20)
                norma = limConvergeLocal / 2.;
        } else {
            (*arqRede.vg1dSP).verificaTramo = 1;
            norma = cicloRedeComp(malha, arqRede, inativo, indativo, normaEvol, bloq);
        }
        cout << "iteracao= " << (*arqRede.vg1dSP).iterRede << " norma= " << norma << " tramos inativos= " << indativo << endl;
        (*arqRede.vg1dSP).iterRede++;
    }
    if ((*arqRede.vg1dSP).iterRede >= 200) {
        cout << "Problema na convergencia da rede, abortando calculo " << endl;
        return -1;
    } else
        return 1;
}

double permanenteSimples(SProd &sistem1, double inichute = -1.) {
    // metodo que gerencia qual forma de marcha e busca de chutes iniciais no processo de convergência do calculo permanente para um tramo simples
    // o solver permanente é baseado em um método marchante, o sistema de convergência é obtido tomando como função objetivo ou a diferenca
    //  entre a pressao do separador e a pressao na celula final obtida a partir da marcha, neste caso, a convergencia e obtida
    // quando esta diferenca é zero, este criterio e usado quando o choke no fim da tubulacao nao esta oferecendo resistencia ao escoamento.
    // Caso a valvula choke do final da tubulacao esteja oferecendo resistencia ao escoamento, "choke ativo", a funcao objetivo e alterado,
    // torna-se a diferenca entre a vazao massica na ultima celula do sistema e a vazao massica que se esta passando pelo choke, obtida
    //  a jusante: a pressao do separador, definida como condicao de contorno
    double ciclo = 1.1e9; // esta variavel sera utilizada como indicacao de que a solucao permanente foi obtida ou ocorreu
    // algum problema de convergencia ou de marcha
    if (sistem1.config().pocinjec == 0) { // nao é um sistema de injecao

        if (sistem1.config().ConContEntrada == 0) { // condicao de contorno na entrada é uma fonte, de gas,
            // de liquido, IPR, fonte de massa, meio poroso radial ou meio poroso 2D
            int reverso = 0;                        // variavel que indica se o escoamento segue o senido montante jusante pre-definido ou se encontra em sentido reverso
            if (sistem1.cell(0).acsr.tipo == 1) { // caso a fonte no inicio do tramop seja uma fonte de gas e a sua vazao seja negativa
                if (sistem1.cell(0).acsr.injg.QGas < 0)
                    reverso = 1;
            } else if (sistem1.cell(0).acsr.tipo == 2) { // caso a fonte no inicio do tramop seja uma fonte de liquido e a sua vazao seja negativa
                if (sistem1.cell(0).acsr.injl.QLiq < 0)
                    reverso = 1;
            } else if (sistem1.cell(0).acsr.tipo == 10) { // caso a fonte no inicio do tramop seja uma fonte de massa e a sua vazao seja negativa
                if ((sistem1.cell(0).acsr.injm.MassC + sistem1.cell(0).acsr.injm.MassG + sistem1.cell(0).acsr.injm.MassP) < 0)
                    reverso = 1;
            }
            // observacao, este teste de escoamento reverso não é feito para IPR, radial poroso ou 2d poroso, pois não se sabe se
            // a fonte de contoirno nestes casos será negativo ou positivo sem se ter uma ideia de qual é a pressao de fundo que o
            // sistema tera de lidar
            if (reverso == 0) { // caso o escoamento siga o sentido montante-jusante pre definido
                if (sistem1.cell(0).acsr.tipo == 1 && sistem1.cell(0).acsr.injg.seco == 1)
                    sistem1.cell(0).flui.RGO =
                        sistem1.cell(0).acsr.injg.FluidoPro.RGO = 1e6 + (*sistem1.vg1dSP).localtiny; // para o caso de fonte
                // de gas seco, neste caso, não existe condensado e, pelas limitações do modelo de propriedade black oil, é necessário admitir
                // um valor de RGO nop sistema de produção, neste caso, impoe-se um valor muito alto de RGO, par que não cause impacto
                // na simulacao de gas seco
                if (sistem1.config().chokep.abertura[0] > 0.6 && sistem1.config().chokep.curvaDinamic == 0) { // choke muito aberto, nao e considerado ativo
                    // a opcao curvaDinamic é uma opcao muito especial, quando é 1 se deseja emular um comportamento dinamico do atuador da valvula,
                    // que respondera aa pressao na linha, aa pressao no atuador e a uma rigidez de mola, alem de considerar a massa do atuador,
                    // neste caso, definir quando a valvula estarar ativa seria uma funcao das condicoes de escoamento
                    if (sistem1.config().lingas == 1 && sistem1.config().gasinj.chuteVaz == 0 &&
                        sistem1.gasCell(0).tipoCC == 0) // caso em que se tem uma linha de servico associada
                        // e a condicao de contorno e pressao de injecao de gas
                        // neste caso, necessita-se de um chute inicial de vazao de injecao, a solucao de pressao de injecao e
                        // mais dificil de se obter no algoritmo utilizado no Marlim 3, por isto, sempre faz-se uma primeira rodada
                        // com uma estimativa inicial obtida de uma solucao com vazao de injecao
                        // caso este chute nao seja fornecido no json
                        // faz-se uma estimativa:
                        sistem1.config().gasinj.vazgas[0] = 150000 * sistem1.gasCell(0).duto.area / (*sistem1.vg1dSP).arearef;
                    // com esta estimativa, faz-se uma primeira solucao
                    // permanente com condicao de contorno vazao de injecao
                    if (sistem1.config().lingas == 1 && sistem1.gasCell(0).tipoCC == 0) {
                        int konta = 0;
                        double multVazGas;
                        sistem1.gasCell(0).tipoCC = 1;
                        sistem1.buscaProdPfundoPerm(inichute);
                        double pref = sistem1.cell(0).pres;
                        double testaPres1 = sistem1.gasCell(0).pres;
                        sistem1.config().gasinj.vazgas[0] *= 1.05;
                        if (inichute < 0)
                            sistem1.buscaProdPfundoPerm(sistem1.cell(0).pres);
                        else
                            sistem1.buscaProdPfundoPerm(inichute);
                        double testaPres2 = sistem1.gasCell(0).pres;
                        if (testaPres1 < testaPres2) {
                            if (sistem1.presiniG > testaPres1)
                                multVazGas = 1.05;
                            else
                                multVazGas = 0.95;
                            sistem1.config().gasinj.vazgas[0] /= 1.05;
                        } else {
                            if (sistem1.presiniG > testaPres1)
                                multVazGas = 0.95;
                            else
                                multVazGas = 1.05;
                            sistem1.config().gasinj.vazgas[0] /= 1.05;
                        }
                        while (fabs(ciclo) > 0.9e9 && konta < 10) {
                            sistem1.gasCell(0).tipoCC = 1; // mudando para condicao vazao de injecao na linha de servico,
                            // para se ter uma primeira estimativa da pressao de fundo
                            // solucao permanente com vazao de injecao
                            if (konta > 0) {
                                ciclo = sistem1.buscaProdPfundoPerm(inichute, konta); // solucao de estimativa com vazao de injecao
                                pref = sistem1.cell(0).pres;                        // utilizando a pressao de fundo obtida da
                                                               // estimativa com vazao de injecao de gas-lift
                            }
                            sistem1.gasCell(0).tipoCC = 0; // voltando para a condicao pressao de injecao
                            // com o valor da pressao de fundo desta solucao, busca-se a solucao
                            // com condicao de contorno original, pressao de injecao
                            ciclo = sistem1.buscaProdPfundoPerm(pref, konta); // nova tentativa de convergencia,agora com um valor de chute inicial
                            // se nao for possivel a convergencia, faz-se uma nova estimativa de vazao
                            // de injecao e tenta-se novamente a solucao com pressao de injecao,
                            // aumenta-se o valor da estimativa de vazao de injecao
                            if (fabs(ciclo) > 0.9e9) {
                                // caso o valor retornado em ciclo seja maior que 0.9e9, isto é um indicativo de que a
                                // convergencia nao foi atingida por se ter pouco injecao de gas, neste caso
                                // retorna-se aa abordagem descrita anteriormente, mas aumenta-se o valor da estimativa de
                                // quanto deve ser a vazao de gas
                                sistem1.config().gasinj.vazgas[0] *= multVazGas;
                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                            }
                            konta++;
                            if (konta > 40) { // o ciclo e repetidoo no maximo 10 vezes, se der errado, termina
                                // a busca permanente
                                ciclo = 1.1e10;
                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                           "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################",
                                           "", "");
                            }
                        }
                    } else
                        ciclo = sistem1.buscaProdPfundoPerm(inichute); // solucao para o caso em nao se tem linha de servico,
                    // ou tendo linha de servico, esta tem uma condicao de contorno=vazao de injecao
                } else { // caso em que o choke esta restrito o suficiente para ser considerado ativo
                    if (sistem1.config().lingas == 1 && sistem1.config().gasinj.chuteVaz == 0 &&
                        sistem1.gasCell(0).tipoCC == 0) // caso em que se tem uma linha de servico associada
                        // e a condicao de contorno e pressao de injecao de gas
                        // neste caso, necessita-se de um chute inicial de vazao de injecao, a solucao de pressao de injecao e
                        // mais dificil de se obter no algoritmo utilizado no Marlim 3, por isto, sempre faz-se uma primeira rodada
                        // com uma estimativa inicial obtida de uma solucao com vazao de injecao
                        // caso este chute nao seja fornecido no json
                        // faz-se uma estimativa:
                        sistem1.config().gasinj.vazgas[0] = 150000 * sistem1.gasCell(0).duto.area / (*sistem1.vg1dSP).arearef;
                    // com esta estimativa, faz-se uma primeira solucao
                    // permanente com condicao de contorno vazao de injecao
                    if (sistem1.config().lingas == 1 && sistem1.gasCell(0).tipoCC == 0) {
                        int konta = 0;
                        double multVazGas;
                        sistem1.gasCell(0).tipoCC = 1;
                        sistem1.buscaProdPfundoPerm2(inichute);
                        double testaPres1 = sistem1.gasCell(0).pres;
                        sistem1.config().gasinj.vazgas[0] *= 1.05;
                        if (inichute < 0)
                            sistem1.buscaProdPfundoPerm2(sistem1.cell(0).pres);
                        else
                            sistem1.buscaProdPfundoPerm2(inichute);
                        double testaPres2 = sistem1.gasCell(0).pres;
                        if (testaPres1 < testaPres2) {
                            if (sistem1.presiniG > testaPres1)
                                multVazGas = 1.05;
                            else
                                multVazGas = 0.95;
                            sistem1.config().gasinj.vazgas[0] /= 1.05;
                        } else {
                            if (sistem1.presiniG > testaPres1)
                                multVazGas = 0.95;
                            else
                                multVazGas = 1.05;
                            sistem1.config().gasinj.vazgas[0] /= 1.05;
                        }
                        while (fabs(ciclo) > 0.9e9 && konta < 10) {
                            sistem1.gasCell(0).tipoCC = 1;
                            // solucao permanente com vazao de injecao
                            if (konta > 0) {
                                ciclo = sistem1.buscaProdPfundoPerm2(inichute, konta); // solucao de estimativa com vazao de injecao
                            }
                            sistem1.gasCell(0).tipoCC = 0; // voltando para a condicao pressao de injecao
                                                           // com o valor da pressao de fundo desta solucao, busca-se a solucao
                                                           // com condicao de contorno original, pressao de injecao
                            double pref = sistem1.cell(0).pres; // utilizando a pressao de fundo obtida da
                            // estimativa com vazao de injecao de gas-lift
                            ciclo = sistem1.buscaProdPfundoPerm2(pref, konta);
                            // se nao for possivel a convergencia, faz-se uma nova estimativa de vazao
                            // de injecao e tenta-se novamente a solucao com pressao de injecao,
                            // aumenta-se o valor da estimativa de vazao de injecao
                            if (fabs(ciclo) > 0.9e9) {
                                // caso o valor retornado em ciclo seja maior que 0.9e9, isto é um indicativo de que a
                                // convergencia nao foi atingida por se ter pouco injecao de gas, neste caso
                                // retorna-se aa abordagem descrita anteriormente, mas aumenta-se o valor da estimativa de
                                // quanto deve ser a vazao de gas
                                sistem1.config().gasinj.vazgas[0] *= multVazGas;
                                cout << "#################NOVA TENTATIVA PARA GL COM CONDICAO DE PRESSAO########################" << endl;
                            }
                            konta++;
                            if (konta > 40) { // o ciclo e reptidoo no maximo 10 vezes, se der errado, termina
                                // a busca permanente
                                ciclo = 1.1e10;
                                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                                logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                           "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################",
                                           "", "");
                            }
                        }
                    } else
                        ciclo = sistem1.buscaProdPfundoPerm2(inichute); // solucao para o caso em nao se tem linha de servico,
                    // ou tendo linha de servico, esta tem uma condicao de contorno=vazao de injecao
                }
            } else
                ciclo = sistem1.buscaProdPfundoPermRev(inichute);
        } else if (sistem1.config().ConContEntrada == 1) { // caso em que a condicao de contorno e pressao na entrada do duto
            // o resultado que se deseja obter, oua  "informaco faltante" é
            // a vazao massica na entrada do duto
            double chutemass = 0.;
            double tit = sistem1.config().CCPres.tit[0]; // o titulo na entrada, junto com a pressão e a razao beta são
            ProFlu fluCP;
            if (sistem1.config().CCPres.indFluido == -1)
                fluCP = sistem1.cell(0).flui;
            else
                fluCP = sistem1.config().flup[sistem1.config().CCPres.indFluido];
            // valçores de condicao de contorno
            if (tit < (1. - (*sistem1.vg1dSP).localtiny) && sistem1.config().tipoFluido == 0) { // caso em que se injeta uma corrente bifasica,
                // preparacao das condicoes de entrada
                // RGO de entrada para a solucao permanente neste caso, observar
                // que nao se indica a RGO de entrada e sim o título de entrada

                sistem1.tempE = sistem1.config().CCPres.temperatura[0]; // temperatura de entrada
                sistem1.presE = sistem1.config().CCPres.pres[0];        // pressão de entrada
                double rcst = sistem1.cell(0).fluicol.MasEspFlu(1.01, 15);
                double rcis = sistem1.cell(0).fluicol.MasEspFlu(sistem1.presE, sistem1.tempE);
                if ((1 - sistem1.betaE) > 1e-10)
                    sistem1.betaE = sistem1.config().CCPres.bet[0] * rcst / rcis; // beta na condicao in-situ
                else
                    sistem1.betaE = 1.;
                double rgoentrada = 0.;
                if (sistem1.config().flashCompleto != 2) { // obtencao da RGO  de entrada a partir do titulo, para casos que nao
                    // envolva modelo composicional, no caso composicional, isto é obtido internamente pelo proprio modelo composicional
                    // propriedades standard
                    double rgST = fluCP.Deng * 1.225;
                    double roST = 141.5 * 1000. / (131.5 + fluCP.API);
                    // propriedades in-situ
                    double rl = fluCP.MasEspoleo(sistem1.presE, sistem1.tempE);
                    double rw = fluCP.Denag * 1000.;
                    double RSent = fluCP.RS(sistem1.presE, sistem1.tempE) * 6.29 / 35.31467;
                    double bsw = fluCP.BSW;
                    double bo = fluCP.BOFunc(sistem1.presE, sistem1.tempE);
                    double ba = fluCP.BAFunc(sistem1.presE, sistem1.tempE);
                    // fracao de agua in-situ
                    double fw = bsw * ba / (bo + ba * bsw - bsw * bo);

                    // volume de liquidos in-situ para uma unidade de volume de oleo em condicao standard (oleo morto)
                    double Ost = (1 - sistem1.betaE) * bo * (1 - fw);
                    double Wst = fw * (1 - sistem1.betaE);
                    double Cst = sistem1.betaE;
                    sistem1.titE = sistem1.config().CCPres.tit[0];
                    if (sistem1.titE < 0) {
                        double titH = fluCP.FracMassHidra(sistem1.presE, sistem1.tempE);
                        double rc = sistem1.cell(0).fluicol.MasEspFlu(sistem1.presE, sistem1.tempE);
                        double rlMix = sistem1.betaE * rc + (1. - sistem1.betaE) * rl;
                        double val1 = ((1. - sistem1.betaE) * rl * titH / (1. - titH));
                        sistem1.titE = val1 / (rlMix + val1);
                    }

                    if ((1 - sistem1.betaE) > 1e-10) {
                        rgoentrada = RSent * rgST * Ost + sistem1.titE * roST * Ost + sistem1.titE * rw * Wst +
                                     sistem1.titE * rcst * Cst;
                        rgoentrada /= ((1 - sistem1.titE) * rgST * Ost);
                    } else
                        rgoentrada = sistem1.cell(0).flui.RGO;
                }

                if (rgoentrada > 1600.)
                    rgoentrada = 1600.; // limite de rgo para o modelo black-oil

                // no caso da solução permanente com condicao de pressao, a vazao massica no inicio da tubulacao é a informacao a ser obtida
                // para tanto, no solver permanente, uma fonte e incluida na celula inicial, como um acessorio, esta fonte e um artificialismo para auxiliar
                // o processo de convergencia, nela, será armazenada nesta fonte a vazao que se deseja obter e armazenada, aqui,
                // se utiliza uma fonte de líquido, pois o titulo de entrada não é 1, caso seja 1, será utilizada uma fonte de gas

                sistem1.cell(0).pres = sistem1.config().CCPres.pres[0];
                sistem1.cell(0).acsr.tipo = 2;
                sistem1.cell(0).acsr.injl.bet = sistem1.betaE;
                sistem1.cell(0).acsr.injl.temp = sistem1.config().CCPres.temperatura[0];
                sistem1.cell(0).acsr.injl.FluidoPro.TempL = fluCP.TempL;
                sistem1.cell(0).acsr.injl.FluidoPro.TempH = fluCP.TempH;
                sistem1.cell(0).acsr.injl.FluidoPro.LVisL = fluCP.LVisL;
                sistem1.cell(0).acsr.injl.FluidoPro.LVisH = fluCP.LVisH;
                if (sistem1.config().flashCompleto != 2) { // se não for composicional
                    sistem1.cell(0).acsr.injl.FluidoPro.RGO = rgoentrada;
                    sistem1.cell(0).acsr.injl.FluidoPro.Deng = fluCP.Deng;
                    sistem1.cell(0).acsr.injl.FluidoPro.Denag = fluCP.Denag;
                    sistem1.cell(0).acsr.injl.FluidoPro.yco2 = fluCP.yco2;
                    sistem1.cell(0).acsr.injl.FluidoPro.API = fluCP.API;
                    sistem1.cell(0).acsr.injl.FluidoPro.BSW = fluCP.BSW;
                    sistem1.cell(0).acsr.injl.FluidoPro.RenovaFluido();
                } else { // caso composicional
                    fluCP.atualizaPropCompStandard();
                    fluCP.atualizaPropComp(sistem1.presE, sistem1.tempE);
                    double rl = fluCP.MasEspLiq(sistem1.presE, sistem1.tempE);
                    double titH = fluCP.FracMassHidra(sistem1.presE, sistem1.tempE);
                    double rc = sistem1.cell(0).fluicol.MasEspFlu(sistem1.presE, sistem1.tempE);
                    double rlMix = sistem1.betaE * rc + (1. - sistem1.betaE) * rl;
                    double val1 = ((1. - sistem1.betaE) * rl * titH / (1. - titH));
                    sistem1.titE = val1 / (rlMix + val1);
                    sistem1.cell(0).acsr.injl.FluidoPro = fluCP;
                }
                sistem1.cell(0).fluiL = &sistem1.cell(0).acsr.injl.FluidoPro;
                chutemass = 1000. * pow(sistem1.cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.); // estimativa
                                                                                                    // de vazao de entrada na tubulkacao
            } else { // caso em que so se tem gas na entrada
                double rgoentrada = 0.;
                sistem1.cell(0).pres = sistem1.config().CCPres.pres[0];
                sistem1.cell(0).acsr.tipo = 1;
                sistem1.cell(0).acsr.injg.temp = sistem1.config().CCPres.temperatura[0];
                sistem1.cell(0).acsr.injg.FluidoPro.TempL = fluCP.TempL;
                sistem1.cell(0).acsr.injg.FluidoPro.TempH = fluCP.TempH;
                sistem1.cell(0).acsr.injg.FluidoPro.LVisL = fluCP.LVisL;
                sistem1.cell(0).acsr.injg.FluidoPro.LVisH = fluCP.LVisH;
                if (sistem1.config().flashCompleto != 2) {
                    double rgST = fluCP.Deng * 1.225;
                    double rg = fluCP.MasEspGas(sistem1.presE, sistem1.tempE);
                    double rl = fluCP.MasEspoleo(sistem1.presE, sistem1.tempE);
                    double rcST = sistem1.cell(0).fluicol.MasEspFlu(1.01, 20.);
                    double rc = sistem1.cell(0).fluicol.MasEspFlu(sistem1.presE, sistem1.tempE);
                    double val1;
                    double val2;
                    double titT;
                    val1 = (rcST / rc) * (rg / rgST) * sistem1.config().CCPres.bet[0] / tit;
                    val2 = (rg / rl) * (1 - tit) / tit;
                    titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                    sistem1.titE = titT;
                    sistem1.betaE = val1 / (val2 + val1);
                    rgoentrada = 1e6 + (*sistem1.vg1dSP).localtiny;
                    if (rgoentrada > 1400.)
                        rgoentrada = 1e6 + (*sistem1.vg1dSP).localtiny;
                    fluCP.RGO = sistem1.cell(0).flui.RGO = sistem1.cell(0).acsr.injg.FluidoPro.RGO = rgoentrada;
                    sistem1.cell(0).acsr.injg.FluidoPro.Deng = fluCP.Deng;
                    sistem1.cell(0).acsr.injg.FluidoPro.Denag = fluCP.Denag;
                    sistem1.cell(0).acsr.injg.FluidoPro.yco2 = fluCP.yco2;
                    sistem1.cell(0).acsr.injg.FluidoPro.API = fluCP.API;
                    sistem1.cell(0).acsr.injg.FluidoPro.BSW = fluCP.BSW;
                    sistem1.cell(0).acsr.injg.FluidoPro.RenovaFluido();
                } else {
                    fluCP.atualizaPropCompStandard();
                    fluCP.atualizaPropComp(sistem1.presE, sistem1.tempE);
                    double rgST = fluCP.Deng * 1.225;
                    double rg = fluCP.MasEspGas(sistem1.presE, sistem1.tempE);
                    double rl = fluCP.MasEspoleo(sistem1.presE, sistem1.tempE);
                    double rcST = sistem1.cell(0).fluicol.MasEspFlu(1.01, 20.);
                    double rc = sistem1.cell(0).fluicol.MasEspFlu(sistem1.presE, sistem1.tempE);
                    double val1 = (rcST / rc) * (rg / rgST) * sistem1.config().CCPres.bet[0] / tit;
                    double val2 = (rg / rl) * (1 - tit) / tit;
                    double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                    sistem1.titE = titT;
                    sistem1.betaE = val1 / (val2 + val1);
                    sistem1.cell(0).acsr.injg.FluidoPro = fluCP;
                }
                if (tit >= (1. - (*sistem1.vg1dSP).localtiny)) {
                    sistem1.cell(0).acsr.injg.seco = 0;
                    sistem1.cell(0).acsr.injg.razCompGas = sistem1.config().CCPres.bet[0];
                } else {
                    sistem1.cell(0).acsr.injg.seco = 0;
                    sistem1.cell(0).acsr.injg.razCompGas = sistem1.config().CCPres.bet[0];
                }
                chutemass = 10000. * pow(sistem1.cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.); // estimativa
                // de vazao de entrada na tubulkacao
                sistem1.cell(0).fluiL = &sistem1.cell(0).acsr.injg.FluidoPro;
            }
            if (inichute > 0)
                chutemass = inichute;
            // solucao permanente para condicao de pressao na entrada sem choke ativo
            if (sistem1.config().chokep.abertura[0] > 0.6)
                ciclo = sistem1.buscaProdPresPresPerm(chutemass);
            // solucao permanente para condicao de pressao na entrada com choke ativo
            else {
                if (sistem1.config().chokep.abertura[0] > 1e-15)
                    ciclo = sistem1.buscaProdPresPresPerm2(chutemass);
                else
                    ciclo = sistem1.buscaProdPresPresPerm3(inichute);
            }
        } else if (sistem1.config().ConContEntrada == 2) { // caso em que a condicao de contorno e pressao e vazao de injecao
            ciclo = sistem1.buscaProdPfundoPerm3(sistem1.config().CCVPres.pres[0]);
        }

    } else { // solucao para o caso em que o sistema é de poco injetor
        sistem1.config().imprimeProfile(sistem1.celula, sistem1.flut, 0, sistem1.indTramo);
        sistem1.config().resumoPermanente(sistem1.celula, sistem1.celulaG, sistem1.pGSup, sistem1.presiniG, sistem1.indTramo);
        if(sistem1.config().nintermi>0)sistem1.config().resumoIntermitencia(sistem1.celula, sistem1.indTramo);
        if(sistem1.config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<sistem1.config().nCelUnit; iCelUni++)
        	sistem1.config().relatorioCelulaUnitaria(sistem1.celula,sistem1.config().celUnit[iCelUni].posicP, sistem1.indTramo);
        }
        // enterramento
        for (int j = 0; j <= sistem1.ncel; j++) {
            if (sistem1.cell(j).calor.difus2D == 1) {
                sistem1.cell(j).calor.poisson2D.imprimePermanente(sistem1.indTramo);
            }
        }
        // ver a definicao das condicoes de contorno de poco injetor na classe leitura
        if (sistem1.config().condpocinj.CC == 1 || sistem1.config().condpocinj.CC == 3)
            ciclo = sistem1.buscaInjPfundoPerm1(inichute);
        else if (sistem1.config().condpocinj.CC == 0)
            ciclo = sistem1.buscaInjPfundoPerm2(inichute);
        else if (sistem1.config().condpocinj.CC == 2)
            ciclo = sistem1.buscaInjPfundoPerm3(inichute);
        else if (sistem1.config().condpocinj.CC == 4)
            ciclo = sistem1.buscaInjPfundoPerm4();
        else
            ciclo = sistem1.buscaInjPfundoPerm5(inichute);
    }
    return ciclo;
}

double SolveTramoSolteiro(SProd &sistem1, double chute0 = -1.) {
    // este metodo age em conjunto com o metodo permanenteSimples, em verdade, e um metodo que se tornou necessario principalmente
    // pelo uso do modelo composicional. O modelo composicional, a despeito da melhor qualidade no cálculo de proipŕiedades,
    // principalmente quando se tem a disputa de várias correntes de hidrocarboneto com composições distintas,
    // tem um tributo de perda de desempenho no simulador, mesmo no calculo permanente. No caso permanente, para agilizar a convergência,
    // um primeiro processo de convergência é feito aplicando um modelo black-oil para depois usar o resultado do calculo black-oil
    // como chute para o modelo composicional, esta abordagem levou a uma melhora no desempenho do calculo permanente composicional
    sistem1.modoPerm = 1;
    double saida = -1.e10;
    // caso sem necessidade de fazer uma pre-simulacao black-oil para aumento de desempenho
    if (sistem1.config().flashCompleto != 2 ||
        (sistem1.config().pocinjec == 1 && sistem1.config().condpocinj.CC == 4))
        saida = permanenteSimples(sistem1, chute0);
    else {
        if (sistem1.config().pocinjec == 0) {
            // alterando as "flags" que indicam que a simulacao e composicional para flags indicando que a simulacao sera black-oil
            sistem1.config().flashCompleto = 0;
            int calclat0 = sistem1.CalcLat;
            sistem1.CalcLat = 0;
            int tipoFluido0 = sistem1.config().tipoFluido;
            sistem1.config().tipoFluido = 0;
            for (int i = 0; i <= sistem1.ncel; i++) {
                sistem1.cell(i).flui.flashCompleto = 0;
                if (sistem1.cell(i).acsr.tipo == 1)
                    sistem1.cell(i).acsr.injg.FluidoPro.flashCompleto = 0;
                else if (sistem1.cell(i).acsr.tipo == 2)
                    sistem1.cell(i).acsr.injl.FluidoPro.flashCompleto = 0;
                else if (sistem1.cell(i).acsr.tipo == 3)
                    sistem1.cell(i).acsr.ipr.FluidoPro.flashCompleto = 0;
                else if (sistem1.cell(i).acsr.tipo == 15) {
                    sistem1.cell(i).acsr.radialPoro.flup.flashCompleto = 0;
                    int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.radialPoro.celula[j].flup.flashCompleto = 0;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 16) {
                    sistem1.cell(i).acsr.poroso2D.dados.flup.flashCompleto = 0;
                    int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[j].flup.flashCompleto = 0;
                    }
                    int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
                    for (int j = 0; j < ncelEle; j++) {
                        sistem1.cell(i).acsr.poroso2D.malha.mlh2d[j].flup.flashCompleto = 0;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 10)
                    sistem1.cell(i).acsr.injm.FluidoPro.flashCompleto = 0;
            }
            // calculo de solucao como black oil
            saida = permanenteSimples(sistem1, chute0);

            // preparando os valores de pressao ou vazao de fundo para serem utilizadas como chutes iniciais para
            // a solucao black-oil
            double inichute;
            if (sistem1.config().ConContEntrada == 0)
                inichute = sistem1.cell(0).pres;
            else {
                double pl = sistem1.cell(1).presaux;
                double tl = sistem1.cell(1).temp;
                if (sistem1.cell(0).acsr.tipo == 1)
                    inichute = (sistem1.cell(1).MC - sistem1.cell(1).Mliqini) * 86400 /
                                   (sistem1.cell(1).flui.Deng * 1.225) +
                               (1. - sistem1.cell(1).flui.BSW) *
                                   sistem1.cell(1).QL * sistem1.cell(1).flui.RS(tl, pl) * 86400;
                else if (sistem1.cell(0).acsr.tipo == 2)
                    inichute = sistem1.cell(1).QL * 86400;
            }
            sistem1.config().flashCompleto = 2;
            sistem1.CalcLat = calclat0;
            sistem1.config().tipoFluido = tipoFluido0;
            for (int i = 0; i <= sistem1.ncel; i++) {
                // retornando as "flags" que indicam qual modelo de calculo de propriedade de fluidos sera utilizada
                // para a condicao composicional
                sistem1.cell(i).flui.flashCompleto = 2;
                if (sistem1.cell(i).acsr.tipo == 1)
                    sistem1.cell(i).acsr.injg.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 2)
                    sistem1.cell(i).acsr.injl.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 3)
                    sistem1.cell(i).acsr.ipr.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 15) {
                    sistem1.cell(i).acsr.radialPoro.flup.flashCompleto = 2;
                    int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.radialPoro.celula[j].flup.flashCompleto = 2;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 16) {
                    sistem1.cell(i).acsr.poroso2D.dados.flup.flashCompleto = 2;
                    int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[j].flup.flashCompleto = 2;
                    }
                    int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
                    for (int j = 0; j < ncelEle; j++) {
                        sistem1.cell(i).acsr.poroso2D.malha.mlh2d[j].flup.flashCompleto = 2;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 10)
                    sistem1.cell(i).acsr.injm.FluidoPro.flashCompleto = 2;
            }
            sistem1.buscaIni = 0;
            // existem duas maneiras de fazer a solucao composicional, uma completa, todas as propriedades sao calculadas a
            // cada ciclo de convergencia, e uma abordagem que exige que, a cada ciclo iterativo, um novo flash seja calculado, celula a celula
            // uma outra abordagem e mais aproximada, lancando mao da solucao black-oil, pode-se ter uma estimativa dos patamares de pressao
            // e de como a mistura de correntes se da e a sua consequente alteracaoi de composicao. Com este cenario inicial, obtido da solucao
            // black-oil, pode-se fazer um conjuntoi de tabelas flashs com o range de pressao e temperatura definido a partir doi resultado
            // desta simulacao, alem da estimativa de qual a composicao das correntes trecho a trecho, nbesta caso,
            // usa-se o metodo preparaTabelaDin
            if (sistem1.config().tabelaDinamica == 1)
                preparaTabDin(sistem1);
            // calculo de solucao composicional
            saida = permanenteSimples(sistem1, inichute);
        } else {
            // neste namesapace faz-se algo similar ao que foi feito anteriormente para sistemas de producao ou exportacao,
            // para sistemas de injecao. Nestre caso, altera-se a "flag" que indica que o calcvulo não e composicional para black oil
            sistem1.config().flashCompleto = 3;
            int calclat0 = sistem1.CalcLat;
            sistem1.CalcLat = 0;
            int tipoFluido0 = sistem1.config().tipoFluido;
            sistem1.config().tipoFluido = 0;
            for (int i = 0; i <= sistem1.ncel; i++) {
                sistem1.cell(i).flui.flashCompleto = 3;
                if (sistem1.cell(i).acsr.tipo == 1)
                    sistem1.cell(i).acsr.injg.FluidoPro.flashCompleto = 3;
                else if (sistem1.cell(i).acsr.tipo == 2)
                    sistem1.cell(i).acsr.injl.FluidoPro.flashCompleto = 3;
                else if (sistem1.cell(i).acsr.tipo == 3)
                    sistem1.cell(i).acsr.ipr.FluidoPro.flashCompleto = 3;
                else if (sistem1.cell(i).acsr.tipo == 15) {
                    sistem1.cell(i).acsr.radialPoro.flup.flashCompleto = 3;
                    int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.radialPoro.celula[j].flup.flashCompleto = 3;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 16) {
                    sistem1.cell(i).acsr.poroso2D.dados.flup.flashCompleto = 3;
                    int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[j].flup.flashCompleto = 3;
                    }
                    int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
                    for (int j = 0; j < ncelEle; j++) {
                        sistem1.cell(i).acsr.poroso2D.malha.mlh2d[j].flup.flashCompleto = 3;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 10)
                    sistem1.cell(i).acsr.injm.FluidoPro.flashCompleto = 3;
            }
            // solucao de poco injetor black oil
            saida = permanenteSimples(sistem1, chute0);

            // retorno para a "flag" que indica moidelo composicional
            sistem1.config().flashCompleto = 2;
            sistem1.CalcLat = calclat0;
            sistem1.config().tipoFluido = tipoFluido0;
            for (int i = 0; i <= sistem1.ncel; i++) {
                sistem1.cell(i).flui.flashCompleto = 2;
                if (sistem1.cell(i).acsr.tipo == 1)
                    sistem1.cell(i).acsr.injg.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 2)
                    sistem1.cell(i).acsr.injl.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 3)
                    sistem1.cell(i).acsr.ipr.FluidoPro.flashCompleto = 2;
                else if (sistem1.cell(i).acsr.tipo == 15) {
                    sistem1.cell(i).acsr.radialPoro.flup.flashCompleto = 2;
                    int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.radialPoro.celula[j].flup.flashCompleto = 2;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 16) {
                    sistem1.cell(i).acsr.poroso2D.dados.flup.flashCompleto = 2;
                    int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[j].flup.flashCompleto = 2;
                    }
                    int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
                    for (int j = 0; j < ncelEle; j++) {
                        sistem1.cell(i).acsr.poroso2D.malha.mlh2d[j].flup.flashCompleto = 2;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 10)
                    sistem1.cell(i).acsr.injm.FluidoPro.flashCompleto = 2;
            }
            // preparacao do chute para poço injetor
            double chute;
            if (sistem1.config().condpocinj.CC == 1 || sistem1.config().condpocinj.CC == 2 || sistem1.config().condpocinj.CC == 3)
                chute = sistem1.cell(0).acsr.injg.QGas;
            else
                chute = sistem1.pGSup;
            // prepara tabela dinamica, caso tenha sido solicitada
            if (sistem1.config().tabelaDinamica == 1)
                preparaTabDin(sistem1);
            // solucao composicional
            saida = permanenteSimples(sistem1, chute);
        }
    }
    // impressão de perfis apos a obtencao de solucao permanente
    if ((*sistem1.vg1dSP).chaverede == 0 && sistem1.config().transiente == 0 && sistem1.config().AP == 0) { // impŕessao de perfis e tendencias, quando a opcao transiente
        // nao esta ativa
        sistem1.config().imprimeProfile(sistem1.celula, sistem1.flut, 0, sistem1.indTramo);
        sistem1.config().resumoPermanente(sistem1.celula, sistem1.celulaG, sistem1.pGSup, sistem1.presiniG, sistem1.indTramo);
        if(sistem1.config().nintermi>0)sistem1.config().resumoIntermitencia(sistem1.celula, sistem1.indTramo);
        if(sistem1.config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<sistem1.config().nCelUnit; iCelUni++)
        	sistem1.config().relatorioCelulaUnitaria(sistem1.celula,sistem1.config().celUnit[iCelUni].posicP, sistem1.indTramo);
        }
        // enterramento
        for (int j = 0; j <= sistem1.ncel; j++) {
            if (sistem1.cell(j).calor.difus2D == 1) {
                sistem1.cell(j).calor.poisson2D.imprimePermanente(sistem1.indTramo);
            }
            if (sistem1.cell(j).acsr.tipo == 16) {

                sistem1.cell(j).acsr.poroso2D.imprimePseudo();
            }
            if (sistem1.cell(j).acsr.tipo == 15) {
                FullMtx<double> matrizsaida(sistem1.cell(j).acsr.radialPoro.nglobal, 11);
                matrizsaida = sistem1.cell(j).acsr.radialPoro.perfil();
                ostringstream saida;
                // saida << pathPrefixoArqSaida << "PerfisPocoRadial" << "-" << kontaTempoImp
                saida << pathPrefixoArqSaida << "PerfisPocoRadial" << "-" << 0
                      << "-" << sistem1.cell(j).acsr.radialPoro.posicMarlim << ".dat";
                string tmp = saida.str();
                ofstream escreveIni(tmp.c_str(), ios_base::out);
                escreveIni << "tempo = " << 0 << endl;
                escreveIni << " raio (m) ;" << " raio (pol.) ;" << " pressao (kgf/cm2) ;" << " vazao total (sm3/d) ;" << " vazao de oleo (sm3/d) ;" << " vazao de gas (sm3/d) ;" << " vazao de agua (sm3/d) ;" << "saturacao de liquido (-) ;" << " saturacao de agua (-) ;" << " fracao volumetrica de gas homogeneo (-) ;" << endl;
                escreveIni << matrizsaida;
                escreveIni.close();
            }
        }
        sistem1.kimpT = 1;
        for (int i = 0; i < sistem1.config().ntendp; i++) {
            sistem1.ImprimeTrendPCab(i);
            sistem1.config().imprimeTrend(sistem1.celula, sistem1.MatTrendP[i], 0, i, 0);
            sistem1.ImprimeTrendP(i);
        }
        if (sistem1.config().lingas == 1) {
            sistem1.config().imprimeProfileG(sistem1.celulaG, sistem1.flutG, 0, sistem1.indTramo);
            for (int i = 0; i < sistem1.config().ntendg; i++) {
                sistem1.ImprimeTrendGCab(i);
                sistem1.config().imprimeTrendG(sistem1.celulaG, sistem1.MatTrendG[i], 0, i, 0, 0);
                sistem1.ImprimeTrendG(i);
            }
        }
    }
    if ((*sistem1.vg1dSP).chaverede == 0 && sistem1.config().transiente == 0) { // impressao do arquivo de log, quando a
        // opcao transiente nao esta ativa
        ostringstream saidaT;
        if (sistem1.indTramo < 0) {
            saidaT << sistem1.tmpLog;
        } else {
            saidaT << "Tramo" << sistem1.indTramo << "-" << sistem1.tmpLog;
        }
        string tmp = saidaT.str();
        ofstream escreveIni(tmp.c_str(), ios_base::out);
        if (fabs(saida) > 1e9) {
            escreveIni << "    PROBLEMAS NA CONVERGENCIA DO PERMANENTE     " << endl;
            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                       "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################",
                       "", "");
        }
        if (sistem1.config().saidaClassica == 1) {
            srand(time(NULL));
            int frase = rand() % 15;
            escreveIni << "*******************************************************************************" << endl;
            escreveIni << "                                  UFA!!!!!!!!                                  " << endl;
            escreveIni << saidaTexto[frase] << endl;
            escreveIni << saidaSubTexto[frase] << endl;
            escreveIni << "*******************************************************************************" << endl;
        } else
            escreveIni << "                                 FIM                                  " << endl;

        time_t now = time(0);
        tm *ltm = localtime(&now); ///////////Retirado de https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
        escreveIni << "datahora = ";
        escreveIni << ltm->tm_mday << "/";
        escreveIni << 1 + ltm->tm_mon << "/";
        escreveIni << 1900 + ltm->tm_year << " ";
        escreveIni << 0 + ltm->tm_hour << ":";
        escreveIni << 0 + ltm->tm_min << ":";
        escreveIni << 0 + ltm->tm_sec;
        escreveIni << endl;

        int diaFim = (ltm->tm_mday);
        int horaFim;
        if (diaFim == diaIni)
            horaFim = ltm->tm_hour;
        else
            horaFim = ltm->tm_hour + 24;
        horaFim *= 3600;
        int minutoFim = 60 * ltm->tm_min;
        int segundoFim = ltm->tm_sec;
        int totalFim = horaFim + minutoFim + segundoFim;
        int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
        escreveIni << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
        escreveIni << "     Versao    " << versao << endl;

        escreveIni.close();
    }
    if (sistem1.config().transiente == 1) {

        for (int j = 0; j <= sistem1.ncel; j++) {
            if (sistem1.cell(j).calor.difus2D == 1) {
                sistem1.cell(j).calor.poisson2D.imprimePermanente(sistem1.indTramo);
            }
            if (sistem1.cell(j).acsr.tipo == 16) {

                sistem1.cell(j).acsr.poroso2D.imprimePseudo();
                if (sistem1.cell(j).acsr.poroso2D.dados.temp.tempoImp[0] <= 1e-15)
                    sistem1.cell(j).acsr.poroso2D.kontaTempo++;
            }
            if (sistem1.cell(j).acsr.tipo == 15) {
                FullMtx<double> matrizsaida(sistem1.cell(j).acsr.radialPoro.nglobal, 11);
                matrizsaida = sistem1.cell(j).acsr.radialPoro.perfil();
                ostringstream saida;
                // saida << pathPrefixoArqSaida << "PerfisPocoRadial" << "-" << kontaTempoImp
                saida << pathPrefixoArqSaida << "PerfisPocoRadial" << "-" << 0
                      << "-" << sistem1.cell(j).acsr.radialPoro.posicMarlim << ".dat";
                string tmp = saida.str();
                ofstream escreveIni(tmp.c_str(), ios_base::out);
                escreveIni << "tempo = " << 0 << endl;
                escreveIni << " raio (m) ;" << " raio (pol.) ;" << " pressao (kgf/cm2) ;" << " vazao total (sm3/d) ;" << " vazao de oleo (sm3/d) ;" << " vazao de gas (sm3/d) ;" << " vazao de agua (sm3/d) ;" << "saturacao de liquido (-) ;" << " saturacao de agua (-) ;" << " fracao volumetrica de gas homogeneo (-) ;" << endl;
                escreveIni << matrizsaida;
                escreveIni.close();
                if (sistem1.cell(j).acsr.radialPoro.temp.tempoImp[0] <= 1e-15)
                    sistem1.cell(j).acsr.radialPoro.kontaTempoImp++;
            }
        }

        if (sistem1.config().tabelaDinamica == 1) {
            sistem1.config().tabelaDinamica = 0;
            for (int i = 0; i <= sistem1.ncel; i++) {
                sistem1.cell(i).flui.tabelaDinamica = 0;
                if (sistem1.cell(i).acsr.tipo == 1)
                    sistem1.cell(i).acsr.injg.FluidoPro.tabelaDinamica = 0;
                else if (sistem1.cell(i).acsr.tipo == 2)
                    sistem1.cell(i).acsr.injl.FluidoPro.tabelaDinamica = 0;
                else if (sistem1.cell(i).acsr.tipo == 3)
                    sistem1.cell(i).acsr.ipr.FluidoPro.tabelaDinamica = 0;
                else if (sistem1.cell(i).acsr.tipo == 15) {
                    sistem1.cell(i).acsr.radialPoro.flup.tabelaDinamica = 0;
                    int ncelRad = sistem1.cell(i).acsr.radialPoro.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.radialPoro.celula[j].flup.tabelaDinamica = 0;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 16) {
                    sistem1.cell(i).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
                    int ncelRad = sistem1.cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int j = 0; j < ncelRad; j++) {
                        sistem1.cell(i).acsr.poroso2D.dados.transfer.celula[j].flup.tabelaDinamica = 0;
                    }
                    int ncelEle = sistem1.cell(i).acsr.poroso2D.malha.nele;
                    for (int j = 0; j < ncelEle; j++) {
                        sistem1.cell(i).acsr.poroso2D.malha.mlh2d[j].flup.tabelaDinamica = 0;
                    }
                } else if (sistem1.cell(i).acsr.tipo == 10)
                    sistem1.cell(i).acsr.injm.FluidoPro.tabelaDinamica = 0;
            }
        }
    }
    return saida;
}

void leituraAPparalelo(string nomeArquivoAP, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson, SProd &sistem1) {
    // metodo em que e feita a sequencia de simulacoes de analise de sensibilidade, observando que a analise de sensibilidade e
    // feita apenas para tramos simples.
    APara analisePara(sistem1.vg1dSP, nomeArquivoAP, sistem1.ncel, sistem1.config().celp,
                      sistem1.config().flup, sistem1.config().bcs, sistem1.config().multiBcs, sistem1.config().fonteg); // construtor do sistema de analise de sensibilidade, ver a classe LerAP. Neste cosntrutor é lido um arquivo
    // json onde se esta registradoas variaveis que estao envolvidas na analise de sensibilidade e os valores que se
    // testara destas variaveis na analise de sensibilidade
    int imprime = 0;
    struct varSaida {
        double presIni;
        double tempFim;
        int falha;
    };
    std::vector<std::pair<int, varSaida>> dadosAP;

    if (analisePara.vfp == 1 || analisePara.vfp == 3)
        analisePara.cabecalhoAP(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,sistem1.config().multiBcs,
								sistem1.config().bvol, sistem1.config().dpreq);
    else if (analisePara.vfp == 0 || analisePara.vfp == 2)
        analisePara.cabecalhoAPImex(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                    sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                    sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,sistem1.config().multiBcs,
									sistem1.config().bvol, sistem1.config().dpreq);

    ostringstream saidaAP;
    saidaAP << pathPrefixoArqSaida << "sucessoAP.dat"; // criacao de arquivo com relatorio de quais pontos de analise de sensibilidade tiveram sucesso e quais falharam

    string tmpAP = saidaAP.str();
    ofstream escreveAP(tmpAP.c_str(), ios_base::out);
    escreveAP << "relatorio de falhas da Analise Parametrica para um Tramo " << endl;
    escreveAP.close();
    int indfalha[analisePara.nVariaveis];
    sistem1.kimpT = 1;
    for (int i = 0; i < sistem1.config().ntendp; i++) {
        sistem1.ImprimeTrendPCab(i);
    }
    if (sistem1.config().lingas == 1) {
        for (int i = 0; i < sistem1.config().ntendg; i++) {
            sistem1.ImprimeTrendGCab(i);
        }
    }

    Ler *vecArq;
    vecArq = new Ler[analisePara.nVariaveis];
    int vecnoextremo[analisePara.nVariaveis];
    int vecnoinicial[analisePara.nVariaveis];
    int vecderivaAnel[analisePara.nVariaveis];
    int vecbloq[analisePara.nVariaveis];
    double vecbetaRev[analisePara.nVariaveis];
    double vecbetaRevini[analisePara.nVariaveis];
    double titRev[analisePara.nVariaveis];
    double titRevini[analisePara.nVariaveis];
    double vecdtCicMin[analisePara.nVariaveis];
    varGlob1D *vg1dTramo;
    vg1dTramo = new varGlob1D[analisePara.nVariaveis];
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++)vg1dTramo[iSeq] = (*sistem1.vg1dSP);
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) {
        vg1dTramo[iSeq].sequenciaAP = iSeq;
        string nomeArquivoLogAP = nomeArquivoLog;
        nomeArquivoLogAP.erase(nomeArquivoLogAP.size() - 4);
        std::string sIseq = std::to_string(iSeq);
        nomeArquivoLogAP = nomeArquivoLogAP + "_AP_" + sIseq + ".log";
        vecnoextremo[iSeq] = sistem1.noextremo;
        vecnoinicial[iSeq] = sistem1.noinicial;
        vecderivaAnel[iSeq] = sistem1.derivaAnel;
        vecbloq[iSeq] = sistem1.bloq;
        vecbetaRev[iSeq] = sistem1.betaRev;
        vecbetaRevini[iSeq] = sistem1.betaRevini;
        titRev[iSeq] = sistem1.titRev;
        titRevini[iSeq] = sistem1.titRevini;
        vecdtCicMin[iSeq] = sistem1.dtCicMin;
        vecArq[iSeq].copiaSemJson(sistem1.arq);
        vecArq[iSeq].vg1dSP = &vg1dTramo[iSeq];
    }
    double chute=-1;
    if(sistem1.config().pocinjec==0){
    	double falha0=SolveTramoSolteiro(sistem1);
    	if(fabs(falha0)<1e9){
    		if(sistem1.config().ConContEntrada==0)chute=sistem1.cell(0).pres;
    		else{
    			double pres=sistem1.cell(0).pres;
    			double temp=sistem1.cell(0).temp;
    			if(sistem1.cell(0).acsr.tipo==1)
    				chute=sistem1.cell(0).QG*86400*sistem1.cell(0).flui.MasEspGas(pres, temp)/(sistem1.cell(0).flui.Deng*1.229);
    			else if(sistem1.cell(0).acsr.tipo==2)
    				chute=sistem1.cell(0).QL*86400*sistem1.cell(0).flui.MasEspLiq(pres, temp)/(sistem1.cell(0).flui.MasEspLiq(1.0, 20));
    		}
    	}
    }
    #pragma omp parallel for num_threads(analisePara.nthrdAP)
	for(int iSeq=0;iSeq<analisePara.nVariaveis;iSeq++){//laco onde e calculado o permanente de cada combinacao
		//algumas variaveis ja se encontram na condicao correta para se apolicar no sistema de simulacao,
		//outras precisam de um "pos procesaamento"
		double vazgasG;
		double presE;
		double tempE;
		double titE;
		double betaE;
		double vazE;
		int indChk;
		double falha=0.;
		SProd sistem2;
		sistem2.copiaSemJson(vecArq[iSeq], vecnoextremo[iSeq],  vecnoinicial[iSeq], vecderivaAnel[iSeq], vecbloq[iSeq],
				vecbetaRev[iSeq],vecbetaRevini[iSeq],titRev[iSeq], titRevini[iSeq], vecdtCicMin[iSeq]);

        // criar objeto de simulacao
        // construtor do objeto que representa o tramo

        if (analisePara.vfp == 1 || analisePara.vfp == 3){
            analisePara.selecaoAPsemImpre(sistem2.ncelGas, sistem2.chokeSup, sistem2.celula, sistem2.celulaG,
                                          sistem2.config().flup,
                                          sistem2.config().IPRS, sistem2.config().valv, sistem2.config().fonteg,
                                          sistem2.config().fontel, sistem2.config().fontem, sistem2.config().furo, sistem2.config().bcs,
										  sistem2.config().multiBcs,sistem2.config().bvol, sistem2.config().dpreq,
                                          sistem2.pGSup, sistem2.temperatura, sistem2.presiniG, sistem2.tempiniG, vazgasG,
                                          presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem2.config().correcao.dPdLHidro, sistem2.config().correcao.dPdLFric,
                                          sistem2.config().correcao.dTdL);

        }
        else if (analisePara.vfp == 0 || analisePara.vfp == 2){
            analisePara.selecaoAPImexsemImpre(sistem2.ncelGas, sistem2.chokeSup, sistem2.celula, sistem2.celulaG,
                                              sistem2.config().flup,
                                              sistem2.config().IPRS, sistem2.config().valv, sistem2.config().fonteg,
                                              sistem2.config().fontel, sistem2.config().fontem, sistem2.config().furo, sistem2.config().bcs,
											  sistem2.config().multiBcs,sistem2.config().bvol, sistem2.config().dpreq,
                                              sistem2.pGSup, sistem2.temperatura, sistem2.presiniG, sistem2.tempiniG, vazgasG,
                                              presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem2.config().correcao.dPdLHidro, sistem2.config().correcao.dPdLFric,
                                              sistem2.config().correcao.dTdL);

        }
        // variavei que precisam de um pos processamento para se encaixar na condicao de simulacao
        if (sistem2.config().lingas > 0 && analisePara.listaV.vgasinj == 1 && analisePara.APGasInj.parserieVazGas > 0) {
            // caso tenha linha de gas e analise de sensibilidade para vazao de injecao, a entrada no json e em stdM3,
            // mas no simulador, deve ser atualizada para a vazao massica
            if (sistem2.config().gasinj.tipoCC == 1) {
                sistem2.gasCell(0).massfonteCH = vazgasG * (sistem2.config().flug.Deng * 1.225) / 86400;
            }
        }
        if (analisePara.listaV.vchk == 1 && analisePara.APCHK.parserieAbre > 0) {
            // carregando o valo de abertuda do choke na variavel do sistema
            sistem2.config().chokep.abertura[0] = analisePara.APCHK.abertura[indChk];
        }
        if (sistem2.config().ConContEntrada == 1 && analisePara.listaV.vpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de pressao
            if (analisePara.APpEntrada.parseriePres > 0)
                sistem2.presE = presE;
            if (analisePara.APpEntrada.parserieBet > 0)
                sistem2.betaE = betaE;
            if (analisePara.APpEntrada.parserieTemp > 0)
                sistem2.tempE = tempE;
            if (analisePara.APpEntrada.parserieTit > 0)
                sistem2.titE = titE;
        } else if (sistem2.config().ConContEntrada == 2 && analisePara.listaV.vvazpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de vazao e pressao
            if (analisePara.APvpEntrada.parseriePres > 0)
                sistem2.presE = presE;
            if (analisePara.APvpEntrada.parserieBet > 0)
                sistem2.betaE = betaE;
            if (analisePara.APvpEntrada.parserieTemp > 0)
                sistem2.tempE = tempE;
            if (analisePara.APvpEntrada.parserieMass > 0) {
                double tit;
                tit = sistem2.config().flup[0].FracMassHidra(sistem2.presE, sistem2.tempE);
                double vazC = vazE * betaE;
                sistem2.cell(0).acsr.injm.MassC = vazC;
                sistem2.cell(0).acsr.injm.MassG = (vazE - vazC) * tit;
                sistem2.cell(0).acsr.injm.MassP = vazE - vazC - sistem2.cell(0).acsr.injm.MassG;
            }
        } else if (sistem2.config().pocinjec == 1 && analisePara.listaV.vpocinj == 1) {
            // carregando os valores de analise de senisbilidade para condicao de contorno a montante de vazao, poco injetor
            sistem2.config().condpocinj.presfundo = presE;
            if (sistem2.config().flashCompleto == 0) {
                sistem2.cell(0).acsr.injl.temp = tempE;
                sistem2.cell(0).acsr.injl.QLiq = vazE;
            } else if (sistem2.config().flashCompleto == 1) {
                sistem2.cell(0).acsr.injg.temp = tempE;
                sistem2.cell(0).acsr.injg.QGas = vazE;
            }
        }
		#pragma omp critical
			cout<<"Resolvendo Sequencia "<<iSeq<<" da analise de sensibilidade"<< endl;
		falha=SolveTramoSolteiro(sistem2,chute);
		if(fabs(falha)>(0.9e10) && chute!=-1){
			falha=SolveTramoSolteiro(sistem2);
		}
		if(fabs(falha)>1e9){
			//caso ocorre uma falha nma busca da solucao para um caso da analise de sensibilidade
			#pragma omp critical
				cout<<"----------------------------- falha no caso -------------------------------"<<endl;
		}

        if (fabs(falha) < 1e9)
            indfalha[iSeq] = 1;
        else
            indfalha[iSeq] = -1;

        varSaida tempSaida;
        tempSaida.presIni = sistem2.cell(0).pres;
        tempSaida.tempFim = sistem2.cell(sistem2.ncel).temp;
        tempSaida.falha = indfalha[iSeq];
		#pragma omp critical
			dadosAP.emplace_back(iSeq,tempSaida);
        // impressão dos perfis e tendencias da analise de sensibilidade, caso sem construcao de tabela de pressao de fundo
        if(analisePara.tipoAP == 0 && analisePara.imprimePerfil==1){
        	sistem2.config().imprimeProfile(sistem2.celula, sistem2.flut, 0, sistem2.indTramo);
        	sistem2.config().resumoPermanente(sistem2.celula, sistem2.celulaG, sistem2.pGSup,
                                     sistem2.presiniG, sistem2.indTramo);
        	if(sistem2.config().nintermi>0)sistem2.config().resumoIntermitencia(sistem2.celula, sistem2.indTramo);
        	if(sistem2.config().nCelUnit>0){
        		for(int iCelUni=0; iCelUni<sistem2.config().nCelUnit; iCelUni++)
        			sistem2.config().relatorioCelulaUnitaria(sistem2.celula,sistem2.config().celUnit[iCelUni].posicP, sistem2.indTramo);
        	}
        	// enterramento
        	for (int j = 0; j <= sistem2.ncel; j++) {
            	if (sistem2.cell(j).calor.difus2D == 1) {
                	sistem2.cell(j).calor.poisson2D.imprimePermanente(sistem2.indTramo);
            	}
        	}
        	sistem2.kimpT = 1;
        	for (int i = 0; i < sistem2.config().ntendp; i++) {
            	sistem2.config().imprimeTrend(sistem2.celula, sistem2.MatTrendP[i], 0, i, 0);
            	sistem2.ImprimeTrendP(i);
        	}
        	if (sistem2.config().lingas == 1) {
            	sistem2.config().imprimeProfileG(sistem2.celulaG, sistem2.flutG, 0, sistem2.indTramo);
            	for (int i = 0; i < sistem2.config().ntendg; i++) {
                	sistem2.config().imprimeTrendG(sistem2.celulaG, sistem2.MatTrendG[i], 0, i, 0, 0);
                	sistem2.ImprimeTrendG(i);
            	}
        	}
        }
        (*sistem2.vg1dSP).contaExit = 0;
    }

    // saidaAP << "sucessoAP.dat";//arquivo com o relatorio de sucessos e falhas da AP

    ofstream escreveAP2(tmpAP.c_str(), ios_base::out);

    ostringstream saidaVarInt;
    saidaVarInt << "variaveisInteresseAP.dat";
    string varInt = saidaVarInt.str();
    std::sort(dadosAP.begin(), dadosAP.end(),
              [](const std::pair<int, varSaida> &a, const std::pair<int, varSaida> &b) -> bool {
                  return a.first < b.first;
              });

    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) {
        double vazgasG;
        double presE;
        double tempE;
        double titE;
        double betaE;
        double vazE;
        int indChk;
        imprime = 1;
        if (analisePara.vfp == 1) // analise de sensibilidade para problemas padrao ou para a curva de fundo Eclipse
            analisePara.selecaoAP(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                  sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                  sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
								  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                  sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                  presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                  sistem1.config().correcao.dTdL, imprime);
        else if (analisePara.vfp == 0) // anaslise de sensibilidade para curva de pressao de fundo Imex
            analisePara.selecaoAPImex(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                      sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                      sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
									  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                      sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                      presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                      sistem1.config().correcao.dTdL, imprime);
        if (indfalha[iSeq] > 0)
            escreveAP2 << iSeq << " : " << " Resultado = " << "sucesso" << endl;
        else
            escreveAP2 << iSeq << " : " << " Resultado = " << "falha" << endl;

        ofstream escreveVarInt(varInt.c_str(), ios_base::app);
        escreveVarInt << iSeq << " ; " << dadosAP[iSeq].second.falha << " ; " << dadosAP[iSeq].second.presIni << " ; " << dadosAP[iSeq].second.tempFim << " ; ";
        escreveVarInt.close();
        if (analisePara.vfp == 1 || analisePara.vfp == 3)
            analisePara.imprimeVarInteresseAP(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG,
                                              sistem1.config().flup,
                                              sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                              sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
											  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq, iSeq);
        else if (analisePara.vfp == 0 || analisePara.vfp == 2)
            analisePara.imprimeVarInteresseAPImex(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG,
                                                  sistem1.config().flup,
                                                  sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                                  sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
												  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq, iSeq);

        if (analisePara.tipoAP != 0) {
            // caso de construcao de tabela VFP:
            double BHP;
            if (indfalha[iSeq] >0)
                BHP = dadosAP[iSeq].second.presIni;
            else
                BHP = -1e10;
            analisePara.tabelaGenerica(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                       sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                       sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
									   sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                       sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                       presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                       sistem1.config().correcao.dTdL, BHP);
        }

    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    int diaFim = (ltm->tm_mday);
    int horaFim;
    if (diaFim == diaIni)
        horaFim = ltm->tm_hour;
    else
        horaFim = ltm->tm_hour + 24;
    horaFim *= 3600;
    int minutoFim = 60 * ltm->tm_min;
    int segundoFim = ltm->tm_sec;
    int totalFim = horaFim + minutoFim + segundoFim;
    int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << "datahora = ";
    escreveAP2 << ltm->tm_mday << "/";
    escreveAP2 << 1 + ltm->tm_mon << "/";
    escreveAP2 << 1900 + ltm->tm_year << " ";
    escreveAP2 << 0 + ltm->tm_hour << ":";
    escreveAP2 << 0 + ltm->tm_min << ":";
    escreveAP2 << 0 + ltm->tm_sec;
    escreveAP2 << endl;
    escreveAP2 << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
    escreveAP2 << "     Versao    " << versao << endl;


	escreveAP2.close();
	delete [] vg1dTramo;
	delete [] vecArq;
	dadosAP.clear();
}

void leituraAPparaleloReserva(string nomeArquivoAP, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson, SProd &sistem1) {
    // metodo em que e feita a sequencia de simulacoes de analise de sensibilidade, observando que a analise de sensibilidade e
    // feita apenas para tramos simples.
    APara analisePara(sistem1.vg1dSP, nomeArquivoAP, sistem1.ncel, sistem1.config().celp,
                      sistem1.config().flup, sistem1.config().bcs, sistem1.config().multiBcs, sistem1.config().fonteg); // construtor do sistema de analise de sensibilidade, ver a classe LerAP. Neste cosntrutor é lido um arquivo
    // json onde se esta registradoas variaveis que estao envolvidas na analise de sensibilidade e os valores que se
    // testara destas variaveis na analise de sensibilidade
    int imprime = 0;

    ostringstream saidaAP;
    saidaAP << "sucessoAP.dat"; // criacao de arquivo com relatorio de quais pontos de analise de sensibilidade tiveram sucesso e quais falharam

    string tmpAP = saidaAP.str();
    ofstream escreveAP(tmpAP.c_str(), ios_base::out);
    escreveAP << "relatório de falhas da Analise Parametrica para um Tramo " << endl;
    escreveAP.close();

    sistem1.kimpT = 1;
    for (int i = 0; i < sistem1.config().ntendp; i++) {
        sistem1.ImprimeTrendPCab(i);
    }
    if (sistem1.config().lingas == 1) {
        for (int i = 0; i < sistem1.config().ntendg; i++) {
            sistem1.ImprimeTrendGCab(i);
        }
    }

    int indfalha[analisePara.nVariaveis];
    if (analisePara.vfp != 1 || analisePara.tipoAP != 0) {
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Parâmetros inadequados para a analise de sensibilidade em paralelo",
                   "", "");
    }
    SProd sistem2[analisePara.nVariaveis];
    varGlob1D vg1dTramo[analisePara.nVariaveis];
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) {
        string nomeArquivoLogAP = nomeArquivoLog;
        sistem2[iSeq] = SProd(sistem1.config().impfile, nomeArquivoLogAP, validacaoJson, sistem1.config().tipoSimulacao, &vg1dTramo[iSeq]);
    }
#pragma omp parallel for num_threads(analisePara.nthrdAP)
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) { // laco onde e calculado o permanente de cada combinacao
        // algumas variaveis ja se encontram na condicao correta para se apolicar no sistema de simulacao,
        // outras precisam de um "pos procesaamento"
        double vazgasG;
        double presE;
        double tempE;
        double titE;
        double betaE;
        double vazE;
        int indChk;
        double falha = 0.;

        // criar objeto de simulacao
        // construtor do objeto que representa o tramo

        analisePara.selecaoAPsemImpre(sistem2[iSeq].ncelGas, sistem2[iSeq].chokeSup, sistem2[iSeq].celula, sistem2[iSeq].celulaG,
                                      sistem2[iSeq].config().flup,
                                      sistem2[iSeq].config().IPRS, sistem2[iSeq].config().valv, sistem2[iSeq].config().fonteg,
                                      sistem2[iSeq].config().fontel, sistem2[iSeq].config().fontem, sistem2[iSeq].config().furo, sistem2[iSeq].config().bcs,
									  sistem2[iSeq].config().multiBcs,sistem2[iSeq].config().bvol, sistem2[iSeq].config().dpreq,
                                      sistem2[iSeq].pGSup, sistem2[iSeq].temperatura, sistem2[iSeq].presiniG, sistem2[iSeq].tempiniG, vazgasG,
                                      presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem2[iSeq].config().correcao.dPdLHidro, sistem2[iSeq].config().correcao.dPdLFric,
                                      sistem2[iSeq].config().correcao.dTdL);
        // variavei que precisam de um pos processamento para se encaixar na condicao de simulacao
        if (sistem2[iSeq].config().lingas > 0 && analisePara.listaV.vgasinj == 1 && analisePara.APGasInj.parserieVazGas > 0) {
            // caso tenha linha de gas e analise de sensibilidade para vazao de injecao, a entrada no json e em stdM3,
            // mas no simulador, deve ser atualizada para a vazao massica
            if (sistem2[iSeq].config().gasinj.tipoCC == 1) {
                sistem2[iSeq].gasCell(0).massfonteCH = vazgasG * (sistem2[iSeq].config().flug.Deng * 1.225) / 86400;
            }
        }
        if (analisePara.listaV.vchk == 1 && analisePara.APCHK.parserieAbre > 0) {
            // carregando o valo de abertuda do choke na variavel do sistema
            sistem2[iSeq].config().chokep.abertura[0] = analisePara.APCHK.abertura[indChk];
        }
        if (sistem2[iSeq].config().ConContEntrada == 1 && analisePara.listaV.vpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de pressao
            if (analisePara.APpEntrada.parseriePres > 0)
                sistem2[iSeq].presE = presE;
            if (analisePara.APpEntrada.parserieBet > 0)
                sistem2[iSeq].betaE = betaE;
            if (analisePara.APpEntrada.parserieTemp > 0)
                sistem2[iSeq].tempE = tempE;
            if (analisePara.APpEntrada.parserieTit > 0)
                sistem2[iSeq].titE = titE;
        } else if (sistem2[iSeq].config().ConContEntrada == 2 && analisePara.listaV.vvazpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de vazao e pressao
            if (analisePara.APvpEntrada.parseriePres > 0)
                sistem2[iSeq].presE = presE;
            if (analisePara.APvpEntrada.parserieBet > 0)
                sistem2[iSeq].betaE = betaE;
            if (analisePara.APvpEntrada.parserieTemp > 0)
                sistem2[iSeq].tempE = tempE;
            if (analisePara.APvpEntrada.parserieMass > 0) {
                double tit;
                tit = sistem2[iSeq].config().flup[0].FracMassHidra(sistem2[iSeq].presE, sistem2[iSeq].tempE);
                double vazC = vazE * betaE;
                sistem2[iSeq].cell(0).acsr.injm.MassC = vazC;
                sistem2[iSeq].cell(0).acsr.injm.MassG = (vazE - vazC) * tit;
                sistem2[iSeq].cell(0).acsr.injm.MassP = vazE - vazC - sistem2[iSeq].cell(0).acsr.injm.MassG;
            }
        } else if (sistem2[iSeq].config().pocinjec == 1 && analisePara.listaV.vpocinj == 1) {
            // carregando os valores de analise de senisbilidade para condicao de contorno a montante de vazao, poco injetor
            sistem2[iSeq].config().condpocinj.presfundo = presE;
            if (sistem2[iSeq].config().flashCompleto == 0) {
                sistem2[iSeq].cell(0).acsr.injl.temp = tempE;
                sistem2[iSeq].cell(0).acsr.injl.QLiq = vazE;
            } else if (sistem2[iSeq].config().flashCompleto == 1) {
                sistem2[iSeq].cell(0).acsr.injg.temp = tempE;
                sistem2[iSeq].cell(0).acsr.injg.QGas = vazE;
            }
        }
        cout << "Resolvendo Sequencia " << iSeq << " da analise de sensibilidade" << endl;
        (*sistem2[iSeq].vg1dSP).sequenciaAP = iSeq;
        falha = SolveTramoSolteiro(sistem2[iSeq]);
        if (fabs(falha) > 0.9e10) {
            // caso ocorre uma falha nma busca da solucao para um caso da analise de sensibilidade
            cout << "----------------------------- falha no caso -------------------------------" << endl;
        }

        if (fabs(falha) < 1e9)
            indfalha[iSeq] = 1;
        else
            indfalha[iSeq] = -1;

        // impressão dos perfis e tendencias da analise de sensibilidade, caso sem construcao de tabela de pressao de fundo
        sistem2[iSeq].config().imprimeProfile(sistem2[iSeq].celula, sistem2[iSeq].flut, 0, sistem2[iSeq].indTramo);
        sistem2[iSeq].config().resumoPermanente(sistem2[iSeq].celula, sistem2[iSeq].celulaG, sistem2[iSeq].pGSup,
                                           sistem2[iSeq].presiniG, sistem2[iSeq].indTramo);
        if(sistem2[iSeq].config().nintermi>0)sistem2[iSeq].config().resumoIntermitencia(sistem2[iSeq].celula, sistem2[iSeq].indTramo);
        if(sistem2[iSeq].config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<sistem2[iSeq].config().nCelUnit; iCelUni++)
        	sistem2[iSeq].config().relatorioCelulaUnitaria(sistem2[iSeq].celula,sistem2[iSeq].config().celUnit[iCelUni].posicP, sistem2[iSeq].indTramo);
        }
        // enterramento
        for (int j = 0; j <= sistem2[iSeq].ncel; j++) {
            if (sistem2[iSeq].cell(j).calor.difus2D == 1) {
                sistem2[iSeq].cell(j).calor.poisson2D.imprimePermanente(sistem2[iSeq].indTramo);
            }
        }
        sistem2[iSeq].kimpT = 1;
        for (int i = 0; i < sistem2[iSeq].config().ntendp; i++) {
            sistem2[iSeq].config().imprimeTrend(sistem2[iSeq].celula, sistem2[iSeq].MatTrendP[i], 0, i, 0);
            sistem2[iSeq].ImprimeTrendP(i);
        }
        if (sistem2[iSeq].config().lingas == 1) {
            sistem2[iSeq].config().imprimeProfileG(sistem2[iSeq].celulaG, sistem2[iSeq].flutG, 0, sistem2[iSeq].indTramo);
            for (int i = 0; i < sistem2[iSeq].config().ntendg; i++) {
                sistem2[iSeq].config().imprimeTrendG(sistem2[iSeq].celulaG, sistem2[iSeq].MatTrendG[i], 0, i, 0, 0);
                sistem2[iSeq].ImprimeTrendG(i);
            }
        }
        (*sistem2[iSeq].vg1dSP).contaExit = 0;
    }

    // saidaAP << "sucessoAP.dat";//arquivo com o relatorio de sucessos e falhas da AP

    ofstream escreveAP2(tmpAP.c_str(), ios_base::out);
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) {
        double vazgasG;
        double presE;
        double tempE;
        double titE;
        double betaE;
        double vazE;
        int indChk;
        imprime = 1;
        if (analisePara.vfp == 1) // analise de sensibilidade para problemas padrao ou para a curva de fundo Eclipse
            analisePara.selecaoAP(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                  sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                  sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
								  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                  sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                  presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                  sistem1.config().correcao.dTdL, imprime);
        else if (analisePara.vfp == 0) // anaslise de sensibilidade para curva de pressao de fundo Imex
            analisePara.selecaoAPImex(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                      sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                      sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
									  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                      sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                      presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                      sistem1.config().correcao.dTdL, imprime);
        if (indfalha[iSeq] > 0)
            escreveAP2 << iSeq << " : " << " Resultado = " << "sucesso" << endl;
        else
            escreveAP2 << iSeq << " : " << " Resultado = " << "falha" << endl;
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    int diaFim = (ltm->tm_mday);
    int horaFim;
    if (diaFim == diaIni)
        horaFim = ltm->tm_hour;
    else
        horaFim = ltm->tm_hour + 24;
    horaFim *= 3600;
    int minutoFim = 60 * ltm->tm_min;
    int segundoFim = ltm->tm_sec;
    int totalFim = horaFim + minutoFim + segundoFim;
    int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << "datahora = ";
    escreveAP2 << ltm->tm_mday << "/";
    escreveAP2 << 1 + ltm->tm_mon << "/";
    escreveAP2 << 1900 + ltm->tm_year << " ";
    escreveAP2 << 0 + ltm->tm_hour << ":";
    escreveAP2 << 0 + ltm->tm_min << ":";
    escreveAP2 << 0 + ltm->tm_sec;
    escreveAP2 << endl;
    escreveAP2 << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
    escreveAP2 << "     Versao    " << versao << endl;

    escreveAP2.close();
}

void leituraAP(string nomeArquivoAP, SProd &sistem1) {
    // metodo em que e feita a sequencia de simulacoes de analise de sensibilidade, observando que a analise de sensibilidade e
    // feita apenas para tramos simples.
    APara analisePara(sistem1.vg1dSP, nomeArquivoAP, sistem1.ncel, sistem1.config().celp,
                      sistem1.config().flup, sistem1.config().bcs,sistem1.config().multiBcs, sistem1.config().fonteg); // construtor do sistema de analise de sensibilidade, ver a classe LerAP. Neste cosntrutor é lido um arquivo
    // json onde se esta registradoas variaveis que estao envolvidas na analise de sensibilidade e os valores que se
    // testara destas variaveis na analise de sensibilidade
    double vazgasG;
    double presE;
    double tempE;
    double titE;
    double betaE;
    double vazE;
    int indChk;
    int imprime = 1;

    ostringstream saidaAP;
    saidaAP << "sucessoAP.dat"; // criacao de arquivo com relatorio de quais pontos de analise de sensibilidade tiveram sucesso e quais falharam

    string tmpAP = saidaAP.str();
    ofstream escreveAP(tmpAP.c_str(), ios_base::out);
    escreveAP << "relatório de falhas da Analise Parametrica para um Tramo " << endl;
    escreveAP.close();
    double falha = 0.;
    for (int iSeq = 0; iSeq < analisePara.nVariaveis; iSeq++) { // laco onde e calculado o permanente de cada combinacao
        // algumas variaveis ja se encontram na condicao correta para se apolicar no sistema de simulacao,
        // outras precisam de um "pos procesaamento"
        if (analisePara.vfp == 1 || analisePara.vfp == 3) // analise de sensibilidade para problemas padrao ou para a curva de fundo Eclipse
            analisePara.selecaoAP(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                  sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                  sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
								  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                  sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                  presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                  sistem1.config().correcao.dTdL, imprime);
        else if (analisePara.vfp == 0 || analisePara.vfp == 2) // anaslise de sensibilidade para curva de pressao de fundo Imex
            analisePara.selecaoAPImex(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                      sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                      sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
									  sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                      sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                      presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                      sistem1.config().correcao.dTdL, imprime);
        // variavei que precisam de um pos processamento para se encaixar na condicao de simulacao
        if (sistem1.config().lingas > 0 && analisePara.listaV.vgasinj == 1 && analisePara.APGasInj.parserieVazGas > 0) {
            // caso tenha linha de gas e analise de sensibilidade para vazao de injecao, a entrada no json e em stdM3,
            // mas no simulador, deve ser atualizada para a vazao massica
            if (sistem1.config().gasinj.tipoCC == 1) {
                sistem1.gasCell(0).massfonteCH = vazgasG * (sistem1.config().flug.Deng * 1.225) / 86400;
            }
        }
        if (analisePara.listaV.vchk == 1 && analisePara.APCHK.parserieAbre > 0) {
            // carregando o valo de abertuda do choke na variavel do sistema
            sistem1.config().chokep.abertura[0] = analisePara.APCHK.abertura[indChk];
        }
        if (sistem1.config().ConContEntrada == 1 && analisePara.listaV.vpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de pressao
            if (analisePara.APpEntrada.parseriePres > 0)
                sistem1.presE = presE;
            if (analisePara.APpEntrada.parserieBet > 0)
                sistem1.betaE = betaE;
            if (analisePara.APpEntrada.parserieTemp > 0)
                sistem1.tempE = tempE;
            if (analisePara.APpEntrada.parserieTit > 0)
                sistem1.titE = titE;
        } else if (sistem1.config().ConContEntrada == 2 && analisePara.listaV.vvazpresent == 1) {
            // carregando os valores de analise de senisbilidadepara condicao de contorno a montante de vazao e pressao
            if (analisePara.APvpEntrada.parseriePres > 0)
                sistem1.presE = presE;
            if (analisePara.APvpEntrada.parserieBet > 0)
                sistem1.betaE = betaE;
            if (analisePara.APvpEntrada.parserieTemp > 0)
                sistem1.tempE = tempE;
            if (analisePara.APvpEntrada.parserieMass > 0) {
                double tit;
                tit = sistem1.config().flup[0].FracMassHidra(sistem1.presE, sistem1.tempE);
                double vazC = vazE * betaE;
                sistem1.cell(0).acsr.injm.MassC = vazC;
                sistem1.cell(0).acsr.injm.MassG = (vazE - vazC) * tit;
                sistem1.cell(0).acsr.injm.MassP = vazE - vazC - sistem1.cell(0).acsr.injm.MassG;
            }
        } else if (sistem1.config().pocinjec == 1 && analisePara.listaV.vpocinj == 1) {
            // carregando os valores de analise de senisbilidade para condicao de contorno a montante de vazao, poco injetor
            sistem1.config().condpocinj.presfundo = presE;
            if (sistem1.config().flashCompleto == 0) {
                sistem1.cell(0).acsr.injl.temp = tempE;
                sistem1.cell(0).acsr.injl.QLiq = vazE;
            } else if (sistem1.config().flashCompleto == 1) {
                sistem1.cell(0).acsr.injg.temp = tempE;
                sistem1.cell(0).acsr.injg.QGas = vazE;
            }
        }
        double chute = -1.;
        cout << "Resolvendo Sequencia " << iSeq << " da analise de sensibilidade" << endl;
        (*sistem1.vg1dSP).sequenciaAP = iSeq;
        if (iSeq == 0 || sistem1.config().ConContEntrada == 2 || sistem1.config().pocinjec == 1) {
            // se for o primeiro caso da analise de sensibilidade, não se usa nenhum valor para o chute, chute=-1
            chute = -1.;
            falha = SolveTramoSolteiro(sistem1, chute);
        } else {
            // para os casos de analise de sensibilidade subsequentes, se usa o valor da solução da analise de sensibilidade anterior como chute
            int tipoChute = sistem1.config().ConContEntrada;
            if (tipoChute == 0) {
                chute = sistem1.cell(0).pres;
            } else {
                if (sistem1.cell(0).acsr.tipo == 2)
                    chute = sistem1.cell(0).acsr.injl.QLiq;
                else if (sistem1.cell(0).acsr.tipo == 1)
                    chute = sistem1.cell(0).acsr.injg.QGas;
            }
            if (sistem1.config().lingas == 1 && sistem1.gasCell(0).tipoCC == 0) {
                sistem1.config().gasinj.chuteVaz = 1;
                sistem1.config().gasinj.vazgas[0] = sistem1.gasCell(0).VGasR * 86400. / (sistem1.gasCell(0).flui.MasEspGas(1., 15.6));
            }
            if (iSeq == 153 || iSeq == 60) {
            }
            if (fabs(falha) > 1e9) {
                chute = -1;
            }
            falha = SolveTramoSolteiro(sistem1, chute);
            // caso a solucao com chute tenha falhado, tenta-se uma nova solucao sem o chute
            if (fabs(falha) > 1e9 && chute > 0) {
                falha = SolveTramoSolteiro(sistem1);
            }
        }
        if (fabs(falha) > 0.9e10) {
            // caso ocorre uma falha nma busca da solucao para um caso da analise de sensibilidade
            cout << "----------------------------- falha no caso -------------------------------" << endl;
        }

        ostringstream saidaAP;
        saidaAP << "sucessoAP.dat"; // arquivo com o relatorio de sucessos e falhas da AP

        string tmpAP = saidaAP.str();
        ofstream escreveAP(tmpAP.c_str(), ios_base::app);
        if (fabs(falha) < 1e9)
            escreveAP << iSeq << " : " << " Resultado = " << "sucesso" << endl;
        else
            escreveAP << iSeq << " : " << " Resultado = " << "falha" << endl;

        escreveAP.close();

        if (analisePara.tipoAP == 0 && analisePara.imprimePerfil==1) {
            // impressão dos perfis e tendencias da analise de sensibilidade, caso sem construcao de tabela de pressao de fundo
            sistem1.config().imprimeProfile(sistem1.celula, sistem1.flut, 0, sistem1.indTramo);
            sistem1.config().resumoPermanente(sistem1.celula, sistem1.celulaG, sistem1.pGSup, sistem1.presiniG, sistem1.indTramo);
            if(sistem1.config().nintermi>0)sistem1.config().resumoIntermitencia(sistem1.celula, sistem1.indTramo);
            if(sistem1.config().nCelUnit>0){
            	for(int iCelUni=0; iCelUni<sistem1.config().nCelUnit; iCelUni++)
            	sistem1.config().relatorioCelulaUnitaria(sistem1.celula,sistem1.config().celUnit[iCelUni].posicP, sistem1.indTramo);
            }
            // enterramento
            for (int j = 0; j <= sistem1.ncel; j++) {
                if (sistem1.cell(j).calor.difus2D == 1) {
                    sistem1.cell(j).calor.poisson2D.imprimePermanente(sistem1.indTramo);
                }
            }
            sistem1.kimpT = 1;
            for (int i = 0; i < sistem1.config().ntendp; i++) {
                if (iSeq == 0)
                    sistem1.ImprimeTrendPCab(i);
                sistem1.config().imprimeTrend(sistem1.celula, sistem1.MatTrendP[i], 0, i, 0);
                sistem1.ImprimeTrendP(i);
            }
            if (sistem1.config().lingas == 1) {
                sistem1.config().imprimeProfileG(sistem1.celulaG, sistem1.flutG, 0, sistem1.indTramo);
                for (int i = 0; i < sistem1.config().ntendg; i++) {
                    if (iSeq == 0)
                        sistem1.ImprimeTrendGCab(i);
                    sistem1.config().imprimeTrendG(sistem1.celulaG, sistem1.MatTrendG[i], 0, i, 0, 0);
                    sistem1.ImprimeTrendG(i);
                }
            }
        } else {
            // caso de construcao de tabela VFP:
            double BHP;
            if (fabs(falha) < 1e9)
                BHP = sistem1.cell(0).pres;
            else
                BHP = -1e10;
            analisePara.tabelaGenerica(sistem1.ncelGas, sistem1.chokeSup, sistem1.celula, sistem1.celulaG, sistem1.config().flup,
                                       sistem1.config().IPRS, sistem1.config().valv, sistem1.config().fonteg,
                                       sistem1.config().fontel, sistem1.config().fontem, sistem1.config().furo, sistem1.config().bcs,
									   sistem1.config().multiBcs,sistem1.config().bvol, sistem1.config().dpreq,
                                       sistem1.pGSup, sistem1.temperatura, sistem1.presiniG, sistem1.tempiniG, vazgasG,
                                       presE, tempE, titE, betaE, vazE, iSeq, indChk, sistem1.config().correcao.dPdLHidro, sistem1.config().correcao.dPdLFric,
                                       sistem1.config().correcao.dTdL, BHP);
        }
        (*sistem1.vg1dSP).contaExit = 0;
    }

    // saidaAP << "sucessoAP.dat";//arquivo com o relatorio de sucessos e falhas da AP

    string tmpAP2 = saidaAP.str();
    ofstream escreveAP2(tmpAP2.c_str(), ios_base::app);

    time_t now = time(0);
    tm *ltm = localtime(&now);
    int diaFim = (ltm->tm_mday);
    int horaFim;
    if (diaFim == diaIni)
        horaFim = ltm->tm_hour;
    else
        horaFim = ltm->tm_hour + 24;
    horaFim *= 3600;
    int minutoFim = 60 * ltm->tm_min;
    int segundoFim = ltm->tm_sec;
    int totalFim = horaFim + minutoFim + segundoFim;
    int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << endl;
    escreveAP2 << "datahora = ";
    escreveAP2 << ltm->tm_mday << "/";
    escreveAP2 << 1 + ltm->tm_mon << "/";
    escreveAP2 << 1900 + ltm->tm_year << " ";
    escreveAP2 << 0 + ltm->tm_hour << ":";
    escreveAP2 << 0 + ltm->tm_min << ":";
    escreveAP2 << 0 + ltm->tm_sec;
    escreveAP2 << endl;
    escreveAP2 << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
    escreveAP2 << "     Versao    " << versao << endl;

    escreveAP2.close();
}

int testaBloqueio(Rede &arqRede, int i, vector<tramoPart> &bloq) {
    int fontePart = 0;
    int tot = 0;
    tramoPart temp;
    for (int j = 0; j < 20; j++)
        temp.part[j] = -1;
    if (arqRede.malha[i].nafluente > 0) {
        if (arqRede.malha[arqRede.malha[i].afluente[0]].ncoleta == 2) {
            for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
                int iaflu = arqRede.malha[i].afluente[j];
                if (arqRede.malha[iaflu].bloqueio >= 0 && arqRede.malha[iaflu].bloqueio != i) {
                    fontePart = 1;
                    temp.part[j] = iaflu;
                    tot++;
                } else if (arqRede.malha[iaflu].bloqueio < 0) {
                    temp.part[j] = iaflu;
                    tot++;
                }
            }
        } else {
            for (int j = 0; j < arqRede.malha[i].nafluente; j++) {
                int iaflu = arqRede.malha[i].afluente[j];
                temp.part[tot] = iaflu;
                tot++;
            }
        }
    }
    temp.ntram = tot;
    bloq.push_back(temp);
    return fontePart;
}

void avaliaPerm(SProd *malha, Rede &arqRede, int narq, int &semPerm) {
    Vcr<int> resolvPerm(narq, 0);
    int totPerm = 0;
    while (totPerm < narq) {
        for (int i = 0; i < narq; i++) {
            if (arqRede.malha[i].perm == 0 && resolvPerm[i] == 0) {
                resolvPerm[i] = 1;
                totPerm++;
                semPerm++;
            } else {
                if (arqRede.malha[i].nafluente == 0 && resolvPerm[i] == 0) {
                    if (malha[i].cell(0).acsr.tipo == 1 && fabs(malha[i].cell(0).acsr.injg.QGas) < 1e-15) {
                        if (malha[i].config().ConContEntrada != 2 || fabs(malha[i].config().CCVPres.mass[0]) < 1e-15) {
                            arqRede.malha[i].perm = 0;
                            resolvPerm[i] = 1;
                            totPerm++;
                            semPerm++;
                        }
                    } else if (malha[i].cell(0).acsr.tipo == 2 && fabs(malha[i].cell(0).acsr.injl.QLiq) < 1e-15) {
                        if (malha[i].config().ConContEntrada != 2 || fabs(malha[i].config().CCVPres.mass[0]) < 1e-15) {
                            arqRede.malha[i].perm = 0;
                            resolvPerm[i] = 1;
                            totPerm++;
                            semPerm++;
                        }
                    } else if (malha[i].cell(0).acsr.tipo == 10 &&
                               fabs(malha[i].cell(0).acsr.injm.MassC + malha[i].cell(0).acsr.injm.MassP +
                                    malha[i].cell(0).acsr.injm.MassG) < 1e-15) {
                        if (malha[i].config().ConContEntrada != 2 || fabs(malha[i].config().CCVPres.mass[0]) < 1e-15) {
                            arqRede.malha[i].perm = 0;
                            resolvPerm[i] = 1;
                            totPerm++;
                            semPerm++;
                        }
                    } else {
                        resolvPerm[i] = 1;
                        totPerm++;
                    }
                } else if (resolvPerm[i] == 0) {
                    int afluPerm = 0;
                    int ativoAflu = 0;
                    for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                        int indaflu = arqRede.malha[i].afluente[iaflu];
                        afluPerm += resolvPerm[indaflu];
                        ativoAflu += arqRede.malha[indaflu].perm;
                    }
                    if (afluPerm == arqRede.malha[i].nafluente) {
                        if (ativoAflu > 0) {
                            if (malha[i].config().ConContEntrada == 2 && fabs(malha[i].config().CCVPres.mass[0]) <= 1e-15) {
                                arqRede.malha[i].perm = 0;
                                resolvPerm[i] = 1;
                                totPerm++;
                                semPerm++;
                            } else {
                                resolvPerm[i] = 1;
                                totPerm++;
                            }
                        } else {

                            if (malha[i].cell(0).acsr.tipo == 1 && fabs(malha[i].cell(0).acsr.injg.QGas) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].cell(1).acsr.tipo == 1 && fabs(malha[i].cell(1).acsr.injg.QGas) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                                if ((malha[i].cell(1).acsr.injg.QGas) > 1e-15) {
                                    malha[i].cell(0).acsr.injg = malha[i].cell(1).acsr.injg;
                                    malha[i].cell(1).acsr.tipo = 0;
                                    malha[i].cell(0).acsr.tipo = 1;
                                    if (malha[i].config().flashCompleto == 2 && malha[i].config().tabelaDinamica == 1) {
                                        malha[i].tabDin[0].celFim = malha[i].tabDin[1].celFim;
                                        malha[i].tabDin.erase(malha[i].tabDin.begin() + 1);
                                        malha[i].ntabDin--;
                                    }
                                }
                            } else if (malha[i].cell(0).acsr.tipo == 2 && fabs(malha[i].cell(0).acsr.injl.QLiq) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].cell(1).acsr.tipo == 2 && fabs(malha[i].cell(1).acsr.injl.QLiq) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                                if ((malha[i].cell(1).acsr.injl.QLiq) > 1e-15) {
                                    malha[i].cell(0).acsr.injl = malha[i].cell(1).acsr.injl;
                                    malha[i].cell(1).acsr.tipo = 0;
                                    malha[i].cell(0).acsr.tipo = 2;
                                    if (malha[i].config().flashCompleto == 2 && malha[i].config().tabelaDinamica == 1) {
                                        malha[i].tabDin[0].celFim = malha[i].tabDin[1].celFim;
                                        malha[i].tabDin.erase(malha[i].tabDin.begin() + 1);
                                        malha[i].ntabDin--;
                                    }
                                }
                            } else if (malha[i].cell(0).acsr.tipo == 10 &&
                                       fabs(malha[i].cell(0).acsr.injm.MassC + malha[i].cell(0).acsr.injm.MassP +
                                            malha[i].cell(0).acsr.injm.MassG) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].cell(1).acsr.tipo == 10 &&
                                       fabs(malha[i].cell(1).acsr.injm.MassC + malha[i].cell(1).acsr.injm.MassP +
                                            malha[i].cell(1).acsr.injm.MassG) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                                if ((malha[i].cell(1).acsr.injm.MassC + malha[i].cell(1).acsr.injm.MassP +
                                     malha[i].cell(1).acsr.injm.MassG) > 1e-15) {
                                    malha[i].cell(0).acsr.injm = malha[i].cell(1).acsr.injm;
                                    malha[i].cell(0).acsr.tipo = 10;
                                    malha[i].cell(1).acsr.tipo = 0;
                                    if (malha[i].config().flashCompleto == 2 && malha[i].config().tabelaDinamica == 1) {
                                        malha[i].tabDin[0].celFim = malha[i].tabDin[1].celFim;
                                        malha[i].tabDin.erase(malha[i].tabDin.begin() + 1);
                                        malha[i].ntabDin--;
                                    }
                                }
                            } else if (malha[i].cell(0).acsr.tipo == 3) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].cell(0).acsr.tipo == 15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].cell(0).acsr.tipo == 16) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else if (malha[i].config().ConContEntrada == 2 && fabs(malha[i].config().CCVPres.mass[0]) > 1e-15) {
                                resolvPerm[i] = 1;
                                totPerm++;
                            } else {
                                arqRede.malha[i].perm = 0;
                                resolvPerm[i] = 1;
                                totPerm++;
                                semPerm++;
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].nafluente > 0 && arqRede.malha[i].perm == 1) {
            int indaflu = arqRede.malha[i].afluente[0];
            if (arqRede.malha[indaflu].ncoleta > 1) {
                int afluPerm = 0;
                int colPerm = 0;
                for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                    int indaflu2 = arqRede.malha[i].afluente[iaflu];
                    afluPerm += arqRede.malha[indaflu2].perm;
                }
                for (int icol = 0; icol < arqRede.malha[indaflu].ncoleta; icol++) {
                    int indacol = arqRede.malha[indaflu].coleta[icol];
                    colPerm += arqRede.malha[indacol].perm;
                    if (malha[indacol].config().ConContEntrada == 2 && fabs(malha[indacol].config().CCVPres.mass[0]) > 1e-15) {
                    }
                }

                if (afluPerm == 0 && colPerm > 1) {
                    arqRede.malha[indaflu].perm = 1;
                    malha[indaflu].cell(0).acsr.tipo = malha[i].cell(0).acsr.tipo;
                    if (malha[i].cell(0).acsr.tipo == 1) {
                        malha[indaflu].cell(0).acsr.injg = malha[i].cell(0).acsr.injg;
                        malha[indaflu].cell(0).acsr.injg.QGas = 500.;
                        malha[indaflu].cell(0).acsr.injg.razCompGas = 0.;
                        malha[indaflu].cell(0).acsr.injg.temp = malha[i].cell(0).temp;
                    } else if (malha[i].cell(0).acsr.tipo == 2) {
                        malha[indaflu].cell(0).acsr.injl = malha[i].cell(0).acsr.injl;
                        malha[indaflu].cell(0).acsr.injl.QLiq = 5.;
                        malha[indaflu].cell(0).acsr.injl.bet = 0.;
                        malha[indaflu].cell(0).acsr.injl.temp = malha[i].cell(0).temp;
                    } else if (malha[i].cell(0).acsr.tipo == 10) {
                        malha[indaflu].cell(0).acsr.injm = malha[i].cell(0).acsr.injm;
                        malha[indaflu].cell(0).acsr.injm.MassP = 0.01;
                        malha[indaflu].cell(0).acsr.injm.MassC = 0.0;
                        malha[indaflu].cell(0).acsr.injm.MassG = 0.004;
                        malha[indaflu].cell(0).acsr.injm.temp = malha[i].cell(0).temp;
                    }
                }
            }
            for (int iaflu = 0; iaflu < arqRede.malha[i].nafluente; iaflu++) {
                int indaflu = arqRede.malha[i].afluente[iaflu];
                int nQG = 0;
                int nQL = 0;
                int nQM = 0;
                double qTotG = 0.;
                double qTotL = 0.;
                double qTotM = 0.;
                if (arqRede.malha[indaflu].perm == 1 && malha[indaflu].cell(0).acsr.tipo == 1) {
                    nQG++;
                    qTotG += malha[indaflu].cell(0).acsr.injg.QGas;
                } else if (arqRede.malha[indaflu].perm == 1 && malha[indaflu].cell(0).acsr.tipo == 2) {
                    nQL++;
                    qTotL += malha[indaflu].cell(0).acsr.injl.QLiq;
                } else if (arqRede.malha[indaflu].perm == 1 && malha[indaflu].cell(0).acsr.tipo == 10) {
                    nQM++;
                    qTotM += (malha[indaflu].cell(0).acsr.injm.MassP + malha[indaflu].cell(0).acsr.injm.MassC +
                              malha[indaflu].cell(0).acsr.injm.MassG);
                }
                if (nQG == arqRede.malha[i].nafluente && malha[i].cell(1).acsr.tipo == 1) {
                    if ((malha[i].cell(1).acsr.injg.QGas + qTotG) < 1e-15)
                        malha[i].cell(1).acsr.injg.QGas += 500.;
                } else if (nQL == arqRede.malha[i].nafluente && malha[i].cell(1).acsr.tipo == 2) {
                    if ((malha[i].cell(1).acsr.injl.QLiq + qTotL) < 1e-15)
                        malha[i].cell(1).acsr.injl.QLiq += 5.;
                } else if (nQM == arqRede.malha[i].nafluente && malha[i].cell(1).acsr.tipo == 10) {
                    double totM = (malha[i].cell(0).acsr.injm.MassP + malha[i].cell(0).acsr.injm.MassC +
                                   malha[i].cell(0).acsr.injm.MassG);
                    if ((totM + qTotL) < 1e-15) {
                        malha[i].cell(1).acsr.injm.MassP += 0.01;
                        malha[i].cell(1).acsr.injm.MassG += 0.004;
                    }
                }
            }
        }
    }
}

void preparaRedeProd(SProd *malha, Rede &arqRede, int narq, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
                     tipoSimulacao_t tipoSimulacaoMrt, int nrede, int &contapermRede, varGlob1D *vG1d) {
    if (narq > 1) { // narq=numero de tramos ativos da rede
        (*arqRede.vg1dSP).simulaTransiente = arqRede.chaveredeT;
        for (int i = 0; i < narq; i++) { // laco em que se constroe os tramos da rede
            vector<tramoPart> bloq;
            int bloqR = testaBloqueio(arqRede, i, bloq); // verifica se em uma rede de manifold existe uma configuracao de derivacao especial
            // em que parte dos tramos coletados no manifold e desviado para a linha de teste
            int tD = arqRede.tabelaDinamica; // rede com modelo composicional e uso de tabela dinamica
            if ((*arqRede.vg1dSP).fluidoRede == 0)
                (*arqRede.vg1dSP).tipoFluidoRedeGlob = 1; // tipoFluidoRede==1 implica em que a rede será uma rede dominada pela fase líquida
            // construtor do tramo temporario com as caracteristicas basicas do tramo i da rede
            SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                             tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD, bloqR, 0, arqRede.malha[i].reverso, 0, 0, 0, arqRede.malha[i].perm);

            if (arqRede.malha[i].ncoleta > 0)
                temporario.noextremo = 0; // caso o tramo não esteja no fim da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            else
                temporario.noextremo = 1; // caso o tramo esteja no fim da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            if (arqRede.malha[i].nafluente > 0)
                temporario.noinicial = 0; // caso o tramo não esteja no inicio da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            else
                temporario.noinicial = 1; // caso o tramo esteja no inicio da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            malha[i] = temporario; // carregando o tramo temporario no tramo i da malha
            malha[i].indTramo = i;
            if (malha[i].config().tabelaDinamica == 1) {
                // caso a opcao de tabela dinamica esteja ligada, como, para o modelo composicional, sempre, primeiramente, se faz um
                // calculo inicial com modelo black oil, a opcao tabela dinamica e desligada temporariamente enquanto se esta no modo blackoil
                for (int j = 0; j <= malha[i].ncel; j++) {
                    malha[i].cell(j).flui.tabelaDinamica = 0;
                    if (malha[i].cell(j).acsr.tipo == 1) {
                        malha[i].cell(j).acsr.injg.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 2) {
                        malha[i].cell(j).acsr.injl.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 3) {
                        malha[i].cell(j).acsr.ipr.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 15) {
                        malha[i].cell(j).acsr.radialPoro.flup.tabelaDinamica = 0;
                        int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                        for (int k = 0; k < ncelRad; k++) {
                            malha[j].cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
                        }
                    } else if (malha[i].cell(j).acsr.tipo == 16) {
                        malha[i].cell(j).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
                        int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                        for (int k = 0; k < ncelRad; k++) {
                            malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
                        }
                        int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                        for (int k = 0; k < ncelEle; k++) {
                            malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
                        }
                    } else if (malha[i].cell(j).acsr.tipo == 10) {
                        malha[i].cell(j).acsr.injm.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 9 && malha[i].cell(j).acsr.fontechk.abertura > 1e-6) {
                        malha[i].cell(j).acsr.fontechk.fluidoP.tabelaDinamica = 0;
                        malha[i].cell(j).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
                    }
                }
            }
            if (arqRede.malha[i].nafluente == 0 && malha[i].config().ConContEntrada == 1) {
                // preparando o tramo para o caso em que este esta no inicio da rede e tem
                // como condicao de contorno e pressao na entrada do tramo, neste caso, tem de se ter uma estimativa inicial de vazao
                // neste tramo, alem disto, e necessario definir detalhes como o beta in-situ e o RGO
                // observar tambem que no solver permanente, e necessario para este tipo de condicao de contorno
                // e necessario colocar um acessorio de fonte, onde e armazendao a estimativa de vazao no inicio do tramo
                // caso o fluido seja dominado por liquido, a fointe de armazenamento sera do tipo liquido
                //  se gas, a fonte sera de gas
                double tit = malha[i].config().CCPres.tit[0];
                /// atencao, verificar se aposicao do if esta correta.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ((*arqRede.vg1dSP).fluidoRede == 0) {
                    if (arqRede.chutHol > 0.1)
                        arqRede.chutHol = 0.1;
                    if (arqRede.relax > 0.5)
                        arqRede.relax = 0.5;
                }
                ProFlu fluCP;
                if (malha[i].config().CCPres.indFluido == -1)
                    fluCP = malha[i].cell(0).flui;
                else
                    fluCP = malha[i].config().flup[malha[i].config().CCPres.indFluido];

                if (tit < (1. - (*arqRede.vg1dSP).localtiny)) {

                    malha[i].tempE = malha[i].config().CCPres.temperatura[0];
                    malha[i].presE = malha[i].config().CCPres.pres[0];
                    double rgoentrada = 0.;
                    double rcst = malha[i].cell(0).fluicol.MasEspFlu(1.01, 15);
                    double rcis = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                    malha[i].betaE = malha[i].config().CCPres.bet[0] * rcst / rcis;

                    if (malha[i].config().flashCompleto != 2) {
                        double rgST = fluCP.Deng * 1.225;
                        double roST = 141.5 * 1000. / (131.5 + fluCP.API);
                        double rl = fluCP.MasEspoleo(malha[i].presE, malha[i].tempE);
                        double rw = fluCP.Denag * 1000.;
                        double RSent = fluCP.RS(malha[i].presE, malha[i].tempE) * 6.29 / 35.31467;
                        double bsw = fluCP.BSW;
                        double bo = fluCP.BOFunc(malha[i].presE, malha[i].tempE);
                        double ba = fluCP.BAFunc(malha[i].presE, malha[i].tempE);
                        double fw = bsw * ba / (bo + ba * bsw - bsw * bo);

                        double Ost = (1 - malha[i].betaE) * bo * (1 - fw);
                        double Wst = fw * (1 - malha[i].betaE);
                        double Cst = malha[i].betaE;
                        malha[i].titE = malha[i].config().CCPres.tit[0];
                        if (malha[i].titE < 0) {
                            double titH = fluCP.FracMassHidra(malha[i].presE, malha[i].tempE);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double rlMix = malha[i].betaE * rc + (1. - malha[i].betaE) * rl;
                            double val1 = ((1. - malha[i].betaE) * rl * titH / (1. - titH));
                            malha[i].titE = val1 / (rlMix + val1);
                        }

                        rgoentrada = RSent * rgST * Ost + malha[i].titE * roST * Ost + malha[i].titE * rw * Wst +
                                     malha[i].titE * rcst * Cst;
                        rgoentrada /= ((1 - malha[i].titE) * rgST * Ost);
                    }

                    if ((*arqRede.vg1dSP).fluidoRede == 1) { // dominado por liquido
                        if (rgoentrada > 1400.)
                            rgoentrada = 1400.;
                        malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                        malha[i].cell(0).acsr.injl.bet = malha[i].betaE;
                        malha[i].cell(0).acsr.injl.temp = malha[i].config().CCPres.temperatura[0];
                        malha[i].cell(0).acsr.injl.FluidoPro.BSW = fluCP.BSW;
                        malha[i].cell(0).acsr.injl.FluidoPro.TempL = fluCP.TempL;
                        malha[i].cell(0).acsr.injl.FluidoPro.TempH = fluCP.TempH;
                        malha[i].cell(0).acsr.injl.FluidoPro.LVisL = fluCP.LVisL;
                        malha[i].cell(0).acsr.injl.FluidoPro.LVisH = fluCP.LVisH;
                        if (malha[i].config().flashCompleto != 2) {
                            malha[i].cell(0).acsr.injl.FluidoPro.RGO = rgoentrada;
                            malha[i].cell(0).acsr.injl.FluidoPro.Deng = fluCP.Deng;
                            malha[i].cell(0).acsr.injl.FluidoPro.Denag = fluCP.Denag;
                            malha[i].cell(0).acsr.injl.FluidoPro.yco2 = fluCP.yco2;
                            malha[i].cell(0).acsr.injl.FluidoPro.API = fluCP.API;
                            malha[i].cell(0).acsr.injl.FluidoPro.RenovaFluido();
                        } else {
                            fluCP.atualizaPropCompStandard();
                            fluCP.atualizaPropComp(malha[i].presE, malha[i].tempE);
                            double rl = fluCP.MasEspLiq(malha[i].presE, malha[i].tempE);
                            double titH = fluCP.FracMassHidra(malha[i].presE, malha[i].tempE);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double rlMix = malha[i].betaE * rc + (1. - malha[i].betaE) * rl;
                            double val1 = ((1. - malha[i].betaE) * rl * titH / (1. - titH));
                            malha[i].titE = val1 / (rlMix + val1);
                            malha[i].cell(0).acsr.injl.FluidoPro = fluCP;
                        }
                        if (fabs(malha[i].cell(0).acsr.injl.QLiq) < 1e-15)
                            malha[i].cell(0).acsr.injl.QLiq = 1000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                        malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injl.FluidoPro;
                    } else {
                        malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                        malha[i].cell(0).acsr.injg.seco = 0;
                        malha[i].cell(0).acsr.injg.temp = malha[i].config().CCPres.temperatura[0];
                        malha[i].cell(0).acsr.injg.FluidoPro.BSW = fluCP.BSW;
                        malha[i].cell(0).acsr.injg.FluidoPro.TempL = fluCP.TempL;
                        malha[i].cell(0).acsr.injg.FluidoPro.TempH = fluCP.TempH;
                        malha[i].cell(0).acsr.injg.FluidoPro.LVisL = fluCP.LVisL;
                        malha[i].cell(0).acsr.injg.FluidoPro.LVisH = fluCP.LVisH;
                        if (malha[i].config().flashCompleto != 2) { // caso nao composicional
                            double rgST = fluCP.Deng * 1.225;
                            double rg = fluCP.MasEspGas(malha[i].presE, malha[i].tempE);
                            double rl = fluCP.MasEspoleo(malha[i].presE, malha[i].tempE);
                            double tit = fluCP.FracMass(malha[i].presE, malha[i].tempE);
                            double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                            double val2 = (rg / rl) * (1 - tit) / tit;
                            double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                            malha[i].titE = titT;
                            malha[i].betaE = val1 / (val2 + val1);
                            rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                            if (rgoentrada > 1400.)
                                rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                            fluCP.RGO = malha[i].cell(0).acsr.injg.FluidoPro.RGO = rgoentrada;
                            malha[i].cell(0).acsr.injg.FluidoPro.Deng = fluCP.Deng;
                            malha[i].cell(0).acsr.injg.FluidoPro.Denag = fluCP.Denag;
                            malha[i].cell(0).acsr.injg.FluidoPro.yco2 = fluCP.yco2;
                            malha[i].cell(0).acsr.injg.FluidoPro.API = fluCP.API;
                            malha[i].cell(0).acsr.injg.FluidoPro.BSW = fluCP.BSW;
                            malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                        } else { // caso composicional
                            fluCP.atualizaPropCompStandard();
                            fluCP.atualizaPropComp(malha[i].presE, malha[i].tempE);
                            double rgST = fluCP.Deng * 1.225;
                            double rg = fluCP.MasEspGas(malha[i].presE, malha[i].tempE);
                            double rl = fluCP.MasEspoleo(malha[i].presE, malha[i].tempE);
                            double tit = fluCP.FracMass(malha[i].presE, malha[i].tempE);
                            double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                            double val2 = (rg / rl) * (1 - tit) / tit;
                            double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                            malha[i].titE = titT;
                            malha[i].betaE = val1 / (val2 + val1);
                            malha[i].cell(0).acsr.injg.FluidoPro = fluCP;
                        }
                        malha[i].cell(0).acsr.injg.fluidocol = malha[i].config().fluc;
                        malha[i].cell(0).acsr.injg.razCompGas = malha[i].config().CCPres.bet[0];
                        if (malha[i].cell(0).acsr.injg.QGas < 1e-15)
                            malha[i].cell(0).acsr.injg.QGas = 1000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                        malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injg.FluidoPro;
                    }
                } else { // dominado por gas
                    double rgoentrada = 0.;
                    malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                    malha[i].cell(0).acsr.injg.temp = malha[i].config().CCPres.temperatura[0];
                    malha[i].cell(0).acsr.injg.seco = 0;
                    malha[i].cell(0).acsr.injg.FluidoPro.TempL = fluCP.TempL;
                    malha[i].cell(0).acsr.injg.FluidoPro.TempH = fluCP.TempH;
                    malha[i].cell(0).acsr.injg.FluidoPro.LVisL = fluCP.LVisL;
                    malha[i].cell(0).acsr.injg.FluidoPro.LVisH = fluCP.LVisH;
                    if (malha[i].config().flashCompleto != 2) { // caso nao composicional
                        double rgST = fluCP.Deng * 1.225;
                        double rg = fluCP.MasEspGas(malha[i].presE, malha[i].tempE);
                        double rl = fluCP.MasEspoleo(malha[i].presE, malha[i].tempE);
                        double tit = fluCP.FracMass(malha[i].presE, malha[i].tempE);
                        double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                        double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                        double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                        double val2 = (rg / rl) * (1 - tit) / tit;
                        double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                        malha[i].titE = titT;
                        malha[i].betaE = val1 / (val2 + val1);
                        rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                        if (rgoentrada > 1400.)
                            rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                        fluCP.RGO = malha[i].cell(0).flui.RGO = malha[i].cell(0).acsr.injg.FluidoPro.RGO = rgoentrada;
                        malha[i].cell(0).acsr.injg.FluidoPro.Deng = fluCP.Deng;
                        malha[i].cell(0).acsr.injg.FluidoPro.Denag = fluCP.Denag;
                        malha[i].cell(0).acsr.injg.FluidoPro.yco2 = fluCP.yco2;
                        malha[i].cell(0).acsr.injg.FluidoPro.API = fluCP.API;
                        malha[i].cell(0).acsr.injg.FluidoPro.BSW = fluCP.BSW;
                        malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                    } else { // caso composicional
                        fluCP.atualizaPropCompStandard();
                        fluCP.atualizaPropComp(malha[i].presE, malha[i].tempE);
                        double rgST = fluCP.Deng * 1.225;
                        double rg = fluCP.MasEspGas(malha[i].presE, malha[i].tempE);
                        double rl = fluCP.MasEspoleo(malha[i].presE, malha[i].tempE);
                        double tit = fluCP.FracMass(malha[i].presE, malha[i].tempE);
                        double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                        double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                        double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                        double val2 = (rg / rl) * (1 - tit) / tit;
                        double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                        malha[i].titE = titT;
                        malha[i].betaE = val1 / (val2 + val1);
                        malha[i].cell(0).acsr.injg.FluidoPro = fluCP;
                    }
                    malha[i].cell(0).acsr.injg.fluidocol = malha[i].config().fluc;
                    malha[i].cell(0).acsr.injg.razCompGas = malha[i].config().CCPres.bet[0];
                    if (malha[i].cell(0).acsr.injg.QGas < 1e-15)
                        malha[i].cell(0).acsr.injg.QGas = 10000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                    malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injl.FluidoPro;
                }
            }
            // fazendo uma estimativa inicial da vazao total de liquido e gas que escoara na rede, para tanto,
            // se buscara as fontes de liquido e gas que estao computadas no inicio do tramo, observar que e uma estimativa imperfeita, pois
            // no caso de condicao de contorno pressao, esta vazao e apenas um chuet, alem de nao contabilizar a vazao em celulas que não a
            // inicial do tramo
            if (malha[i].cell(0).acsr.tipo == 2 && arqRede.malha[i].perm == 1)
                (*arqRede.vg1dSP).somavaz += malha[i].cell(0).acsr.injl.QLiq / 86400.;
            else if (malha[i].cell(0).acsr.tipo == 1 && arqRede.malha[i].perm == 1)
                (*arqRede.vg1dSP).somavazG += malha[i].cell(0).acsr.injg.QGas / 86400.;
            else if (malha[i].cell(0).acsr.tipo == 10 && arqRede.malha[i].perm == 1 && malha[i].config().ConContEntrada != 2) {
                double massTot = malha[i].cell(0).acsr.injm.MassP + malha[i].cell(0).acsr.injm.MassG;
                double rhopstd = 141.5 * 1000. / (malha[i].cell(0).acsr.injm.FluidoPro.API + 131.5);
                double bsw = malha[i].cell(0).acsr.injm.FluidoPro.BSW;
                double rhogstd = malha[i].cell(0).acsr.injm.FluidoPro.Deng * 1.225;
                double tit = malha[i].cell(0).acsr.injm.FluidoPro.FracMassHidra(1., 20.);
                (*arqRede.vg1dSP).somavaz += massTot * (1 - tit) * ((1 - bsw) / rhopstd + bsw / 1000.) +
                                             malha[i].cell(0).acsr.injm.MassC / malha[i].cell(0).acsr.injm.fluidocol.rholStd;
                (*arqRede.vg1dSP).somavazG += massTot * tit / rhogstd;
            }
            (*arqRede.vg1dSP).somaarea += malha[i].cell(0).duto.area; // somaarea sera utilizado para fazer uma distribuicao inicial da vazao em cada tramo, tomando como
            // referencia a razao entre a area do tramo e somaarea, esta razao sera aplicada para distribuir o chute inicial da vazao de cada tramo,
            // tomando somavaz e somavazG e aplicando a razao de area
            if (arqRede.malha[i].perm == 0)
                contapermRede++; // contabilizando se o tramo esta com o permanente ligado, o numero de tramops dan rede
            // pode ser diferente do numero de tramso que efetivamente poarticipara do solver permanente da rede
            temporario.config().tabelaDinamica = 0;
        }
        for (int i = 0; i < narq; i++) {
            malha[i].chuteHol = arqRede.chutHol; // o chutehol serve como uma estimativa inicial do holdup na rede para o calculo inicial
            // da variacao de pressao em cada tramo da rede devido aa hidrostatica
            if (arqRede.malha[i].ncoleta > 0) {
                int icol = arqRede.malha[i].coleta[0];
                malha[i].fluiRevRede = malha[icol].cell(0).flui; // computando o fluido do tramo coletor, caso no afluente
                // ocorra reversao do escoamento pre defindio, o fluido que escoarah no tramo afluente vira do fluido que se encontra no coletor
                if (malha[i].config().flashCompleto == 2)
                    malha[i].fluiRevRede.atualizaPropCompStandard();
            }
        }
    } else {
        int i = 0;
        int tD = arqRede.tabelaDinamica;
        SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                         tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD);
        temporario.noextremo = 1;
        temporario.noinicial = 1;
        malha[i] = temporario;
        malha[i].indTramo = i;

        if (malha[i].config().tabelaDinamica == 1) {
            for (int j = 0; j <= malha[i].ncel; j++) {
                malha[i].cell(j).flui.tabelaDinamica = 0;
                if (malha[i].cell(j).acsr.tipo == 1) {
                    malha[i].cell(j).acsr.injg.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 2) {
                    malha[i].cell(j).acsr.injl.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 3) {
                    malha[i].cell(j).acsr.ipr.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 15) {
                    malha[i].cell(j).acsr.radialPoro.flup.tabelaDinamica = 0;
                    int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                    for (int k = 0; k < ncelRad; k++) {
                        malha[j].cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
                    }
                } else if (malha[i].cell(j).acsr.tipo == 16) {
                    malha[i].cell(j).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
                    int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int k = 0; k < ncelRad; k++) {
                        malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
                    }
                    int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                    for (int k = 0; k < ncelEle; k++) {
                        malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
                    }
                } else if (malha[i].cell(j).acsr.tipo == 10) {
                    malha[i].cell(j).acsr.injm.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 9 && malha[i].cell(j).acsr.fontechk.abertura > 1e-6) {
                    malha[i].cell(j).acsr.fontechk.fluidoP.tabelaDinamica = 0;
                    malha[i].cell(j).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
                }
            }
        }
    }
}

void solveRedeProd(SProd *malha, Rede &arqRede, int narq,
                   Vcr<int> &inativo, int indativo, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
                   tipoSimulacao_t tipoSimulacaoMrt, int nrede, int &contapermRede, varGlob1D *vG1d, vector<noRede> &normaEvol, vector<tramoPart> &bloq) {
    int convergido = 0;
    if (narq > 1) {
        for (int i = 0; i < narq; i++)
            (void)testaBloqueio(arqRede, i, bloq);
        int semPerm = 0;
        // preprocessamento final, em avaliaPerm faz-se um "pente fino" em tramos sem vazao, estes são descartados da rede
        // para o caluculo permanente, semperm indica o numero de tramos que se observou não ter condicao de participar
        // do calculo permanente da rede por não ter vazao nele
        avaliaPerm(malha, arqRede, narq, semPerm);
        // Vcr<int> inativo(narq,0);//alteracao4
        if (semPerm < narq) { // numero de tramos sem vazao eh inferior ao numero de tramos da rede, pode-se resolver a rede
            for (int i = 0; i < narq; i++) {
                if (arqRede.malha[i].ncoleta > 0 && arqRede.malha[i].perm == 1) {
                    int indicadorVP = i;
                    verificaTramoVazPres(indicadorVP, arqRede, malha);
                }
            }
            for (int i = 0; i < narq; i++) {
                if (arqRede.malha[i].nafluente == 0 && malha[i].cell(0).acsr.tipo == 1 && malha[i].cell(0).flui.flashCompleto == 1) {
                    malha[i].cell(0).acsr.injg.FluidoPro = malha[i].cell(0).flui;
                    malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                }
            }

            int chequeModFlui = 1;
            for (int i = 1; i < narq; i++) {
                if (malha[i].config().flashCompleto != malha[i - 1].config().flashCompleto &&
                    arqRede.malha[i].perm * arqRede.malha[i - 1].perm > 0) {
                    cout << "Malha com tramos com modelo de fluidos diferentes, tramos " << i << " e " << i - 1 << endl;
                    chequeModFlui = 0;
                }
                if (chequeModFlui == 0)
                    NumError("Malha com tramos com modelos de fluido diferentes");

                if (malha[i].config().tabelaDinamica != malha[i - 1].config().tabelaDinamica &&
                    arqRede.malha[i].perm * arqRede.malha[i - 1].perm > 0) {
                    cout << "Malha com tramos com modelo de fluidos diferentes, tramos - Tabela Dinamica" << i << " e " << i - 1 << endl;
                    chequeModFlui = 0;
                }
                if (chequeModFlui == 0)
                    NumError("Malha com tramos com modelos de fluido diferentes-Tabela Dinamica");
            }

            if (contapermRede < narq) {
                while ((*arqRede.vg1dSP).restartRede == 1) {
                    vector<int> coletorfinal;
                    for (int i = 0; i < narq; i++)
                        if (arqRede.malha[i].ncoleta == 0 && inativo[i] == 0 && arqRede.malha[i].perm == 1 && malha[i].config().ConContEntrada != 2)
                            coletorfinal.push_back(i);
                    sort(coletorfinal.begin(), coletorfinal.end());
                    int ncoletfim = coletorfinal.size();
                    Vcr<double> razcolet(narq - indativo - 0 * contapermRede);
                    Vcr<double> prescolet(narq - indativo - 0 * contapermRede);
                    if (arqRede.chute == 0) {
                        for (int i = 0; i < ncoletfim; i++)
                            chutePresRede(coletorfinal[i], malha, arqRede, arqRede.chutHol, razcolet, prescolet);
                        for (int i = 0; i < narq - indativo; i++) {
                            int testaAflu = 1;
                            for (int j = 0; j < ncoletfim; j++) {
                                if (i == coletorfinal[j]) {
                                    testaAflu = 0;
                                }
                            }
                            if (testaAflu == 1 && arqRede.malha[i].presimposta == 0)
                                malha[i].pGSup = prescolet[i] / razcolet[i];
                        }
                    } else {
                        Vcr<int> testeJus(narq);
                        Vcr<int> testeMon(narq);
                        for (int i = 0; i < narq; i++) {
                            testeJus[i] = 0;
                            testeMon[i] = 0;
                        }
                        for (int i = 0; i < narq; i++) {
                            if (arqRede.malha[i].ncoleta > 0) {
                                int col = arqRede.malha[i].coleta[0];
                                int nafl = arqRede.malha[col].nafluente;
                                for (int j = 0; j < nafl; j++) {
                                    int ind = arqRede.malha[col].afluente[j];
                                    if (arqRede.malha[ind].presimposta == 0 && testeJus[ind] == 0) {
                                        malha[ind].pGSup = arqRede.malha[i].presJus;
                                        testeJus[ind] = 1;
                                    }
                                }
                                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                                    int ind = arqRede.malha[i].coleta[j];
                                    if (testeMon[ind] == 0) {
                                        malha[ind].cell(0).pres = arqRede.malha[i].presJus;
                                        testeMon[ind] = 1;
                                    }
                                }
                            }
                            if (arqRede.malha[i].nafluente == 0) {
                                malha[i].cell(0).pres = arqRede.malha[i].presMon;
                                testeMon[i] = 1;
                            }
                        }
                    }

                    for (int i = 0; i < narq; i++) {
                        int nvecNo = normaEvol.size();
                        int testavec = 0;
                        if (arqRede.malha[i].nafluente > 0) {
                            for (int j = 0; j < nvecNo; j++) {
                                int marcaAflu = arqRede.malha[i].afluente[0];
                                if (normaEvol[j].naflu > 0) {
                                    for (int k = 0; k < normaEvol[j].naflu; k++) {
                                        if (marcaAflu == normaEvol[j].aflu[k])
                                            testavec = 1;
                                    }
                                }
                            }

                            if (testavec == 0) {
                                noRede temp;
                                temp.naflu = arqRede.malha[i].nafluente;
                                int indAflu = arqRede.malha[i].afluente[0];
                                temp.ncole = arqRede.malha[indAflu].ncoleta;
                                for (int j = 0; j < temp.naflu; j++)
                                    temp.aflu[j] = arqRede.malha[i].afluente[j];
                                for (int j = 0; j < temp.ncole; j++)
                                    temp.cole[j] = arqRede.malha[indAflu].coleta[j];
                                temp.cadastrado = 1;
                                normaEvol.push_back(temp);
                            }
                        }
                    }
                    (*arqRede.vg1dSP).restartRede = 0;
                    indativo = 0;
                    if (narq - semPerm > 1) {
                        if (malha[0].config().flashCompleto != 2)
                            convergido = convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                        else {
                            cout << "Iniciando modo Black-Oil" << endl;
                            Vcr<int> calclat0(narq);
                            Vcr<int> tipoFluido0(narq);
                            alteraModoFluidoCompBlack(malha, narq, calclat0, tipoFluido0);
                            convergido = convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                            if (convergido > 0) {
                                cout << "Iniciando modo Composicional" << endl;
                                alteraModoFluidoBlackComp(malha, narq, calclat0, tipoFluido0);
                                if (malha[0].config().tabelaDinamica == 1)
                                    cicloRedeCompCego(malha, arqRede, inativo, indativo, bloq);
                                (*arqRede.vg1dSP).iterRede = 1;
                                convergido = convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                            }
                        }
                    } else {
                        int konta = 0;
                        while (arqRede.malha[konta].perm == 0) {
                            konta++;
                        }
                        SolveTramoSolteiro(malha[konta], malha[konta].config().chutePerm);
                    }
                }
                if (convergido > 0)
                    cout << "Convergencia atingida" << endl;
                else
                    cout << "Falha na Convergencia" << endl;
                for (int i = 0; i < narq; i++) {
                    if (inativo[i] == 0 && arqRede.malha[i].perm == 1) {
                        malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
                        malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
                        if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
                        if(malha[i].config().nCelUnit>0){
                        	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
                        		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
                        }
                        // enterramento
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).calor.difus2D == 1) {
                                malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                            }
                        }
                        if (malha[i].config().lingas == 1)
                            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo, nrede);
                        // enterramento
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).calor.difus2D == 1) {
                                malha[i].cell(j).calor.poisson2D.imprimePermanente(i);
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0 && arqRede.malha[i].perm == 1) {
                    malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
                    malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
                    if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
                    if(malha[i].config().nCelUnit>0){
                    	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
                    		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
                    }
                    malha[i].kimpT++;
                    for (int j = 0; j < malha[i].config().ntendp; j++) {
                        malha[i].ImprimeTrendPCab(j, nrede);
                        malha[i].config().imprimeTrend(malha[i].celula, malha[i].MatTrendP[j], 0, j, 0);
                        malha[i].ImprimeTrendP(j, nrede);
                    }
                    if (malha[i].config().lingas == 1) {
                        malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo);
                        for (int j = 0; j < malha[i].config().ntendg; j++) {
                            malha[i].ImprimeTrendGCab(j, nrede);
                            malha[i].config().imprimeTrendG(malha[i].celulaG, malha[i].MatTrendG[j], 0, j, 0, 0);
                            malha[i].ImprimeTrendG(j, nrede);
                        }
                    }
                    // enterramento
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).calor.difus2D == 1) {
                            malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                        }
                    }
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        malha[i].cell(j).tempini = malha[i].cell(j).temp;
                        malha[i].cell(j).presLini = malha[i].cell(j).presL;
                        malha[i].cell(j).presini = malha[i].cell(j).pres;
                        malha[i].cell(j).MRini = malha[i].cell(j).MR;
                        malha[i].cell(j).MRini = malha[i].cell(j).MR;
                        malha[i].cell(j).MliqiniR0 = malha[i].cell(j).MliqiniR;
                        malha[i].cell(j).presLiniBuf = malha[i].cell(j).presL;
                        malha[i].cell(j).MRiniBuf = malha[i].cell(j).MR;
                        malha[i].cell(j).alfLini = malha[i].cell(j).alfL;
                        malha[i].cell(j).alfRini = malha[i].cell(j).alfR;
                        malha[i].cell(j).alfini = malha[i].cell(j).alf;
                        malha[i].cell(j).betLini = malha[i].cell(j).betL;
                        malha[i].cell(j).betRini = malha[i].cell(j).betR;
                        malha[i].cell(j).betini = malha[i].cell(j).bet;
                        malha[i].cell(j).FWini = malha[i].cell(j).FW;
                        malha[i].cell(j).QLRini = malha[i].cell(j).QLR;
                        malha[i].cell(j).alfPigEini = malha[i].cell(j).alfPigE;
                        malha[i].cell(j).alfPigDini = malha[i].cell(j).alfPigD;
                        malha[i].cell(j).betPigEini = malha[i].cell(j).betPigE;
                        malha[i].cell(j).betPigDini = malha[i].cell(j).betPigD;

                        malha[i].cell(j).rpC = malha[i].cell(j).rpCi =
                            malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rgC = malha[i].cell(j).rgCi =
                            malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rcC = malha[i].cell(j).rcCi =
                            malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rpL = malha[i].cell(j).rpLi =
                            malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                        malha[i].cell(j).rgL = malha[i].cell(j).rgLi =
                            malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                        malha[i].cell(j).rcL = malha[i].cell(j).rcLi =
                            malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).presL, malha[i].cell(j).tempL);

                        malha[i].cell(j).mipC = malha[i].cell(j).flui.ViscOleo(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).migC = malha[i].cell(j).flui.ViscGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).micC = malha[i].cell(j).fluicol.VisFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);

                        if (j > 0) {
                            malha[i].cell(j - 1).rpR = malha[i].cell(j - 1).rpRi = malha[i].cell(j).rpC;
                            malha[i].cell(j - 1).rgR = malha[i].cell(j - 1).rgRi = malha[i].cell(j).rgC;
                            malha[i].cell(j - 1).rcR = malha[i].cell(j - 1).rcRi = malha[i].cell(j).rcC;

                            malha[i].cell(j - 1).mipR = malha[i].cell(j).mipC;
                            malha[i].cell(j - 1).migR = malha[i].cell(j).migC;
                            malha[i].cell(j - 1).micR = malha[i].cell(j).micC;
                        }
                        if (j == malha[i].ncel) {
                            malha[i].cell(j).rpR = malha[i].cell(j).rpRi = malha[i].cell(j).rpC;
                            malha[i].cell(j).rgR = malha[i].cell(j).rgRi = malha[i].cell(j).rgC;
                            malha[i].cell(j).rcR = malha[i].cell(j).rcRi = malha[i].cell(j).rcC;

                            malha[i].cell(j).mipR = malha[i].cell(j).mipC;
                            malha[i].cell(j).migR = malha[i].cell(j).migC;
                            malha[i].cell(j).micR = malha[i].cell(j).micC;
                        }
                    }
                }
            }
        } else {
            cout << "ATENCAO REDE SEM VAZAO, POR FAVOR VERIFICAR SEUS DADOS DE ENTRADA" << endl;
            convergido = 1000;
        }
    } else {
        int i = 0;
        SolveTramoSolteiro(malha[0], malha[0].config().chutePerm);
        malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
        malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
        if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
        if(malha[i].config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
        		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
        }
        malha[i].kimpT++;
        for (int j = 0; j < malha[i].config().ntendp; j++) {
            malha[i].ImprimeTrendPCab(j, nrede);
            malha[i].config().imprimeTrend(malha[i].celula, malha[i].MatTrendP[j], 0, j, 0);
            malha[i].ImprimeTrendP(j, nrede);
        }
        if (malha[i].config().lingas == 1) {
            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo);
            for (int j = 0; j < malha[i].config().ntendg; j++) {
                malha[i].ImprimeTrendGCab(j, nrede);
                malha[i].config().imprimeTrendG(malha[i].celulaG, malha[i].MatTrendG[j], 0, j, 0, 0);
                malha[i].ImprimeTrendG(j, nrede);
            }
        }
        // enterramento
        for (int j = 0; j <= malha[i].ncel; j++) {
            if (malha[i].cell(j).calor.difus2D == 1) {
                malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
            }
        }
        if (malha[i].config().lingas == 1)
            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo, nrede);
    }

    ostringstream relatSucesso;
    relatSucesso << pathPrefixoArqSaida << "relatorioSucessoRede" << ".dat";
    string tmp = relatSucesso.str();
    ofstream escreveRelatorioSucesso(tmp.c_str(), ios_base::app);
    if (convergido > 0 && convergido < 1000) {
        for (int i = 0; i < narq; i++) {
            if (narq > 1)
                escreveRelatorioSucesso << nrede << "      ;" << i << "      ;" << arqRede.malha[i].perm << "     ;" << 1 - inativo[i] << "       " << endl;
            else
                escreveRelatorioSucesso << nrede << "      ;" << i << "      ;" << 1 << "     ;" << 1 << "       " << endl;
        }
    } else if (convergido == 1000) {
        escreveRelatorioSucesso << nrede << "      ;" << "   Rede sem Vazão    " << "     ;" << "   Rede sem Vazão    " << "     ;" << "   Rede sem Vazão    " << "       " << endl;
    } else {
        escreveRelatorioSucesso << nrede << "      ;" << "   Falha na Convergência    " << "     ;" << "   Falha na Convergência    " << "     ;" << "   Falha na Convergência    " << "       " << endl;
    }

    escreveRelatorioSucesso.close();

    if (arqRede.chaveredeT == 0) { // impressao do arquivo de log, quando a
        // opcao transiente nao esta ativa
        ostringstream saidaRede;
        saidaRede << pathPrefixoArqSaida << "LogEvento.dat";
        string tmp = saidaRede.str();
        ofstream escreveIni(tmp.c_str(), ios_base::app);
        srand(time(NULL));
        int frase = rand() % 15;
        escreveIni << "                                  Rede interna = " << nrede << endl;
        escreveIni << "*******************************************************************************" << endl;
        escreveIni << "                                  UFA!!!!!!!!                                  " << endl;
        escreveIni << saidaTexto[frase] << endl;
        escreveIni << saidaSubTexto[frase] << endl;
        escreveIni << "*******************************************************************************" << endl;

        time_t now = time(0);
        tm *ltm = localtime(&now); ///////////Retirado de https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
        escreveIni << "datahora = ";
        escreveIni << ltm->tm_mday << "/";
        escreveIni << 1 + ltm->tm_mon << "/";
        escreveIni << 1900 + ltm->tm_year << " ";
        escreveIni << 0 + ltm->tm_hour << ":";
        escreveIni << 0 + ltm->tm_min << ":";
        escreveIni << 0 + ltm->tm_sec;
        escreveIni << endl;

        int diaFim = (ltm->tm_mday);
        int horaFim;
        if (diaFim == diaIni)
            horaFim = ltm->tm_hour;
        else
            horaFim = ltm->tm_hour + 24;
        horaFim *= 3600;
        int minutoFim = 60 * ltm->tm_min;
        int segundoFim = ltm->tm_sec;
        int totalFim = horaFim + minutoFim + segundoFim;
        int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
        escreveIni << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
        escreveIni << "     Versao    " << versao << endl;

        escreveIni.close();
    }
}

void RedeProd(SProd *malha, Rede &arqRede, int narq,
              Vcr<int> &inativo, int indativo, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
              tipoSimulacao_t tipoSimulacaoMrt, int nrede, int contapermRede, varGlob1D *vG1d, vector<noRede> &normaEvol, vector<tramoPart> &bloq) {
    // neste metodo e feita a construcao de um objeto rede e a sua resolucao permanente
    if (narq > 1) { // narq=numero de tramos ativos da rede
        (*arqRede.vg1dSP).simulaTransiente = arqRede.chaveredeT;
        for (int i = 0; i < narq; i++) {                 // laco em que se constroe os tramos da rede
            int bloqR = testaBloqueio(arqRede, i, bloq); // verifica se em uma rede de manifold existe uma configuracao de derivacao especial
            // em que parte dos tramos coletados no manifold e desviado para a linha de teste
            int tD = arqRede.tabelaDinamica; // rede com modelo composicional e uso de tabela dinamica
            if ((*arqRede.vg1dSP).fluidoRede == 0)
                (*arqRede.vg1dSP).tipoFluidoRedeGlob = 1; // tipoFluidoRede==1 implica em que a rede será uma rede dominada pela fase líquida
            // construtor do tramo temporario com as caracteristicas basicas do tramo i da rede
            SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                             tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD, bloqR, 0, arqRede.malha[i].reverso, 0, 0, 0, arqRede.malha[i].perm);

            if (arqRede.malha[i].ncoleta > 0)
                temporario.noextremo = 0; // caso o tramo não esteja no fim da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            else
                temporario.noextremo = 1; // caso o tramo esteja no fim da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            if (arqRede.malha[i].nafluente > 0)
                temporario.noinicial = 0; // caso o tramo não esteja no inicio da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            else
                temporario.noinicial = 1; // caso o tramo esteja no inicio da rede, tomando como referencia o sentido
            // de escoamento assumido para a rede
            malha[i] = temporario; // carregando o tramo temporario no tramo i da malha
            malha[i].indTramo = i;
            if (malha[i].config().tabelaDinamica == 1) {
                // caso a opcao de tabela dinamica esteja ligada, como, para o modelo composicional, sempre, primeiramente, se faz um
                // calculo inicial com modelo black oil, a opcao tabela dinamica e desligada temporariamente enquanto se esta no modo blackoil
                for (int j = 0; j <= malha[i].ncel; j++) {
                    malha[i].cell(j).flui.tabelaDinamica = 0;
                    if (malha[i].cell(j).acsr.tipo == 1) {
                        malha[i].cell(j).acsr.injg.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 2) {
                        malha[i].cell(j).acsr.injl.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 3) {
                        malha[i].cell(j).acsr.ipr.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 15) {
                        malha[i].cell(j).acsr.radialPoro.flup.tabelaDinamica = 0;
                        int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                        for (int k = 0; k < ncelRad; k++) {
                            malha[j].cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
                        }
                    } else if (malha[i].cell(j).acsr.tipo == 16) {
                        malha[i].cell(j).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
                        int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                        for (int k = 0; k < ncelRad; k++) {
                            malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
                        }
                        int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                        for (int k = 0; k < ncelEle; k++) {
                            malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
                        }
                    } else if (malha[i].cell(j).acsr.tipo == 10) {
                        malha[i].cell(j).acsr.injm.FluidoPro.tabelaDinamica = 0;
                    } else if (malha[i].cell(j).acsr.tipo == 9 && malha[i].cell(j).acsr.fontechk.abertura > 1e-6) {
                        malha[i].cell(j).acsr.fontechk.fluidoP.tabelaDinamica = 0;
                        malha[i].cell(j).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
                    }
                }
            }
            if (arqRede.malha[i].nafluente == 0 && malha[i].config().ConContEntrada == 1) {
                // preparando o tramo para o caso em que este esta no inicio da rede e tem
                // como condicao de contorno e pressao na entrada do tramo, neste caso, tem de se ter uma estimativa inicial de vazao
                // neste tramo, alem disto, e necessario definir detalhes como o beta in-situ e o RGO
                // observar tambem que no solver permanente, e necessario para este tipo de condicao de contorno
                // e necessario colocar um acessorio de fonte, onde e armazendao a estimativa de vazao no inicio do tramo
                // caso o fluido seja dominado por liquido, a fointe de armazenamento sera do tipo liquido
                //  se gas, a fonte sera de gas
                double tit = malha[i].config().CCPres.tit[0];
                /// atencao, verificar se aposicao do if esta correta.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ((*arqRede.vg1dSP).fluidoRede == 0) {
                    if (arqRede.chutHol > 0.1)
                        arqRede.chutHol = 0.1;
                    if (arqRede.relax > 0.5)
                        arqRede.relax = 0.5;
                }
                if (tit < (1. - (*arqRede.vg1dSP).localtiny)) {

                    malha[i].tempE = malha[i].config().CCPres.temperatura[0];
                    malha[i].presE = malha[i].config().CCPres.pres[0];
                    double rgoentrada = 0.;
                    double rcst = malha[i].cell(0).fluicol.MasEspFlu(1.01, 15);
                    double rcis = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                    malha[i].betaE = malha[i].config().CCPres.bet[0] * rcst / rcis;

                    if (malha[i].config().flashCompleto != 2) {
                        double rgST = malha[i].cell(0).flui.Deng * 1.225;
                        double roST = 141.5 * 1000. / (131.5 + malha[i].cell(0).flui.API);
                        double rw = malha[i].cell(0).flui.Denag * 1000.;
                        double RSent = malha[i].cell(0).flui.RS(malha[i].presE, malha[i].tempE) * 6.29 / 35.31467;
                        double bsw = malha[i].cell(0).flui.BSW;
                        double bo = malha[i].cell(0).flui.BOFunc(malha[i].presE, malha[i].tempE);
                        double ba = malha[i].cell(0).flui.BAFunc(malha[i].presE, malha[i].tempE);
                        double fw = bsw * ba / (bo + ba * bsw - bsw * bo);

                        double Ost = (1 - malha[i].betaE) * bo * (1 - fw);
                        double Wst = fw * (1 - malha[i].betaE);
                        double Cst = malha[i].betaE;
                        malha[i].titE = malha[i].config().CCPres.tit[0];

                        rgoentrada = RSent * rgST * Ost + malha[i].titE * roST * Ost + malha[i].titE * rw * Wst +
                                     malha[i].titE * rcst * Cst;
                        rgoentrada /= ((1 - malha[i].titE) * rgST * Ost);
                    }

                    if ((*arqRede.vg1dSP).fluidoRede == 1) { // dominado por liquido
                        if (rgoentrada > 1400.)
                            rgoentrada = 1400.;
                        malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                        malha[i].cell(0).acsr.injl.bet = malha[i].betaE;
                        malha[i].cell(0).acsr.injl.temp = malha[i].config().CCPres.temperatura[0];
                        malha[i].cell(0).acsr.injl.FluidoPro.BSW = malha[i].cell(0).flui.BSW;
                        malha[i].cell(0).acsr.injl.FluidoPro.TempL = malha[i].cell(0).flui.TempL;
                        malha[i].cell(0).acsr.injl.FluidoPro.TempH = malha[i].cell(0).flui.TempH;
                        malha[i].cell(0).acsr.injl.FluidoPro.LVisL = malha[i].cell(0).flui.LVisL;
                        malha[i].cell(0).acsr.injl.FluidoPro.LVisH = malha[i].cell(0).flui.LVisH;
                        if (malha[i].config().flashCompleto != 2) {
                            malha[i].cell(0).acsr.injl.FluidoPro.RGO = rgoentrada;
                            malha[i].cell(0).acsr.injl.FluidoPro.Deng = malha[i].cell(0).flui.Deng;
                            malha[i].cell(0).acsr.injl.FluidoPro.Denag = malha[i].cell(0).flui.Denag;
                            malha[i].cell(0).acsr.injl.FluidoPro.yco2 = malha[i].cell(0).flui.yco2;
                            malha[i].cell(0).acsr.injl.FluidoPro.API = malha[i].cell(0).flui.API;
                            malha[i].cell(0).acsr.injl.FluidoPro.RenovaFluido();
                        } else {
                            malha[i].cell(0).flui.atualizaPropCompStandard();
                            malha[i].cell(0).flui.atualizaPropComp(malha[i].presE, malha[i].tempE);
                            double rl = malha[i].cell(0).flui.MasEspLiq(malha[i].presE, malha[i].tempE);
                            double titH = malha[i].cell(0).flui.FracMassHidra(malha[i].presE, malha[i].tempE);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double rlMix = malha[i].betaE * rc + (1. - malha[i].betaE) * rl;
                            double val1 = ((1. - malha[i].betaE) * rl * titH / (1. - titH));
                            malha[i].titE = val1 / (rlMix + val1);
                            malha[i].cell(0).acsr.injl.FluidoPro = malha[i].cell(0).flui;
                        }
                        if (fabs(malha[i].cell(0).acsr.injl.QLiq) < 1e-15)
                            malha[i].cell(0).acsr.injl.QLiq = 1000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                        malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injl.FluidoPro;
                    } else {
                        malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                        malha[i].cell(0).acsr.injg.seco = 0;
                        malha[i].cell(0).acsr.injg.temp = malha[i].config().CCPres.temperatura[0];
                        malha[i].cell(0).acsr.injg.FluidoPro.BSW = malha[i].cell(0).flui.BSW;
                        malha[i].cell(0).acsr.injg.FluidoPro.TempL = malha[i].cell(0).flui.TempL;
                        malha[i].cell(0).acsr.injg.FluidoPro.TempH = malha[i].cell(0).flui.TempH;
                        malha[i].cell(0).acsr.injg.FluidoPro.LVisL = malha[i].cell(0).flui.LVisL;
                        malha[i].cell(0).acsr.injg.FluidoPro.LVisH = malha[i].cell(0).flui.LVisH;
                        if (malha[i].config().flashCompleto != 2) { // caso nao composicional
                            double rgST = malha[i].cell(0).flui.Deng * 1.225;
                            double rg = malha[i].cell(0).flui.MasEspGas(malha[i].presE, malha[i].tempE);
                            double rl = malha[i].cell(0).flui.MasEspoleo(malha[i].presE, malha[i].tempE);
                            double tit = malha[i].cell(0).flui.FracMass(malha[i].presE, malha[i].tempE);
                            double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                            double val2 = (rg / rl) * (1 - tit) / tit;
                            double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                            malha[i].titE = titT;
                            malha[i].betaE = val1 / (val2 + val1);
                            rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                            if (rgoentrada > 1400.)
                                rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                            malha[i].cell(0).flui.RGO = malha[i].cell(0).acsr.injg.FluidoPro.RGO = rgoentrada;
                            malha[i].cell(0).acsr.injg.FluidoPro.Deng = malha[i].cell(0).flui.Deng;
                            malha[i].cell(0).acsr.injg.FluidoPro.Denag = malha[i].cell(0).flui.Denag;
                            malha[i].cell(0).acsr.injg.FluidoPro.yco2 = malha[i].cell(0).flui.yco2;
                            malha[i].cell(0).acsr.injg.FluidoPro.API = malha[i].cell(0).flui.API;
                            malha[i].cell(0).acsr.injg.FluidoPro.BSW = malha[i].cell(0).flui.BSW;
                            malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                        } else { // caso composicional
                            malha[i].cell(0).flui.atualizaPropCompStandard();
                            malha[i].cell(0).flui.atualizaPropComp(malha[i].presE, malha[i].tempE);
                            double rgST = malha[i].cell(0).flui.Deng * 1.225;
                            double rg = malha[i].cell(0).flui.MasEspGas(malha[i].presE, malha[i].tempE);
                            double rl = malha[i].cell(0).flui.MasEspoleo(malha[i].presE, malha[i].tempE);
                            double tit = malha[i].cell(0).flui.FracMass(malha[i].presE, malha[i].tempE);
                            double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                            double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                            double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                            double val2 = (rg / rl) * (1 - tit) / tit;
                            double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                            malha[i].titE = titT;
                            malha[i].betaE = val1 / (val2 + val1);
                            malha[i].cell(0).acsr.injg.FluidoPro = malha[i].cell(0).flui;
                        }
                        malha[i].cell(0).acsr.injg.fluidocol = malha[i].config().fluc;
                        malha[i].cell(0).acsr.injg.razCompGas = malha[i].config().CCPres.bet[0];
                        if (malha[i].cell(0).acsr.injg.QGas < 1e-15)
                            malha[i].cell(0).acsr.injg.QGas = 1000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                        malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injg.FluidoPro;
                    }
                } else { // dominado por gas
                    double rgoentrada = 0.;
                    malha[i].cell(0).pres = malha[i].config().CCPres.pres[0];
                    malha[i].cell(0).acsr.injg.temp = malha[i].config().CCPres.temperatura[0];
                    malha[i].cell(0).acsr.injg.seco = 0;
                    malha[i].cell(0).acsr.injg.FluidoPro.TempL = malha[i].cell(0).flui.TempL;
                    malha[i].cell(0).acsr.injg.FluidoPro.TempH = malha[i].cell(0).flui.TempH;
                    malha[i].cell(0).acsr.injg.FluidoPro.LVisL = malha[i].cell(0).flui.LVisL;
                    malha[i].cell(0).acsr.injg.FluidoPro.LVisH = malha[i].cell(0).flui.LVisH;
                    if (malha[i].config().flashCompleto != 2) { // caso nao composicional
                        double rgST = malha[i].cell(0).flui.Deng * 1.225;
                        double rg = malha[i].cell(0).flui.MasEspGas(malha[i].presE, malha[i].tempE);
                        double rl = malha[i].cell(0).flui.MasEspoleo(malha[i].presE, malha[i].tempE);
                        double tit = malha[i].cell(0).flui.FracMass(malha[i].presE, malha[i].tempE);
                        double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                        double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                        double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                        double val2 = (rg / rl) * (1 - tit) / tit;
                        double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                        malha[i].titE = titT;
                        malha[i].betaE = val1 / (val2 + val1);
                        rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                        if (rgoentrada > 1400.)
                            rgoentrada = 1e6 + (*arqRede.vg1dSP).localtiny;
                        malha[i].cell(0).flui.RGO = malha[i].cell(0).acsr.injg.FluidoPro.RGO = rgoentrada;
                        malha[i].cell(0).acsr.injg.FluidoPro.Deng = malha[i].cell(0).flui.Deng;
                        malha[i].cell(0).acsr.injg.FluidoPro.Denag = malha[i].cell(0).flui.Denag;
                        malha[i].cell(0).acsr.injg.FluidoPro.yco2 = malha[i].cell(0).flui.yco2;
                        malha[i].cell(0).acsr.injg.FluidoPro.API = malha[i].cell(0).flui.API;
                        malha[i].cell(0).acsr.injg.FluidoPro.BSW = malha[i].cell(0).flui.BSW;
                        malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                    } else { // caso composicional
                        malha[i].cell(0).flui.atualizaPropCompStandard();
                        malha[i].cell(0).flui.atualizaPropComp(malha[i].presE, malha[i].tempE);
                        double rgST = malha[i].cell(0).flui.Deng * 1.225;
                        double rg = malha[i].cell(0).flui.MasEspGas(malha[i].presE, malha[i].tempE);
                        double rl = malha[i].cell(0).flui.MasEspoleo(malha[i].presE, malha[i].tempE);
                        double tit = malha[i].cell(0).flui.FracMass(malha[i].presE, malha[i].tempE);
                        double rcST = malha[i].cell(0).fluicol.MasEspFlu(1.01, 20.);
                        double rc = malha[i].cell(0).fluicol.MasEspFlu(malha[i].presE, malha[i].tempE);
                        double val1 = (rcST / rc) * (rg / rgST) * malha[i].config().CCPres.bet[0] / tit;
                        double val2 = (rg / rl) * (1 - tit) / tit;
                        double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
                        malha[i].titE = titT;
                        malha[i].betaE = val1 / (val2 + val1);
                        malha[i].cell(0).acsr.injg.FluidoPro = malha[i].cell(0).flui;
                    }
                    malha[i].cell(0).acsr.injg.fluidocol = malha[i].config().fluc;
                    malha[i].cell(0).acsr.injg.razCompGas = malha[i].config().CCPres.bet[0];
                    if (malha[i].cell(0).acsr.injg.QGas < 1e-15)
                        malha[i].cell(0).acsr.injg.QGas = 10000. * pow(malha[i].cell(0).duto.dia, 2.) / pow(6. * 2.54 / 100, 2.);
                    malha[i].cell(0).fluiL = &malha[i].cell(0).acsr.injl.FluidoPro;
                }
            }
            // fazendo uma estimativa inicial da vazao total de liquido e gas que escoara na rede, para tanto,
            // se buscara as fontes de liquido e gas que estao computadas no inicio do tramo, observar que e uma estimativa imperfeita, pois
            // no caso de condicao de contorno pressao, esta vazao e apenas um chuet, alem de nao contabilizar a vazao em celulas que não a
            // inicial do tramo
            if (malha[i].cell(0).acsr.tipo == 2 && arqRede.malha[i].perm == 1)
                (*arqRede.vg1dSP).somavaz += malha[i].cell(0).acsr.injl.QLiq / 86400.;
            else if (malha[i].cell(0).acsr.tipo == 1 && arqRede.malha[i].perm == 1)
                (*arqRede.vg1dSP).somavazG += malha[i].cell(0).acsr.injg.QGas / 86400.;
            else if (malha[i].cell(0).acsr.tipo == 10 && arqRede.malha[i].perm == 1 && malha[i].config().ConContEntrada != 2) {
                double massTot = malha[i].cell(0).acsr.injm.MassP + malha[i].cell(0).acsr.injm.MassG;
                double rhopstd = 141.5 * 1000. / (malha[i].cell(0).acsr.injm.FluidoPro.API + 131.5);
                double bsw = malha[i].cell(0).acsr.injm.FluidoPro.BSW;
                double rhogstd = malha[i].cell(0).acsr.injm.FluidoPro.Deng * 1.225;
                double tit = malha[i].cell(0).acsr.injm.FluidoPro.FracMassHidra(1., 20.);
                (*arqRede.vg1dSP).somavaz += massTot * (1 - tit) * ((1 - bsw) / rhopstd + bsw / 1000.) +
                                             malha[i].cell(0).acsr.injm.MassC / malha[i].cell(0).acsr.injm.fluidocol.rholStd;
                (*arqRede.vg1dSP).somavazG += massTot * tit / rhogstd;
            }
            (*arqRede.vg1dSP).somaarea += malha[i].cell(0).duto.area; // somaarea sera utilizado para fazer uma distribuicao inicial da vazao em cada tramo, tomando como
            // referencia a razao entre a area do tramo e somaarea, esta razao sera aplicada para distribuir o chute inicial da vazao de cada tramo,
            // tomando somavaz e somavazG e aplicando a razao de area
            if (arqRede.malha[i].perm == 0)
                contapermRede++; // contabilizando se o tramo esta com o permanente ligado, o numero de tramops dan rede
            // pode ser diferente do numero de tramso que efetivamente poarticipara do solver permanente da rede
            temporario.config().tabelaDinamica = 0;
        }
        for (int i = 0; i < narq; i++) {
            malha[i].chuteHol = arqRede.chutHol; // o chutehol serve como uma estimativa inicial do holdup na rede para o calculo inicial
            // da variacao de pressao em cada tramo da rede devido aa hidrostatica
            if (arqRede.malha[i].ncoleta > 0) {
                int icol = arqRede.malha[i].coleta[0];
                malha[i].fluiRevRede = malha[icol].cell(0).flui; // computando o fluido do tramo coletor, caso no afluente
                // ocorra reversao do escoamento pre defindio, o fluido que escoarah no tramo afluente vira do fluido que se encontra no coletor
                if (malha[i].config().flashCompleto == 2)
                    malha[i].fluiRevRede.atualizaPropCompStandard();
            }
        }
        int semPerm = 0;
        // preprocessamento final, em avaliaPerm faz-se um "pente fino" em tramos sem vazao, estes são descartados da rede
        // para o caluculo permanente, semperm indica o numero de tramos que se observou não ter condicao de participar
        // do calculo permanente da rede por não ter vazao nele
        avaliaPerm(malha, arqRede, narq, semPerm);
        // Vcr<int> inativo(narq,0);//alteracao4
        if (semPerm < narq) { // numero de tramos sem vazao eh inferior ao numero de tramos da rede, pode-se resolver a rede
            for (int i = 0; i < narq; i++) {
                if (arqRede.malha[i].nafluente == 0 && malha[i].cell(0).acsr.tipo == 1 && malha[i].cell(0).flui.flashCompleto == 1) {
                    malha[i].cell(0).acsr.injg.FluidoPro = malha[i].cell(0).flui;
                    malha[i].cell(0).acsr.injg.FluidoPro.RenovaFluido();
                }
            }

            int chequeModFlui = 1;
            for (int i = 1; i < narq; i++) {
                if (malha[i].config().flashCompleto != malha[i - 1].config().flashCompleto &&
                    arqRede.malha[i].perm * arqRede.malha[i - 1].perm > 0) {
                    cout << "Malha com tramos com modelo de fluidos diferentes, tramos " << i << " e " << i - 1 << endl;
                    chequeModFlui = 0;
                }
                if (chequeModFlui == 0)
                    NumError("Malha com tramos com modelos de fluido diferentes");

                if (malha[i].config().tabelaDinamica != malha[i - 1].config().tabelaDinamica &&
                    arqRede.malha[i].perm * arqRede.malha[i - 1].perm > 0) {
                    cout << "Malha com tramos com modelo de fluidos diferentes, tramos - Tabela Dinamica" << i << " e " << i - 1 << endl;
                    chequeModFlui = 0;
                }
                if (chequeModFlui == 0)
                    NumError("Malha com tramos com modelos de fluido diferentes-Tabela Dinamica");
            }

            if (contapermRede < narq) {
                while ((*arqRede.vg1dSP).restartRede == 1) {
                    vector<int> coletorfinal;
                    for (int i = 0; i < narq; i++)
                        if (arqRede.malha[i].ncoleta == 0 && inativo[i] == 0 && arqRede.malha[i].perm == 1 && malha[i].config().ConContEntrada != 2)
                            coletorfinal.push_back(i);
                    sort(coletorfinal.begin(), coletorfinal.end());
                    int ncoletfim = coletorfinal.size();
                    Vcr<double> razcolet(narq - indativo - 0 * contapermRede);
                    Vcr<double> prescolet(narq - indativo - 0 * contapermRede);
                    if (arqRede.chute == 0) {
                        for (int i = 0; i < ncoletfim; i++)
                            chutePresRede(coletorfinal[i], malha, arqRede, arqRede.chutHol, razcolet, prescolet);
                        for (int i = 0; i < narq - indativo; i++) {
                            int testaAflu = 1;
                            for (int j = 0; j < ncoletfim; j++) {
                                if (i == coletorfinal[j]) {
                                    testaAflu = 0;
                                }
                            }
                            if (testaAflu == 1 && arqRede.malha[i].presimposta == 0)
                                malha[i].pGSup = prescolet[i] / razcolet[i];
                        }
                    } else {
                        Vcr<int> testeJus(narq);
                        Vcr<int> testeMon(narq);
                        for (int i = 0; i < narq; i++) {
                            testeJus[i] = 0;
                            testeMon[i] = 0;
                        }
                        for (int i = 0; i < narq; i++) {
                            if (arqRede.malha[i].ncoleta > 0) {
                                int col = arqRede.malha[i].coleta[0];
                                int nafl = arqRede.malha[col].nafluente;
                                for (int j = 0; j < nafl; j++) {
                                    int ind = arqRede.malha[col].afluente[j];
                                    if (arqRede.malha[ind].presimposta == 0 && testeJus[ind] == 0) {
                                        malha[ind].pGSup = arqRede.malha[i].presJus;
                                        testeJus[ind] = 1;
                                    }
                                }
                                for (int j = 0; j < arqRede.malha[i].ncoleta; j++) {
                                    int ind = arqRede.malha[i].coleta[j];
                                    if (testeMon[ind] == 0) {
                                        malha[ind].cell(0).pres = arqRede.malha[i].presJus;
                                        testeMon[ind] = 1;
                                    }
                                }
                            }
                            if (arqRede.malha[i].nafluente == 0) {
                                malha[i].cell(0).pres = arqRede.malha[i].presMon;
                                testeMon[i] = 1;
                            }
                        }
                    }

                    for (int i = 0; i < narq; i++) {
                        int nvecNo = normaEvol.size();
                        int testavec = 0;
                        if (arqRede.malha[i].nafluente > 0) {
                            for (int j = 0; j < nvecNo; j++) {
                                int marcaAflu = arqRede.malha[i].afluente[0];
                                if (normaEvol[j].naflu > 0) {
                                    for (int k = 0; k < normaEvol[j].naflu; k++) {
                                        if (marcaAflu == normaEvol[j].aflu[k])
                                            testavec = 1;
                                    }
                                }
                            }

                            if (testavec == 0) {
                                noRede temp;
                                temp.naflu = arqRede.malha[i].nafluente;
                                int indAflu = arqRede.malha[i].afluente[0];
                                temp.ncole = arqRede.malha[indAflu].ncoleta;
                                for (int j = 0; j < temp.naflu; j++)
                                    temp.aflu[j] = arqRede.malha[i].afluente[j];
                                for (int j = 0; j < temp.ncole; j++)
                                    temp.cole[j] = arqRede.malha[indAflu].coleta[j];
                                temp.cadastrado = 1;
                                normaEvol.push_back(temp);
                            }
                        }
                    }
                    (*arqRede.vg1dSP).restartRede = 0;
                    indativo = 0;
                    if (narq - semPerm > 1) {
                        if (malha[0].config().flashCompleto != 2)
                            convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                        else {
                            cout << "Iniciando modo Black-Oil" << endl;
                            Vcr<int> calclat0(narq);
                            Vcr<int> tipoFluido0(narq);
                            alteraModoFluidoCompBlack(malha, narq, calclat0, tipoFluido0);
                            convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                            cout << "Iniciando modo Composicional" << endl;
                            alteraModoFluidoBlackComp(malha, narq, calclat0, tipoFluido0);
                            if (malha[0].config().tabelaDinamica == 1)
                                cicloRedeCompCego(malha, arqRede, inativo, indativo, bloq);
                            (*arqRede.vg1dSP).iterRede = 1;
                            convergeRede(malha, arqRede, narq, inativo, indativo, normaEvol, bloq);
                        }
                    } else {
                        int konta = 0;
                        while (arqRede.malha[konta].perm == 0) {
                            konta++;
                        }
                        SolveTramoSolteiro(malha[konta], malha[konta].config().chutePerm);
                    }
                }
                cout << "Convergencia atingida" << endl;
                for (int i = 0; i < narq; i++) {
                    if (inativo[i] == 0 && arqRede.malha[i].perm == 1) {
                        malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
                        malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
                        if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
                        if(malha[i].config().nCelUnit>0){
                        	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
                        		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
                        }
                        // enterramento
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).calor.difus2D == 1) {
                                malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                            }
                        }
                        if (malha[i].config().lingas == 1)
                            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo, nrede);
                        // enterramento
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].cell(j).calor.difus2D == 1) {
                                malha[i].cell(j).calor.poisson2D.imprimePermanente(i);
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < narq; i++) {
                if (inativo[i] == 0 && arqRede.malha[i].perm == 1) {
                    malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
                    malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
                    if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
                    if(malha[i].config().nCelUnit>0){
                    	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
                    		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
                    }
                    // enterramento
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).calor.difus2D == 1) {
                            malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                        }
                    }
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        malha[i].cell(j).tempini = malha[i].cell(j).temp;
                        malha[i].cell(j).presLini = malha[i].cell(j).presL;
                        malha[i].cell(j).presini = malha[i].cell(j).pres;
                        malha[i].cell(j).MRini = malha[i].cell(j).MR;
                        malha[i].cell(j).MRini = malha[i].cell(j).MR;
                        malha[i].cell(j).MliqiniR0 = malha[i].cell(j).MliqiniR;
                        malha[i].cell(j).presLiniBuf = malha[i].cell(j).presL;
                        malha[i].cell(j).MRiniBuf = malha[i].cell(j).MR;
                        malha[i].cell(j).alfLini = malha[i].cell(j).alfL;
                        malha[i].cell(j).alfRini = malha[i].cell(j).alfR;
                        malha[i].cell(j).alfini = malha[i].cell(j).alf;
                        malha[i].cell(j).betLini = malha[i].cell(j).betL;
                        malha[i].cell(j).betRini = malha[i].cell(j).betR;
                        malha[i].cell(j).betini = malha[i].cell(j).bet;
                        malha[i].cell(j).FWini = malha[i].cell(j).FW;
                        malha[i].cell(j).QLRini = malha[i].cell(j).QLR;
                        malha[i].cell(j).alfPigEini = malha[i].cell(j).alfPigE;
                        malha[i].cell(j).alfPigDini = malha[i].cell(j).alfPigD;
                        malha[i].cell(j).betPigEini = malha[i].cell(j).betPigE;
                        malha[i].cell(j).betPigDini = malha[i].cell(j).betPigD;

                        malha[i].cell(j).rpC = malha[i].cell(j).rpCi =
                            malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rgC = malha[i].cell(j).rgCi =
                            malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rcC = malha[i].cell(j).rcCi =
                            malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).rpL = malha[i].cell(j).rpLi =
                            malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                        malha[i].cell(j).rgL = malha[i].cell(j).rgLi =
                            malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                        malha[i].cell(j).rcL = malha[i].cell(j).rcLi =
                            malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).presL, malha[i].cell(j).tempL);

                        malha[i].cell(j).mipC = malha[i].cell(j).flui.ViscOleo(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).migC = malha[i].cell(j).flui.ViscGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                        malha[i].cell(j).micC = malha[i].cell(j).fluicol.VisFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);

                        if (j > 0) {
                            malha[i].cell(j - 1).rpR = malha[i].cell(j - 1).rpRi = malha[i].cell(j).rpC;
                            malha[i].cell(j - 1).rgR = malha[i].cell(j - 1).rgRi = malha[i].cell(j).rgC;
                            malha[i].cell(j - 1).rcR = malha[i].cell(j - 1).rcRi = malha[i].cell(j).rcC;

                            malha[i].cell(j - 1).mipR = malha[i].cell(j).mipC;
                            malha[i].cell(j - 1).migR = malha[i].cell(j).migC;
                            malha[i].cell(j - 1).micR = malha[i].cell(j).micC;
                        }
                        if (j == malha[i].ncel) {
                            malha[i].cell(j).rpR = malha[i].cell(j).rpRi = malha[i].cell(j).rpC;
                            malha[i].cell(j).rgR = malha[i].cell(j).rgRi = malha[i].cell(j).rgC;
                            malha[i].cell(j).rcR = malha[i].cell(j).rcRi = malha[i].cell(j).rcC;

                            malha[i].cell(j).mipR = malha[i].cell(j).mipC;
                            malha[i].cell(j).migR = malha[i].cell(j).migC;
                            malha[i].cell(j).micR = malha[i].cell(j).micC;
                        }
                    }
                }
            }
        } else {
            cout << "ATENCAO REDE SEM VAZAO, POR FAVOR VERIFICAR SEUS DADOS DE ENTRADA" << endl;
        }
    } else {

        int i = 0;
        int tD = arqRede.tabelaDinamica;
        SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                         tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD);
        temporario.noextremo = 1;
        temporario.noinicial = 1;
        malha[i] = temporario;
        malha[i].indTramo = i;

        if (malha[i].config().tabelaDinamica == 1) {
            for (int j = 0; j <= malha[i].ncel; j++) {
                malha[i].cell(j).flui.tabelaDinamica = 0;
                if (malha[i].cell(j).acsr.tipo == 1) {
                    malha[i].cell(j).acsr.injg.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 2) {
                    malha[i].cell(j).acsr.injl.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 3) {
                    malha[i].cell(j).acsr.ipr.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 15) {
                    malha[i].cell(j).acsr.radialPoro.flup.tabelaDinamica = 0;
                    int ncelRad = malha[j].cell(i).acsr.radialPoro.ncel;
                    for (int k = 0; k < ncelRad; k++) {
                        malha[j].cell(i).acsr.radialPoro.celula[k].flup.tabelaDinamica = 0;
                    }
                } else if (malha[i].cell(j).acsr.tipo == 16) {
                    malha[i].cell(j).acsr.poroso2D.dados.flup.tabelaDinamica = 0;
                    int ncelRad = malha[j].cell(i).acsr.poroso2D.dados.transfer.ncel;
                    for (int k = 0; k < ncelRad; k++) {
                        malha[j].cell(i).acsr.poroso2D.dados.transfer.celula[k].flup.tabelaDinamica = 0;
                    }
                    int ncelEle = malha[j].cell(i).acsr.poroso2D.malha.nele;
                    for (int k = 0; k < ncelEle; k++) {
                        malha[j].cell(i).acsr.poroso2D.malha.mlh2d[k].flup.tabelaDinamica = 0;
                    }
                } else if (malha[i].cell(j).acsr.tipo == 10) {
                    malha[i].cell(j).acsr.injm.FluidoPro.tabelaDinamica = 0;
                } else if (malha[i].cell(j).acsr.tipo == 9 && malha[i].cell(j).acsr.fontechk.abertura > 1e-6) {
                    malha[i].cell(j).acsr.fontechk.fluidoP.tabelaDinamica = 0;
                    malha[i].cell(j).acsr.fontechk.fluidoPamb.tabelaDinamica = 0;
                }
            }
        }
        SolveTramoSolteiro(malha[0], malha[0].config().chutePerm);
        malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
        malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
        if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
        if(malha[i].config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
        		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
        }
        // enterramento
        for (int j = 0; j <= malha[i].ncel; j++) {
            if (malha[i].cell(j).calor.difus2D == 1) {
                malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
            }
        }
        if (malha[i].config().lingas == 1)
            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo, nrede);
    }

    ostringstream relatSucesso;
    relatSucesso << pathPrefixoArqSaida << "relatorioSucessoRede" << ".dat";
    string tmp = relatSucesso.str();
    ofstream escreveRelatorioSucesso(tmp.c_str(), ios_base::app);
    for (int i = 0; i < narq; i++) {
        if (narq > 1)
            escreveRelatorioSucesso << nrede << "      ;" << i << "      ;" << arqRede.malha[i].perm << "     ;" << 1 - inativo[i] << "       " << endl;
        else
            escreveRelatorioSucesso << nrede << "      ;" << i << "      ;" << 1 << "     ;" << 1 << "       " << endl;
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    int diaFim = (ltm->tm_mday);
    int horaFim;
    if (diaFim == diaIni)
        horaFim = ltm->tm_hour;
    else
        horaFim = ltm->tm_hour + 24;
    horaFim *= 3600;
    int minutoFim = 60 * ltm->tm_min;
    int segundoFim = ltm->tm_sec;
    int totalFim = horaFim + minutoFim + segundoFim;
    int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
    escreveRelatorioSucesso << "datahora = ";
    escreveRelatorioSucesso << ltm->tm_mday << "/";
    escreveRelatorioSucesso << 1 + ltm->tm_mon << "/";
    escreveRelatorioSucesso << 1900 + ltm->tm_year << " ";
    escreveRelatorioSucesso << 0 + ltm->tm_hour << ":";
    escreveRelatorioSucesso << 0 + ltm->tm_min << ":";
    escreveRelatorioSucesso << 0 + ltm->tm_sec;
    escreveRelatorioSucesso << endl;
    escreveRelatorioSucesso << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
    escreveRelatorioSucesso << "     Versao    " << versao << endl;

    escreveRelatorioSucesso.close();
}

void TransAnel(int narq, int nfontes, int *indfonte, int *indtramo, int *posicfonte, int indAnel,
               vector<fonteposic> dreno, SProd *malha, Rede &arqRede) {

    Vcr<double> razMast0(narq);
    Vcr<double> razMast(narq);
    Vcr<int> celpos(narq);
    Vcr<int> kontasnp(narq, 1);
    if (malha[indAnel].config().ConContEntrada == 1) {
        malha[indAnel].cell(0).acsr.injg.QGas = 0.;
        malha[indAnel].cell(0).MC = malha[indAnel].cell(0).fontemassGR;
        malha[indAnel].cell(0).Mliqini = 0.;
        malha[indAnel].cell(0).Mliqini0 = malha[indAnel].cell(0).Mliqini;
        malha[indAnel].cell(0).QL = 0;
        malha[indAnel].cell(0).QLini = malha[indAnel].cell(0).QL;
        double rhog = malha[indAnel].cell(0).flui.MasEspGas(
            malha[indAnel].cell(0).pres, malha[indAnel].cell(0).temp);
        malha[indAnel].cell(0).QG = (malha[indAnel].cell(0).MC - malha[indAnel].cell(0).Mliqini) / rhog;
        malha[indAnel].cell(0).QGini = malha[indAnel].cell(0).QG;
        malha[indAnel].cell(0).fontemassLR = 0.;
        malha[indAnel].cell(0).fontemassCR = 0.;
        malha[indAnel].cell(0).fontemassGR = 0.;
        malha[indAnel].cell(0).fluiL = &malha[indAnel].cell(0).flui;
    }

    for (int i = 0; i < narq; i++) {
        malha[i].modoPerm = 0;
        malha[i].config().tfinal = malha[i].tfinal = (*arqRede.vg1dSP).TmaxR;
    }
    double dtteste;
    (*arqRede.vg1dSP).lixo5R = 0.;
    while ((*arqRede.vg1dSP).lixo5R < (*arqRede.vg1dSP).TmaxR) {
        dtteste = 1e6;
        for (int i = 0; i < narq; i++) {
            malha[i].determinaDT();
            malha[i].dtauxCFL = malha[i].dt;
            malha[i].dtauxFinal = malha[i].dt;
            if (malha[i].dt < dtteste)
                dtteste = malha[i].dt;
        }
        dtteste = dtteste / 1.;
        for (int i = 0; i < narq; i++) {
            malha[i].dt = dtteste;
            for (int j = 0; j <= malha[i].ncel; j++) {
                malha[i].cell(j).dt = malha[i].dt;
                malha[i].cell(j).dt2 = malha[i].dt;
                malha[i].cell(j).dtPig = malha[i].dt;
            }
        }
        if ((*arqRede.vg1dSP).lixo5 < 1e-5) {
            for (int i = 0; i < narq; i++)
                malha[i].renovaTemp();
        }

        if ((*arqRede.vg1dSP).lixo5R >= 600) {
        }

        for (int i = 0; i < narq; i++) {
            malha[i].reinicia = 0;
            celpos[i] = malha[i].config().master1.posic;
            razMast0[i] = malha[i].cell(celpos[i]).acsr.chk.AreaGarg / malha[i].cell(celpos[i]).duto.area;
            malha[i].solveLinGas();

            malha[i].config().atualiza(malha[i].noinicial, malha[i].noextremo, malha[i].derivaAnel,
                                  malha[i].chokeSup,
                                  malha[i].chokeInj, malha[i].celula,
                                  malha[i].celulaG, malha[i].pGSup, malha[i].temperatura,
                                  malha[i].presiniG, malha[i].tempiniG,
                                  malha[i].presE, malha[i].tempE, malha[i].titE, malha[i].betaE, (*arqRede.vg1dSP).lixo5R);
            razMast[i] = malha[i].cell(celpos[i]).acsr.chk.AreaGarg / malha[i].cell(celpos[i]).duto.area;
            malha[i].modeloCompleto = malha[i].config().correcaoMassaEspLiq;
            if (razMast[i] != razMast0[i])
                malha[i].modeloCompleto = 0;
            if (malha[i].modeloCompleto == 1)
                malha[i].avaliaVariaDpDt(razMast[i], razMast0[i]);
        }
        for (int i = 0; i < dreno.size(); i++) {
            int ifposic = posicfonte[i];
            malha[indAnel].cell(ifposic).acsr.injg.QGas = 0.;
            for (int j = 0; j < dreno[i].nmani; j++) {
                int itramo = dreno[i].mani[j];
                malha[indAnel].cell(ifposic).acsr.injg.QGas -= malha[itramo].gasCell(0).VGasR * 86400 /
                                                                 malha[itramo].gasCell(0).flui.MasEspGas(1.0, 15.6);
            }
        }

        int modeloCompletoGlob = 0;
        for (int i = 0; i < narq; i++) {
            if (malha[i].modeloCompleto == 1)
                modeloCompletoGlob = 1;
        }

        for (int kontaAcop = 0; kontaAcop <= 1 * modeloCompletoGlob; kontaAcop++) {
            int reiniGlob = 0;
            for (int i = 0; i < narq; i++) {
                if (i != indAnel) {
                    if (malha[i].modeloCompleto == 0) {
                        for (int j = 0; j <= malha[i].ncel; j++)
                            malha[i].cell(j).m2d = 0.;
                    } else {
                        for (int j = 0; j <= malha[i].ncel; j++) {
                            if (malha[i].TransMassModel == 0)
                                malha[i].cell(j).m2d = 1.;
                            else
                                malha[i].cell(j).m2d = 0.;
                        }
                    }
                    malha[i].EvoluiFrac(1., 1., kontaAcop);
                }
            }
            for (int i = 0; i < narq; i++)
                if (malha[i].reinicia == -1)
                    reiniGlob = -1;
            if (reiniGlob < 0) {
                for (int i = 0; i < narq; i++) {
                    if (i != indAnel)
                        malha[i].ReiniEvolFrac0();
                    if (malha[i].dt < dtteste) {
                        dtteste = malha[i].dt;
                        malha[i].dtauxFinal = malha[i].dt;
                    }
                }
                for (int i = 0; i < narq; i++) {
                    malha[i].dt = dtteste;
                    if (i != indAnel) {
                        malha[i].ReiniEvolFrac();
                        malha[i].EvoluiFrac(1., 1., kontaAcop);
                    }
                    malha[i].reinicia = 0;
                }
            }
            if ((*arqRede.vg1dSP).lixo5R <= 1e-15)
                (*arqRede.vg1dSP).dtrede = dtteste;

            for (int i = 0; i < narq; i++) {

                malha[i].AtualizaPig();
                malha[i].atenuaDtMax();
                malha[i].calcCCpres(1., 1.);
                if (malha[i].cell(malha[i].ncel).alf < 0.05 && malha[i].masChkSup == 1)
                    malha[i].cell(malha[i].ncel).alf = 0.05;
                if (malha[i].cell(celpos[i]).alf < 0.05 && razMast[i] <= malha[i].config().master1.razareaativ)
                    malha[i].cell(celpos[i]).alf = 0.05;
                if (i != indAnel)
                    malha[i].renovaterm();
                malha[i].SolveAcopPV();
                if (kontaAcop == 1 * modeloCompletoGlob)
                    malha[i].renova();
                else {
                    for (int j = 0; j <= malha[i].ncel; j++)
                        malha[i].cell(j).dpdt = (malha[i].termolivreP[2 * j + 1] -
                                                   malha[i].cell(j).pres) /
                                                  malha[i].cell(j).dt;
                }
                if (kontaAcop != 1 * modeloCompletoGlob) {
                    malha[i].indpigPini = malha[i].indpigP;
                    for (int j = 0; j <= malha[i].ncel; j++)
                        malha[i].cell(j).FeiticoDoTempo2();
                }
            }
        }

        (*arqRede.vg1dSP).iterRedeT = 0;
        while ((*arqRede.vg1dSP).iterRedeT < 1) {

            for (int i = 0; i < narq; i++) {
                malha[i].SolveTrans();
                for (int kontasnp = 0; kontasnp < malha[i].config().nsnp; kontasnp++) {
                    if (fabs((*arqRede.vg1dSP).lixo5 - malha[i].config().tempsnp[kontasnp]) < malha[i].dt) {
                        WriteSnapShot(malha[i], kontasnp, i); // registro do arquivo SNP
                    }
                }
            }
            (*arqRede.vg1dSP).iterRedeT++;
        }

        for (int i = 0; i < dreno.size(); i++) {
            int ifposic = posicfonte[i];
            for (int j = 0; j < dreno[i].nmani; j++) {
                int itramo = dreno[i].mani[j];
                malha[itramo].presiniG = malha[indAnel].cell(ifposic).pres;
                malha[itramo].tempiniG = malha[indAnel].cell(ifposic).temp;
            }
        }

        (*arqRede.vg1dSP).lixo5R += dtteste;
        (*arqRede.vg1dSP).dtrede = dtteste;
    }

    for (int i = 0; i < narq; i++)
        WriteSnapShot(malha[i], malha[i].config().nsnp, i);
}

void calcPeriAnelGL(SProd *malha, int narq, int nfontes, int *indfonte, int *indtramo, int *posicfonte, int indAnel,
                    int iter, vector<fonteposic> dreno) {

    for (int i = malha[indAnel].config().ninjgas; i < malha[indAnel].config().ninjgas + dreno.size(); i++) {
        int ifposic = posicfonte[i - malha[indAnel].config().ninjgas];
        for (int j = 0; j < dreno[i - malha[indAnel].config().ninjgas].nmani; j++) {
            int itramo = dreno[i - malha[indAnel].config().ninjgas].mani[j];
            int ndirve = dreno[i - malha[indAnel].config().ninjgas].nmani;
            malha[itramo].config().gasinj.vazgas[0] = -malha[indAnel].cell(ifposic).acsr.injg.QGas / ndirve;
            malha[itramo].chokeInj.tempEstag = malha[itramo].tempiniG =
                malha[itramo].config().gasinj.temperatura[0] = malha[indAnel].cell(ifposic).temp;
            malha[itramo].gasCell(0).tipoCC = 1;
            if (malha[itramo].config().chokep.abertura[0] > 0.6) {
                if (iter == 0)
                    malha[itramo].buscaProdPfundoPerm();
                malha[itramo].gasCell(0).tipoCC = 0;
                malha[itramo].chokeInj.presEstag = malha[itramo].presiniG = malha[itramo].config().gasinj.presinj[0] =
                    malha[indAnel].cell(ifposic).pres;
                malha[itramo].buscaProdPfundoPerm(malha[itramo].cell(0).pres);
            } else {
                if (iter == 0)
                    malha[itramo].buscaProdPfundoPerm2();
                malha[itramo].gasCell(0).tipoCC = 0;
                malha[itramo].chokeInj.presEstag = malha[itramo].presiniG = malha[itramo].config().gasinj.presinj[0] = malha[indAnel].cell(ifposic).pres;
                malha[itramo].buscaProdPfundoPerm2(malha[itramo].cell(0).pres);
            }
        }
    }
    for (int i = malha[indAnel].config().ninjgas; i < malha[indAnel].config().ninjgas + dreno.size(); i++) {
        int ifposic = posicfonte[i - malha[indAnel].config().ninjgas];
        malha[indAnel].cell(ifposic).acsr.injg.QGas = 0.;
        for (int j = 0; j < dreno[i - malha[indAnel].config().ninjgas].nmani; j++) {
            int itramo = dreno[i - malha[indAnel].config().ninjgas].mani[j];
            malha[indAnel].cell(ifposic).acsr.injg.QGas -= malha[itramo].gasCell(0).VGasR * 86400 /
                                                             malha[itramo].gasCell(0).flui.MasEspGas(1.0, 15.6);
        }
    }
}

double calcErroGL(SProd *malha, int narq, int nfontes, int *posicfonte, int indAnel,
                  vector<fonteposic> dreno) {
    double erro = 0;
    for (int i = 0; i < malha[indAnel].config().ninjgas + dreno.size(); i++) {
        int ifposic;
        if (i < malha[indAnel].config().ninjgas)
            ifposic = malha[indAnel].config().fonteg[i].posicP;
        else
            ifposic = posicfonte[i - malha[indAnel].config().ninjgas];
        erro += malha[indAnel].cell(ifposic).acsr.injg.QGas;
    }
    if (malha[indAnel].config().ConContEntrada == 1)
        erro += malha[indAnel].cell(0).acsr.injg.QGas;
    double divisor = 1;
    if (fabs(malha[indAnel].cell(0).acsr.injg.QGas) > 0)
        divisor = fabs(malha[indAnel].cell(0).acsr.injg.QGas);
    return erro / divisor;
}

double objetCC0(SProd *malha, int narq, int nfontes, int *indfonte, int *indtramo, int *posicfonte, int indAnel,
                double pchute, int iter, vector<fonteposic> dreno) {
    malha[indAnel].marchaProdPerm1(pchute);
    calcPeriAnelGL(malha, narq, nfontes, indfonte, indtramo, posicfonte, indAnel, iter, dreno);
    return calcErroGL(malha, narq, nfontes, posicfonte, indAnel, dreno);
}

double SIGN(double a, double b) {
    return (b >= 0 ? 1.0 : -1.0) * fabs(a);
}
double zriddr(SProd *malha, double x1, double x2,
              int narq, int nfontes, int *indfonte, int *indtramo, int *posicfonte, int indAnel,
              vector<fonteposic> dreno) {
    double xacc = 1e-6;
    int maxit = 100;
    double fl = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicfonte, indAnel, x1, 1, dreno);
    double fh = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicfonte, indAnel, x2, 1, dreno);
    if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
        double xl = x1;
        double xh = x2;
        double ans = -1.e20;
        for (int j = 0; j < maxit; j++) {
            double xm = 0.5 * (xl + xh);
            double fm = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicfonte, indAnel, xm, 1, dreno);
            double s = sqrt(fm * fm - fl * fh);
            if (s == 0.0)
                return ans;
            double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
            if (fabs(xnew - ans) <= xacc) {
                return ans;
            }
            ans = xnew;
            double fnew = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicfonte, indAnel, ans, 1, dreno);
            if (fabs(fnew) <= xacc)
                return ans;
            if (SIGN(fm, fnew) != fm) {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fnew;
            } else if (SIGN(fl, fnew) != fl) {
                xh = ans;
                fh = fnew;
            } else if (SIGN(fh, fnew) != fh) {
                xl = ans;
                fl = fnew;
            } else
                return -100000.;
            if (fabs(xh - xl) <= xacc)
                return ans;
        }
        return 100000.;
    } else {
        if (fabs(fl) <= xacc)
            return x1;
        if (fabs(fh) <= xacc)
            return x2;
        return -10000.;
    }
}

void conectaPrincipal(SProd *malha, int iP, int iS, int itera = 1) {
    for (int i = malha[iP].PrimSecIniRedeP; i >= malha[iP].PrimSecFimRedeP; i--) {
        int j = malha[iP].PrimSecIniRedeP + malha[iP].SecPrimIniRedeP - i;
        if (itera == 1)
            malha[iP].cell(i).calor.Textern2 = malha[iS].cell(j).calor.Tcamada[0][0];
        else
            malha[iP].cell(i).calor.Textern2 = malha[iS].cell(j).calor.Textern1;
        malha[iP].cell(i).calor.betext = malha[iS].cell(j).calor.betint;

        int icam = malha[iP].cell(i).calor.geom.ncamadas - 1;
        int idisc = malha[iP].cell(i).calor.ncamada[icam] - 1;
        if (itera == 1)
            malha[iS].cell(j).calor.Tint2 = malha[iP].cell(i).calor.Tcamada[icam][idisc];
        else
            malha[iS].cell(j).calor.Tint2 = malha[iS].cell(j).calor.Textern1;
        malha[iP].cell(i).calor.colunaDia = malha[iS].cell(j).duto.dia;
        malha[iP].cell(i).calor.geom.b = malha[iS].cell(j).calor.geom.a;

        if (itera == 1)
            malha[iP].cell(i).calor.Textern1 = malha[iS].cell(j).temp;
        else
            malha[iP].cell(i).calor.Textern1 = malha[iS].cell(j).calor.Textern1;
        double cpg = malha[iS].cell(j).flui.CalorGas(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double rhog = malha[iS].cell(j).flui.MasEspGas(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double condg = malha[iS].cell(j).flui.CondGas(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double viscg = malha[iS].cell(j).flui.ViscGas(malha[iS].cell(j).pres, malha[iS].cell(j).temp) * 1.e-3;
        double cplp = malha[iS].cell(j).flui.CalorLiq(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double rholp = malha[iS].cell(j).flui.MasEspLiq(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double condlp = malha[iS].cell(j).fluicol.CondLiq(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double visclp = malha[iS].cell(j).flui.ViscOleo(malha[iS].cell(j).pres, malha[iS].cell(j).temp) * 1.e-3;
        double cplc = malha[iS].cell(j).fluicol.CalorLiq(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double rholc = malha[iS].cell(j).fluicol.MasEspFlu(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double condlc = malha[iS].cell(j).fluicol.CondLiq(malha[iS].cell(j).pres, malha[iS].cell(j).temp);
        double visclc = malha[iS].cell(j).fluicol.VisFlu(malha[iS].cell(j).pres, malha[iS].cell(j).temp) * 1.e-3;
        double alfExt = malha[iS].cell(j).alf;
        double betExt = malha[iS].cell(j).bet;
        double cplMix = (1. - betExt) * cplp + betExt * cplc;
        double rholMix = (1. - betExt) * rholp + betExt * rholc;
        double condlMix = (1. - betExt) * condlp + betExt * condlc;
        double visclMix = (1. - betExt) * visclp + betExt * visclc;
        if (itera == 1)
            malha[iP].cell(i).calor.Vextern1 = (malha[iS].cell(j).QG) / (malha[iS].cell(j).duto.area) +
                                                 malha[iS].cell(j).QL / (malha[iS].cell(j).duto.area);
        else
            malha[iP].cell(i).calor.Vextern1 = 100.;

        malha[iP].cell(i).calor.kextern1 = alfExt * condg + (1. - alfExt) * condlMix;
        malha[iP].cell(i).calor.cpextern1 = alfExt * cpg + (1. - alfExt) * cplMix;
        malha[iP].cell(i).calor.rhoextern1 = alfExt * rhog + (1. - alfExt) * rholMix;
        malha[iP].cell(i).calor.viscextern1 = alfExt * viscg + (1. - alfExt) * visclMix;
    }
}

double chutePresRedeParalelaSec(SProd *malha, int iP, int iS, double pchute, int indPartida) {
    double betaChute = 0.;
    double taux;
    if (malha[iS].cell(0).acsr.tipo == 2)
        betaChute = malha[iS].cell(0).acsr.injl.bet;
    for (int i = indPartida; i > 0; i--) {
        taux = malha[iS].config().celp[i].textern;
        if (taux < malha[iS].config().tmin)
            taux = malha[iS].config().tmin;
        double rhol = malha[iS].cell(i).flui.MasEspLiq(pchute, taux);
        double rhog = malha[iS].cell(i).flui.MasEspGas(pchute, taux);
        // se o sistema nao for um anel principal de GL
        // para a estimativa de pressao, se utilizara so a hidrostatica de liquido, o que dara um
        // chute inicial de pressao muito alta
        double alfa = 0.;
        if (betaChute < 0.5) {
            double tit = malha[iS].cell(i).flui.FracMassHidra(pchute, taux);
            alfa = tit * rhol / (rhog - tit * rhog + tit * rhol);
        }
        if (malha[iS].cell(0).acsr.tipo == 1)
            alfa = 1.;
        double rhomix = (1. - alfa) * rhol + alfa * rhog;
        double dxmed = 0.5 * (malha[iS].cell(i).dx + malha[iS].cell(i - 1).dx);
        // avanco da pressao por meio da hidrostatica e da perda por friccao estimada
        pchute += (rhomix * 9.81 * sin(malha[iS].cell(i).duto.teta) * dxmed) / 98066.5;
        if (malha[iS].config().usaTabela == 1 && (malha[iS].config().tabent.pmax - pchute) < (*malha[iS].vg1dSP).localtiny)
            pchute = 0.9 * malha[iS].config().tabent.pmax;
    }
    if (pchute > 1000)
        pchute = 1000.; // pressao maxima de chute
    return pchute;
}

void SolveRedeParalelaTrans(SProd *malha, Rede &arqRede, int nrede) {
    int narq = arqRede.nsisprod;
    Vcr<int> celpos(narq);
    Vcr<int> kontasnp(narq, 1);

    Vcr<double> fonteG(narq);
    Vcr<double> fonteP(narq);
    Vcr<double> fonteC(narq);

    int nfontes = arqRede.nfonteR;
    int iS = arqRede.tramoSecundario;
    int iP = arqRede.tramoPrimario;

    for (int i = 0; i < narq; i++) {
        malha[i].modoPerm = 0;
        malha[i].config().tfinal = malha[i].tfinal = (*arqRede.vg1dSP).TmaxR;
    }

    double dtteste;
    (*arqRede.vg1dSP).lixo5R = 0.;
    while ((*arqRede.vg1dSP).lixo5R < (*arqRede.vg1dSP).TmaxR) {
        for (int i = 0; i < narq; i++) {
            if ((*arqRede.vg1dSP).lixo5R >= 7480) {
            }
        }

        if (malha[iP].verificaAcopRedeP == 1 && malha[iS].verificaAcopRedeS == 1) {
            conectaPrincipal(malha, iP, iS);
            int nacopP = malha[iP].PrimSecIniRedeP - malha[iP].PrimSecFimRedeP;
            for (int iacop = 0; iacop <= nacopP; iacop++) {
                malha[iS].cell(malha[iS].SecPrimIniRedeP + iacop).fluxcalAcopRedeP =
                    malha[iP].cell(malha[iP].PrimSecIniRedeP - iacop).fluxcalmed;
            }
        }

        for (int i = 0; i < narq; i++) {
            if (malha[i].config().calculaEnvelope == 1 && (*arqRede.vg1dSP).lixo5 <= (*arqRede.vg1dSP).TmaxR) { //*vg1dSP).lixo5>0 && //chris - Hidratos
                FA_Hidrato solverHidrato(malha[i]);
                solverHidrato.solverHidrato();
            }

            if (malha[i].config().lingas > 0 && malha[i].config().calculaEnvelope == 1 && (*arqRede.vg1dSP).lixo5 <= (*arqRede.vg1dSP).TmaxR) { //*vg1dSP).lixo5>0 && //chris - Hidratos
                FA_Hidrato_Servico solverHidratoG(malha[i]);
                solverHidratoG.solverHidratoG();
            }

            if (malha[i].config().flashCompleto == 2 && (*arqRede.vg1dSP).lixo5 < 1e-15) {
                malha[i].atualizaMiniTab();
            }

            if ((*arqRede.vg1dSP).lixo5 < 1e-15) {
                for (int j = 0; j < malha[i].config().ntendp; j++) {
                    malha[i].config().imprimeTrend(malha[i].celula, malha[i].MatTrendP[j], (*arqRede.vg1dSP).lixo5, i, malha[i].ntrend[j]);
                }
                malha[i].renovaTemp();
            }
        }

        dtteste = 1e6;
        for (int i = 0; i < narq; i++) {
            malha[i].modeloCompleto = malha[i].config().correcaoMassaEspLiq;
            malha[i].determinaDT();
            malha[i].dtauxCFL = malha[i].dt;
            malha[i].dtauxFinal = malha[i].dt;
            if (malha[i].dt < dtteste)
                dtteste = malha[i].dt;
        }
        for (int i = 0; i < narq; i++) {
            malha[i].dt = dtteste;
            for (int j = 0; j <= malha[i].ncel; j++) {
                malha[i].cell(j).dt = malha[i].dt;
                malha[i].cell(j).dt2 = malha[i].dt;
                malha[i].cell(j).dtPig = malha[i].dt;
            }
        }
        for (int i = 0; i < narq; i++) {
            malha[i].reinicia = 0;
            celpos[i] = malha[i].config().master1.posic;
            // razMast0[i]=malha[i].cell(celpos[i]).acsr.chk.AreaGarg/malha[i].cell(celpos[i]).duto.area;//caso so Master
            malha[i].aberturaVal0(); // caso varias valvulas
            malha[i].solveLinGas();
            int presinterna = 1;
            malha[i].config().atualiza(malha[i].noinicial, presinterna, malha[i].derivaAnel,
                                  malha[i].chokeSup,
                                  malha[i].chokeInj, malha[i].celula,
                                  malha[i].celulaG, malha[i].pGSup, malha[i].temperatura,
                                  malha[i].presiniG, malha[i].tempiniG,
                                  malha[i].presE, malha[i].tempE, malha[i].titE, malha[i].betaE, (*arqRede.vg1dSP).lixo5R);
            malha[i].atualizaCC1();
            for (int j = 0; j <= malha[i].config().nvalv; j++)
                malha[i].vRazMastCrit[j] = 0.5; // caso varias valvulas
            malha[i].aberturaVal1();            // caso varias valvulas
            for (int j = 0; j <= malha[i].config().nvalv; j++)
                if (malha[i].vRazMast1[j] != malha[i].vRazMast0[j])
                    malha[i].modeloCompleto = 0; // caso varias valvulas
            if (malha[i].modeloCompleto == 1)
                malha[i].avaliaVariaDpDt(); // caso varias valvulas
            if (malha[i].modeloCompleto == 0)
                malha[i].config().cicloAcopTerm = 0;
            else
                malha[i].config().cicloAcopTerm = 1;
        }
        int modeloCompletoGlob = 1;
        for (int i = 0; i < narq; i++) {
            if (malha[i].modeloCompleto == 0)
                modeloCompletoGlob = 0;
        }
        int ciclomax = modeloCompletoGlob;

        for (int i = 0; i < narq; i++) {
            malha[i].abertoini = malha[i].aberto;
            malha[i].tempoabertoini = malha[i].tempoaberto;
            malha[i].masChkSupini = malha[i].masChkSup;
            malha[i].mudaModoChkini = malha[i].mudaModoChk;
        }

        for (int kontaAcop = 0; kontaAcop <= 1 * modeloCompletoGlob; kontaAcop++) {
            int reiniGlob = 0;

            for (int i = 0; i < narq; i++) {

                if (malha[i].modeloCompleto == 0) {
                    for (int j = 0; j <= malha[i].ncel; j++)
                        malha[i].cell(j).m2d = 0.;
                } else {
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        double razDp = 0.1;
                        double razDT = 1;
                        if (j < celpos[i] && malha[i].cell(celpos[i]).acsr.chk.AreaGarg < 1e-15 * malha[i].cell(celpos[i]).acsr.chk.AreaTub) {
                            razDT = 1;
                        } else if (j == celpos[i] + 1 &&
                                   malha[i].cell(celpos[i]).acsr.chk.AreaGarg < 1e-15 * malha[i].cell(celpos[i]).acsr.chk.AreaTub) {
                            razDT = 1;
                        }
                        if ((fabs(malha[i].cell(j).dpdtIni) / malha[i].cell(j).pres < razDp) && fabs(malha[i].cell(j).dTdtIni) < razDT) {
                            if (malha[i].TransMassModel == 0)
                                malha[i].cell(j).m2d = 1.;
                            else
                                malha[i].cell(j).m2d = 0.;
                            malha[i].cell(j).mudaDT = 1.;
                        } else {
                            malha[i].cell(j).m2d = 0.;
                            malha[i].cell(j).mudaDT = 0.;
                        }
                    }
                }
                if (malha[i].config().estabCol == 1) {
                    for (int j = 0; j <= celpos[i]; j++) {
                        malha[i].cell(j).m2d = 0.;
                        malha[i].cell(j).mudaDT = 0.;
                        malha[i].cell(j).estabCol = 1;
                    }
                }

                malha[i].EvoluiFrac(1., 0., kontaAcop);
                for (int j = 0; j <= malha[i].ncel; j++) {
                    if (malha[i].cell(j).acsr.tipo == 15) {
                        malha[i].cell(j).acsr.radialPoro.avancoSW(malha[i].dt);
                        if (malha[i].cell(j).acsr.radialPoro.reinicia == -1) {
                            if (malha[i].reinicia > -1)
                                malha[i].reinicia = -1;
                        }
                    } else if (malha[i].cell(j).acsr.tipo == 16) {
                        malha[i].cell(j).acsr.poroso2D.avancoSW(malha[i].dt);
                        if (malha[i].cell(j).acsr.poroso2D.reinicia == -1) {
                            if (malha[i].reinicia > -1)
                                malha[i].reinicia = -1;
                        }
                    }
                }
                if (malha[i].config().correcaoMassaEspLiq == 1) {
                    for (int j = 0; j < malha[i].ncel; j++)
                        malha[i].cell(j + 1).mudaDTL = malha[i].cell(j).mudaDT;
                }

                if (kontaAcop == 0 && malha[i].config().controleDTvalv == 1) {
                    if (malha[i].reinicia > -1)
                        malha[i].restringeDTporValv(); // caso varias valvulas
                }
            }
            for (int i = 0; i < narq; i++) {
                if (malha[i].reinicia == -1)
                    reiniGlob = -1;
            }
            if (reiniGlob < 0) {
                for (int i = 0; i < narq; i++) {
                    malha[i].ReiniEvolFrac0();
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).acsr.tipo == 15) {
                            malha[i].cell(j).acsr.radialPoro.reavaliaDT(malha[i].dt);
                        }
                        if (malha[i].cell(j).acsr.tipo == 16) {
                            malha[i].cell(j).acsr.poroso2D.reavaliaDT(malha[i].dt);
                        }
                    }
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).acsr.tipo == 15) {
                            malha[i].cell(j).acsr.radialPoro.reiniciaEvoluiSW(malha[i].dt);
                        }
                        if (malha[i].cell(j).acsr.tipo == 16) {
                            malha[i].cell(j).acsr.poroso2D.reiniciaEvoluiSW(malha[i].dt);
                        }
                    }
                    if (malha[i].dt < dtteste) {
                        dtteste = malha[i].dt;
                    }
                }
                for (int i = 0; i < narq; i++) {
                    malha[i].dt = dtteste;
                    malha[i].dtauxFinal = malha[i].dt;
                    malha[i].ReiniEvolFrac();
                    malha[i].EvoluiFrac(1., 0., kontaAcop);
                    malha[i].reinicia = 0;
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        if (malha[i].cell(j).acsr.tipo == 15) {
                            malha[i].cell(j).acsr.radialPoro.avancoSWcorrec();
                        }
                        if (malha[i].cell(j).acsr.tipo == 16) {
                            malha[i].cell(j).acsr.poroso2D.avancoSWcorrec();
                        }
                    }
                }
            }
            if ((*arqRede.vg1dSP).lixo5R <= 1e-15)
                (*arqRede.vg1dSP).dtrede = dtteste;
            for (int i = 0; i < narq; i++) {
                malha[i].AtualizaPig();

                if (kontaAcop == 0)
                    malha[i].dtCicMin = malha[i].dt;

                if (kontaAcop == 1 * modeloCompletoGlob)
                    malha[i].atenuaDtMax();

                if (modeloCompletoGlob == 1) {
                    fonteC[i] = malha[i].cell(malha[i].ncel).fontemassCR;
                    fonteP[i] = malha[i].cell(malha[i].ncel).fontemassLR;
                    fonteG[i] = malha[i].cell(malha[i].ncel).fontemassGR;
                }

                malha[i].calcCCpres();
                malha[i].renovaterm();
                if (malha[i].cell(malha[i].ncel).alf < 0.05 && malha[i].masChkSup == 1)
                    malha[i].cell(malha[i].ncel).alf = 0.05;
                // caso varias valvulas
                for (int j = 0; j <= malha[i].config().nvalv; j++) {
                    int celposAux;
                    if (j > 0)
                        celposAux = malha[i].config().valv[j - 1].posicP;
                    else
                        celposAux = celpos[i];
                    if (malha[i].cell(celposAux).alf < 0.05 && malha[i].vRazMast1[j] <= malha[i].config().master1.razareaativ)
                        malha[i].cell(celposAux).alf = 0.05;
                }
                malha[i].SolveAcopPV();

                if (kontaAcop < 1 * modeloCompletoGlob) {
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        malha[i].cell(j).dpdt = 1 * (malha[i].termolivreP[2 * j + 1] - malha[i].cell(j).pres) / malha[i].cell(j).dt;
                        malha[i].cell(j).dpdtIni = malha[i].cell(j).dpdt;
                    }
                }
                if (kontaAcop == 1 * modeloCompletoGlob || malha[i].config().cicloAcopTerm == 1) {
                    malha[i].renova();
                }
                if (malha[i].config().cicloAcopTerm == 1 && modeloCompletoGlob == 1) {
                    if (kontaAcop < 1 * modeloCompletoGlob)
                        for (int j = 0; j <= malha[i].ncel; j++)
                            malha[i].cell(j).dpdt = malha[i].cell(j).d2pdt2;
                    malha[i].marchaEnergTrans(kontaAcop, ciclomax);
                }
                if (kontaAcop != 1 * modeloCompletoGlob) {
                    for (int j = 0; j <= malha[i].ncel; j++) {
                        malha[i].cell(j).FeiticoDoTempo2();
                        if (malha[i].cell(j).acsr.tipo == 15) {
                            malha[i].cell(j).acsr.radialPoro.FeiticoDoTempoSW();
                        } else if (malha[i].cell(j).acsr.tipo == 16) {
                            malha[i].cell(j).acsr.poroso2D.FeiticoDoTempoSW();
                        }
                    }

                    malha[i].cell(malha[i].ncel).fontemassCR = fonteC[i];
                    malha[i].cell(malha[i].ncel).fontemassLR = fonteP[i];
                    malha[i].cell(malha[i].ncel).fontemassGR = fonteG[i];

                    malha[i].aberto = malha[i].abertoini;
                    malha[i].tempoaberto = malha[i].tempoabertoini;
                }
            }
        }
        for (int i = 0; i < narq; i++) {
            if (modeloCompletoGlob == 0 || malha[i].config().cicloAcopTerm == 0) {
                for (int ciclo = 0; ciclo <= ciclomax; ciclo++) {
                    malha[i].marchaEnergTrans(ciclo, ciclomax);
                }
            }
        }
        for (int i = 0; i < nfontes; i++) {
            int celS = arqRede.conexFR[i].noS;
            int celP = arqRede.conexFR[i].noP;
            malha[iP].cell(celP).acsr.fontechk.tamb = malha[iS].cell(celS).temp;
            malha[iP].cell(celP).acsr.fontechk.pamb = malha[iS].cell(celS).pres;
            malha[iP].cell(celP).acsr.fontechk.betISamb = malha[iS].cell(celS).bet;
            malha[iP].cell(celP).acsr.fontechk.fluidoPamb = malha[iS].cell(celS).flui;
            double alf = malha[iP].cell(celP).alf;
            double bet = malha[iP].cell(celP).bet;
            double pres = malha[iP].cell(celP).pres;
            double temp = malha[iP].cell(celP).temp;
            malha[iP].cell(celP).acsr.fontechk.fluidoP = malha[iP].cell(celP).flui;
            malha[iP].cell(celP).acsr.fontechk.fluidocol = malha[iP].cell(celP).fluicol;
            double rhog = malha[iP].cell(celP).flui.MasEspGas(pres, temp);
            double rhoP = malha[iP].cell(celP).flui.MasEspLiq(pres, temp);
            double rhoC = malha[iP].cell(celP).fluicol.MasEspFlu(pres, temp);
            malha[iP].cell(celP).acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
            malha[iP].cell(celP).acsr.fontechk.titamb = malha[iS].cell(celS).acsr.fontechk.titT;
        }
        for (int i = 0; i < nfontes; i++) {
            int celS = arqRede.conexFR[i].noS;
            int celP = arqRede.conexFR[i].noP;
            malha[iS].cell(celS).acsr.fontechk.tamb = malha[iP].cell(celP).temp;
            malha[iS].cell(celS).acsr.fontechk.pamb = malha[iP].cell(celP).pres;
            malha[iS].cell(celS).acsr.fontechk.betISamb = malha[iP].cell(celP).bet;
            malha[iS].cell(celS).acsr.fontechk.fluidoPamb = malha[iP].cell(celP).flui;
            double alf = malha[iS].cell(celS).alf;
            double bet = malha[iS].cell(celS).bet;
            double pres = malha[iS].cell(celS).pres;
            double temp = malha[iS].cell(celS).temp;
            malha[iS].cell(celS).acsr.fontechk.fluidoP = malha[iS].cell(celS).flui;
            malha[iS].cell(celS).acsr.fontechk.fluidocol = malha[iS].cell(celS).fluicol;
            double rhog = malha[iS].cell(celS).flui.MasEspGas(pres, temp);
            double rhoP = malha[iS].cell(celS).flui.MasEspLiq(pres, temp);
            double rhoC = malha[iS].cell(celS).fluicol.MasEspFlu(pres, temp);
            malha[iS].cell(celS).acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
            malha[iS].cell(celS).acsr.fontechk.titamb = malha[iP].cell(celP).acsr.fontechk.titT;
        }
        for (int i = 0; i < narq; i++) {
            malha[i].SolveTrans(1., 1., 0., nrede);
            for (int kontasnp = 0; kontasnp < malha[i].config().nsnp; kontasnp++) {
                if (fabs((*arqRede.vg1dSP).lixo5 - malha[i].config().tempsnp[kontasnp]) < malha[i].dt) {
                    WriteSnapShot(malha[i], kontasnp, i); // registro do arquivo SNP
                }
            }
        }

        (*arqRede.vg1dSP).lixo5R += dtteste;
        (*arqRede.vg1dSP).dtrede = dtteste;
    }

    for (int i = 0; i < narq; i++) {
        WriteSnapShot(malha[i], malha[i].config().nsnp, i);
    }
}

void RedeParalela(SProd *malha, Rede &arqRede, int narq,
                  Vcr<int> &inativo, int indativo, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
                  int transiente, tipoSimulacao_t tipoSimulacaoMrt, varGlob1D *vG1d) {

    int nfontes = arqRede.nfonteR;
    int iS = arqRede.tramoSecundario;
    int iP = arqRede.tramoPrimario;
    int tD = arqRede.tabelaDinamica;
    for (int i = 0; i < narq; i++) {
        SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                         tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD, 0, 0, 0, 0, 0, 0, arqRede.malha[i].perm);
        malha[i] = temporario;
        malha[i].indTramo = i;
    }
    for (int i = 0; i < nfontes; i++) {
        arqRede.conexFR[i].noP = malha[iP].config().furo[arqRede.conexFR[i].noP].posicP;
        arqRede.conexFR[i].noS = malha[iS].config().furo[arqRede.conexFR[i].noS].posicP;
        int celS = arqRede.conexFR[i].noS;
        int celP = arqRede.conexFR[i].noP;
        if (malha[iS].cell(celS).acsr.tipo == 9 && malha[iP].cell(celP).acsr.tipo == 9) {
            malha[iS].cell(celS).fonteCompartilhada = 1;
            malha[iS].cell(celS).multiFComp = -1.;
            malha[iP].cell(celP).fonteCompartilhada = 1;
            malha[iP].cell(celP).multiFComp = 1.;
        } else {
            NumError("Rede paralela com fontes compartilhadas que não são do tipo 'fonteChoke'");
        }
    }
    if (malha[iS].config().ConContEntrada != 1 &&
        (malha[iS].cell(0).acsr.tipo != 1 && malha[iS].cell(0).acsr.tipo != 2 && malha[iS].cell(0).acsr.tipo != 10))
        NumError("Tramo secundario da Rede paralela deve ter condicao de contorno de pressao, ou fonte de líquido, gás ou massa no início do tramo");

    malha[iP].redeParalelaP = 1;
    malha[iS].redeParalelaS = 1;
    malha[iP].redeParalelaCCsecundario = malha[iS].config().ConContEntrada;
    malha[iS].redeParalelaCCsecundario = malha[iS].config().ConContEntrada;
    double massGas = 0.;
    double massLiqP = 0.;
    double massLiqC = 0.;
    if (malha[iS].config().ConContEntrada == 0) {
        if (malha[iS].cell(0).acsr.tipo == 1) {
            massGas = malha[iS].cell(0).acsr.injg.QGas * malha[iS].cell(0).acsr.injg.FluidoPro.Deng * 1.225 / 86400;
            massLiqP = 0.;
            massLiqC = 0.;
        } else if (malha[iS].cell(0).acsr.tipo == 2) {
            massGas = 0.;
            double beta = malha[iS].cell(0).acsr.injg.QGas * malha[iS].cell(0).acsr.injl.bet;
            double rhoP = 141.5 / (131.5 + malha[iS].cell(0).acsr.injl.FluidoPro.API);
            double rhoC = malha[iS].cell(0).acsr.injl.fluidocol.rholStd;
            massLiqP = (1. - beta) * malha[iS].cell(0).acsr.injl.QLiq * rhoP / 86400;
            massLiqC = beta * malha[iS].cell(0).acsr.injl.QLiq * rhoC / 86400;
        } else if (malha[iS].cell(0).acsr.tipo == 10) {
            massGas = malha[iS].cell(0).acsr.injm.MassG;
            massLiqP = malha[iS].cell(0).acsr.injm.MassP;
            massLiqC = malha[iS].cell(0).acsr.injm.MassC;
        }
        massGas /= nfontes;
        massLiqP /= nfontes;
        massLiqC /= nfontes;
        for (int i = 0; i < nfontes; i++) {
            malha[iP].indFonteRedeParalelaIni.emplace_back(arqRede.conexFR[i].noP);
            malha[iP].fonteMgRedeParalelaIni.emplace_back(massGas);
            malha[iP].fonteMpRedeParalelaIni.emplace_back(massLiqP);
            malha[iP].fonteMcRedeParalelaIni.emplace_back(massLiqC);
        }
    }

    malha[iP].verificaAcopRedeP = 0;
    malha[iP].verificaAcopRedeS = 0;
    malha[iS].verificaAcopRedeP = 0;
    malha[iS].verificaAcopRedeS = 0;

    malha[iP].PrimSecIniRedeP = malha[iP].config().acopPriRedeParalelaini();
    malha[iP].PrimSecFimRedeP = malha[iP].config().acopPriRedeParalelafim();
    malha[iS].PrimSecIniRedeP = malha[iP].PrimSecIniRedeP;
    malha[iS].PrimSecFimRedeP = malha[iP].PrimSecFimRedeP;
    malha[iS].SecPrimIniRedeP = malha[iS].config().acopSecRedeParalelaini();
    malha[iS].SecPrimFimRedeP = malha[iS].config().acopSecRedeParalelafim();
    malha[iP].SecPrimIniRedeP = malha[iS].SecPrimIniRedeP;
    malha[iP].SecPrimFimRedeP = malha[iS].SecPrimFimRedeP;
    if (malha[iP].PrimSecIniRedeP >= 0 && malha[iP].PrimSecFimRedeP >= 0 && malha[iS].PrimSecIniRedeP >= 0 && malha[iS].PrimSecFimRedeP >= 0) {
        malha[iP].verificaAcopRedeP = 1;
        malha[iS].verificaAcopRedeS = 1;
    }
    if (malha[iP].verificaAcopRedeP == 1 && malha[iS].verificaAcopRedeS == 1) {
        conectaPrincipal(malha, iP, iS, 0);
        int nacopP = malha[iP].PrimSecIniRedeP - malha[iP].PrimSecFimRedeP;
        int nacopS = malha[iS].SecPrimFimRedeP - malha[iS].SecPrimIniRedeP;
        if (nacopP != nacopS)
            NumError("Acoplamento termico entre tramo primario e secundario com numero de celulas diferente");
    }

    double titSec;
    if (malha[iS].config().ConContEntrada == 1) {
        titSec = malha[iS].config().CCPres.tit[0];
        malha[iS].hidroTramoSecundario(titSec);
    }

    double erro = 1000.;
    double resulP = 1000;
    double resulP0 = 0.;
    double resulS = 1000;
    double resulS0 = 0.;
    int konta = 0;
    (*vG1d).iterRede = konta;
    double resulmed;
    while (erro > 1.e-5) {
        for (int i = 0; i < nfontes; i++) {
            int celS = arqRede.conexFR[i].noS;
            int celP = arqRede.conexFR[i].noP;
            malha[iP].cell(celP).acsr.fontechk.tamb = malha[iS].cell(celS).temp;
            malha[iP].cell(celP).acsr.fontechk.pamb = malha[iS].cell(celS).pres;
            malha[iP].cell(celP).acsr.fontechk.betISamb = malha[iS].cell(celS).bet;
            malha[iP].cell(celP).acsr.fontechk.fluidoPamb = malha[iS].cell(celS).flui;
            double alf;
            double bet;
            double pres;
            double temp;
            double rhog;
            double rhoP;
            double rhoC;
            if (konta > 0) {
                alf = malha[iP].cell(celP).alf;
                bet = malha[iP].cell(celP).bet;
                pres = malha[iP].cell(celP).pres;
                temp = malha[iP].cell(celP).temp;
                malha[iP].cell(celP).acsr.fontechk.fluidoP = malha[iP].cell(celP).flui;
                malha[iP].cell(celP).acsr.fontechk.fluidocol = malha[iP].cell(celP).fluicol;
                rhog = malha[iP].cell(celP).flui.MasEspGas(pres, temp);
                rhoP = malha[iP].cell(celP).flui.MasEspLiq(pres, temp);
                rhoC = malha[iP].cell(celP).fluicol.MasEspFlu(pres, temp);
                malha[iP].cell(celP).acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
                malha[iP].cell(celP).acsr.fontechk.titamb = malha[iS].cell(celS).acsr.fontechk.titT;
            }
        }
        if (malha[iP].verificaAcopRedeP == 1 && malha[iS].verificaAcopRedeS == 1) {
            int nacopP = malha[iP].PrimSecIniRedeP - malha[iP].PrimSecFimRedeP;
            for (int iacop = 0; iacop <= nacopP; iacop++) {
                malha[iP].cell(malha[iP].PrimSecIniRedeP - iacop).resAcopRedeP =
                    malha[iS].cell(malha[iS].SecPrimIniRedeP + iacop).calor.resGlob;
            }
        }
        resulmed = (0.5 * resulP0 + 0.5 * resulP);
        resulP0 = resulP;
        if (konta > 0 && resulP0 < 1e6 && resulP < 1e6)
            resulP = SolveTramoSolteiro(malha[iP], resulmed);
        else
            resulP = SolveTramoSolteiro(malha[iP]);
        if (malha[iP].verificaAcopRedeP == 1 && malha[iS].verificaAcopRedeS == 1) {
            int nacopP = malha[iP].PrimSecIniRedeP - malha[iP].PrimSecFimRedeP;
            for (int iacop = 0; iacop <= nacopP; iacop++) {
                malha[iS].cell(malha[iS].SecPrimIniRedeP + iacop).fluxcalAcopRedeP =
                    malha[iP].cell(malha[iP].PrimSecIniRedeP - iacop).fluxcalmed;
            }
        }
        for (int i = 0; i < nfontes; i++) {
            int celS = arqRede.conexFR[i].noS;
            int celP = arqRede.conexFR[i].noP;
            malha[iS].cell(celS).acsr.fontechk.tamb = malha[iP].cell(celP).temp;
            malha[iS].cell(celS).acsr.fontechk.pamb = malha[iP].cell(celP).pres;
            malha[iS].cell(celS).acsr.fontechk.betISamb = malha[iP].cell(celP).bet;
            malha[iS].cell(celS).acsr.fontechk.fluidoPamb = malha[iP].cell(celP).flui;
            double alf = malha[iS].cell(celS).alf;
            double bet = malha[iS].cell(celS).bet;
            double pres = malha[iS].cell(celS).pres;
            double temp = malha[iS].cell(celS).temp;
            malha[iS].cell(celS).acsr.fontechk.fluidoP = malha[iS].cell(celS).flui;
            malha[iS].cell(celS).acsr.fontechk.fluidocol = malha[iS].cell(celS).fluicol;
            double rhog = malha[iS].cell(celS).flui.MasEspGas(pres, temp);
            double rhoP = malha[iS].cell(celS).flui.MasEspLiq(pres, temp);
            double rhoC = malha[iS].cell(celS).fluicol.MasEspFlu(pres, temp);
            malha[iS].cell(celS).acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
            malha[iS].cell(celS).acsr.fontechk.titamb = malha[iP].cell(celP).acsr.fontechk.titT;
        }
        resulmed = (0.5 * resulS0 + 0.5 * resulS);
        resulS0 = resulS;
        if (malha[iS].config().ConContEntrada == 1) {
            if (konta > 0 && resulS0 < 1e10 && resulS < 1e10)
                resulS = SolveTramoSolteiro(malha[iS], resulmed);
            else
                resulS = SolveTramoSolteiro(malha[iS]);
        } else {
            int indPartida = arqRede.conexFR[0].noS;
            int indPartidaP = arqRede.conexFR[0].noP;
            double pchute = chutePresRedeParalelaSec(malha, iP, iS, malha[iP].cell(indPartidaP).pres, indPartida);
            resulS = SolveTramoSolteiro(malha[iS], pchute);
        }
        if (konta > 0) {
            if (fabs(resulS0) > 1e-15 && fabs(resulP0) > 1e-15)
                erro = fabs(resulP - resulP0) / resulP0 + fabs(resulS - resulS0) / resulS0;
            else if (fabs(resulP0) > 1e-15)
                erro = fabs(resulP - resulP0) / resulP0 + fabs(resulS - resulS0);
            else if (fabs(resulS0) > 1e-15)
                erro = fabs(resulP - resulP0) + fabs(resulS - resulS0) / resulS0;
            else
                erro = fabs(resulP - resulP0) + fabs(resulS - resulS0);
        }
        konta++;
        (*vG1d).iterRede = konta;
    }

    int nrede = 0;
    for (int i = 0; i < 2; i++) {
        malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo, nrede);
        malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo, nrede);
        if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo,nrede);
        if(malha[i].config().nCelUnit>0){
        	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
        		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo,nrede);
        }
        malha[i].kimpT++;
        for (int j = 0; j < malha[i].config().ntendp; j++) {
            malha[i].ImprimeTrendPCab(j, nrede);
            malha[i].config().imprimeTrend(malha[i].celula, malha[i].MatTrendP[j], 0, j, 0);
            malha[i].ImprimeTrendP(j, nrede);
        }
        if (malha[i].config().lingas == 1) {
            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo);
            for (int j = 0; j < malha[i].config().ntendg; j++) {
                malha[i].ImprimeTrendGCab(j, nrede);
                malha[i].config().imprimeTrendG(malha[i].celulaG, malha[i].MatTrendG[j], 0, j, 0, 0);
                malha[i].ImprimeTrendG(j, nrede);
            }
        }
        // enterramento
        for (int j = 0; j <= malha[i].ncel; j++) {
            if (malha[i].cell(j).calor.difus2D == 1) {
                malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
            }
        }
        if (malha[i].config().lingas == 1)
            malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo, nrede);
    }

    if (arqRede.chaveredeT == 1) {
        for (int i = 0; i < 2; i++) {
            // preparacao de algumas variaveis para a partida na simulacao transiente:
            if (malha[i].config().ConContEntrada == 1) {
                if (malha[i].cell(0).acsr.tipo == 2) {
                    malha[i].cell(0).acsr.injl.QLiq = 0.;
                    malha[i].cell(0).MC = malha[i].cell(0).fontemassLR + malha[i].cell(0).fontemassCR + malha[i].cell(0).fontemassGR;
                    malha[i].cell(0).MCini = malha[i].cell(0).MC;
                    malha[i].cell(0).Mliqini = malha[i].cell(0).fontemassLR + malha[i].cell(0).fontemassCR;
                    malha[i].cell(0).Mliqini0 = malha[i].cell(0).Mliqini;
                    double rholMix = malha[i].cell(0).bet * malha[i].cell(0).fluicol.MasEspFlu(
                                                                  malha[i].cell(0).pres, malha[i].cell(0).temp);
                    rholMix += (1. - malha[i].cell(0).bet) * malha[i].cell(0).flui.MasEspLiq(
                                                                   malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QL = malha[i].cell(0).Mliqini / rholMix;
                    malha[i].cell(0).QLini = malha[i].cell(0).QL;
                    double rhog = malha[i].cell(0).flui.MasEspGas(
                        malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QG = (malha[i].cell(0).MC - malha[i].cell(0).Mliqini) / rhog;
                    malha[i].cell(0).QGini = malha[i].cell(0).QG;
                    malha[i].cell(0).fontemassLR = 0.;
                    malha[i].cell(0).fontemassCR = 0.;
                    malha[i].cell(0).fontemassGR = 0.;
                } else if (malha[i].cell(0).acsr.tipo == 1) {
                    malha[i].cell(0).acsr.injg.QGas = 0.;
                    malha[i].cell(0).MC = malha[i].cell(0).fontemassGR;
                    malha[i].cell(0).MCini = malha[i].cell(0).MC;
                    malha[i].cell(0).Mliqini = 0.;
                    malha[i].cell(0).Mliqini0 = malha[i].cell(0).Mliqini;
                    malha[i].cell(0).QL = 0.;
                    malha[i].cell(0).QLini = malha[i].cell(0).QL;
                    double rhog = malha[i].cell(0).flui.MasEspGas(
                        malha[i].cell(0).pres, malha[i].cell(0).temp);
                    malha[i].cell(0).QG = (malha[i].cell(0).MC - malha[i].cell(0).Mliqini) / rhog;
                    malha[i].cell(0).QGini = malha[i].cell(0).QG;
                    malha[i].cell(0).fontemassLR = 0.;
                    malha[i].cell(0).fontemassCR = 0.;
                    malha[i].cell(0).fontemassGR = 0.;
                }
            }
            for (int j = 0; j <= malha[i].ncel; j++) {
                malha[i].cell(j).rpC = malha[i].cell(j).rpCi =
                    malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).pres, malha[i].cell(j).temp);
                malha[i].cell(j).rgC = malha[i].cell(i).rgCi =
                    malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                malha[i].cell(j).rcC = malha[i].cell(j).rcCi =
                    malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);
                malha[i].cell(j).rpL = malha[i].cell(j).rpLi =
                    malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                malha[i].cell(j).rgL = malha[i].cell(j).rgLi =
                    malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).presL, malha[i].cell(j).tempL);
                malha[i].cell(j).rcL = malha[i].cell(j).rcLi =
                    malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).presL, malha[i].cell(j).tempL);

                malha[i].cell(j).mipC = malha[i].cell(j).flui.ViscOleo(malha[i].cell(j).pres, malha[i].cell(j).temp);
                malha[i].cell(j).migC = malha[i].cell(j).flui.ViscGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
                malha[i].cell(j).micC = malha[i].cell(j).fluicol.VisFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);

                if (j > 0) {
                    malha[i].cell(j - 1).rpR = malha[i].cell(j - 1).rpRi = malha[i].cell(j).rpC;
                    malha[i].cell(j - 1).rgR = malha[i].cell(j - 1).rgRi = malha[i].cell(j).rgC;
                    malha[i].cell(j - 1).rcR = malha[i].cell(j - 1).rcRi = malha[i].cell(j).rcC;

                    malha[i].cell(j - 1).mipR = malha[i].cell(j).mipC;
                    malha[i].cell(j - 1).migR = malha[i].cell(j).migC;
                    malha[i].cell(j - 1).micR = malha[i].cell(j).micC;
                }
                if (j == malha[i].ncel) {
                    malha[i].cell(j).rpR = malha[i].cell(j).rpRi = malha[i].cell(j).rpC;
                    malha[i].cell(j).rgR = malha[i].cell(j).rgRi = malha[i].cell(j).rgC;
                    malha[i].cell(j).rcR = malha[i].cell(j).rcRi = malha[i].cell(j).rcC;

                    malha[i].cell(j).mipR = malha[i].cell(j).mipC;
                    malha[i].cell(j).migR = malha[i].cell(j).migC;
                    malha[i].cell(j).micR = malha[i].cell(j).micC;
                }
            }
        }

        SolveRedeParalelaTrans(malha, arqRede, nrede);
    }
}

void RedeAnelGL(SProd *malha, Rede &arqRede, int narq,
                Vcr<int> &inativo, int indativo, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
                int transiente, tipoSimulacao_t tipoSimulacaoMrt, varGlob1D *vG1d) {
    int nfontes = narq - 1;
    vector<fonteposic> dreno;
    double *compfonte;
    compfonte = new double[nfontes];
    int *indfonte;
    indfonte = new int[nfontes];
    int *indtramo;
    indtramo = new int[nfontes + 1];
    int *posicfonte;
    posicfonte = new int[nfontes];
    int indAnel = 0;
    int konta = 0;
    int tD = arqRede.tabelaDinamica;
    for (int i = 0; i < narq; i++) {
        if (arqRede.malha[i].tipoanel == 0) {
            SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog,
                             tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD, 0, 0, 0, 0, 0, 0, arqRede.malha[i].perm);
            compfonte[konta] = arqRede.malha[i].compfonte;
            indfonte[konta] = i;
            indtramo[i] = konta;
            konta++;
            temporario.noextremo = 1;
            temporario.noinicial = 1;
            temporario.derivaAnel = 1;
            malha[i] = temporario;
            malha[i].indTramo = i;
            if (dreno.size() == 0) {
                fonteposic temp;
                temp.nmani = 1;
                temp.mani[temp.nmani - 1] = i;
                temp.posicmani = arqRede.malha[i].compfonte;
                dreno.push_back(temp);
            } else {
                int dimen = dreno.size();
                for (int j = 0; j < dimen; j++) {
                    if (fabs(dreno[j].posicmani - arqRede.malha[i].compfonte) < (*arqRede.vg1dSP).localtiny) {
                        dreno[j].nmani++;
                        dreno[j].mani[dreno[j].nmani - 1] = i;
                        j = dimen;
                    } else {
                        fonteposic temp;
                        temp.nmani = 1;
                        temp.mani[temp.nmani - 1] = i;
                        temp.posicmani = arqRede.malha[i].compfonte;
                        dreno.push_back(temp);
                        j = dimen;
                    }
                }
            }
        } else
            indAnel = i;
    }
    double *compdreno;
    compdreno = new double[dreno.size()];
    int *inddreno;
    inddreno = new int[dreno.size()];
    int *posicdreno;
    posicdreno = new int[dreno.size()];
    for (int i = 0; i < dreno.size(); i++)
        compdreno[i] = dreno[i].posicmani;
    SProd temporario(pathArqEntrada + arqRede.impfiles[indAnel], nomeArquivoLog,
                     tipoValidacaoMrt, tipoSimulacaoMrt, vG1d, tD, 0, 0, 0, compdreno, posicdreno, dreno.size(), arqRede.malha[indAnel].perm);
    for (int i = 0; i < dreno.size(); i++)
        for (int j = 0; j < dreno[i].nmani; j++)
            posicfonte[dreno[i].mani[j]] = posicdreno[i];
    temporario.noextremo = 1;
    temporario.noinicial = 1;
    temporario.derivaAnel = 0;
    malha[indAnel] = temporario;
    for (int i = 0; i < dreno.size(); i++) {
        int iposp = posicdreno[i];
        InjGas injgasMRT(0, 0, malha[indAnel].config().flug);
        malha[indAnel].cell(iposp).acsr.tipo = 1;
        malha[indAnel].cell(iposp).acsr.injg = injgasMRT;
        malha[indAnel].cell(iposp).acsr.injg.FluidoPro = malha[indAnel].cell(iposp).flui;
    }
    malha[indAnel].indTramo = indAnel;
    malha[indAnel].chokeSup.AreaGarg = 0.;
    if (malha[indAnel].config().ConContEntrada == 1) {
        malha[indAnel].tempE = malha[indAnel].config().CCPres.temperatura[0];
        malha[indAnel].presE = malha[indAnel].config().CCPres.pres[0];
        malha[indAnel].titE = malha[indAnel].config().CCPres.tit[0];
        malha[indAnel].betaE = malha[indAnel].config().CCPres.bet[0];
    }

    for (int i = 0; i < narq; i++) {
        malha[i].modoPerm = 1;
        if (arqRede.malha[i].tipoanel == 0) {
            int ncelg = malha[i].ncelGas;
            for (int j = 0; j <= ncelg; j++) {
                malha[i].gasCell(j).flui.Deng = malha[indAnel].cell(0).flui.Deng;
                malha[i].gasCell(j).flui.yco2 = malha[indAnel].cell(0).flui.yco2;
                malha[i].gasCell(j).flui.corrC = malha[indAnel].cell(0).flui.corrC;
                malha[i].gasCell(j).flui.RenovaFluido();
            }
        }
    }

    for (int i = 0; i < dreno.size(); i++)
        malha[indAnel].cell(posicdreno[i]).acsr.injg.QGas = 0.;
    for (int i = 0; i < narq; i++) {
        if (i != indAnel) {
            if (malha[i].config().gasinj.chuteVaz == 0 && malha[indAnel].config().ConContEntrada == 1) {
                malha[i].config().gasinj.vazgas[0] = 150000. * malha[i].gasCell(0).duto.area / (*arqRede.vg1dSP).arearef;
                int ifposic = posicfonte[indtramo[i]];
                malha[indAnel].cell(ifposic).acsr.injg.QGas -= malha[i].config().gasinj.vazgas[0];
            } else if (malha[i].config().gasinj.chuteVaz == 1) {
                int ifposic = posicfonte[indtramo[i]];
                malha[indAnel].cell(ifposic).acsr.injg.QGas -= malha[i].config().gasinj.vazgas[0];
            }
        }
    }
    if (malha[indAnel].config().ConContEntrada == 1) {
        malha[indAnel].cell(0).pres = malha[indAnel].config().CCPres.pres[0];
        malha[indAnel].cell(0).acsr.injg.temp = malha[indAnel].config().CCPres.temperatura[0];
        malha[indAnel].cell(0).flui.RGO = malha[indAnel].cell(0).acsr.injg.FluidoPro.RGO =
            1e6 + (*arqRede.vg1dSP).localtiny;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.Deng = malha[indAnel].cell(0).flui.Deng;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.yco2 = malha[indAnel].cell(0).flui.yco2;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.RenovaFluido();
        malha[indAnel].cell(0).fluiL = &malha[indAnel].cell(0).acsr.injg.FluidoPro;
        double pres = malha[indAnel].cell(0).pres;
        double temp = malha[indAnel].cell(0).temp = malha[indAnel].cell(0).calor.Textern1;
        for (int i = 1; i <= malha[indAnel].ncel; i++) {
            double rhog = malha[indAnel].cell(i).flui.MasEspGas(pres, temp);
            double tet = malha[indAnel].cell(i).duto.teta;
            pres = malha[indAnel].cell(i).pres = pres - rhog * 9.81 * malha[indAnel].cell(i).dx * sin(tet) / 98066.5;
            temp = malha[indAnel].cell(i).temp = malha[indAnel].cell(i).calor.Textern1;
        }
    } else {
        malha[indAnel].cell(0).acsr.injg.FluidoPro.Deng = malha[indAnel].cell(0).flui.Deng;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.yco2 = malha[indAnel].cell(0).flui.yco2;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.corrC = malha[indAnel].cell(0).flui.corrC;
        malha[indAnel].cell(0).acsr.injg.FluidoPro.RenovaFluido();
        double qori = 0.;
        for (int i = 0; i < malha[indAnel].config().ninjgas; i++) {
            int iinjg = malha[indAnel].config().fonteg[i].posicP;
            qori += malha[indAnel].cell(iinjg).acsr.injg.QGas;
        }
        double sumid = -qori / dreno.size();
        int fontemax = 0;
        int pocodist = -1;
        for (int i = malha[indAnel].config().ninjgas; i < malha[indAnel].config().ninjgas + dreno.size(); i++) {
            int ifposic = posicdreno[i - malha[indAnel].config().ninjgas];
            if (ifposic > fontemax) {
                fontemax = ifposic;
                pocodist = i - malha[indAnel].config().ninjgas;
            }
            malha[indAnel].cell(ifposic).acsr.injg.QGas = sumid;
            malha[indAnel].cell(ifposic).acsr.injg.temp =
                malha[indAnel].cell(ifposic).calor.Textern1;
        }
        malha[indAnel].cell(0).flui.RGO =
            malha[indAnel].cell(0).acsr.injg.FluidoPro.RGO = 1e6 + (*arqRede.vg1dSP).localtiny; // verificar
        int itramo = dreno[pocodist].mani[0];
        malha[itramo].config().gasinj.vazgas[0] = -sumid / dreno[pocodist].nmani;
        malha[itramo].gasCell(0).tipoCC = 1;
        if (malha[itramo].config().chokep.abertura[0] > 0.6)
            malha[itramo].buscaProdPfundoPerm();
        else
            malha[itramo].buscaProdPfundoPerm2();

        int ncel = malha[indAnel].ncel;
        malha[indAnel].cell(ncel).pres = malha[itramo].gasCell(0).pres;
        double pres = malha[indAnel].cell(ncel).pres;
        double temp = malha[indAnel].cell(ncel).temp = malha[indAnel].cell(ncel).calor.Textern1;
        for (int i = ncel - 1; i >= 0; i--) {
            double rhog = malha[indAnel].cell(i + 1).flui.MasEspGas(pres, temp);
            double tet = malha[indAnel].cell(i + 1).duto.teta;
            pres = malha[indAnel].cell(i).pres = pres + rhog * 9.81 * malha[indAnel].cell(i + 1).dx * sin(tet) / 98066.5;
            temp = malha[indAnel].cell(i).temp = malha[indAnel].cell(i).calor.Textern1;
        }
    }
    double erro = 1000.;
    int iter = 0;

    if (malha[indAnel].config().ConContEntrada == 0) {
        double pchute = malha[indAnel].cell(0).pres;
        double pchute2 = pchute;
        double pbuffer = pchute;
        erro = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicdreno, indAnel, pchute, iter, dreno);
        cout << "iteracao= " << iter << " norma= " << erro << endl;
        if (erro < 0.) {
            while (erro < 0.) {
                pchute2 = 0.9 * pchute2;
                iter++;
                erro = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicdreno, indAnel, pchute2, iter,
                                dreno);
                if (erro < 0.)
                    pbuffer = pchute2;
                cout << "iteracao= " << iter << " norma= " << erro << endl;
            }
        } else {
            while (erro > 0.) {
                pchute2 = 1.1 * pchute2;
                iter++;
                erro = objetCC0(malha, narq, nfontes, indfonte, indtramo, posicdreno, indAnel, pchute2, iter,
                                dreno);
                if (erro > 0.)
                    pbuffer = pchute2;
                cout << "iteracao= " << iter << " norma= " << erro << endl;
            }
        }
        pchute = pbuffer;
        zriddr(malha, pchute, pchute2, narq, nfontes, indfonte, indtramo, posicdreno, indAnel,
               dreno);

    } else {
        while (fabs(erro) > 0.0001) {
            calcPeriAnelGL(malha, narq, nfontes, indfonte, indtramo, posicdreno, indAnel, iter, dreno);
            erro = calcErroGL(malha, narq, nfontes, posicdreno, indAnel, dreno);
            double chutemass = 0.;
            for (int i = 0; i < malha[indAnel].config().ninjgas + dreno.size(); i++) {
                int ifposic;
                if (i < malha[indAnel].config().ninjgas)
                    ifposic = malha[indAnel].config().fonteg[i].posicP;
                else
                    ifposic = posicdreno[i - malha[indAnel].config().ninjgas];
                chutemass -= malha[indAnel].cell(ifposic).acsr.injg.QGas;
            }
            malha[indAnel].marchaProdPresPres1(chutemass);
            iter++;
            cout << "iteracao= " << iter << " norma= " << erro << endl;
        }
    }

    cout << "Convergencia atingida" << endl;
    for (int i = 0; i < narq; i++)
        if (inativo[i] == 0) {
            malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo);
            malha[i].config().resumoPermanente(malha[i].celula, malha[i].celulaG, malha[i].pGSup, malha[i].presiniG, malha[i].indTramo);
            if(malha[i].config().nintermi>0)malha[i].config().resumoIntermitencia(malha[i].celula, malha[i].indTramo);
            if(malha[i].config().nCelUnit>0){
            	for(int iCelUni=0; iCelUni<malha[i].config().nCelUnit; iCelUni++)
            		malha[i].config().relatorioCelulaUnitaria(malha[i].celula,malha[i].config().celUnit[iCelUni].posicP, malha[i].indTramo);
            }
            // enterramento
            for (int j = 0; j <= malha[i].ncel; j++) {
                if (malha[i].cell(j).calor.difus2D == 1) {
                    malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                }
            }
            if (malha[i].config().lingas == 1)
                malha[i].config().imprimeProfileG(malha[i].celulaG, malha[i].flutG, 0, malha[i].indTramo);
        }

    for (int i = 0; i <= malha[indAnel].ncel; i++) {
        malha[indAnel].cell(i).alf = 1.;
        malha[indAnel].cell(i).alfL = 1.;
    }

    for (int i = 0; i < narq; i++) {

        for (int j = 0; j <= malha[i].ncel; j++) {
            malha[i].cell(j).tempini = malha[i].cell(j).temp;
            malha[i].cell(j).presLini = malha[i].cell(j).presL;
            malha[i].cell(j).presini = malha[i].cell(j).pres;
            malha[i].cell(j).MRini = malha[i].cell(j).MR;
            malha[i].cell(j).MRini = malha[i].cell(j).MR;
            malha[i].cell(j).MliqiniR0 = malha[i].cell(j).MliqiniR;
            malha[i].cell(j).presLiniBuf = malha[i].cell(j).presL;
            malha[i].cell(j).MRiniBuf = malha[i].cell(j).MR;
            malha[i].cell(j).alfLini = malha[i].cell(j).alfL;
            malha[i].cell(j).alfRini = malha[i].cell(j).alfR;
            malha[i].cell(j).alfini = malha[i].cell(j).alf;
            malha[i].cell(j).betLini = malha[i].cell(j).betL;
            malha[i].cell(j).betRini = malha[i].cell(j).betR;
            malha[i].cell(j).betini = malha[i].cell(j).bet;
            malha[i].cell(j).FWini = malha[i].cell(j).FW;
            malha[i].cell(j).QLRini = malha[i].cell(j).QLR;
            malha[i].cell(j).alfPigEini = malha[i].cell(j).alfPigE;
            malha[i].cell(j).alfPigDini = malha[i].cell(j).alfPigD;
            malha[i].cell(j).betPigEini = malha[i].cell(j).betPigE;
            malha[i].cell(j).betPigDini = malha[i].cell(j).betPigD;

            malha[i].cell(j).rpC = malha[i].cell(j).rpCi =
                malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).pres, malha[i].cell(j).temp);
            malha[i].cell(j).rgC = malha[i].cell(j).rgCi =
                malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
            malha[i].cell(j).rcC = malha[i].cell(j).rcCi =
                malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);
            malha[i].cell(j).rpL = malha[i].cell(j).rpLi =
                malha[i].cell(j).flui.MasEspLiq(malha[i].cell(j).presL, malha[i].cell(j).tempL);
            malha[i].cell(j).rgL = malha[i].cell(j).rgLi =
                malha[i].cell(j).flui.MasEspGas(malha[i].cell(j).presL, malha[i].cell(j).tempL);
            malha[i].cell(j).rcL = malha[i].cell(j).rcLi =
                malha[i].cell(j).fluicol.MasEspFlu(malha[i].cell(j).presL, malha[i].cell(j).tempL);

            malha[i].cell(j).mipC = malha[i].cell(j).flui.ViscOleo(malha[i].cell(j).pres, malha[i].cell(j).temp);
            malha[i].cell(j).migC = malha[i].cell(j).flui.ViscGas(malha[i].cell(j).pres, malha[i].cell(j).temp);
            malha[i].cell(j).micC = malha[i].cell(j).fluicol.VisFlu(malha[i].cell(j).pres, malha[i].cell(j).temp);

            if (j > 0) {
                malha[i].cell(j - 1).rpR = malha[i].cell(j - 1).rpRi = malha[i].cell(j).rpC;
                malha[i].cell(j - 1).rgR = malha[i].cell(j - 1).rgRi = malha[i].cell(j).rgC;
                malha[i].cell(j - 1).rcR = malha[i].cell(j - 1).rcRi = malha[i].cell(j).rcC;

                malha[i].cell(j - 1).mipR = malha[i].cell(j).mipC;
                malha[i].cell(j - 1).migR = malha[i].cell(j).migC;
                malha[i].cell(j - 1).micR = malha[i].cell(j).micC;
            }
            if (j == malha[i].ncel) {
                malha[i].cell(j).rpR = malha[i].cell(j).rpRi = malha[i].cell(j).rpC;
                malha[i].cell(j).rgR = malha[i].cell(j).rgRi = malha[i].cell(j).rgC;
                malha[i].cell(j).rcR = malha[i].cell(j).rcRi = malha[i].cell(j).rcC;

                malha[i].cell(j).mipR = malha[i].cell(j).mipC;
                malha[i].cell(j).migR = malha[i].cell(j).migC;
                malha[i].cell(j).micR = malha[i].cell(j).micC;
            }
        }
    }
    malha[indAnel].chokeSup.AreaGarg = 0.;
    for (int i = 0; i < malha[indAnel].config().chokep.parserie; i++)
        malha[indAnel].config().chokep.abertura[i] = 0.;
    (*arqRede.vg1dSP).modoTransiente = 1;
    (*arqRede.vg1dSP).RGOMax = 14000.;
    if (arqRede.chaveredeT == 1)
        TransAnel(narq, nfontes, indfonte, indtramo, posicdreno, indAnel, dreno, malha, arqRede);

    delete[] indfonte;
    delete[] indtramo;
    delete[] posicfonte;
    delete[] compfonte;
    delete[] inddreno;
    delete[] posicdreno;
    delete[] compdreno;
}

int chutePresRedeInj(int indprod, SProd *malha, Rede &arqRede, double chutehol,
                     Vcr<double> &razcolet, Vcr<double> &prescolet) {
    double vaz = (*arqRede.vg1dSP).somavaz * malha[indprod].cell(0).duto.area / (*arqRede.vg1dSP).somaarea;
    double presno = malha[indprod].hidroreversoInj(chutehol, vaz);
    int contarecur = 1;
    for (int i = 0; i < arqRede.malha[indprod].nafluente; i++) {
        int aflu = arqRede.malha[indprod].afluente[i];
        malha[aflu].config().condpocinj.presfundo = presno;
        razcolet[aflu] += 1.;
        prescolet[aflu] += presno;
        if (arqRede.malha[aflu].nafluente > 0)
            contarecur += chutePresRedeInj(aflu, malha, arqRede, chutehol, razcolet, prescolet);
    }
    return contarecur;
}

double cicloRedeInj(SProd *malha, Rede &arqRede, Vcr<int> &inativo, int indativo) {
    double norma = 0.;
    int nnos = 0;
    int narq = arqRede.nsisprod - 0 * indativo;
    (*arqRede.vg1dSP).relax = arqRede.relax;
    Vcr<int> Resolv(narq, 0);
    Vcr<int> IndNorma(narq, 0);
    int ciclomalha = 0;
    double valor;
    cout << "#################Iniciando ciclo iterativo: " << (*arqRede.vg1dSP).iterRede << endl;
    while (ciclomalha < narq - 1) {
        for (int i = 0; i < narq; i++) {
            if (arqRede.malha[i].nafluente == 0 && Resolv[i] == 0 && inativo[i] == 0) {
                malha[i].modoPerm = 1;
                if (malha[i].config().condpocinj.CC == 3)
                    valor = malha[i].buscaInjPfundoPerm1();
                else
                    valor = malha[i].buscaInjPfundoPerm5();
                if (valor < -1e9 || valor > 1e9) {
                    inativo[i] = 1;
                    inativoColetor(i, arqRede, inativo, narq);
                    if (valor < -1e9) {
                        aviso(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    } else {
                        aviso2(i);
                        (*arqRede.vg1dSP).iterRede = 200;
                    }
                    (*arqRede.vg1dSP).restartRede = 1;
                    return 9000.;
                }
                Resolv[i] = 1;
                ciclomalha++;
                cout << "convergencia tramo: " << i << endl;
            }
        }
        for (int i = 0; i < narq; i++) {
            if (i >= 10) {
            }
            int naflu = arqRede.malha[i].nafluente;
            int cicloaflu = 0;
            if (arqRede.malha[i].nafluente > 0 && Resolv[i] == 0) {
                for (int j = 0; j < narq; j++) {
                    for (int k = 0; k < naflu; k++) {
                        if (arqRede.malha[i].afluente[k] == j && Resolv[j] == 1)
                            cicloaflu++;
                    }
                }
                if (cicloaflu == naflu && inativo[i] == 0) {
                    Vcr<double> temp(naflu);
                    Vcr<double> cpl(naflu);
                    Vcr<double> Mliq(naflu);
                    Vcr<double> Qliq(naflu);
                    Vcr<double> Qcomp(naflu);
                    Vcr<double> Mgas(naflu);
                    Vcr<double> Qgas(naflu);

                    double cpmist = 0;
                    double tempmist = 0;
                    double qlmistStd = 0.;
                    double qgmistStd = 0.;
                    if (malha[i].config().flashCompleto == 0)
                        malha[i].cell(0).acsr.injl.fluidocol = malha[i].cell(0).fluicol;
                    else
                        malha[i].cell(0).acsr.injg.FluidoPro = malha[i].cell(0).flui;

                    for (int k = 0; k < naflu; k++) {
                        int ind = arqRede.malha[i].afluente[k];
                        if (inativo[ind] == 0) {
                            int fim = malha[ind].ncel - 1;
                            temp[k] = (malha[ind].cell(fim + 1).dx * malha[ind].cell(fim + 1).temp + malha[ind].cell(fim + 1).dxL * malha[ind].cell(fim + 1).tempL) / (malha[ind].cell(fim + 1).dx + malha[ind].cell(fim + 1).dxL);
                            double pres = malha[ind].cell(fim + 1).presaux;
                            cpl[k] = malha[ind].cell(fim).fluicol.CalorLiq(pres, temp[k]);
                            if (malha[ind].config().flashCompleto == 0) {
                                Mliq[k] = malha[ind].cell(fim).Mliqini;
                                Qliq[k] = malha[ind].cell(fim + 1).QL;
                                Qcomp[k] = Qliq[k];
                                cpmist += (Mliq[k] * cpl[k]);
                                tempmist += (Mliq[k] * cpl[k]) * temp[k];
                                qlmistStd += (Qcomp[k] * malha[i].cell(0).fluicol.MasEspFlu(pres, temp[k]) / malha[i].cell(0).fluicol.MasEspFlu(1.001, 15.));
                            } else {
                                Mgas[k] = malha[ind].cell(fim).MC - malha[ind].cell(fim).Mliqini;
                                Qliq[k] = malha[ind].cell(fim + 1).QG;
                                Qcomp[k] = 0.;
                                cpmist += (Mgas[k] * cpl[k]);
                                tempmist += (Mgas[k] * cpl[k]) * temp[k];
                                qgmistStd += (Qgas[k] * malha[i].cell(0).flui.MasEspGas(pres, temp[k]) / malha[i].cell(0).flui.MasEspGas(1.001, 15.));
                            }
                        }
                    }
                    if (cpmist > (*arqRede.vg1dSP).localtiny)
                        tempmist = tempmist / cpmist;
                    else
                        tempmist = 0.;
                    if (malha[i].config().flashCompleto == 0)
                        qlmistStd *= 86400;
                    else
                        qgmistStd *= 86400;

                    int ind = arqRede.malha[i].afluente[0];
                    int nderiva = arqRede.malha[ind].ncoleta;
                    Vcr<int> ordCol(nderiva);
                    vector<double> dcol;
                    for (int icol = 0; icol < nderiva; icol++) {
                        int aux = arqRede.malha[ind].coleta[icol];
                        dcol.push_back(malha[aux].cell(0).duto.dia);
                    }
                    sort(dcol.begin(), dcol.end());
                    Vcr<int> carregado(nderiva, 0);
                    for (int icol = 0; icol < nderiva; icol++) {
                        int aux = arqRede.malha[ind].coleta[icol];
                        for (int icol2 = 0; icol2 < nderiva; icol2++) {
                            if (dcol[icol2] == malha[aux].cell(0).duto.dia && carregado[icol2] == 0) {
                                ordCol[icol2] = aux;
                                carregado[icol2] = 1;
                                icol2 = nderiva;
                            }
                        }
                    }
                    for (int icol = 0; icol < nderiva; icol++) {
                        if (icol == nderiva - 1) {
                            int aux = ordCol[icol];
                            if (malha[aux].config().flashCompleto == 0) {
                                malha[aux].cell(0).acsr.injl.bet = 1.;
                                malha[aux].cell(0).acsr.injl.temp = tempmist;
                                malha[aux].cell(0).acsr.injl.QLiq = qlmistStd;
                            } else {
                                malha[aux].cell(0).acsr.injg.temp = tempmist;
                                malha[aux].cell(0).acsr.injg.QGas = qgmistStd;
                            }

                            double pini = malha[aux].cell(0).pres;
                            malha[aux].modoPerm = 1;
                            if (malha[aux].cell(malha[aux].ncel).acsr.tipo == 3) {
                                malha[aux].config().condpocinj.CC = 0;
                                valor = malha[aux].buscaInjPfundoPerm2();
                            } else {
                                malha[aux].config().condpocinj.CC = 5;
                                valor = malha[aux].buscaInjPfundoPerm5();
                            }
                            if (valor < -1e9 || valor > 1e9) {
                                inativo[aux] = 1;
                                inativoColetor(aux, arqRede, inativo, narq);
                                inativoAfluente(aux, arqRede, inativo);
                                if (valor < -1e9) {
                                    aviso(aux);
                                    (*arqRede.vg1dSP).iterRede = 200;
                                } else {
                                    aviso2(aux);
                                    (*arqRede.vg1dSP).iterRede = 200;
                                }
                                (*arqRede.vg1dSP).restartRede = 1;
                                return 9000.;
                            }
                            Resolv[aux] = 1;
                            cout << "convergencia tramo: " << aux << endl;
                            ciclomalha++;
                            if (inativo[aux] == 0) {
                                if (IndNorma[aux] == 0) {
                                    norma += pow(malha[aux].cell(0).pres - pini, 2.);
                                    IndNorma[aux] = 1;
                                    nnos++;
                                }
                                for (int iaflu = 0; iaflu < arqRede.malha[aux].nafluente; iaflu++) {
                                    int indaflu = arqRede.malha[aux].afluente[iaflu];
                                    pini = malha[indaflu].config().condpocinj.presfundo;
                                    if (IndNorma[indaflu] == 0) {
                                        IndNorma[indaflu] = 1;
                                        norma += pow(malha[indaflu].config().condpocinj.presfundo - pini, 2.);
                                        nnos++;
                                    }
                                    malha[indaflu].config().condpocinj.presfundo = ((*arqRede.vg1dSP).relax) * malha[aux].cell(0).pres + (1. - (*arqRede.vg1dSP).relax) * malha[indaflu].config().condpocinj.presfundo;
                                }
                            }
                            int indaflu = arqRede.malha[aux].afluente[0];
                            for (int icol2 = 0; icol2 < nderiva - 1; icol2++) {
                                malha[ordCol[icol2]].cell(0).pres = malha[indaflu].config().condpocinj.presfundo;
                                malha[ordCol[icol2]].pGSup = malha[indaflu].config().condpocinj.presfundo;
                                malha[ordCol[icol2]].config().condpocinj.presinj = malha[ordCol[icol2]].pGSup;
                            }
                        } else {
                            double chutemass = 0.;
                            double area = 0.;
                            int aux = ordCol[icol];
                            if (malha[aux].config().flashCompleto == 0) {
                                malha[aux].cell(0).acsr.injl.bet = 1.;
                                malha[aux].cell(0).acsr.injl.temp = tempmist;
                            } else
                                malha[aux].cell(0).acsr.injg.temp = tempmist;
                            malha[aux].config().condpocinj.CC = 1;

                            for (int icol2 = icol; icol2 < nderiva; icol2++)
                                area += malha[ordCol[icol2]].cell(0).duto.area;
                            if (malha[aux].config().flashCompleto == 0) {
                                chutemass = qlmistStd * malha[aux].cell(0).duto.area / area;
                                malha[aux].cell(0).acsr.injl.QLiq = chutemass;
                            } else {
                                chutemass = qgmistStd * malha[aux].cell(0).duto.area / area;
                                malha[aux].cell(0).acsr.injg.QGas = chutemass;
                            }
                            malha[aux].modoPerm = 1;
                            if (malha[aux].cell(malha[aux].ncel).acsr.tipo == 3)
                                malha[aux].config().condpocinj.CC = 1;
                            else
                                malha[aux].config().condpocinj.CC = 3;
                            valor = malha[aux].buscaInjPfundoPerm1();
                            cout << "convergencia tramo: " << aux << endl;
                            if (valor < -1e9 || valor > 1e9) {
                                inativo[aux] = 1;
                                inativoColetor(aux, arqRede, inativo, narq);
                                inativoAfluente(aux, arqRede, inativo);
                                if (valor < -1e9) {
                                    aviso(aux);
                                    (*arqRede.vg1dSP).iterRede = 200;
                                } else {
                                    aviso2(aux);
                                    (*arqRede.vg1dSP).iterRede = 200;
                                }
                                (*arqRede.vg1dSP).restartRede = 1;
                                return 9000.;
                            }
                            Resolv[aux] = 1;
                            ciclomalha++;
                            if (malha[aux].config().flashCompleto == 0)
                                qlmistStd -= malha[aux].cell(0).acsr.injl.QLiq;
                            else
                                qgmistStd -= malha[aux].cell(0).acsr.injg.QGas;
                        }
                    }
                }
            }
        }
    }
    return sqrt(norma) / nnos;
}

void RedeInj(SProd *malha, Rede &arqRede, int narq, string nomeArquivoLog, tipoValidacaoJson_t tipoValidacaoMrt,
             tipoSimulacao_t tipoSimulacaoMrt, varGlob1D *vG1d) {
    int testaFlashComp = 0;
    for (int i = 0; i < narq; i++) {
        SProd temporario(pathArqEntrada + arqRede.impfiles[i], nomeArquivoLog, tipoValidacaoMrt, tipoSimulacaoMrt, vG1d);

        if (arqRede.malha[i].ncoleta > 0)
            temporario.noextremo = 0; // caso o tramo não esteja no fim da rede, tomando como referencia o sentido
        // de escoamento assumido para a rede
        else
            temporario.noextremo = 1; // caso o tramo esteja no fim da rede, tomando como referencia o sentido
        // de escoamento assumido para a rede
        if (arqRede.malha[i].nafluente > 0)
            temporario.noinicial = 0; // caso o tramo não esteja no inicio da rede, tomando como referencia o sentido
        // de escoamento assumido para a rede
        else
            temporario.noinicial = 1; // caso o tramo esteja no inicio da rede, tomando como referencia o sentido
        // de escoamento assumido para a rede

        malha[i] = temporario;
        malha[i].indTramo = i;
        if (malha[i].cell(0).acsr.tipo == 2)
            (*arqRede.vg1dSP).somavaz += malha[i].cell(0).acsr.injl.QLiq / 86400.;
        (*arqRede.vg1dSP).somaarea += malha[i].cell(0).duto.area;
        testaFlashComp += malha[i].config().flashCompleto;
    }
    if (testaFlashComp != 0 && testaFlashComp != narq && testaFlashComp != 2 * narq)
        NumError("Rede com modelos de fluido diferentes");
    Vcr<int> inativo(narq, 0);
    int indativo = 0;
    while ((*arqRede.vg1dSP).restartRede == 1) {
        vector<int> coletorfinal;
        for (int i = 0; i < narq; i++)
            if (arqRede.malha[i].ncoleta == 0 && inativo[i] == 0)
                coletorfinal.push_back(i);
        sort(coletorfinal.begin(), coletorfinal.end());
        int ncoletfim = coletorfinal.size();
        Vcr<double> razcolet(narq - indativo);
        Vcr<double> prescolet(narq - indativo);
        for (int i = 0; i < ncoletfim; i++)
            chutePresRedeInj(coletorfinal[i], malha, arqRede, arqRede.chutHol, razcolet, prescolet);
        for (int i = 0; i < narq - indativo; i++) {
            int testaAflu = 1;
            for (int j = 0; j < ncoletfim; j++) {
                if (i == coletorfinal[j]) {
                    testaAflu = 0;
                }
            }
            if (testaAflu == 1) {
                malha[i].config().condpocinj.presfundo = prescolet[i] / razcolet[i];
                if (arqRede.malha[i].ncoleta > 0) {
                    int ncol = arqRede.malha[i].ncoleta;
                    for (int j = 0; j < ncol; j++) {
                        int qcol = arqRede.malha[i].coleta[j];
                        malha[qcol].pGSup = malha[i].config().condpocinj.presfundo;
                        malha[qcol].config().condpocinj.presinj = malha[i].config().condpocinj.presfundo;
                    }
                }
                if (malha[i].config().condpocinj.vazinj > 1e-5) {
                    malha[i].config().condpocinj.CC = 5;
                } else if (malha[i].config().condpocinj.presinj > 1e-5) {
                    malha[i].config().condpocinj.CC = 3;
                }
            }
        }
        double norma = 1000.;
        (*arqRede.vg1dSP).restartRede = 0;
        indativo = 0;
        while (norma > 0.001 * (*arqRede.vg1dSP).relax && indativo < narq && (*arqRede.vg1dSP).restartRede == 0) {
            indativo = 0;
            for (int i = 0; i < narq; i++)
                indativo += inativo[i];
            norma = cicloRedeInj(malha, arqRede, inativo, indativo);
            cout << "iteracao= " << (*arqRede.vg1dSP).iterRede << " norma= " << norma << " tramos inativos= " << indativo << endl;
            for (int i = 0; i < narq - indativo; i++) {
                for (int j = 0; j < ncoletfim; j++) {
                    if (i == coletorfinal[j]) {
                    }
                }
            }
            (*arqRede.vg1dSP).iterRede++;
        }
    }
    cout << "Convergencia atingida" << endl;
    for (int i = 0; i < narq; i++)
        if (inativo[i] == 0) {
            malha[i].config().imprimeProfile(malha[i].celula, malha[i].flut, 0, malha[i].indTramo);
            // enterramento
            for (int j = 0; j <= malha[i].ncel; j++) {
                if (malha[i].cell(j).calor.difus2D == 1) {
                    malha[i].cell(j).calor.poisson2D.imprimePermanente(malha[i].indTramo);
                }
            }
        }
    if (arqRede.nsisprod > 0)
        delete[] malha;
}

/*
 * showUsage.
 */
static void showUsage(string name) {
    cerr << "Usage: " << name << " <opcoes(s)> " << endl;
    cerr << "Marlim transiente versao " << versao << endl;
    cerr << "Opcoes:\n"
         << "\t-h,--help\t\tExibir esta mensagem de execucao\n"
         << "\t-i,--input ARQUIVO\tRedefinir o arquivo de entrada\n"
         << "\t-o,--output ARQUIVO\tRedefinir o arquivo LOG de resultado\n"
         << "\t-p,--path CAMINHO\tRedefinir o caminho de outros arquivos de entrada (PVTSIM)\n"
         << "\t-d,--dir PASTA\t\tRedefinir a pasta de saida de resultados\n"
         << "\t-v,--validate <JSON>  \tValidar a estrutura/formato JSON da entrada\n"
         << "\t              <SCHEMA>\tValidar o item anterior e a sintaxe de preenchimento da entrada\n"
         << "\t              <RULES> \tValidar os itens anteriores e as regras de negocio dos parametros da entrada\n"
         << "\t              <OFF>   \tDesligar a validacao da entrada e nao interromper a execucao\n"
         << "\t-s,--simulate <TRANSIENTE>\tSimular um sistema transiente\n"
         << "\t              <INJETOR>   \tSimular um sistema de poco injetor\n"
         << "\t              <REDE>      \tSimular um sistema de rede\n"
         << endl;
}

/*
 * Funcao:
 * Especificar o path de todos os arquivos de entrada e saida a partir dos paths dos arquivos
 * de entrada e de log do simulador.
 * Parametros:
 * pathComArquivoEntrada - path completo de arquivo de entrada da simulacao
 * pathComArquivoLog - path completo de arquivo de log da simulacao
 * tipoSimulacao - tipo da simulacao: rede, transiente, poco injetor
 */
void determinarPathArqEntSai(string pathComNomeArqEntrada, string pathComNomeArqLog,
                             tipoSimulacao_t tipoSimulacao, string diretorioSaida) {
    // posicao do ponto que separa o nome e a extensao do arquivo de entrada
    int posPontoExtensaoArqEntrada = pathComNomeArqEntrada.size() - 1;
    while (posPontoExtensaoArqEntrada > 0 && pathComNomeArqEntrada.at(posPontoExtensaoArqEntrada) != '.') {
        posPontoExtensaoArqEntrada--;
    }
    // posicao da ultima barra do path do arquivo de entrada, caso exista
    int posBarraSeparacaoPathArqEntrada = posPontoExtensaoArqEntrada;
    while (posBarraSeparacaoPathArqEntrada >= 0 && pathComNomeArqEntrada.at(posBarraSeparacaoPathArqEntrada) != '/' && pathComNomeArqEntrada.at(posBarraSeparacaoPathArqEntrada) != '\\') {
        posBarraSeparacaoPathArqEntrada--;
    }
    // path de arquivo de entrada incluindo a barra
    pathArqEntrada = pathComNomeArqEntrada.substr(0, posBarraSeparacaoPathArqEntrada + 1);
    // caso o path ainda nao tenha sido definido pelo parametro de execucao do simulador
    if (pathArqExtEntrada.length() == 0) {
        // define o path do arquivo de entrada com a barra
        pathArqExtEntrada = pathArqEntrada;
    }
    // nome do arquivo de entrada sem a extensao
    string nomeArqEntrada = pathComNomeArqEntrada.substr(posBarraSeparacaoPathArqEntrada + 1, posPontoExtensaoArqEntrada - posBarraSeparacaoPathArqEntrada - 1);
    // caso nao seja simulacao POCO_INJETOR
    if (tipoSimulacao == tipoSimulacao_t::poco_injetor) {
        // posicao da ultima barra do path do arquivo de entrada, caso exista
        int posBarraSeparacaoPathArqLog = pathComNomeArqLog.size() - 1;
        while (posBarraSeparacaoPathArqLog >= 0 && pathComNomeArqLog.at(posBarraSeparacaoPathArqLog) != '/' && pathComNomeArqLog.at(posBarraSeparacaoPathArqLog) != '\\') {
            posBarraSeparacaoPathArqLog--;
        }
        if (diretorioSaida == "")
            // define o path dos arquivos de saida como o mesmo do arquivo de log, com a barra
            pathPrefixoArqSaida = pathComNomeArqLog.substr(0, posBarraSeparacaoPathArqLog + 1);
        else
            pathPrefixoArqSaida = diretorioSaida;
        if (pathPrefixoArqSaida[pathPrefixoArqSaida.size() - 1] != '/') {
            pathPrefixoArqSaida = pathPrefixoArqSaida + "/";
        }
        // define o path como o mesmo do arquivo de log
        arqSaidaSnapShot = nomeArqEntrada + ".snp";

    } else {
        // posicao da ultima barra do path do arquivo de entrada, caso exista
        int posBarraSeparacaoPathArqLog = pathComNomeArqLog.size() - 1;
        while (posBarraSeparacaoPathArqLog >= 0 && pathComNomeArqLog.at(posBarraSeparacaoPathArqLog) != '/' && pathComNomeArqLog.at(posBarraSeparacaoPathArqLog) != '\\') {
            posBarraSeparacaoPathArqLog--;
        }
        if (diretorioSaida == "")
            // define o path dos arquivos de saida como o mesmo do arquivo de log, com a barra
            pathPrefixoArqSaida = pathComNomeArqLog.substr(0, posBarraSeparacaoPathArqLog + 1);
        else
            pathPrefixoArqSaida = diretorioSaida;
        if (pathPrefixoArqSaida[pathPrefixoArqSaida.size() - 1] != '/') {
            pathPrefixoArqSaida = pathPrefixoArqSaida + "/";
        }
        // define o path como o mesmo do arquivo de log
        arqSaidaSnapShot = nomeArqEntrada + ".snp";
    }
}

/*
 * Funcao:
 * Especificar o path de todos os arquivos de entrada do PVTSim e do SnapShot a partir do parametro de entrada.
 * Parametros:
 * path - path completo de entrada dos arquivos PVTSim e SnapShot para a simulacao
 */
void determinarPathArqExtEntrada(string path) {
    // caso o path possua a barra no ultimo caracter
    if (path.at(path.length() - 1) == '/' || path.at(path.length() - 1) == '\\') {
        // define o path como o original
        pathArqExtEntrada = path;
    } else {
        // define o path com a barra
        pathArqExtEntrada = path + BARRA;
    }
}

/*
 * Funcao:
 * Especificar o path de todos os arquivos de saida.
 * Parametros:
 * path - path completo de saida dos arquivos
 */
void determinarPathSaida(string path) {
    diretorioSaida = path;
}

/*
 * Funcao:
 * Tratar sinais para encerramento do programa gerados pelo usuário ou oriundos de operações de erro, inválidas, S.O.
 * Parametros:
 * sigNumber - numero que indica o tipo de sinal de interrupcao
 */
void signalHandler(int sigNumber) {
    string signalDescription = "";
    int exitCode = EXIT_FAILURE;
    switch (sigNumber) {
    case SIGABRT:
        signalDescription = "SIGABRT - Encerramento anormal da simulacao, por uma chamada 'abort, provavelmente devido a sintaxe do arquivo json', por exemplo";
        break;
    case SIGFPE:
        signalDescription = "SIGFPE - Operacao aritmetica erronea, 'divisao por zero' ou 'overflow em operacao', por exemplo";
        break;
    case SIGILL:
        signalDescription = "SIGILL - Instrucao ilegal detectada";
        break;
    case SIGINT:
        signalDescription = "SIGINT - Encerramento solicitado pelo usuario";
        exitCode = EXIT_SUCCESS;
        break;
    case SIGSEGV:
        signalDescription = "SIGSEGV - Acesso invalido a memoria";
        break;
    case SIGTERM:
        signalDescription = "SIGTERM - Requisicao de termino enviada ao programa";
        exitCode = EXIT_SUCCESS;
        break;
    }

    if (tipoSimulacao == tipoSimulacao_t::transiente) {
        WriteSnapShot(*ptrSistemaProducao);
    }
    // incluir log INFO de cancelamento da simulacao
    logger.log(LOGGER_FALHA, LOG_INFO_SIMULATION_CANCELLED, to_string(sigNumber), "", signalDescription);
    // gravar arquivo de log
    logger.writeOutputLog();
    // caso nao seja simulacao POCO_INJETOR
    if (tipoSimulacao != tipoSimulacao_t::poco_injetor && tipoSimulacao != tipoSimulacao_t::convecNatural) {
        // fechar arquivo de registro dos arquivos de saida de dados da simulacao
        arqRelatorioPerfis.close();
    }
    exit(exitCode);
}

int main(int argc, char **argv) {

    nowGlobIni = time(0);
    ltmGlobIni = localtime(&nowGlobIni);

    diaIni = ltmGlobIni->tm_mday;
    horaIni = ltmGlobIni->tm_hour;
    minutoIni = ltmGlobIni->tm_min;
    segundoIni = ltmGlobIni->tm_sec;

    // configurar nome do arquivo padrao do simulador simples em formato json do simulador
    string nomeArquivoEntrada = "teste1.json"; // nome padrao para arquivo de entrada json
    // configurar nome do arquivo de log em formato json do simulador
    string nomeArquivoLog = "simulacao.log";
    // configurar para processar o arquivo MRT e executar o simulador
    tipoValidacaoJson_t validacaoJson = tipoValidacaoJson_t::none;

    string nomeArquivoAS = "AS.json"; // nome padrao para arquivo de entrada de analise de sensibilidade

    int saidaClassica = 0;

    /*
     * Processamento dos argumentos de entrada
     */
    srand(time(NULL));
    int frase = rand() % 15;

    // percorrer a lista de parametros de entrada
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        int ultimo = arg.size() - 1;
        if ((arg.at(0) == '\'') && (arg.at(ultimo) == '\'')) {
            arg = arg.substr(1, ultimo - 1);
        }
        if ((arg == "-h") || (arg == "--help")) {
            showUsage(argv[0]);
            return EXIT_SUCCESS;
        } else if ((arg == "-v") || (arg == "--validate")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar a opcao de validacao a partir do parametro de entrada
                string opcao = argv[++i];
                // confirmar qual opcao foi escolhida
                if (opcao.compare("JSON") == 0) {
                    validacaoJson = tipoValidacaoJson_t::json;
                } else if (opcao.compare("SCHEMA") == 0) {
                    validacaoJson = tipoValidacaoJson_t::schema;
                } else if (opcao.compare("RULES") == 0) {
                    validacaoJson = tipoValidacaoJson_t::rules;
                } else if (opcao.compare("OFF") == 0) {
                    validacaoJson = tipoValidacaoJson_t::off;
                } else {
                    showUsage(argv[0]);
                    cerr << "opcao --validate requer uma das subopcoes <JSON>, <SCHEMA>, <RULES> ou <OFF>." << endl;
                    return EXIT_FAILURE;
                }
                // Forcando com que NUNCA haja validacao com o schema
                validacaoJson = tipoValidacaoJson_t::off;
            } else {
                showUsage(argv[0]);
                cerr << "opcao --validate requer uma subopcao." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-s") || (arg == "--simulate")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar a opcao de simulacao a partir do parametro de entrada
                string opcao = argv[++i];
                // confirmar qual opcao foi escolhida
                if (opcao.compare("TRANSIENTE") == 0) {
                    tipoSimulacao = tipoSimulacao_t::transiente;
                } else if (opcao.compare("INJETOR") == 0) {
                    tipoSimulacao = tipoSimulacao_t::poco_injetor;
                } else if (opcao.compare("CONVECNAT") == 0) {
                    tipoSimulacao = tipoSimulacao_t::convecNatural;
                    saidaClassica = 1;
                    nomeArquivoEntrada = "parametros.json";
                    if ((arg == "-i") || (arg == "--input")) {
                        // caso exista outro parametro apos esta opcao
                        if (i + 1 < argc) {
                            // recuperar o nome do arquivo a partir do parametro de entrada
                            nomeArquivoEntrada = argv[++i];
                        } else {
                            showUsage(argv[0]);
                            cerr << "opcao --input requer um arquivo." << endl;
                            return EXIT_FAILURE;
                        }
                    }
                } else if (opcao.compare("REDE") == 0) {
                    tipoSimulacao = tipoSimulacao_t::rede;
                    // caso nome do arquivo padrao do simulador ainda inalterado
                    if (nomeArquivoEntrada.compare("teste1.json") == 0) {
                        // configurar nome do arquivo padrao do simulador rede
                        nomeArquivoEntrada = "rede.json";
                    }
                } else {
                    showUsage(argv[0]);
                    cerr << "opcao --simulate requer uma das subopcoes <TRANSIENTE>, <INJETOR> ou <REDE>." << endl;
                    return EXIT_FAILURE;
                }
            } else {
                showUsage(argv[0]);
                cerr << "opcao --simulate requer uma subopcao." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-i") || (arg == "--input")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar o nome do arquivo a partir do parametro de entrada
                nomeArquivoEntrada = argv[++i];
            } else {
                showUsage(argv[0]);
                cerr << "opcao --input requer um arquivo." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-o") || (arg == "--output")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar o nome do arquivo a partir do parametro de entrada
                nomeArquivoLog = argv[++i];
            } else {
                showUsage(argv[0]);
                cerr << "opcao --output requer um arquivo." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-p") || (arg == "--path")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar o caminho de outros arquivos de entrada a partir do parametro
                determinarPathArqExtEntrada(argv[++i]);
            } else {
                showUsage(argv[0]);
                cerr << "opcao --path requer um caminho." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-d") || (arg == "--dir")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                // recuperar o caminho de outros arquivos de entrada a partir do parametro
                determinarPathSaida(argv[++i]);
            } else {
                showUsage(argv[0]);
                cerr << "opcao --dir requer um caminho." << endl;
                return EXIT_FAILURE;
            }
        } else if ((arg == "-t") || (arg == "--threads")) {
            // caso exista outro parametro apos esta opcao
            if (i + 1 < argc) {
                string nthr = argv[++i];
            } else {
                showUsage(argv[0]);
                cerr << "opcao --threads requer um valor." << endl;
                return EXIT_FAILURE;
            }
        } else {
            showUsage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    // register signal SIGABRT and signal handler (Abnormal termination of the program, such as a call to abort)
    signal(SIGABRT, signalHandler);
    // register signal SIGFPE and signal handler (An erroneous arithmetic operation, such as a divide by zero or an operation resulting in overflow)
    signal(SIGFPE, signalHandler);
    // register signal SIGILL and signal handler (Detection of an illegal instruction)
    signal(SIGILL, signalHandler);
    // register signal SIGINT and signal handler (Receipt of an interactive attention signal)
    signal(SIGINT, signalHandler);
    // register signal SIGSEGV and signal handler (An invalid access to storage)
    signal(SIGSEGV, signalHandler);
    // register signal SIGTERM and signal handler (A termination request sent to the program)
    signal(SIGTERM, signalHandler);

    // Apagar arquivo de "debug do modulo composicional FORTRAN" caso exista:

        const std::string debugFileName = "Composicional_Debug.txt";

        std::remove(debugFileName.c_str());

    /*
     * Processamento da simulacao
     */

    try {
        logger.changeResultadoSimulacao(nomeArquivoEntrada, nomeArquivoLog);

        // obter path dos arquivos de saida e o nome do arquivo de SnapShot
        determinarPathArqEntSai(nomeArquivoEntrada, nomeArquivoLog, tipoSimulacao, diretorioSaida);

        // caso NAO simulacao POCO_INJETOR
        if (tipoSimulacao != tipoSimulacao_t::poco_injetor && tipoSimulacao != tipoSimulacao_t::convecNatural) {
            // arquivo de registro dos arquivos de saida de dados da simulacao
            arqRelatorioPerfis.open(pathPrefixoArqSaida + ARQUIVO_RELATORIO_PERFIS, ios_base::out);
        }

        // caso simulacao de sistema de rede
        if (tipoSimulacao == tipoSimulacao_t::rede) { // simulacao do tipo rede
            logRede = 1;
            varGlob1D vg1dRedeSimples = varGlob1D();
            vg1dRedeSimples.chaverede = 1;
            nomeRedePrincipal = nomeArquivoEntrada;
            Rede arqRede(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, &vg1dRedeSimples); // leitura
            // do arquivo rede
            int narq = arqRede.nsisprod; // numero de tramos
            vg1dRedeSimples.narq = narq;
            vg1dRedeSimples.chaveredeT = arqRede.chaveredeT; // indica se vai haver simulacao de rede transiente
            vg1dRedeSimples.chaveAnelGL = arqRede.anelGL;    // indica se a rede e do tipo anel de gas-lift
            vg1dRedeSimples.chaveRedeParalela = arqRede.redeParalela;
            vg1dRedeSimples.fluidoRede = arqRede.fluidoRede; // para o caso de rede, indicacao se a rede vai ser de com um fluido dominado por gas ou liquido
            vg1dRedeSimples.ntrd = arqRede.nthrRede;
            if (vg1dRedeSimples.chaveredeT > 0)
                vg1dRedeSimples.TmaxR = arqRede.TmaxR; // tempo maximo de simulacao transiente para rede
            else
                vg1dRedeSimples.TmaxR = 1;

            // malha=new SProd [narq];//vetor com os objetos tramos
            // Vcr<int> inativo(narq,0);//indica se algum tramo falhou na sua busca pela soluçao transiente
            // falha de convergencia
            if ((*arqRede.vg1dSP).chaveAnelGL == 0 && (*arqRede.vg1dSP).chaveRedeParalela == 0) { // nao é anel de GL
                // caso nao seja simulacao de POCO INJETOR
                // SProd* malha;//alteracao7
                if (arqRede.injec == 0) {
                    // metodo que faz o preprocessamento da rede, avalia quais tramos da rede se encontram conectados, reorganizando a rede
                    // em subredes, chamadas de "redes internas", neste método, alem de se verificar as conexoes de redes internas,
                    // e construido arquivos json de entrada de cada rede interna

                    descarteTramo(arqRede, narq);
                    preProcRede(arqRede, narq);
                    ostringstream relatSucesso;
                    relatSucesso << pathPrefixoArqSaida << "relatorioSucessoRede" << ".dat";
                    string tmp = relatSucesso.str();
                    ofstream escreveRelatorioSucesso(tmp.c_str(), ios_base::out);
                    escreveRelatorioSucesso << "# RedeInterna      ;" << "# Tramo      ;" << "Permanente     ;" << "Ativo       " << endl;
                    escreveRelatorioSucesso.close();
                    if (arqRede.apenasPreProc == 0) {     // sucesso no pre processamento
                        SProd **malha;                    // ponteiro para armazenamentodas redes internas que compoe a rede original
                        malha = new SProd *[redeLeitura]; // alocamento de memoria para o vetor de redes internas
                        Rede *arqRedeTemp;
                        arqRedeTemp = new Rede[redeLeitura];
                        varGlob1D *vg1dRede;
                        vg1dRede = new varGlob1D[redeLeitura];
                        int *contapermRede;
                        contapermRede = new int[redeLeitura];
                        varGlob1D vg1dTemp = varGlob1D();
                        for (int i = 0; i < redeLeitura; i++)
                            vg1dRede[i] = vg1dTemp;
                        vg1dRedeSimples.ntrdGlob = 1;
                        if (redeLeitura > 1) {
                            vg1dRedeSimples.ntrdGlob = arqRede.nthrRede;
                            arqRede.nthrRede = 1;
                        }
                        for (int i = 0; i < redeLeitura; i++) {
                            vg1dRede[i].chaverede = 1;
                            vg1dRede[i].chaveredeT = arqRede.chaveredeT; // indica se vai haver simulacao de rede transiente
                            vg1dRede[i].chaveAnelGL = arqRede.anelGL;    // indica se a rede e do tipo anel de gas-lift
                            vg1dRede[i].fluidoRede = arqRede.fluidoRede; // para o caso de rede, indicacao se a rede vai ser de com um fluido dominado por gas ou liquido
                            vg1dRede[i].ntrd = arqRede.nthrRede;
                            vg1dRede[i].ntrdGlob = vg1dRedeSimples.ntrdGlob;
                            if (vg1dRede[i].chaveredeT > 0)
                                vg1dRede[i].TmaxR = arqRede.TmaxR; // tempo maximo de simulacao transiente para rede
                            else
                                vg1dRede[i].TmaxR = 1;
                            vg1dRede[i].qualRede = i;

                            ostringstream entrada;
                            entrada << pathPrefixoArqSaida << "RedeInterna-" << i << ".json";
                            string tmp = entrada.str();
                            // leitura do arquivo json da rede interna i
                            Rede arqRedeTemp2(tmp, nomeArquivoLog, validacaoJson, tipoSimulacao, &vg1dRede[i]); // leitura
                            arqRedeTemp[i] = arqRedeTemp2;
                            vg1dRede[i].narq = arqRedeTemp[i].nsisprod; // numero de tramos da rede interna i
                            vg1dRede[i].narq = arqRedeTemp[i].nsisprod;
                            vg1dRede[i].fluidoRede = arqRedeTemp[i].fluidoRede;
                        }
                        for (int i = 0; i < redeLeitura; i++) {
                            contapermRede[i] = 0;
                            malha[i] = new SProd[vg1dRede[i].narq];
                            preparaRedeProd(malha[i], arqRedeTemp[i], vg1dRede[i].narq, nomeArquivoLog, validacaoJson,
                                            tipoSimulacao_t::transiente, i, contapermRede[i], &vg1dRede[i]);
                        }
#pragma omp parallel for num_threads(vg1dRedeSimples.ntrdGlob)
                        for (int i = 0; i < redeLeitura; i++) {
                            Vcr<int> inativo(vg1dRede[i].narq, 0); // indica se algum tramo falhou na sua busca pela solu�ao transiente
                            int indativo = 0;                      // indica o numero de tramso que foram inativados na solucao permanente por
                            // falha de convergencia
                            cout << "!!!!!! Resolvendo REDE INTERNA " << i << "!!!!!!" << "\n";
                            vector<noRede> normaEvol;
                            vector<tramoPart> bloq;
                            solveRedeProd(malha[i], arqRedeTemp[i], vg1dRede[i].narq, inativo, indativo, nomeArquivoLog, validacaoJson,
                                          tipoSimulacao_t::transiente, i, contapermRede[i], &vg1dRede[i], normaEvol, bloq); // construcao dos objetos tramos de uma rede
                            // neste metodo e criado o vetor de tramos da rede, avaliado se algum tramo encontra-se em uma condicao sem vazao,
                            // neste caso o tramo e retirado da rede interna, e e feita a resolucao permenente da rede interna, apos esta resolucao
                            // os perfis permanentes sao gravados.
                            // de producao e busca de sua solucao permanente
                            saidaClassica = malha[i][0].config().saidaClassica;
                            // solucao transiente de uma rede de producao classica, sem anel.
                            if (vg1dRede[i].chaveredeT == 1 && arqRede.injec == 0 && (*arqRede.vg1dSP).chaveAnelGL == 0) {
                                vg1dRede[i].RGOMax = 14000.;
                                vg1dRede[i].modoTransiente = 1;
                                SolveRedeTrans(malha[i], arqRedeTemp[i], inativo, indativo, i); // metodo em que se faz a resolucao transiente da rede interna
                                vg1dRede[i].modoTransiente = 0;
                            }
                            vg1dRede[i].restartRede = 1;
                            // reiniciando os parametros globais de solucao de rede para a resolucao da proxima rede interna
                            vg1dRede[i].restartRede = 1;
                            vg1dRede[i].somaarea = 0;
                            vg1dRede[i].somavaz = 0;
                            vg1dRede[i].somavazG = 0;
                            vg1dRede[i].relax = 0.5;
                            vg1dRede[i].iterRedeT = 0;
                            vg1dRede[i].erroRede = 1000;
                            vg1dRede[i].arearef = 0.25 * M_PI * (7 * 7 - 5 * 5) * 2.54 * 2.54 / 10000;

                            normaEvol.clear();
                            bloq.clear();
                        }

                        if (redeLeitura > 0) {
                            if (saidaClassica == 1) {
                                cout << "*******************************************************************************" << endl;
                                cout << "                                  UFA!!!!!!!!                                  " << endl;
                                cout << saidaTexto[frase] << endl;
                                cout << saidaSubTexto[frase] << endl;
                                cout << "*******************************************************************************" << endl;
                            } else
                                cout << "                                  FIM                                  " << endl;
                            // deletando as redes internas
                            for (int kRede = 0; kRede < redeLeitura; kRede++)
                                delete[] malha[kRede];
                            delete[] malha; // deletando os objetos tramos
                            delete[] vg1dRede;
                            delete[] arqRedeTemp;
                            delete[] contapermRede;
                        }

                        ostringstream relatSucesso;
                        relatSucesso << pathPrefixoArqSaida << "relatorioSucessoRede" << ".dat";
                        string tmp = relatSucesso.str();
                        ofstream escreveRelatorioSucesso(tmp.c_str(), ios_base::app);
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        int diaFim = (ltm->tm_mday);
                        int horaFim;
                        if (diaFim == diaIni)
                            horaFim = ltm->tm_hour;
                        else
                            horaFim = ltm->tm_hour + 24;
                        horaFim *= 3600;
                        int minutoFim = 60 * ltm->tm_min;
                        int segundoFim = ltm->tm_sec;
                        int totalFim = horaFim + minutoFim + segundoFim;
                        int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
                        escreveRelatorioSucesso << endl;
                        escreveRelatorioSucesso << endl;
                        escreveRelatorioSucesso << endl;
                        escreveRelatorioSucesso << "datahora = ";
                        escreveRelatorioSucesso << ltm->tm_mday << "/";
                        escreveRelatorioSucesso << 1 + ltm->tm_mon << "/";
                        escreveRelatorioSucesso << 1900 + ltm->tm_year << " ";
                        escreveRelatorioSucesso << 0 + ltm->tm_hour << ":";
                        escreveRelatorioSucesso << 0 + ltm->tm_min << ":";
                        escreveRelatorioSucesso << 0 + ltm->tm_sec;
                        escreveRelatorioSucesso << endl;
                        escreveRelatorioSucesso << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
                        escreveRelatorioSucesso << "     Versao    " << versao << endl;

                        escreveRelatorioSucesso.close();
                    }

                } else { // caso rede de injecao, resolucao permanente da rede de injecao
                    // OBS: sistemas de injecao so tem solucao permanente
                    SProd *malha;
                    vg1dRedeSimples.narq = narq;
                    malha = new SProd[narq];                                                                                       // vetor com os objetos tramos
                    RedeInj(malha, arqRede, narq, nomeArquivoLog, validacaoJson, tipoSimulacao_t::poco_injetor, &vg1dRedeSimples); // construcao
                    // dos objetos tramos de uma rede de injecao e busca de sua solucao permanente
                    if (saidaClassica == 1) {
                        cout << "*******************************************************************************" << endl;
                        cout << "                                  UFA!!!!!!!!                                  " << endl;
                        cout << saidaTexto[frase] << endl;
                        cout << saidaSubTexto[frase] << endl;
                        cout << "*******************************************************************************" << endl;
                    } else
                        cout << "                                  FIM                                  " << endl;
                    if (arqRede.nsisprod > 0)
                        delete[] malha; // deletando os objetos tramos
                }
            } else {
                if ((*arqRede.vg1dSP).chaveAnelGL == 1) {
                    // solucao de rede de anel de Gas Lift
                    // SProd* malha;//alteracao7
                    SProd *malha;
                    vg1dRedeSimples.narq = narq;
                    malha = new SProd[narq];   // vetor com os objetos tramos
                    Vcr<int> inativo(narq, 0); // indica se algum tramo falhou na sua busca pela soluçao transiente
                    int indativo = 0;          // indica o numero de tramso que foram inativados na solucao permanente por
                    // falha de convergencia
                    RedeAnelGL(malha, arqRede, narq, inativo, indativo, nomeArquivoLog,
                               validacaoJson, arqRede.chaveredeT, tipoSimulacao_t::transiente, &vg1dRedeSimples); // cosntrucao dos objetos tramos
                    // de uma rede de anel de GL, dentro deste metodo e feita a solucao permanente e
                    // eventualmente a solucao transiente tambem

                    saidaClassica = malha[0].config().saidaClassica;
                    if (saidaClassica == 1) {
                        cout << "*******************************************************************************" << endl;
                        cout << "                                  UFA!!!!!!!!                                  " << endl;
                        cout << saidaTexto[frase] << endl;
                        cout << saidaSubTexto[frase] << endl;
                        cout << "*******************************************************************************" << endl;
                    } else
                        cout << "                                  FIM                                  " << endl;
                    if (arqRede.nsisprod > 0)
                        delete[] malha; // deletando os objetos tramos
                } else {
                    SProd *malha;
                    vg1dRedeSimples.narq = narq;
                    malha = new SProd[narq];   // vetor com os objetos tramos
                    Vcr<int> inativo(narq, 0); // indica se algum tramo falhou na sua busca pela soluçao transiente
                    int indativo = 0;          // indica o numero de tramso que foram inativados na solucao permanente por
                    // falha de convergencia

                    RedeParalela(malha, arqRede, narq,
                                 inativo, indativo, nomeArquivoLog, validacaoJson, arqRede.chaveredeT,
                                 tipoSimulacao_t::transiente, &vg1dRedeSimples);

                    saidaClassica = malha[0].config().saidaClassica;
                    if (saidaClassica == 1) {
                        cout << "*******************************************************************************" << endl;
                        cout << "                                  UFA!!!!!!!!                                  " << endl;
                        cout << saidaTexto[frase] << endl;
                        cout << saidaSubTexto[frase] << endl;
                        cout << "*******************************************************************************" << endl;
                    } else
                        cout << "                                  FIM                                  " << endl;
                    if (arqRede.nsisprod > 0)
                        delete[] malha; // deletando os objetos tramos
                }
            }
            // solucao transiente de uma rede de producao classica, sem anel.

        } else if (tipoSimulacao != tipoSimulacao_t::convecNatural) { // solucao de tramos solitarios, sem rede
            varGlob1D vg1dTramo = varGlob1D();
            // criar objeto de simulacao
            // construtor do objeto que representa o tramo
            SProd sistem1(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, &vg1dTramo);
            ptrSistemaProducao = &sistem1;

            // atualiza objeto de logger com registros anteriores
            logger.setNomeArqEntrada(nomeArquivoEntrada);
            saidaClassica = sistem1.config().saidaClassica;

            // caso habilitado o calculo do permanente
            if (sistem1.config().perm == 1) {
                if (sistem1.config().AP == 0)
                    SolveTramoSolteiro(sistem1, sistem1.config().chutePerm); // solucao simples de permanente para um tramo simples
                else {
    			ostringstream arquivoAPtemp;
    			arquivoAPtemp<<pathArqEntrada<<sistem1.config().arquivoAP;
    			string arquivoAP = arquivoAPtemp.str();
    			if(sistem1.config().paralelAP==0)
    				leituraAP(arquivoAP, sistem1);//solucao de analise de sensibilidade para tramo simples
    			else
    				leituraAPparalelo(arquivoAP, nomeArquivoLog, validacaoJson, sistem1);//solucao de analise de sensibilidade para tramo simples
                }

            } else if (sistem1.config().perm == 2) { // quando a condicao de contorno nao e o resultado de um regime permanente mas a "foto" em algum momento de
                // uma simulacao transiente
                // caso definido o restart da simulacao
                ReadSnapShot(sistem1);
            }
            if (sistem1.config().transiente == 1) { // opcao transiente ativa, inicio do processo transiente
                // preparacao de algumas variaveis para a partida na simulacao transiente:
                if (sistem1.config().ConContEntrada == 1) {
                    if (sistem1.cell(0).acsr.tipo == 2) {
                        sistem1.cell(0).acsr.injl.QLiq = 0.;
                        sistem1.cell(0).MC = sistem1.cell(0).fontemassLR + sistem1.cell(0).fontemassCR + sistem1.cell(0).fontemassGR;
                        sistem1.cell(0).MCini = sistem1.cell(0).MC;
                        sistem1.cell(0).Mliqini = sistem1.cell(0).fontemassLR + sistem1.cell(0).fontemassCR;
                        sistem1.cell(0).Mliqini0 = sistem1.cell(0).Mliqini;
                        double rholMix = sistem1.cell(0).bet * sistem1.cell(0).fluicol.MasEspFlu(
                                                                     sistem1.cell(0).pres, sistem1.cell(0).temp);
                        rholMix += (1. - sistem1.cell(0).bet) * sistem1.cell(0).flui.MasEspLiq(
                                                                      sistem1.cell(0).pres, sistem1.cell(0).temp);
                        sistem1.cell(0).QL = sistem1.cell(0).Mliqini / rholMix;
                        sistem1.cell(0).QLini = sistem1.cell(0).QL;
                        double rhog = sistem1.cell(0).flui.MasEspGas(
                            sistem1.cell(0).pres, sistem1.cell(0).temp);
                        sistem1.cell(0).QG = (sistem1.cell(0).MC - sistem1.cell(0).Mliqini) / rhog;
                        sistem1.cell(0).QGini = sistem1.cell(0).QG;
                        sistem1.cell(0).fontemassLR = 0.;
                        sistem1.cell(0).fontemassCR = 0.;
                        sistem1.cell(0).fontemassGR = 0.;
                    } else if (sistem1.cell(0).acsr.tipo == 1) {
                        sistem1.cell(0).acsr.injg.QGas = 0.;
                        sistem1.cell(0).MC = sistem1.cell(0).fontemassGR;
                        sistem1.cell(0).MCini = sistem1.cell(0).MC;
                        sistem1.cell(0).Mliqini = 0.;
                        sistem1.cell(0).Mliqini0 = sistem1.cell(0).Mliqini;
                        sistem1.cell(0).QL = 0.;
                        sistem1.cell(0).QLini = sistem1.cell(0).QL;
                        double rhog = sistem1.cell(0).flui.MasEspGas(
                            sistem1.cell(0).pres, sistem1.cell(0).temp);
                        sistem1.cell(0).QG = (sistem1.cell(0).MC - sistem1.cell(0).Mliqini) / rhog;
                        sistem1.cell(0).QGini = sistem1.cell(0).QG;
                        sistem1.cell(0).fontemassLR = 0.;
                        sistem1.cell(0).fontemassCR = 0.;
                        sistem1.cell(0).fontemassGR = 0.;
                    }
                }
                for (int i = 0; i <= sistem1.ncel; i++) {
                    sistem1.cell(i).rpC = sistem1.cell(i).rpCi =
                        sistem1.cell(i).flui.MasEspLiq(sistem1.cell(i).pres, sistem1.cell(i).temp);
                    sistem1.cell(i).rgC = sistem1.cell(i).rgCi =
                        sistem1.cell(i).flui.MasEspGas(sistem1.cell(i).pres, sistem1.cell(i).temp);
                    sistem1.cell(i).rcC = sistem1.cell(i).rcCi =
                        sistem1.cell(i).fluicol.MasEspFlu(sistem1.cell(i).pres, sistem1.cell(i).temp);
                    sistem1.cell(i).rpL = sistem1.cell(i).rpLi =
                        sistem1.cell(i).flui.MasEspLiq(sistem1.cell(i).presL, sistem1.cell(i).tempL);
                    sistem1.cell(i).rgL = sistem1.cell(i).rgLi =
                        sistem1.cell(i).flui.MasEspGas(sistem1.cell(i).presL, sistem1.cell(i).tempL);
                    sistem1.cell(i).rcL = sistem1.cell(i).rcLi =
                        sistem1.cell(i).fluicol.MasEspFlu(sistem1.cell(i).presL, sistem1.cell(i).tempL);

                    sistem1.cell(i).mipC = sistem1.cell(i).flui.ViscOleo(sistem1.cell(i).pres, sistem1.cell(i).temp);
                    sistem1.cell(i).migC = sistem1.cell(i).flui.ViscGas(sistem1.cell(i).pres, sistem1.cell(i).temp);
                    sistem1.cell(i).micC = sistem1.cell(i).fluicol.VisFlu(sistem1.cell(i).pres, sistem1.cell(i).temp);

                    if (i > 0) {
                        sistem1.cell(i - 1).rpR = sistem1.cell(i - 1).rpRi = sistem1.cell(i).rpC;
                        sistem1.cell(i - 1).rgR = sistem1.cell(i - 1).rgRi = sistem1.cell(i).rgC;
                        sistem1.cell(i - 1).rcR = sistem1.cell(i - 1).rcRi = sistem1.cell(i).rcC;

                        sistem1.cell(i - 1).mipR = sistem1.cell(i).mipC;
                        sistem1.cell(i - 1).migR = sistem1.cell(i).migC;
                        sistem1.cell(i - 1).micR = sistem1.cell(i).micC;
                    }
                    if (i == sistem1.ncel) {
                        sistem1.cell(i).rpR = sistem1.cell(i).rpRi = sistem1.cell(i).rpC;
                        sistem1.cell(i).rgR = sistem1.cell(i).rgRi = sistem1.cell(i).rgC;
                        sistem1.cell(i).rcR = sistem1.cell(i).rcRi = sistem1.cell(i).rcC;

                        sistem1.cell(i).mipR = sistem1.cell(i).mipC;
                        sistem1.cell(i).migR = sistem1.cell(i).migC;
                        sistem1.cell(i).micR = sistem1.cell(i).micC;
                    }
                }
                // determinacao do incremento de tempo inicial da simulacao transiente
                sistem1.determinaDT();
                // prepara o modelo de transferência de massa e atualiza os termos fontes existentes em um tramo para o inicio da simulacao transiente
                sistem1.renovaTemp();
                sistem1.modoPerm = 0;
                vg1dTramo.modoTransiente = 1;
                (*sistem1.vg1dSP).modoTransiente = 1;
                if (sistem1.config().modoDifus3D == 1) {
                    int icel = sistem1.config().celAcop[0].indCel;
                    double hiv = sistem1.cell(icel).calor.hi;
                    double hev = sistem1.cell(icel).calor.he; // sistem1.cell(icel).calor.he;
                    double tiv = sistem1.cell(icel).temp;
                    sistem1.poisson3D = solverP3D(sistem1.config().modoDifus3DJson, sistem1.vg1dSP,
                                                  sistem1.config().nacop, sistem1.config().geoAcop, hiv, hev, tiv);
                    if (sistem1.config().nacop != sistem1.poisson3D.dados.CC.nAcop)
                        NumError(
                            "O numero de acoplamentos indicados no json principal é diferente do numero de acoplamentos indicados no parse do solver 3D");
                    else {
                        for (int iacop0 = 0; iacop0 < sistem1.config().nacop; iacop0++) {
                            int inexistente = 1;
                            for (int iacop1 = 0; iacop1 < sistem1.config().nacop; iacop1++) {
                                if (sistem1.config().celAcop[iacop0].rotulo == sistem1.poisson3D.dados.CC.rotuloAcop[iacop1]) {
                                    inexistente = 0;
                                    sistem1.acertaIndAcop.push_back(iacop1);
                                    break;
                                }
                            }
                            if (inexistente == 1)
                                NumError("Rotulo de acoplamento indicado mas não encontrado");
                        }
                    }

                    for (int iacop = 0; iacop < sistem1.config().nacop; iacop++) {
                        int icelAcop = sistem1.config().celAcop[iacop].indCel;
                        int iacop1 = sistem1.acertaIndAcop[iacop];
                        sistem1.poisson3D.dados.tInt[iacop1] = sistem1.cell(icelAcop).temp;
                        double hiCel = sistem1.cell(icelAcop).calor.hInt();
                        sistem1.poisson3D.dados.hI[iacop1] = hiCel;
                    }
                    double tamb[sistem1.poisson3D.dados.CC.nRic];
                    for (int idir = 0; idir < sistem1.poisson3D.dados.CC.nRic; idir++) {
                        tamb[idir] = sistem1.poisson3D.dados.CC.ccRic[idir].valAmb[0];
                        sistem1.poisson3D.dados.CC.ccRic[idir].valAmb[0] = sistem1.cell(icel).temp;
                    }
                    sistem1.poisson3D.inicializaTransientePoisson();
                    int kontaperm = 0;
                    double deltFic = 1.;
                    double qparede0 = 1e15;
                    double erroParede = 1e101;
                    (*sistem1.vg1dSP).tempo = 0.;
                    deltFic = sistem1.poisson3D.defineDeltPoisson();
                    while (erroParede > 1e-1) {
                        sistem1.poisson3D.transientePoissonDummy(deltFic, kontaperm);
                        deltFic *= 2.0;
                        if (deltFic > 100.)
                            deltFic = 100.;
                        kontaperm++;
                        erroParede = fabs(sistem1.poisson3D.dados.qTotal[0] - qparede0);
                        cout << "erroParede= " << (sistem1.poisson3D.dados.qTotal[0] - qparede0);
                        cout << endl;
                        qparede0 = sistem1.poisson3D.dados.qTotal[0];
                    }
                    kontaperm = 0;
                    while ((*sistem1.vg1dSP).tempo < 2000.) {
                        deltFic = sistem1.poisson3D.defineDeltPoisson();
                        if ((*sistem1.vg1dSP).tempo > 100. && (*sistem1.vg1dSP).tempo < 1000.) {
                            for (int iric = 0; iric < sistem1.poisson3D.dados.CC.nRic; iric++) {
                                sistem1.poisson3D.dados.CC.ccRic[iric].valAmb[0] =
                                    sistem1.cell(icel).temp + (tamb[iric] - sistem1.cell(icel).temp) *
                                                                    ((*sistem1.vg1dSP).tempo - 100.) / 900.;
                            }
                        }
                        sistem1.poisson3D.transientePoissonDummy(deltFic, kontaperm);
                        (*sistem1.vg1dSP).tempo += deltFic;
                        kontaperm++;
                        cout << "tempo= " << (*sistem1.vg1dSP).tempo;
                        cout << endl;
                    }
                    sistem1.poisson3D.malha.imprime(0);
                }
                while ((*sistem1.vg1dSP).lixo5 < sistem1.tfinal) {
                    vg1dTramo.RGOMax = 14000.;
                    (*sistem1.vg1dSP).RGOMax = 14000.; // RGO maximo admitido para um eventual modelo black oil, isto e necessario pois como o modelo
                    // foi originalmente baseado em um sistema black-oil, teria problemas nas equacoes de transferencia de massa caso
                    // nao tenha um valor, mesmo que pequeno de oleo como referencia para o oleo morto
                    sistem1.SolveTrans(); // metodo onde se faz o avanco de tempo transiente
                    for (int kontasnp = 0; kontasnp < sistem1.config().nsnp; kontasnp++) {
                        if (fabs((*sistem1.vg1dSP).lixo5 - sistem1.config().tempsnp[kontasnp]) < sistem1.dt) {
                            WriteSnapShot(sistem1, (*sistem1.vg1dSP).lixo5); // registro do arquivo SNP
                        }
                    }
                }
                WriteSnapShot(sistem1, (*sistem1.vg1dSP).lixo5); // registro do SNP do tempo final
            }
            if (saidaClassica == 1) {
                cout << "*******************************************************************************" << endl;
                cout << "                                  UFA!!!!!!!!                                  " << endl;
                cout << saidaTexto[frase] << endl;
                cout << saidaSubTexto[frase] << endl;
                cout << "*******************************************************************************" << endl;
            } else if (sistem1.config().saidaTela == 1)
                cout << "                                  FIM                                  " << endl;

            nowGlobFim = time(0);
            ltmGlobFim = localtime(&nowGlobFim);
        } else {
            varGlob1D vg1dConvec2D = varGlob1D();
            // solucao de caso 2D de conveccao natural em espaco confinado
            // construtor do solver de volumes finitos para uma malha e condicao de contorno especificada em um json
            solv2D resolucao(nomeArquivoEntrada, nomeArquivoLog, &vg1dConvec2D);
            // solver do problema de conveccao natural em espaco confinado 2D
            resolucao.resolve();
        }
    } catch (exception &excInt) {
        cout << "EXCECAO INESPERADA: " << excInt.what() << endl;
        // incluir falha
        logger.log(LOGGER_FALHA, LOG_ERR_UNEXPECTED_EXCEPTION, "", "", excInt.what());
        // gravar arquivo de log
        logger.writeOutputLog();
        // encerrar a aplicacao
        exit(EXIT_SUCCESS);
    }

    // incluir log INFO de conclusao da simulacao
    logger.log(LOGGER_INFO, LOG_INFO_SIMULATION_FINISHED, "", "", "SUCESSO");
    // gravar arquivo de log
    logger.writeOutputLog(saidaClassica);

    // caso nao seja simulacao POCO_INJETOR
    if (tipoSimulacao != tipoSimulacao_t::poco_injetor) {
        // fechar arquivo de registro dos arquivos de saida de dados da simulacao
        arqRelatorioPerfis.close();
    }

    return EXIT_SUCCESS;
}
