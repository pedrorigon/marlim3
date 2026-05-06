// ---

#ifndef _FA_HIDRATO
#define _FA_HIDRATO

#include <vector>
#include <string>
#include <tuple>
#include "SisProd.h"


class FA_Hidrato {
public:
    FA_Hidrato(const SProd& sistema);
    ~FA_Hidrato();

    //std::vector<double> V_w;
    //double V_w;

    void solverHidrato();

    bool checkHidrato(double P, double T);

    double TurnerHidrato(double P, double T, double P_eq, double T_eq,
                              double A_s, double V_h, double V_w,
                              double r_d, double r_p,
                              const std::string& estruturaHidratos,
                              double A_s_input,
                              double &A_s_eff_out);

    double Euler1ordemGas(int i, double j_G, double taxaCinetica,
                               double A, double MM_g);

    double Euler1ordemAgua(int i, double j_W, double taxaCinetica,
                                double A, double eta, double MM_w);

    double Euler1ordemHidrato(double j_H, double taxaCinetica,
                                   double A, double eta, double MM_h, double rho_h);

private:
    const SProd& sistema;

    std::vector<double> temperaturaCurva;
    std::vector<double> pressaoCurva;
    std::vector<double> temperaturaCurvaDeslocada;
    std::vector<double> pressaoCurvaDeslocada;

    double K_Hamm_Etanol, K_Hamm_MEG;
    double MM_H, MM_G, MM_W, W_Hamm;
    double M_Etanol, M_MEG; //não é mais empregada, nem como input de entrada
    double coefEsteq;
    double rhoH;
    double r_d, r_p; //k1_sI, k2_sI, k1_sII, k2_sII,
    string inibidor; //chris - Hidratos
    string estruturaHidratos;
    //double V_w;
    //std::vector<double> Vw;

    void carregarCurvaHidrato(const std::string& nomeArquivo);
    void salvarCurvaDeslocada(const std::string& nomeArquivo);

    //bool checkHidrato(double P_atual, double T_atual);
    double interpolar(double x, const std::vector<double>& xData, const std::vector<double>& yData);


    std::tuple<std::vector<double>, std::vector<double>> gerarCurvaComInibidor(
        const std::vector<double>& tempBase,
        const std::vector<double>& pressBase,
        double K, double M, double w);
};

#endif

