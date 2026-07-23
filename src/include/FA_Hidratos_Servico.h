// ---

#ifndef _FA_HIDRATO_SERVICO
#define _FA_HIDRATO_SERVICO

#include <vector>
#include <string>
#include <tuple>
#include "SisProd_old.h"


class FA_Hidrato_Servico {
public:
    FA_Hidrato_Servico(const SProd& sistema);
    ~FA_Hidrato_Servico();


    void solverHidratoG();

    bool checkHidratoG(double P, double T);

    double TurnerHidratoG(double P, double T, double P_eq, double T_eq,
                              double A_s, double V_h, double V_w,
                              double r_d, double r_p,
                              const std::string& estruturaHidratos,
                              double A_s_input,
                              double &A_s_eff_out);

    double Euler1ordemGasG(int i, double j_G, double taxaCinetica,
                               double A, double MM_g);

    double Euler1ordemAguaG(int i, double j_W, double taxaCinetica,
                                double A, double eta, double MM_w);

    double Euler1ordemHidratoG(double j_H, double taxaCinetica,
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

    double fracFWcarregada; //chris - Hidratos


    void carregarCurvaHidratoG(const std::string& nomeArquivo);
    void salvarCurvaDeslocadaG(const std::string& nomeArquivo);

    //bool checkHidrato(double P_atual, double T_atual);
    double interpolarG(double x, const std::vector<double>& xData, const std::vector<double>& yData);


    std::tuple<std::vector<double>, std::vector<double>> gerarCurvaComInibidorG(
        const std::vector<double>& tempBase,
        const std::vector<double>& pressBase,
        double K, double M, double w);
};

#endif


