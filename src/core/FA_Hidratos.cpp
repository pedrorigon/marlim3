#include "FA_Hidratos.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>

//CHECK DE HIDRATO: pensar na questão de unidades:
//-- geralmente as tabelas vem em [bar] e o Marlim 3 é em [kgf/cm2].
//Conversão: 1bar --> 1.01972 kgf/cm2
//1kgf/cm2 --> 0.98066864218981 bar

            std::string datahora_atual() {
                std::time_t t = std::time(nullptr);
                std::tm* now = std::localtime(&t);

                std::ostringstream ss;
                ss << "datahora = "
                   << std::setfill('0')
                   << std::setw(2) << now->tm_mday << "/"
                   << std::setw(2) << now->tm_mon + 1 << "/"
                   << (now->tm_year + 1900) << " "
                   << std::setw(2) << now->tm_hour << ":"
                   << std::setw(2) << now->tm_min  << ":"
                   << std::setw(2) << now->tm_sec;

                return ss.str();
            }

//std::vector<double> FVH_global; //chris - hidratos
constexpr double G_MOL_TO_KG_MOL = 1.0 / 1000.0; //não precisa mais --> Chris: desde que o user ente com g/mol.

FA_Hidrato::FA_Hidrato(const SProd& sistemaRef) : sistema(sistemaRef) {


}

FA_Hidrato::~FA_Hidrato() {
   // std::cout << "Destruindo FA_Hidrato." << std::endl;
}

void FA_Hidrato::carregarCurvaHidrato(const std::string& nomeArquivo) {
}



void FA_Hidrato::salvarCurvaDeslocada(const std::string& nomeArquivo) {

}


std::tuple<std::vector<double>, std::vector<double>> FA_Hidrato::gerarCurvaComInibidor(const std::vector<double>& tempBase,
	    const std::vector<double>& pressBase,
	    double K, double M, double w) {

}


bool FA_Hidrato::checkHidrato(double P_atual, double T_atual) {
    double P_eq = interpolar(T_atual, temperaturaCurvaDeslocada, pressaoCurvaDeslocada);
    return -1e15;
}

void FA_Hidrato::logModeloHidratoTxt()
{
}


double FA_Hidrato::interpolar(double x,
                                   const std::vector<double>& xData,
                                   const std::vector<double>& yData) {return -1e15;
}

double FA_Hidrato::TurnerHidrato(double P, double T, double P_eq, double T_eq,
        double A_s, double V_h, double V_w,
        double r_d, double r_p, const std::string& estruturaHidratos,   double A_s_input, double &A_s_eff_out) {return -1e15;
}

void FA_Hidrato::salvarMassaHidratoTotal(const std::string& nomeArquivo)
{
}

//*** Euler assume a mesma Massa Molar de Hidrato, pois está é praticamente igual entre as estruturas sI (119.5g/mol)
// e sII (117.9g/mol)
double FA_Hidrato::Euler1ordemHidrato(double j_H, double taxaCinetica,
                                           double A, double eta, double MM_h, double rho_h) {return -1e15;
}

double FA_Hidrato::Euler1ordemGas(int i, double j_G, double taxaCinetica,
                                       double A, double MM_g) {return -1e15;
}

double FA_Hidrato::Euler1ordemAgua(int i, double j_W, double taxaCinetica,
                                        double A, double eta, double MM_w) {return -1e15;
}

void FA_Hidrato::solverHidrato1() {

}

void FA_Hidrato::solverHidrato2() {

}

void FA_Hidrato::solverHidrato3() {


}

void FA_Hidrato::solverHidrato() {
}
