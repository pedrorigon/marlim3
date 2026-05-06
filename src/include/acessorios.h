#ifndef _ACESSORIO
#define _ACESSORIO

#define _USE_MATH_DEFINES // para M_PI
#include <math.h>

#include <algorithm>
#include <fstream>
#include <complex>
using namespace std;
#include "Acidentes2.h"
#include "Bcsm2.h"
#include "multiBCS.h"
#include "BombaVol.h"
//#include "BombaVolVap.h"
#include "FonteMas.h"
#include "FonteMassCHK.h"
#include "PropFlu.h"
#include "PropVapor.h"
#include "FonteMasVap.h"
#include "PorosoRad-Simples.h"
#include "solverPoroso.h"
#include "variaveisGlobais1D.h"
class acessorio{
   public:
   int tipo;//0-> sem acess�rio,1->Fonte de g�s, 2->Fonte de l�quido,3->IPR,
   //4-> BCS, 5-> choke,6-> mudan�a de d�metro, 7-> delp, 8-> bomba volum�trica, 10->fonteMult
   //, 11-> IPR Vapor, 12->massica vapor, 14-> bombavolvap, 15->porosoRadial, 16->poroso2D, 17->multiBCS

   //9-> fontechoke

   BomCentSub bcs;
   multiBomCentSub multibcs;
   BomVol bvol;
   //BomVolVapor bvolvap;
   MudaArea mda;
   choke chk;
   IPR ipr;
   IPRVap iprvap;
   InjGas injg;
   InjLiq injl;
   InjMult injm;//alteracao7
   InjMultVap injmvap;
   fontemaschk fontechk;
   PorosRadSimp radialPoro;
   solverPoro poroso2D;
   ProFlu fluido;
   double delp;
   int tipoCompGas;
   double fatPoli;
   double eficLiq;
   double eficGas;
   acessorio(const int vtipo=0);//construtor
   acessorio(const acessorio& antigo);//construtor de c�pia
   acessorio& operator=(const acessorio& antigo);//sobrecarga do operador =
   //double dp(const double pres, const double temp,const Ccelula<double> cel)const{
   //return 0.;}
   //double dT(const double pres, const double temp,const Ccelula<double> cel)const{
   //return 0.;}
   //double dM(const double pres, const double temp,const Ccelula<double> cel)const{
   //return 0.;}
};

inline void carrega(const BomCentSub& acsr, acessorio& grupo){grupo.bcs=acsr;grupo.tipo=3;}
inline void carrega(const multiBomCentSub& acsr, acessorio& grupo){grupo.multibcs=acsr;grupo.tipo=17;}
inline void carrega(const MudaArea& acsr, acessorio& grupo){grupo.mda=acsr;grupo.tipo=5;}
inline void carrega(const choke& acsr, acessorio& grupo){grupo.chk=acsr;grupo.tipo=4;}
inline void carrega(const IPR& acsr, acessorio& grupo)
{
grupo.ipr=acsr;
grupo.fluido=acsr.FluidoPro;
grupo.tipo=3;}
inline void carrega(const InjGas& acsr, acessorio& grupo)
{
	grupo.injg=acsr;
	grupo.fluido=acsr.FluidoPro;
	grupo.tipo=1;}
inline void carrega(const InjLiq& acsr, acessorio& grupo)
{
	grupo.injl=acsr;
	grupo.fluido=acsr.FluidoPro;
	grupo.tipo=2;}
inline void carrega(const InjMult& acsr, acessorio& grupo)//alteracao7
{
	grupo.injm=acsr;
	grupo.fluido=acsr.FluidoPro;
	grupo.tipo=10;}
inline void carrega(const IPRVap& acsr, acessorio& grupo)
{grupo.iprvap=acsr;grupo.tipo=11;}
inline void carrega(const InjMultVap& acsr, acessorio& grupo)//alteracao7
{grupo.injmvap=acsr; grupo.tipo=12;}
//inline void carrega(const BomVolVapor& acsr, acessorio& grupo)//alteracao7
//{grupo.bvolvap=acsr; grupo.tipo=14;}







#endif
