/*
 * ChokeGas.h
 *
 * R10 Phase 1 — Gas/liquid flow through choke and Venturi restrictions.
 *
 * Port from src/include/chokegas.h and src/core/chokegas.cpp into the
 * marlim::sisprod2 namespace following the modular architecture.
 *
 * Functions / class declared here (implemented in ChokeGas.cpp):
 *   - ChokeGas — public API wrapping throat mass flow, pressure-ratio,
 *     throat-temperature, Venturi flow, and fluid-density calculations.
 *
 * Migration reference: issues/sisprod-migration-plan.md, R10 Phase 1.
 */

#ifndef MARLIM_CHOKEGAS_H_
#define MARLIM_CHOKEGAS_H_

#include <cmath>
#include <algorithm>
#include "PropFlu.h"

namespace marlim {
namespace sisprod2 {

//----------------------------------------------------------------------------
// Constant correlation table: ventCR[5][9][10]
//
// Indexed as [sgg_group][pressure_group][temperature_group].
// sgg_group: maps SGG (gas specific gravity) to one of 5 bins.
// pressure_group: maps upstream pressure (bar) to one of 9 bins (0–400 bar).
// temperature_group: maps upstream temperature (C) to one of 10 bins.
//----------------------------------------------------------------------------

// clang-format off
inline constexpr double kVentCR[5][9][10] = {
    // sgg_group = 0  (SGG ≈ 0.55)
    {
        {0.6697, 0.6685, 0.6671, 0.6655, 0.6639, 0.6623, 0.6606, 0.6590, 0.6574, 0.6559},  // P =   0 bar
        {0.7044, 0.6959, 0.6891, 0.6835, 0.6788, 0.6746, 0.6709, 0.6677, 0.6647, 0.6621},  // P =  50 bar
        {0.7486, 0.7286, 0.7142, 0.7032, 0.6945, 0.6874, 0.6814, 0.6763, 0.6719, 0.6681},  // P = 100 bar
        {0.7971, 0.7633, 0.7401, 0.7232, 0.7102, 0.6998, 0.6915, 0.6845, 0.6786, 0.6736},  // P = 150 bar
        {0.8370, 0.7936, 0.7633, 0.7411, 0.7242, 0.7110, 0.7004, 0.6917, 0.6845, 0.6784},  // P = 200 bar
        {0.8613, 0.8148, 0.7807, 0.7552, 0.7356, 0.7202, 0.7078, 0.6976, 0.6892, 0.6822},  // P = 250 bar
        {0.8719, 0.8265, 0.7918, 0.7648, 0.7437, 0.7268, 0.7132, 0.7020, 0.6928, 0.6850},  // P = 300 bar
        {0.8732, 0.8308, 0.7971, 0.7703, 0.7486, 0.7311, 0.7168, 0.7049, 0.6950, 0.6867},  // P = 350 bar
        {0.8685, 0.8299, 0.7982, 0.7722, 0.7508, 0.7332, 0.7186, 0.7064, 0.6962, 0.6875},  // P = 400 bar
    },
    // sgg_group = 1  (SGG ≈ 0.65)
    {
        {0.6640, 0.6624, 0.6607, 0.6590, 0.6572, 0.6555, 0.6539, 0.6523, 0.6508, 0.6494},
        {0.7094, 0.6976, 0.6887, 0.6816, 0.6758, 0.6709, 0.6666, 0.6630, 0.6598, 0.6570},
        {0.7750, 0.7436, 0.7228, 0.7078, 0.6964, 0.6874, 0.6801, 0.6740, 0.6690, 0.6647},
        {0.8498, 0.7949, 0.7596, 0.7353, 0.7175, 0.7039, 0.6933, 0.6847, 0.6777, 0.6718},
        {0.9016, 0.8362, 0.7914, 0.7597, 0.7365, 0.7189, 0.7052, 0.6943, 0.6855, 0.6781},
        {0.9246, 0.8602, 0.8129, 0.7778, 0.7513, 0.7309, 0.7149, 0.7021, 0.6918, 0.6833},
        {0.9286, 0.8698, 0.8242, 0.7887, 0.7611, 0.7393, 0.7219, 0.7079, 0.6965, 0.6870},
        {0.9222, 0.8699, 0.8277, 0.7936, 0.7662, 0.7441, 0.7262, 0.7115, 0.6994, 0.6894},
        {0.9101, 0.8642, 0.8258, 0.7940, 0.7677, 0.7460, 0.7281, 0.7132, 0.7009, 0.6905},
    },
    // sgg_group = 2  (SGG ≈ 0.75)
    {
        {0.6590, 0.6572, 0.6553, 0.6535, 0.6517, 0.6500, 0.6484, 0.6469, 0.6455, 0.6442},
        {0.7193, 0.7029, 0.6910, 0.6820, 0.6749, 0.6691, 0.6642, 0.6602, 0.6567, 0.6537},
        {0.8257, 0.7710, 0.7389, 0.7175, 0.7022, 0.6906, 0.6816, 0.6743, 0.6683, 0.6634},
        {0.9429, 0.8494, 0.7930, 0.7565, 0.7313, 0.7130, 0.6992, 0.6884, 0.6798, 0.6728},
        {0.9978, 0.9016, 0.8356, 0.7898, 0.7572, 0.7333, 0.7152, 0.7012, 0.6901, 0.6811},
        {1.0083, 0.9230, 0.8592, 0.8115, 0.7758, 0.7487, 0.7278, 0.7115, 0.6984, 0.6880},
        {0.9984, 0.9254, 0.8676, 0.8222, 0.7865, 0.7585, 0.7363, 0.7187, 0.7045, 0.6929},
        {0.9796, 0.9175, 0.8664, 0.8247, 0.7908, 0.7633, 0.7410, 0.7229, 0.7081, 0.6960},
        {0.9573, 0.9044, 0.8596, 0.8219, 0.7904, 0.7642, 0.7426, 0.7246, 0.7097, 0.6973},
    },
    // sgg_group = 3  (SGG ≈ 0.85)
    {
        {0.6548, 0.6528, 0.6509, 0.6490, 0.6472, 0.6456, 0.6440, 0.6426, 0.6413, 0.6401},
        {0.7366, 0.7125, 0.6964, 0.6849, 0.6761, 0.6691, 0.6635, 0.6589, 0.6550, 0.6517},
        {0.9346, 0.8211, 0.7664, 0.7341, 0.7127, 0.6974, 0.6860, 0.6770, 0.6699, 0.6640},
        {1.0950, 0.9402, 0.8476, 0.7907, 0.7537, 0.7281, 0.7096, 0.6957, 0.6849, 0.6763},
        {1.1251, 0.9930, 0.8996, 0.8342, 0.7882, 0.7552, 0.7310, 0.7126, 0.6985, 0.6873},
        {1.1083, 1.0017, 0.9195, 0.8571, 0.8099, 0.7742, 0.7469, 0.7258, 0.7093, 0.6962},
        {1.0769, 0.9905, 0.9204, 0.8644, 0.8198, 0.7846, 0.7567, 0.7345, 0.7168, 0.7025},
        {1.0418, 0.9709, 0.9114, 0.8621, 0.8214, 0.7882, 0.7611, 0.7390, 0.7209, 0.7061},
        {1.0072, 0.9482, 0.8974, 0.8542, 0.8177, 0.7870, 0.7614, 0.7400, 0.7223, 0.7075},
    },
    // sgg_group = 4  (SGG ≈ 0.95)
    {
        {0.6511, 0.6491, 0.6471, 0.6453, 0.6436, 0.6420, 0.6405, 0.6391, 0.6379, 0.6368},
        {0.7677, 0.7285, 0.7058, 0.6905, 0.6794, 0.6710, 0.6644, 0.6591, 0.6547, 0.6510},
        {1.1985, 0.9259, 0.8157, 0.7617, 0.7298, 0.7087, 0.6936, 0.6823, 0.6735, 0.6666},
        {1.3126, 1.0831, 0.9358, 0.8451, 0.7885, 0.7514, 0.7258, 0.7072, 0.6933, 0.6825},
        {1.2785, 1.1115, 0.9866, 0.8966, 0.8326, 0.7870, 0.7540, 0.7296, 0.7112, 0.6969},
        {1.2199, 1.0943, 0.9937, 0.9151, 0.8547, 0.8085, 0.7731, 0.7460, 0.7249, 0.7083},
        {1.1609, 1.0631, 0.9815, 0.9148, 0.8608, 0.8176, 0.7831, 0.7556, 0.7336, 0.7159},
        {1.1069, 1.0285, 0.9614, 0.9048, 0.8575, 0.8183, 0.7860, 0.7595, 0.7378, 0.7200},
        {1.1069, 1.0285, 0.9614, 0.9048, 0.8575, 0.8183, 0.7860, 0.7595, 0.7378, 0.7200},
    },
};
// clang-format on

//----------------------------------------------------------------------------
// ChokeGas
//
// Models gas and liquid flow through a choke or Venturi restriction.
// Mirrors every field and method from the legacy ``chokegas.h`` class.
//----------------------------------------------------------------------------

class ChokeGas {
public:
    ProFlu flui;          // PVT fluid-property model.
    int tipo;             // Choke / flow-model type  (0 = simple gas, 1 = Venturi, ≥2 = liquid).
    double areagarg;      // Choke throat area (m²).
    double cd;            // Gas discharge coefficient.
    double cdliq;         // Liquid discharge coefficient.
    double areafole;      // Venturi or restriction area (m²).
    double pcalib;        // Calibration pressure (kgf/cm²).
    double tcalib;        // Calibration temperature (°C).
    double dextern;       // External diameter (m).
    double presEstag;     // Upstream pressure (kgf/cm²).
    double presGarg;      // Throat pressure (kgf/cm²).
    double tempEstag;     // Upstream temperature (°C).
    double frec;          // Gas correction factor.
    double frecliq;       // Liquid correction factor.
    double qGarg;         // Volumetric flow through throat (m³/day or m³/s depending on caller).
    double tempGarg;      // Throat temperature (°C).

    //--- construction --------------------------------------------------------
    ChokeGas()
        : flui(), tipo(0), areagarg(0.), cd(1.), cdliq(1.), areafole(0.),
          pcalib(1.), tcalib(1.), dextern(1.),
          presEstag(10.), presGarg(10.), tempEstag(0.),
          frec(0.), frecliq(0.), qGarg(0.), tempGarg(0.) {}

    //--- public API (ported verbatim from legacy chokegas.cpp) ----------------
    double calculateMassFlow(int fluido = 0, double salin = 0.);
    double throatTemperature();
    double criticalPressureRatio() const;
    double fluidDensity(double pres, double temp, double salin) const;

    // Pressure-ratio helpers (used internally & by tests):
    double calcPressureRatio(double mass, double rp) const;
    double calcThroatArea() const { return areagarg; }
    double getQgarg() const { return qGarg; }
    double getTempGarg() const { return tempGarg; }
    double getCd() const { return cd; }

    //--- internal helpers (used by other R10 code / future tests) ------------
    int sggBin(double sgg) const;
    double interpolateVentCR(double sgg, double pbar, double tc) const;
    double calcAdiabaticIndex() const;
    double calcZFactor() const { return flui.Zdran(presEstag, tempEstag); }
    double gasDensity(double p, double t) const { return flui.MasEspGas(p, t); }

private:
    //--- private numerics (ported verbatim) -----------------------------------
    double fraiz(double kad, double rt, double mass, double rp) const;
    double derraiz(double kad, double rp) const;
    double newtonRoot(double kad, double rt, double mass, double rp) const;
    double razpresSimples(double mass, double rp) const;
    double fMachVenturi(double Mach, double kad) const;
    double deriFMachVenturi(double Mach, double kad) const;
    double findMachVenturi(double kad) const;
    double criticalPressureRatioVenturi() const;

    //--- unit-conversion helpers (ported) -------------------------------------
    double psia(double p) const { return (p * 0.9678411) * 14.69595; }
    double faren(double t) const { return 1.8 * t + 32; }
};

}  // namespace sisprod2
}  // namespace marlim

#endif // MARLIM_CHOKEGAS_H_
