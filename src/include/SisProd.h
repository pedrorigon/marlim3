/*
 * SisProd.h
 *
 * Public interface of the SisProd architecture.
 *
 * This header keeps declarations, value types and inline accessors for:
 * - Root finding and solvers
 * - Stream mixing and fluid models
 * - Drift-flux correlations (R03)
 * - Black-oil properties (R04)
 * - Pressure-gradient correlations (R05)
 * - Hydraulic friction and equipment (R06/R10)
 *
 * Implementation files:
 *   - SisProd.cpp           : Main implementation
 *   - BlackOilProperties.cpp: R04 black-oil correlations
 *   - DriftFluxCorrelations.cpp: R03 drift-flux (C0/Ud)
 *   - HydraulicFriction.cpp : R06 friction + R10 GLV
 *   - PressureGradientEngine.cpp: R05 pressure-gradient dispatcher
 *
 * Legacy interface preserved in SisProd_old.h (class SProd).
 *
 * Implementation selection controlled by MARLIM_USE_NEW_SISPROD (see CMakeLists.txt).
 */

#ifndef MARLIM_SISPROD_H_
#define MARLIM_SISPROD_H_

#include <cstddef>
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace marlim {
namespace sisprod2 {

// ---------------------------------------------------------------------------
// Physical constants. Kept in the header because they are compile-time values
// shared by several translation units.
// ---------------------------------------------------------------------------
namespace constants {
const double kGravity = 9.80665;            // m / s^2
const double kStandardPressure = 1.01325e5; // Pa (one atmosphere)
const double kPi = 3.14159265358979323846;
const double kPressureFloor = 1.0e3;        // Pa, floor for density evaluation
} // namespace constants

// Relative-tolerance floating-point comparison (defined in the .cpp).
bool almostEqual(double a, double b, double relativeTolerance);

// ---------------------------------------------------------------------------
// Implementation selection flag. Reports which SisProd implementation the build
// selected through the MARLIM_USE_NEW_SISPROD compile option.
// ---------------------------------------------------------------------------
bool usingNewSisProd();
const char *activeSisProdImplementation();

// ---------------------------------------------------------------------------
// 1. Typed results (replace signalling by cout and termination by exit).
// ---------------------------------------------------------------------------
enum class SolveStatus {
    Ok,
    MaxIterations,
    OutOfRange,
    NotBracketed,
    NotANumber
};

const char *toString(SolveStatus status);

struct SolveResult {
    SolveStatus status;
    double value;      // e.g. the bottomhole pressure that was found
    int iterations;
    std::string detail;

    SolveResult() : status(SolveStatus::Ok), value(0.0), iterations(0) {}

    bool ok() const { return status == SolveStatus::Ok; }
};

struct StepResult {
    bool advanced;        // the time step was accepted
    bool needsRollback;   // at least one cell left the physical range
    double timeStepUsed;  // the time step actually applied
    std::string detail;

    StepResult() : advanced(false), needsRollback(false), timeStepUsed(0.0) {}
};

// ---------------------------------------------------------------------------
// 2. Simulation context (encapsulates the shared global state that today is
//    varGlob1D* vg1dSP and the (*vg1dSP).lixo5 clock). Trivial accessors stay
//    inline; there is no behaviour to move out.
// ---------------------------------------------------------------------------
class SimContext {
  public:
    SimContext()
        : simulationTime_(0.0), threadCount_(1), localTiny_(1e-15),
          maxGasOilRatio_(1e7) {}

    double simulationTime() const { return simulationTime_; }
    int threadCount() const { return threadCount_; }
    double localTiny() const { return localTiny_; }
    double maxGasOilRatio() const { return maxGasOilRatio_; }

    void setSimulationTime(double time) { simulationTime_ = time; }
    void setThreadCount(int count) { threadCount_ = count > 0 ? count : 1; }
    void setLocalTiny(double value) { localTiny_ = value; }
    void setMaxGasOilRatio(double value) { maxGasOilRatio_ = value; }

    void advanceTime(double timeStep) { simulationTime_ += timeStep; }

  private:
    double simulationTime_;
    int threadCount_;
    double localTiny_;
    double maxGasOilRatio_;
};

// ---------------------------------------------------------------------------
// 3. Diagnostics and observability (structured logging that replaces the cout
//    calls scattered through the physics).
// ---------------------------------------------------------------------------
enum class LogLevel { Debug, Info, Warning, Error };

const char *toString(LogLevel level);

struct LogRecord {
    LogLevel level;
    std::string operation;
    int sectionIndex;   // tramo index, or -1
    double time;
    std::string message;

    LogRecord() : level(LogLevel::Info), sectionIndex(-1), time(0.0) {}
};

class Diagnostics {
  public:
    virtual ~Diagnostics() {}
    virtual void log(const LogRecord &record) = 0;
    virtual void reportSolve(const std::string &operation,
                             const SolveResult &result) = 0;
};

// Silent sink, useful in tests and in parallel regions.
class NullDiagnostics : public Diagnostics {
  public:
    void log(const LogRecord &) override {}
    void reportSolve(const std::string &, const SolveResult &) override {}
};

// Writes structured lines to a std::ostream. Behaviour defined in the .cpp.
class OstreamDiagnostics : public Diagnostics {
  public:
    explicit OstreamDiagnostics(std::ostream &out);
    void log(const LogRecord &record) override;
    void reportSolve(const std::string &operation,
                     const SolveResult &result) override;
    int failureCount() const { return failureCount_; }

  private:
    std::ostream &out_;
    int failureCount_;
};

// ---------------------------------------------------------------------------
// 4. Generic root finders. Brent as the production solver and bisection as an
//    independent reference used by the automatic comparison. The residual is
//    injected as a std::function, which removes the multMarcha indirection and
//    unifies the duplicated zbrent and zriddr routines.
// ---------------------------------------------------------------------------
typedef std::function<double(double)> ResidualFunction;

class RootFinder {
  public:
    static SolveResult solveBracketed(const ResidualFunction &residual,
                                      double lowerBound, double upperBound,
                                      double tolerance = 1e-9,
                                      int maxIterations = 100);

    static SolveResult solveBisection(const ResidualFunction &residual,
                                      double lowerBound, double upperBound,
                                      double tolerance = 1e-9,
                                      int maxIterations = 200);
};

// ---------------------------------------------------------------------------
// 4b. Flow-correlation selector (used by ProductionColumn and R05 dispatcher).
//     Declared here — before ProductionColumn — so that marchToWellheadPhysical
//     can use it as a default parameter.
// ---------------------------------------------------------------------------
/// Correlation selector. IDs mirror the legacy integer codes used by
/// SProd::executarCorrelacao so the two paths can be cross-checked.
enum class FlowCorrelationId : int {
    PoettmannCarpenter = 0,
    BaxendellThomas    = 1,
    FancherBrown       = 2,
    HagedornBrown      = 3,
    DunsRos            = 4,
    Orkiszewski        = 5,
    BeggsAndBrill      = 6,
    MukherjeeeBrill    = 7,
    Aziz               = 8,
    Gray               = 9,
};

// ---------------------------------------------------------------------------
// 4c. R04 — Black-oil fluid-property correlations (pure, no global state).
//     Defined in SisProd_r04.cpp (requires PropFlu.h / MARLIM_BUILD).
//     Mirrors ProFlu methods for the flashCompleto==0 analytical branch.
//
// Units: pressure in kgf/cm2, temperature in °C, output in SI where noted.
// ---------------------------------------------------------------------------
/// Gas compressibility factor (Dranchuk-Abou-Kassem + Gopal seed).
/// @param Deng  Gas specific gravity (air = 1).
/// @param PC    Critical pressure (psia).
/// @param TC    Critical temperature (Rankine).
double zFactor(double pres, double temp, double Deng, double PC, double TC);

/// In-situ gas density via real-gas equation of state (kg/m3).
double gasDensityBlackOil(double pres, double temp,
                           double Deng, double PC, double TC);

/// Gas dynamic viscosity via Lee-Kesler correlation (cP).
double gasViscosityBlackOil(double pres, double temp,
                             double Deng, double PC, double TC);

/// Solution gas-oil ratio via Standing correlation (ft3/bbl).
/// @param Avb,Bvb,Cvb  Standing/Vasquez-Beggs coefficients from arq.*
double solutionGOR(double pres, double temp, double API,
                   double Deng, double Avb, double Bvb, double Cvb);

/// Oil formation volume factor via Vasquez-Beggs (RB/STB), below bubble point.
/// @param rs  Solution GOR in ft3/bbl (from solutionGOR or measured).
double oilFVF(double pres, double temp, double API, double Deng, double rs);

/// Water density via Meehan correlation (kg/m3).
double waterDensityBlackOil(double pres, double temp, double Denag);

/// Water formation volume factor (RB/STB).
double waterFVFBlackOil(double pres, double temp, double Denag);

/// Oil density for the analytical black-oil branch (kg/m3).
double oilDensityBlackOil(double pres, double temp, double API,
                          double Deng, double rs, double rDgD = 1.0);

/// Oil-water mixture density for the analytical black-oil branch (kg/m3).
double liquidDensityBlackOil(double pres, double temp, double API,
                             double Deng, double BSW, double Denag,
                             double rs, double rDgD = 1.0);

/// Water viscosity correlation used by ProFlu::VisAgua (cP).
double waterViscosityBlackOil(double temp);

/// Dead-oil viscosity via Beggs-Robinson (cP).
double deadOilViscosityBeggsRobinson(double temp, double API);

/// Dead-oil viscosity via ASTM interpolation (cP).
double deadOilViscosityASTM(double temp, double API,
                            double TempL, double LVisL,
                            double TempH, double LVisH);

/// Saturated live-oil viscosity via Beggs-Robinson (cP).
double oilViscosityBlackOil(double rs, double deadOilViscosity);

/// Oil-water mixture viscosity in the analytical branch (cP).
double liquidViscosityBlackOil(double pres, double temp, double API,
                               double Deng, double BSW, double Denag,
                               double rs, double rDgD,
                               double TempL, double LVisL,
                               double TempH, double LVisH);

/// Liquid specific heat (oil-water mixture) in J/(kg K).
double liquidSpecificHeatBlackOil(double pres, double temp, double API,
                                  double Deng, double BSW, double Denag,
                                  double rs, double rDgD = 1.0);

/// Gas specific heat in J/(kg K).
double gasSpecificHeatBlackOil(double pres, double temp, double Deng,
                               double PCis, double TCis,
                               double yco2 = 0.0, double rDgL = 1.0);

/// Derivative of oil density with respect to temperature (kg/m3/°C).
double liquidDensityDerivativeTBlackOil(double pres, double temp, double API,
                                        double Deng, double rs,
                                        double rDgD = 1.0);

/// Liquid thermal conductivity in W/(m K) for the analytical black-oil branch.
double liquidThermalConductivityBlackOil(double pres, double temp,
                                         double API, double Deng,
                                         double BSW, double Denag,
                                         double rs, double rDgD = 1.0);

/// Gas thermal conductivity in W/(m K) for the analytical black-oil branch.
double gasThermalConductivityBlackOil(double pres, double temp);

/// Liquid Joule-Thomson term used in the analytical RenovaTempPerm branch.
double liquidJouleThomsonBlackOil(double pres, double temp, double API,
                                  double Deng, double BSW, double Denag,
                                  double rs, double rDgD = 1.0,
                                  double liquidSimple = 0.0);

/// Gas Joule-Thomson term used in the analytical RenovaTempPerm branch.
double gasJouleThomsonBlackOil(double pres, double temp,
                               double Deng, double PCis, double TCis,
                               double rhog = -1.0);

/// Liquid enthalpy in J/kg for the analytical black-oil branch.
double liquidEnthalpyBlackOil(double pres, double temp,
                              double API, double Deng,
                              double BSW, double Denag,
                              double rs, double rDgD = 1.0);

/// Gas enthalpy in J/kg for the analytical black-oil branch.
double gasEnthalpyBlackOil(double pres, double temp,
                           double Deng, double PCis, double TCis,
                           double yco2 = 0.0, double rDgL = 1.0);

// ---------------------------------------------------------------------------
// 5. Standard-condition streams and stream mixing (Strategy). Deduplicates the
//    repeated core of RenovaMassPerm and its variants.
// ---------------------------------------------------------------------------
struct StandardStream {
    double oilRate;       // standard oil volumetric rate
    double waterRate;     // standard water volumetric rate
    double gasRate;       // standard gas volumetric rate
    double oilApi;        // API gravity of the oil
    double gasDensity;    // gas relative density
    double co2Fraction;   // molar CO2 fraction in the gas
    double waterDensity;  // water relative density

    StandardStream()
        : oilRate(0.0), waterRate(0.0), gasRate(0.0), oilApi(30.0),
          gasDensity(0.8), co2Fraction(0.0), waterDensity(1.0) {}

    double basicSedimentAndWater() const {
        const double liquid = oilRate + waterRate;
        return liquid > 0.0 ? waterRate / liquid : 0.0;
    }
    double gasOilRatio() const {
        return oilRate > 0.0 ? gasRate / oilRate : 0.0;
    }

    static double apiToRelativeDensity(double api) {
        return 141.5 / (131.5 + api);
    }
    static double relativeDensityToApi(double relativeDensity) {
        return 141.5 / relativeDensity - 131.5;
    }
};

// Combines two standard streams (defined in the .cpp).
StandardStream blend(const StandardStream &upstream,
                     const StandardStream &source, double tiny = 1e-15);

// Accessory types, today encoded as celula[i-1].acsr.tipo.
enum class AccessoryType {
    None = 0,
    GasSource = 1,
    LiquidSource = 2,
    Ipr = 3,
    MassSource = 10
};

// Mass-marching strategy. Replaces the giant duplicated switch of
// RenovaMassPerm. Behaviour defined in the .cpp.
class StreamMixing {
  public:
    static StandardStream march(const StandardStream &upstream,
                                AccessoryType type, const StandardStream &source,
                                const SimContext &context);
};

// ---------------------------------------------------------------------------
// 6. Black-oil fluid model (pure, deterministic). Constructors and getters stay
//    inline; the property correlations live in the .cpp.
// ---------------------------------------------------------------------------
class FluidModel {
  public:
    FluidModel()
        : liquidDensity_(800.0), standardGasDensity_(1.2),
          standardPressure_(constants::kStandardPressure), viscosity_(2.0e-3) {}

    FluidModel(double liquidDensity, double standardGasDensity, double viscosity)
        : liquidDensity_(liquidDensity), standardGasDensity_(standardGasDensity),
          standardPressure_(constants::kStandardPressure),
          viscosity_(viscosity) {}

    double liquidDensity() const { return liquidDensity_; }
    double viscosity() const { return viscosity_; }

    double gasDensityAt(double pressure) const;
    double massQualityFromGasOilRatio(double gasOilRatio) const;
    double mixtureDensity(double pressure, double massQuality) const;

  private:
    double liquidDensity_;
    double standardGasDensity_;
    double standardPressure_;
    double viscosity_;
};

/// Lightweight analytical black-oil adapter that maps a `StandardStream`
/// plus local pressure/temperature into the R04 pure-property helpers.
/// This is the first real consumer-side bridge between the new pressure march
/// and the ported analytical `flashCompleto==0` property correlations.
struct BlackOilState {
    double pressureKgfCm2;
    double temperatureC;
    double api;
    double gasSpecificGravity;
    double co2Fraction;
    double waterCut;
    double waterRelativeDensity;
    double solutionGor;
    double oilFvf;
    double waterFvf;
    double gasDensity;
    double oilDensity;
    double liquidDensity;
    double gasViscosity;
    double liquidViscosity;
    double liquidSpecificHeat;
    double gasSpecificHeat;
    double oilDensityDerivativeT;
    double zFactorGas;

    BlackOilState()
        : pressureKgfCm2(1.0), temperatureC(60.0), api(30.0),
          gasSpecificGravity(0.75), co2Fraction(0.0), waterCut(0.0),
          waterRelativeDensity(1.0), solutionGor(0.0), oilFvf(1.0),
          waterFvf(1.0), gasDensity(1.2), oilDensity(850.0),
          liquidDensity(850.0), gasViscosity(0.01), liquidViscosity(1.0),
          liquidSpecificHeat(2000.0), gasSpecificHeat(2000.0),
          oilDensityDerivativeT(0.0), zFactorGas(1.0) {}
};

/// Build a local black-oil property snapshot from a `StandardStream`.
/// Pressure is given in Pa and temperature in K to match the new architecture;
/// the implementation converts to the legacy analytical-correlation units.
BlackOilState makeBlackOilState(const StandardStream &stream,
                                double pressurePa,
                                double temperatureK);

// Darcy friction factor (defined in the .cpp).
double darcyFrictionFactor(double reynolds);

// ---------------------------------------------------------------------------
// 7. Production column and its geometry. Bottomhole is index 0, wellhead is the
//    last index. Each segment may carry an accessory source.
// ---------------------------------------------------------------------------
struct PipeSegment {
    double length;        // m
    double diameter;      // m
    double inclination;   // rad from horizontal (vertical = pi/2)
    AccessoryType accessory;
    StandardStream sourceStream;

    PipeSegment()
        : length(100.0), diameter(0.15), inclination(constants::kPi / 2.0),
          accessory(AccessoryType::None) {}

    double crossSectionArea() const {
        const double radius = 0.5 * diameter;
        return constants::kPi * radius * radius;
    }
    double verticalRise() const; // defined in the .cpp (uses std::sin)
};

struct ProfilePoint {
    double pressure;
    double mixtureDensity;
    double velocity;
    double gasOilRatio;
};

class ProductionColumn {
  public:
    ProductionColumn(const FluidModel &fluid, double massFlowRate)
        : fluid_(fluid), massFlowRate_(massFlowRate) {}

    void addSegment(const PipeSegment &segment) { segments_.push_back(segment); }
    void reserveSegments(std::size_t count) { segments_.reserve(count); }

    const FluidModel &fluid() const { return fluid_; }
    std::size_t segmentCount() const { return segments_.size(); }
    double massFlowRate() const { return massFlowRate_; }

    void setInletStream(const StandardStream &stream) { inletStream_ = stream; }
    const StandardStream &inletStream() const { return inletStream_; }

    double totalVerticalHeight() const;

    // Marches the pressure from bottomhole to wellhead and returns the wellhead
    // pressure. This loop is inherently sequential (each cell depends on the
    // previous one), so it must not be parallelised. Parallelism is applied one
    // level up, across independent columns (see solveBatch).
    double marchToWellhead(double bottomholePressure, const SimContext &context,
                           std::vector<ProfilePoint> *profileOut = nullptr) const;

    /// Same march but using a full two-phase pressure-gradient correlation
    /// (default: Beggs & Brill) instead of the simplified homogeneous model.
    /// This is the R05 physics-complete path.  Results differ from
    /// marchToWellhead because of holdup slip; the two paths converge as
    /// the flow becomes homogeneous (high Re, small drift velocity).
    double marchToWellheadPhysical(
        double bottomholePressure, const SimContext &context,
        FlowCorrelationId correlationId = FlowCorrelationId::BeggsAndBrill,
        std::vector<ProfilePoint> *profileOut = nullptr) const;

  private:
    FluidModel fluid_;
    double massFlowRate_;
    StandardStream inletStream_;
    std::vector<PipeSegment> segments_;
};

// ---------------------------------------------------------------------------
// 7b. First real R04 consumer extracted from the legacy RenovaTransMassPerm
//     block. Keeps the legacy formula but sources RS/Bo/Ba from BlackOilState.
// ---------------------------------------------------------------------------
struct PhaseTransferSideInput {
    double pressurePa;
    double temperatureK;
    double pressureAuxPa;
    double liquidRate;
    double waterCut;
    double gasSpecificGravity;
    double dissolvedGasDensityRatio;
    double pigFraction;
    StandardStream stream;

    PhaseTransferSideInput()
        : pressurePa(constants::kStandardPressure), temperatureK(333.15),
          pressureAuxPa(constants::kStandardPressure), liquidRate(0.0),
          waterCut(0.0), gasSpecificGravity(0.75),
          dissolvedGasDensityRatio(1.0), pigFraction(0.0) {}
};

struct PhaseTransferInput {
    PhaseTransferSideInput center;
    PhaseTransferSideInput left;
    PhaseTransferSideInput right;
    double cellLength;
    bool accessoryIsNone;

    PhaseTransferInput()
        : cellLength(1.0), accessoryIsNone(true) {}
};

/// Output bundle for the first integrated mass-march + phase-transfer step.
struct MassMarchResult {
    StandardStream mixedStream;
    double phaseTransferRate;

    MassMarchResult() : phaseTransferRate(0.0) {}
};

struct ThermalSideInput {
    double pressurePa;
    double temperatureK;
    double gasHoldup;
    double waterCut;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
    StandardStream stream;

    ThermalSideInput()
        : pressurePa(constants::kStandardPressure), temperatureK(333.15),
          gasHoldup(0.0), waterCut(0.0), gasSuperficialVelocity(0.0),
          liquidSuperficialVelocity(0.0) {}
};

struct ThermalFlowSnapshot {
    double liquidDensity;
    double gasDensity;
    double liquidSpecificHeat;
    double gasSpecificHeat;
    double liquidEnthalpy;
    double gasEnthalpy;
    double latentHeat;
    double liquidJouleThomson;
    double gasJouleThomson;
    double liquidViscosityPaS;
    double gasViscosityPaS;
    double mixedConductivity;
    double mixedSpecificHeat;
    double mixedDensity;
    double mixedViscosityPaS;

    ThermalFlowSnapshot()
        : liquidDensity(0.0), gasDensity(0.0), liquidSpecificHeat(0.0),
          gasSpecificHeat(0.0), liquidEnthalpy(0.0),
          gasEnthalpy(0.0), latentHeat(0.0),
          liquidJouleThomson(0.0),
          gasJouleThomson(0.0),
          liquidViscosityPaS(0.0), gasViscosityPaS(0.0),
          mixedConductivity(0.0), mixedSpecificHeat(0.0),
          mixedDensity(0.0), mixedViscosityPaS(0.0) {}
};

/// Pure helper extracted from the liquid/thermal-property block used by
/// `RenovaTempPerm`, sourcing analytical black-oil properties from R04.
ThermalFlowSnapshot computeThermalFlowSnapshot(const ThermalSideInput &in);

/// Pure helper extracted from the legacy `RenovaTransMassPerm` formula.
/// Returns the phase-transfer source term divided by cell length.
double computePhaseTransferRate(const PhaseTransferInput &in);

// ---------------------------------------------------------------------------
// 7c. R07 — Thermal marching updated from RenovaTempPerm.
//     Pure helper that computes temperature update based on energy balance
//     including convection, Joule-Thomson, heat transfer, and potential energy.
// ---------------------------------------------------------------------------
/// Input bundle for thermal update computation.
/// Mirrors the local variables used in RenovaTempPerm for energy balance.
struct ThermalUpdateInput {
    double currentTemperatureC;      // Current cell temperature in Celsius
    double upstreamTemperatureC;       // Temperature from upstream cell (for dT/dx)
    double pressurePa;               // Local pressure in Pascals
    double upstreamPressurePa;       // Upstream pressure (for dp/dx)
    double dx;                       // Cell length in meters
    double diameter;                 // Pipe inner diameter in meters
    double gasSuperficialVelocity;   // Gas superficial velocity (m/s)
    double liquidSuperficialVelocity; // Liquid superficial velocity (m/s)
    double gasHoldup;                // Gas volume fraction
    double waterCut;                 // Water cut in liquid phase
    double thermalConductivity;       // Mixed thermal conductivity (W/m·K)
    double externalTemperatureC;    // External/geothermal temperature (°C)
    double thermalResistance;         // Overall thermal resistance (K·m²/W)
    double massSourceL;              // Liquid mass source term (kg/m³/s)
    double massSourceG;              // Gas mass source term (kg/m³/s)
    double sourceTemperatureC;       // Temperature of mass source (°C)
    double latentHeatTerm;           // Latent heat term (W/m³)
    bool hasMassSource;              // Whether there's an active mass source
    
    ThermalUpdateInput()
        : currentTemperatureC(60.0), upstreamTemperatureC(60.0),
          pressurePa(constants::kStandardPressure), upstreamPressurePa(constants::kStandardPressure),
          dx(100.0), diameter(0.15),
          gasSuperficialVelocity(0.0), liquidSuperficialVelocity(0.0),
          gasHoldup(0.0), waterCut(0.0),
          thermalConductivity(1.0), externalTemperatureC(60.0),
          thermalResistance(1.0), massSourceL(0.0), massSourceG(0.0),
          sourceTemperatureC(60.0), latentHeatTerm(0.0),
          hasMassSource(false) {}
};

/// Output bundle for thermal update computation.
struct ThermalUpdateResult {
    double updatedTemperatureC;   // Updated temperature after march
    double heatFlux;            // Heat flux to external medium (W/m²)
    double dTdx;                // Temperature gradient computed
    bool lowFlowMode;           // True if fallback to external temperature was used
    int substepsUsed;           // Number of substeps for stability
    
    ThermalUpdateResult()
        : updatedTemperatureC(60.0), heatFlux(0.0), dTdx(0.0),
          lowFlowMode(false), substepsUsed(1) {}
};

/// Pure helper extracted from RenovaTempPerm.
/// Computes temperature update based on energy balance equation including:
/// - Convective term (dT/dx)
/// - Joule-Thomson expansion
/// - Heat transfer with external medium
/// - Potential energy change
/// - Mass source terms
/// - Latent heat effects
/// 
/// @param thermalProps Pre-computed thermal properties from computeThermalFlowSnapshot
/// @param in Input parameters for the thermal update
/// @param context Simulation context for tiny values and limits
/// @return ThermalUpdateResult with updated temperature and metadata
ThermalUpdateResult computeThermalUpdate(const ThermalFlowSnapshot &thermalProps,
                                          const ThermalUpdateInput &in,
                                          const SimContext &context);

// ---------------------------------------------------------------------------
// 8. Facade (TramoEngine) and the shared steady-state helpers.
// ---------------------------------------------------------------------------
struct SteadyStateRequest {
    double separatorPressure;   // wellhead boundary condition (Pa)
    double searchMargin;        // upper search margin above the hydrostatic guess
    double tolerance;
    int maxIterations;
    int sectionIndex;

    SteadyStateRequest()
        : separatorPressure(constants::kStandardPressure), searchMargin(2.0),
          tolerance(1e-3), maxIterations(100), sectionIndex(-1) {}

    SteadyStateRequest(double separator, double margin, double tol, int maxIter,
                       int section)
        : separatorPressure(separator), searchMargin(margin), tolerance(tol),
          maxIterations(maxIter), sectionIndex(section) {}
};

// Shared by the solver and the comparison so both paths solve the same
// equation. Defined in the .cpp.
ResidualFunction makeWellheadResidual(const ProductionColumn &column,
                                      const SimContext &context,
                                      double separatorPressure);

void wellheadSearchBounds(const ProductionColumn &column,
                          const SteadyStateRequest &request, double &lowerBound,
                          double &upperBound);

class TramoEngine {
  public:
    TramoEngine(SimContext &context, Diagnostics &diagnostics)
        : context_(context), diagnostics_(diagnostics) {}

    SolveResult solveBracketed(const std::string &operation,
                               const ResidualFunction &residual,
                               double lowerBound, double upperBound,
                               double tolerance, int maxIterations);

    SolveResult solveBottomholePressure(const ProductionColumn &column,
                                        const SteadyStateRequest &request);

    StandardStream marchMass(const StandardStream &upstream, AccessoryType type,
                             const StandardStream &source) const;

    MassMarchResult marchMassWithPhaseTransfer(
        const StandardStream &upstream, AccessoryType type,
        const StandardStream &source,
        const PhaseTransferInput &phaseTransferInput) const;

    ThermalFlowSnapshot buildThermalSnapshot(const ThermalSideInput &in) const;

    /// R07 thermal advancement step.
    /// Computes temperature update using energy balance from RenovaTempPerm logic.
    /// Integrates thermal property snapshot with convective, Joule-Thomson,
    /// heat transfer and source terms.
    ThermalUpdateResult advanceThermalStep(
        double currentTempC, double upstreamTempC,
        double pressurePa, double upstreamPressurePa,
        double dx, double diameter,
        double ugs, double uls,
        double gasHoldup, double waterCut,
        double externalTempC, double thermalResistance,
        const ThermalFlowSnapshot &props) const;

    const SimContext &context() const { return context_; }

  private:
    SimContext &context_;
    Diagnostics &diagnostics_;
};

// ---------------------------------------------------------------------------
// 9. Batch solver with OpenMP. Each column is solved independently, so the loop
//    is embarrassingly parallel and safe. Defined in the .cpp.
// ---------------------------------------------------------------------------
std::vector<SolveResult>
solveBatch(const std::vector<ProductionColumn> &columns,
           const std::vector<SteadyStateRequest> &requests,
           const SimContext &context);

// ---------------------------------------------------------------------------
// 10. Trend buffer — lightweight channel-based ring for observability.
//     Header-only (trivial inline); no behaviour moved to the .cpp.
// ---------------------------------------------------------------------------
class TrendBuffer {
  public:
    TrendBuffer(std::size_t channelCount, std::size_t capacity)
        : channels_(channelCount), capacity_(capacity) {}

    void append(std::size_t channel, double value) {
        if (channel < channels_.size() && channels_[channel].size() < capacity_)
            channels_[channel].push_back(value);
    }

    std::size_t channelCount() const { return channels_.size(); }
    std::size_t size(std::size_t channel) const {
        return channel < channels_.size() ? channels_[channel].size() : 0;
    }
    const std::vector<double> &channel(std::size_t index) const {
        return channels_[index];
    }

    void reset() {
        for (std::size_t c = 0; c < channels_.size(); ++c)
            channels_[c].clear();
    }

  private:
    std::vector<std::vector<double> > channels_;
    std::size_t capacity_;
};

// ---------------------------------------------------------------------------
// 11. R03 — Two-phase drift-flux correlations (C0 / Ud).
//
// All functions are pure (no global state, no virtual calls). They mirror the
// legacy SProd methods verbatim; the only change is that they are free
// functions taking explicit parameters instead of using `arq.*` fields.
//
// Naming convention:
//   c0  — distribution parameter (dimensionless)
//   ud  — drift velocity (m/s)
//   tet — pipe inclination angle (rad, +pi/2 = vertical upward)
//   rug — absolute roughness (m)
//   dia — inner diameter (m)
//   alf — void fraction (gas volume fraction, 0–1)
//   reymix  — mixture Reynolds number
//   reymixL — liquid-only Reynolds number (used by the Mod variant)
//   ug1, ul1 — gas / liquid volumetric flow rates (m³/s)
//   rhol, rhog — liquid / gas densities (kg/m³)
//   tensup — surface tension (N/m)
//   correcHor — horizontal correction factor for Ud
// ---------------------------------------------------------------------------
void bhagwatGhajar(double rhol, double rhog, double tensup, double alf,
                   double reymix, double reymixL, double ug1, double ul1,
                   double dia, double rug, double tet,
                   double &c0, double &ud, double correcHor);

void bhagwatGhajarMod(double rhol, double rhog, double tensup, double alf,
                      double reymix, double reymixL, double ug1, double ul1,
                      double dia, double rug, double tet,
                      double &c0, double &ud, double correcHor);

void choi(double rhol, double rhog, double tensup, double alf,
          double reymix, double reymixL, double ug1, double ul1,
          double dia, double rug, double tet,
          double &c0, double &ud, double correcHor);

void hibikiIshii(double rhol, double rhog, double tensup, double alf,
                 double reymix, double reymixL, double ug1, double ul1,
                 double dia, double rug, double tet,
                 double &c0, double &ud, double correcHor);

void francaLahey(double rhol, double rhog, double tensup, double alf,
                 double reymix, double reymixL, double ug1, double ul1,
                 double dia, double rug, double tet,
                 double &c0, double &ud, double correcHor);

// Dispatcher correlations — select the kernel via an integer mode parameter
// (replaces the arq.CorreDisper / arq.CorreAnular / arq.CorreEstrat switch).
// mode: 0=Choi 1=BhagwatGhajar 2=FrancaLahey 3=HibikiIshii
//       4=BhagwatGhajarMod 5=hybrid(angle-blended BhagwatGhajarMod/Choi)
void c0UdDisperso(double rhol, double rhog, double tensup, double alf,
                  double reymix, double reymixL, double ug1, double ul1,
                  double dia, double rug, double tet,
                  double &c0, double &ud, double correcHor, int mode);

void c0UdAnularChurn(double rhol, double rhog, double tensup, double alf,
                     double reymix, double reymixL, double ug1, double ul1,
                     double dia, double rug, double tet,
                     double &c0, double &ud, double correcHor, int mode);

void c0UdEstratificado(double rhol, double rhog, double tensup, double alf,
                       double reymix, double reymixL, double ug1, double ul1,
                       double dia, double rug, double tet,
                       double &c0, double &ud, double correcHor, int mode);

// ---------------------------------------------------------------------------
// 12. R05 — Two-phase pressure-gradient correlation dispatcher.
//
// FlowCorrelationId is declared in section 4b (above ProductionColumn).
// PressureGradientInput/Result are typed bundles that replace the long
// argument lists of the legacy GradientCorrelations free functions.
// computePressureGradient() calls the requested function and returns
// the unified result.  All physics is in GradientCorrelations.cpp (existing,
// pure C); this wrapper adds type safety and unit-tests the dispatch.
// ---------------------------------------------------------------------------

/// Input bundle for computePressureGradient().
struct PressureGradientInput {
    double inclinationRad;      // pipe inclination (rad, +π/2 = vertical up)
    double diameter;            // inner diameter (m)
    double roughness;           // absolute roughness (m)
    double pressure;            // local pressure (Pa)
    double mixtureVelocity;     // in-situ mixture velocity (m/s)
    double liquidFraction;      // input liquid volume fraction (no-slip)
    double gasDensity;          // in-situ gas density (kg/m³)
    double liquidDensity;       // in-situ liquid density (kg/m³)
    double gasViscosity;        // gas dynamic viscosity (Pa·s)
    double liquidViscosity;     // liquid dynamic viscosity (Pa·s)
    double surfaceTension;      // gas–liquid surface tension (N/m)
    double temperature;         // local temperature (K)
    double compressibilityFactor; // gas Z-factor (dimensionless)
    double waterFraction;       // water cut in liquid phase (0–1)
    double productionRate;      // volumetric production rate (m³/s) — FancherBrown only

    PressureGradientInput()
        : inclinationRad(constants::kPi / 2.0),
          diameter(0.10), roughness(1e-4), pressure(1e6),
          mixtureVelocity(1.0), liquidFraction(0.5),
          gasDensity(1.2), liquidDensity(850.0),
          gasViscosity(1.5e-5), liquidViscosity(2e-3),
          surfaceTension(0.03), temperature(320.0),
          compressibilityFactor(1.0), waterFraction(0.0),
          productionRate(0.0) {}
};

/// Output bundle for computePressureGradient().
struct PressureGradientResult {
    double holdup;       // in-situ liquid holdup (0–1)
    double frictionGrad; // friction pressure gradient (Pa/m)
    double gravityGrad;  // gravity pressure gradient (Pa/m)
    double accelGrad;    // acceleration pressure gradient (Pa/m)
    double totalGrad;    // total pressure gradient (Pa/m)
    double reynolds;     // mixture Reynolds number
    unsigned char flowType; // flow-pattern flag (correlation-specific)
    SolveStatus status;     // Ok unless the correlationId was unknown

    PressureGradientResult()
        : holdup(0.0), frictionGrad(0.0), gravityGrad(0.0),
          accelGrad(0.0), totalGrad(0.0), reynolds(0.0),
          flowType(0), status(SolveStatus::Ok) {}
};

/// Dispatcher: select the correlation by id and compute the pressure gradient.
/// Defined in SisProd_r05.cpp.
PressureGradientResult computePressureGradient(FlowCorrelationId correlationId,
                                               const PressureGradientInput &in);

// ---------------------------------------------------------------------------
// 13. R06/R10 — Hydraulic friction and gas-lift valve kernels.
//
// Pure functions with no global state. Ported verbatim from the SProd
// methods of the same name. References:
//   - colebrookFrictionFactor: Colebrook–White iterative friction factor
//     (identical algorithm used in bhagwatGhajar and in SProd::fric/Rey).
//   - areaValvCali: gas-lift GLV aperture calculation (R10, SProd::areaValvCali).
// ---------------------------------------------------------------------------
/// Darcy friction factor via Colebrook–White equation.
/// Turbulent (Re > 2400): Halland seed + two refinement steps (converges to
/// within 1e-5 of the exact implicit solution for any Re and ε/D).
/// Laminar: 64/Re.
///
/// @param reynolds   Mixture (or single-phase) Reynolds number.
/// @param relRoughness  Relative pipe roughness ε/D (dimensionless).
/// @returns Darcy friction factor f (dimensionless, f > 0).
double colebrookFrictionFactor(double reynolds, double relRoughness);

/// Gas-lift GLV opening fraction (0 = closed, 1 = fully open).
/// Mirrors SProd::areaValvCali verbatim; all inputs are in the legacy
/// unit system (pressures in psi-gauge unless noted, diameter in metres,
/// temperature in degrees Fahrenheit).
///
/// @param PCal     Calibration (dome) pressure (psi).
/// @param TCal     Calibration temperature (°F).
/// @param PVO      Casing pressure at valve depth (psi).
/// @param PT       Tubing pressure at valve depth (psi).
/// @param dextern  Outer diameter of the GLV body (m).
/// @param areagarg Throat area of the GLV (m²).
/// @param Rvalv    Bellows area / port area ratio (dimensionless, 0 < Rvalv < 1).
/// @param Temp     Bottomhole temperature (°F).
/// @returns Valve opening fraction in [0, 1].
double areaValvCali(double PCal, double TCal, double PVO, double PT,
                    double dextern, double areagarg, double Rvalv, double Temp);

// ---------------------------------------------------------------------------
// 14. Scenario builder used by tests and comparison (defined in the .cpp).
// ---------------------------------------------------------------------------
ProductionColumn buildVerticalWell(const FluidModel &fluid, double massFlowRate,
                                   std::size_t segmentCount, double segmentLength,
                                   double diameter, double inletGasOilRatio);

// ---------------------------------------------------------------------------
// 15. Automatic comparison between two implementations.
// ---------------------------------------------------------------------------
struct ComparisonRecord {
    std::string name;
    double reference;
    double candidate;
    double absoluteError;
    bool withinTolerance;
};

ComparisonRecord makeComparison(const std::string &name, double reference,
                                double candidate, double absoluteTolerance);

std::vector<ComparisonRecord> runReferenceComparison(double tolerance);

bool comparisonPassed(const std::vector<ComparisonRecord> &records);

// ---------------------------------------------------------------------------
// 16. Test entry points (defined in the .cpp).
// ---------------------------------------------------------------------------
bool runAllTests(bool verbose);
bool runSelfTest();

} // namespace sisprod2
} // namespace marlim

#endif // MARLIM_SISPROD_H_
