// The runner supplies production bodies; property evaluations are controlled callbacks.
#include "SisProdConstants.h"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace sisprod;
using Word = std::uint64_t;
Word bits(double value) { return std::bit_cast<Word>(value); }
void require(bool condition) {
    if (!condition) throw std::runtime_error("inlet condition contract failed");
}

struct Properties {
    double gasDensity = 12.5, liquidDensity = 830.25, oilDensity = 790.5;
    double mixtureQuality = .17, hydrocarbonQuality = .81;
    double standardComplementaryDensity = 1000.25, complementaryDensity = 920.5;
    double relativeGasDensity = .71;
    double inletQuality = .37, inletFraction = .23, standardRatio = .000013;
    double pressure = 120.25, temperature = 57.5;
    int propertyFlags = 0;
};

struct Probe {
    const Properties &properties;
    const double &inletQuality, &inletFraction;
    std::vector<Word> events;
    void record(Word operation, double pressure = 0., double temperature = 0.) {
        events.insert(events.end(), {operation, bits(pressure), bits(temperature),
                                     bits(inletQuality), bits(inletFraction)});
        if (properties.propertyFlags) std::feraiseexcept(properties.propertyFlags);
    }
};

struct Fluid {
    Probe *probe;
    double Deng;
    double MasEspGas(double p, double t) const {
        probe->record(1, p, t); return probe->properties.gasDensity;
    }
    double MasEspLiq(double p, double t) const {
        probe->record(2, p, t); return probe->properties.liquidDensity;
    }
    double FracMassHidra(double p, double t) const {
        probe->record(3, p, t); return probe->properties.mixtureQuality;
    }
    double MasEspoleo(double p, double t) const {
        probe->record(6, p, t); return probe->properties.oilDensity;
    }
    double FracMass(double p, double t) const {
        probe->record(7, p, t); return probe->properties.hydrocarbonQuality;
    }
};
struct ComplementaryFluid {
    Probe *probe;
    double standardDensity() const {
        probe->record(4); return probe->properties.standardComplementaryDensity;
    }
    double MasEspFlu(double p, double t) const {
        probe->record(5, p, t); return probe->properties.complementaryDensity;
    }
};
struct Cel { Fluid flui; ComplementaryFluid fluicol; };
struct Input {
    int ConContEntrada, tipoFluido, flashCompleto;
    struct { double *bet; } CCPres;
};
struct TransientStepState {
    Cel *cells;
    Input &input;
    double &inletPressure, &inletTemperature, &inletQuality, &inletComplementaryFraction;
};

// @PRODUCTION_FUNCTIONS@

std::vector<Word> output;
Word cases = 0;
void check(const Properties &properties, int boundary, int type, int flash, int rounding) {
    double quality = properties.inletQuality, fraction = properties.inletFraction;
    double pressure = properties.pressure, temperature = properties.temperature;
    Probe probe{properties, quality, fraction, {}};
    Cel cell{{&probe, properties.relativeGasDensity}, {&probe}};
    std::array<double, 3> ratios{properties.standardRatio, 7.25, -31.5};
    const bool compositional = boundary == 1 && type == 0 && flash == 2;
    const bool gas = boundary == 1 && type == 1;
    Input input{boundary, type, flash, {gas ? ratios.data() : nullptr}};
    TransientStepState state{compositional || gas ? &cell : nullptr, input,
                             pressure, temperature, quality, fraction};
    require(std::fesetround(rounding) == 0);
    std::feclearexcept(FE_ALL_EXCEPT);
    refreshInletCondition(state);
    const int flags = std::fetestexcept(FE_ALL_EXCEPT);
    if (!compositional && !gas) {
        require(probe.events.empty());
        require(bits(quality) == bits(properties.inletQuality));
        require(bits(fraction) == bits(properties.inletFraction));
        require(flags == 0);
    } else {
        require(probe.events.size() == 25);
        const std::array<Word, 5> expected = compositional ? std::array<Word, 5>{1,2,3,4,5}
                                                         : std::array<Word, 5>{1,6,7,4,5};
        for (std::size_t index = 0; index < expected.size(); ++index) {
            require(probe.events[5*index] == expected[index]);
            require(probe.events[5*index+3] == bits(properties.inletQuality));
            require(probe.events[5*index+4] == bits(properties.inletFraction));
        }
        if (compositional) require(bits(fraction) == bits(properties.inletFraction));
    }
    require(bits(pressure) == bits(properties.pressure));
    require(bits(temperature) == bits(properties.temperature));
    output.insert(output.end(), {cases++, static_cast<Word>(boundary), static_cast<Word>(type),
                                static_cast<Word>(flash), static_cast<Word>(rounding),
                                bits(quality), bits(fraction), static_cast<Word>(flags), probe.events.size()});
    output.insert(output.end(), probe.events.begin(), probe.events.end());
}

int main(int argc, char **argv) {
    require(argc == 2);
    std::array<Properties, 20> properties;
    properties[1].mixtureQuality = properties[1].hydrocarbonQuality = 0.;
    properties[2].mixtureQuality = properties[2].hydrocarbonQuality = 1.;
    properties[3].mixtureQuality = properties[3].hydrocarbonQuality = std::nextafter(1., 0.);
    properties[4].mixtureQuality = properties[4].hydrocarbonQuality = std::nextafter(0., 1.);
    properties[5].inletFraction = properties[5].standardRatio = 0.;
    properties[6].inletFraction = 1.;
    properties[7].inletFraction = properties[7].standardRatio = -0.;
    properties[8].gasDensity = 0.;
    properties[9].liquidDensity = properties[9].oilDensity = 0.;
    properties[10].complementaryDensity = 0.;
    properties[11].relativeGasDensity = 0.;
    properties[12].mixtureQuality = properties[12].hydrocarbonQuality = std::bit_cast<double>(0x7ff8000000000001ull);
    properties[13].mixtureQuality = properties[13].hydrocarbonQuality = std::bit_cast<double>(0x7ff0000000000001ull);
    properties[14].gasDensity = std::numeric_limits<double>::infinity();
    properties[15].inletFraction = std::bit_cast<double>(0x7ff8000000000042ull);
    properties[16].gasDensity = 1e200; properties[16].oilDensity = 1e-100;
    properties[17].pressure = -0.; properties[17].temperature = -100.;
    properties[18].standardComplementaryDensity = 0.;
    properties[19].propertyFlags = FE_UNDERFLOW;
    const std::array<int, 4> roundings{FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    for (const auto &data : properties)
        for (int boundary : {-3, 0, 1, 2})
            for (int type : {-1, 0, 1, 2})
                for (int flash : {-1, 0, 1, 2, 3})
                    for (int rounding : roundings) check(data, boundary, type, flash, rounding);
    std::mt19937_64 random(481028);
    for (int index = 0; index < 256; ++index) {
        auto unit = [&]() { return (random() % 1000000 + 1) / 1000001.; };
        Properties data;
        data.gasDensity = 1. + 100. * unit();
        data.liquidDensity = 500. + 600. * unit();
        data.oilDensity = 400. + 500. * unit();
        data.mixtureQuality = unit(); data.hydrocarbonQuality = unit();
        data.inletFraction = unit(); data.standardRatio = .003 * unit();
        data.relativeGasDensity = .4 + unit();
        data.pressure = 1. + 500. * unit(); data.temperature = -20. + 200. * unit();
        check(data, 1, index % 2, index % 2 ? 3 : 2, roundings[index % 4]);
    }
    require(std::fesetround(FE_TONEAREST) == 0);
    auto *file = std::fopen(argv[1], "wb"); require(file != nullptr);
    require(std::fwrite(output.data(), sizeof(Word), output.size(), file) == output.size());
    require(std::fclose(file) == 0);
    std::printf("%llu inlet scenarios; exact output/event/fenv serialization\n", static_cast<unsigned long long>(cases));
}
