#ifndef TRANSIENTKINETICRESPONSE_H_
#define TRANSIENTKINETICRESPONSE_H_

#include <cmath>

namespace sisprod::thermal::detail {

/// Applies the 10 K/s kinetic-response safeguard only when kinetic work
/// reinforces an excessive total temperature change. Compensating work and
/// unrelated heat sources retain their contributions. Nonfinite values remain
/// visible to the caller's state validation. This is a numerical safeguard,
/// not a physical bound on temperature or a conservative energy correction.
inline double stabilizeKineticTemperatureResponse(double initialTemperature, double temperature,
                                            double kineticResponse, double timeStep) {
    constexpr double maximumTemperatureRate = 10.;
    if (!std::isfinite(temperature) || !std::isfinite(kineticResponse) ||
        !std::isfinite(timeStep) || timeStep <= 0.)
        return temperature;
    const double maximumChange = maximumTemperatureRate * timeStep;
    const double change = temperature - initialTemperature;
    if (std::abs(kineticResponse) <= maximumChange || std::abs(change) <= maximumChange ||
        std::signbit(change) != std::signbit(kineticResponse))
        return temperature;
    const double remainingChange = change - kineticResponse;
    const double stableChange = std::signbit(remainingChange) == std::signbit(kineticResponse)
        ? remainingChange + std::copysign(maximumChange, kineticResponse)
        : std::copysign(maximumChange, change);
    return initialTemperature + stableChange;
}

} // namespace sisprod::thermal::detail

#endif // TRANSIENTKINETICRESPONSE_H_
