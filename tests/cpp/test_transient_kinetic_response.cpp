#include "TransientKineticResponse.h"

#include <cmath>
#include <cstdio>
#include <limits>

int main() {
    using sisprod::thermal::detail::stabilizeKineticTemperatureResponse;
    int failures = 0;
    const auto check = [&](double initial, double temperature, double response,
                           double dt, double expected, const char *message) {
        const double actual = stabilizeKineticTemperatureResponse(initial, temperature, response, dt);
        if (actual != expected) {
            std::fprintf(stderr, "%s: got %.17g, expected %.17g\n", message, actual, expected);
            ++failures;
        }
    };
    check(40., 42., 100., 1., 42., "Preserve compensated kinetic work");
    check(40., 32., 100., .5, 32., "Preserve compensating cooling");
    check(40., 48., -100., .5, 48., "Preserve compensating heating");
    check(40., 160., 100., 1., 70., "Retain other positive heat");
    check(40., -80., -100., 1., 10., "Retain other negative heat");
    check(40., 120., 100., 1., 50., "Bound excessive compensated heating");
    check(40., -40., -100., 1., 30., "Bound excessive compensated cooling");
    check(40., 140., 0., 1., 140., "Preserve nonkinetic heating");
    check(40., 50., 100., 1., 50., "Preserve threshold equality");
    check(40., 140., 100., 0., 140., "Leave invalid time step for validation");
    const double infinity = std::numeric_limits<double>::infinity();
    check(40., infinity, 100., 1., infinity, "Preserve invalid temperature");
    check(40., 42., infinity, 1., 42., "Preserve finite solution with unavailable response");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!std::isnan(stabilizeKineticTemperatureResponse(40., nan, 100., 1.)))
        ++failures;
    return failures ? 1 : 0;
}
