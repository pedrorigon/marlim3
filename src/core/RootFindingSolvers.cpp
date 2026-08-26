#include "RootFindingSolvers.h"

// Only reached from here. It declares `using namespace std;` at file scope, and
// it declares templates named zbrent and zriddr in the global namespace -- the
// two the header goes out of its way to avoid colliding with. Confining it to
// this translation unit is the whole reason reportIterationLimit exists.
#include "FerramentasNumericas.h"

namespace rootfinding {

int sign(double value) {
    if (value <= 0.)
        return -1;
    else
        return 1;
}

void reportIterationLimit(const char *message) {
    NumError(message);
}

}  // namespace rootfinding
