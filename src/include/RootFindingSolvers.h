#ifndef ROOTFINDINGSOLVERS_H_
#define ROOTFINDINGSOLVERS_H_

/// Generic root-finding algorithms.
///
/// These three solvers were members of SProd and took `prod` and `tipoCC`,
/// two production-domain selectors, purely to hand them back to
/// SProd::multMarcha on every evaluation. Neither means anything to Brent's
/// method. They are gone from the signatures: the caller builds a callable
/// that has already captured them, and the solver sees only `objective(x)`.
///
/// Everything here lives in namespace rootfinding, and that is load-bearing
/// rather than tidy. FerramentasNumericas.h already declares templates named
/// zbrent and zriddr at global scope, and SisProd.h includes it. Free
/// functions of those names would be ambiguous at best; at worst a call would
/// resolve silently to the other overload, which takes `const T *const` and
/// would compile.

// fabs and sqrt. This is the include SisProd.cpp used when the solvers lived
// there, so overload resolution is the one they were written against.
#include <math.h>

namespace rootfinding {

/// Magnitude of the first argument carrying the sign of the second.
///
/// inline, and defined here rather than in the .cpp, because zriddr calls it
/// three times per iteration. Out of line it would become a cross-TU call in
/// the one solver that is actually hot.
inline double SIGN(double a, double b) {
    return (b >= 0 ? 1.0 : -1.0) * fabs(a);
}

/// Sign of a value, as -1 or 1, with zero counting as negative.
///
/// It has no caller anywhere in the project and never had one; it is preserved
/// because removing dead code is a behaviour change this programme is not
/// authorised to make. Out of line precisely because nothing calls it.
int sign(double var);

/// Reports that a solver exhausted its iteration budget.
///
/// A one-line forward to NumError, which is declared in FerramentasNumericas.h
/// together with `using namespace std;` at file scope and the colliding
/// zbrent/zriddr templates. Reaching it through the .cpp keeps all of that out
/// of every translation unit that includes this header.
void reportIterationLimit(const char *message);

/// Finds a root by the false-position method, bisecting the bracket.
///
/// Reachable only from zbrent, which nothing calls, so it never executes. Its
/// verification is structural: refactor-harness/solver-move.py inverts the move
/// and compares the token stream against the pristine baseline.
template <typename Objective>
double falsacorda(double a, double b, Objective &&objective) {
    double u = objective(a);
    double e = b - a;
    double c;
    int maxit = 100;
    double delta = 0.001;
    double epsn = 0.001;
    double multFC = 0.5;

    for (int k = 1; k <= maxit; k++) { // this block treats the 'falsacorda' properly
        e = b - a;
        e *= 0.5;
        c = a + e;
        double w = objective(c);
        if (fabs(e) < delta || fabs(w) < epsn)
            return c;
        ((u > 0 && w < 0) || (u < 0 && w > 0)) ? (b = c) : (a = c, u = w);
    }
    return c;
}

/// Finds a root by Brent's method: bracketing with inverse quadratic
/// interpolation, falling back to false position when the interval does not
/// bracket a sign change.
///
/// Measured over the whole tree, this has NO call site. It is moved as it
/// stands, and never runs. Two consequences worth stating where they will be
/// read: L2 and L3 cannot see a defect introduced here, and because an
/// uninstantiated template is only parsed, neither can the compiler. What
/// covers it is solver-move.py for the move and verify-solvers.sh, which
/// instantiates and exercises it, for everything after.
template <typename Objective>
double zbrent(double x1, double x2, Objective &&objective, double tol, double epsn, int maxit) {
    double EPS = epsn;
    double a = x1;
    double b = x2;
    double c = x2;
    double fa = objective(a);
    double fb = objective(b);
    if (fabs(fa) > 1e9 || fabs(fb) > 1e9)
        return 1e10;
    double e = 0.;
    double d, fc, p, q, r, s, tol1, xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        double val;
        val = falsacorda(x1, x2, objective);
        return val;
    } else {
        fc = fb;
        for (int iter = 0; iter < maxit; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a;
                fc = fa;
                e = d = b - a;
            }
            if (fabs(fc) < fabs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);
            if (fabs(xm) <= tol1 || fb == 0.0)
                return b;
            if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0)
                    q = -q;
                p = fabs(p);
                double min1 = 3.0 * xm * q - fabs(tol1 * q);
                double min2 = fabs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (fabs(d) > tol1)
                b += d;
            else
                b += SIGN(tol1, xm);
            fb = objective(b);
        }
        reportIterationLimit("Metodo Van Winjngaarden-Dekker-Brent para calcular zero de funcaoo atingiu maximo de iteracoes");
        return 0.0;
    }
}

/// Finds a root by Ridders' method. The only solver here that executes.
///
/// Two callables, not one, because the original evaluates the objective in two
/// different ways. Twelve of its fourteen evaluations are raw; two are divided
/// by a convergence-monitor base and recorded in a member the outer pressure
/// loops read back. That scaling and that write are domain feedback, so they
/// travel in `monitor` and the solver keeps only the composition
/// `monitor(objective(x))` -- the same order of operations the original had.
///
/// `revPerm` selects nothing today: all three branches that test it have
/// identical arms. They are kept verbatim rather than collapsed. Collapsing
/// would be behaviour-preserving -- reading an int member has no side effect --
/// but the branch is the only surviving evidence that someone meant to treat
/// the reverse march differently here, and erasing it would make that
/// unrecoverable from the code. See A2-02 to A2-04 in evidencia/anomalias.md.
///
/// `minit` is a minimum iteration count derived from the input deck at the
/// binding site. It also gates three early returns that would otherwise be
/// unconditional, which is how the division at A2-05 becomes reachable.
template <typename Objective, typename Monitor>
double zriddr(double x1, double x2, Objective &&objective, Monitor &&monitor,
              int revPerm, int minit) {
    double xacc = 1e-5;
    int maxit = 100;
    double fmin;
    double xmin;
    double fl;
    double fh;
    if (revPerm == 0) {
        fl = objective(x1);
        fh = objective(x2);
    } else {
        fl = objective(x1);
        fh = objective(x2);
    }
    if (fabs(fl) > 1e9 || fabs(fh) > 1e9)
        return 1e10;
    if (fl >= 0.) {
        if (fl > 0.9e10) {
            x1 *= 0.9999;
            fl = objective(x1);
        } else {
            int konta = 0;
            while (konta < 100 && fl > 0.) {
                if (revPerm == 0) {
                    if (x2 < x1)
                        x1 *= 1.0001;
                    else
                        x1 *= 0.999;
                } else {
                    if (x2 < x1)
                        x1 *= 1.0001;
                    else
                        x1 *= 0.999;
                }
                fl = objective(x1);
                konta++;
            }
        }
    } else if (fh <= 0.) {
        if (fh < -0.9e10) {
            x2 *= 1.00001;
            fh = objective(x2);
        } else {
            int konta = 0;
            while (konta < 100 && fh < 0.) {
                if (revPerm == 0) {
                    if (x1 < x2)
                        x2 *= 1.0001;
                    else
                        x2 *= 0.999;
                } else {
                    if (x1 < x2)
                        x2 *= 1.0001;
                    else
                        x2 *= 0.999;
                }
                fh = objective(x2);
                konta++;
            }
        }
    }
    if (fabs(fh) < fabs(fl)) {
        fmin = fh;
        xmin = x2;
    } else {
        fmin = fl;
        xmin = x1;
    }
    if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
        double xl = x1;
        double xh = x2;
        double ans = -1.e20;
        for (int j = 0; j < maxit; j++) {
            double xm = 0.5 * (xl + xh);
            double fm = monitor(objective(xm));
            if (fabs(fm) < fabs(fmin)) {
                fmin = fm;
                xmin = xm;
            }
            double s = sqrt(fm * fm - fl * fh);
            if (s == 0.0) {
                fmin = objective(xmin);
                if(j>minit)return xmin;
            }
            double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
            if (fabs(xnew - ans) <= xacc) {
                fmin = objective(xmin);
                if(j>minit)return xmin;
            }
            ans = xnew;
            double fnew = monitor(objective(ans));
            if (fabs(fnew) < fabs(fmin)) {
                fmin = fnew;
                xmin = ans;
            }
            if (fabs(fnew) <= xacc) {
                fmin = objective(xmin);
                if(j>minit)return xmin;
            }
            if (SIGN(fm, fnew) != fm) {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fnew;
            } else if (SIGN(fl, fnew) != fl) {
                xh = ans;
                fh = fnew;
            } else if (SIGN(fh, fnew) != fh) {
                xl = ans;
                fl = fnew;
            } else
                return -1.e10;
            if (fabs(xh - xl) <= xacc) {
                fmin = objective(xmin);
                return xmin;
            }
        }
        return 1.e10;
    } else {
        if (fabs(fl) <= xacc) {
            return x1;
        }
        if (fabs(fh) <= xacc) {
            return x2;
        }
        return -1e10;
    }
}

}  // namespace rootfinding

#endif  // ROOTFINDINGSOLVERS_H_
