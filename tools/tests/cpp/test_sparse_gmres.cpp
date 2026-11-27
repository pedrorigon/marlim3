// GMRES/FGMRES: an exact breakdown converges; a non-finite residual or a singular system ends the run as a failure.
#include "Matriz.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#define MARLIM_TEST_CAN_FORK 1
#endif

// Globals defined by the engine entry point.
string versao = "test_sparse_gmres";
int logRede = 0;
string nomeRedePrincipal;
Logger logger("");
int nthrdMatriz = 1;

namespace {

int failures = 0;
bool finished = false;

void check(bool condition, const char *what) {
    if (!condition) {
        std::printf("FAIL %s\n", what);
        ++failures;
    }
}

void failOnEarlyExit() {
    if (!finished) {
        std::printf("FAIL the solver ended the process instead of returning\n");
        std::fflush(stdout);
        std::_Exit(3);
    }
}

// 2I: A r = 2r for any residual, so the first Arnoldi step breaks down exactly.
SparseMtx<double> twiceIdentity() {
    double values[3] = {2., 2., 2.};
    int columns[3] = {0, 1, 2};
    int rowStart[4] = {0, 1, 2, 3};
    return SparseMtx<double>(3, 3, values, columns, rowStart);
}

// Rows 0 and 1 are equal, so A is singular and b = (1, 0, 0) has no solution.
SparseMtx<double> singularMatrix() {
    double values[5] = {1., 1., 1., 1., 2.};
    int columns[5] = {0, 1, 0, 1, 2};
    int rowStart[4] = {0, 2, 4, 5};
    return SparseMtx<double>(3, 5, values, columns, rowStart);
}

int solve(SparseMtx<double> &matrix, Vcr<double> &x, const Vcr<double> &b, bool flexible) {
    double eps = 1e-10;
    int iterations = 10;
    return flexible ? matrix.FGMRES(x, b, eps, iterations, 1, 2) : matrix.GMRES(x, b, eps, iterations, 1, 2);
}

void checkExactBreakdown(bool flexible) {
    SparseMtx<double> matrix = twiceIdentity();
    double rhs[3] = {2., 4., 6.};
    Vcr<double> b(3, rhs), x(3);
    check(solve(matrix, x, b, flexible) == 0, flexible ? "FGMRES status" : "GMRES status");
    for (int i = 0; i < 3; i++)
        check(std::fabs(x[i] - rhs[i] / 2.) < 1e-12, flexible ? "FGMRES solution" : "GMRES solution");
}

#ifdef MARLIM_TEST_CAN_FORK
enum class Failure { NonFiniteRightHandSide, SingularSystem };

// Exit status of a child that solves a system without a valid solution.
int childExitStatus(Failure failure, bool flexible, const std::string &logPath) {
    std::fflush(stdout);
    const pid_t child = fork();
    if (child == 0) {
        finished = true;
        logger.changeResultadoSimulacao("test_sparse_gmres", logPath);
        const bool singular = failure == Failure::SingularSystem;
        SparseMtx<double> matrix = singular ? singularMatrix() : twiceIdentity();
        double rhs[3] = {singular ? 1. : std::numeric_limits<double>::infinity(), 0., 0.};
        Vcr<double> b(3, rhs), x(3);
        solve(matrix, x, b, flexible);
        std::_Exit(0);
    }
    int status = 0;
    waitpid(child, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

void checkFailure(Failure failure, bool flexible, const char *diagnostic) {
    const std::string logPath = std::string("test_sparse_gmres_") + (flexible ? "fgmres" : "gmres") + ".log";
    std::remove(logPath.c_str());
    check(childExitStatus(failure, flexible, logPath) == EXIT_FAILURE, flexible ? "FGMRES exit status" : "GMRES exit status");
    std::ifstream in(logPath);
    const std::string log((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    check(log.find("\"status\": \"falha\"") != std::string::npos, "log status");
    check(log.find(diagnostic) != std::string::npos, diagnostic);
    std::remove(logPath.c_str());
}
#endif

} // namespace

int main() {
    logger.changeResultadoSimulacao("test_sparse_gmres", "test_sparse_gmres_main.log");
    std::atexit(failOnEarlyExit);
    checkExactBreakdown(false);
    checkExactBreakdown(true);
#ifdef MARLIM_TEST_CAN_FORK
    checkFailure(Failure::NonFiniteRightHandSide, false, "resíduo inicial não finito");
    checkFailure(Failure::NonFiniteRightHandSide, true, "resíduo inicial não finito");
    checkFailure(Failure::SingularSystem, false, "sistema singular");
    checkFailure(Failure::SingularSystem, true, "sistema singular");
#else
    std::printf("SKIP failure checks: fork is not available\n");
#endif
    finished = true;
    std::remove("test_sparse_gmres_main.log");
    if (failures == 0)
        std::printf("test_sparse_gmres: all checks passed\n");
    return failures == 0 ? 0 : 1;
}
