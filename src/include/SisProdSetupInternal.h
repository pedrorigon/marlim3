#ifndef SISPRODSETUPINTERNAL_H_
#define SISPRODSETUPINTERNAL_H_

#include <iosfwd>
#include <string>
#include <vector>

class Cel;
template <class T> class FullMtx;

namespace {

/// Applies an operation to liquid, IPR, multiple-source and porous-reservoir
/// fluids, in source order. Gas-injection fluids are excluded. A contextual
/// operation receives the source cell's or reservoir cell's pressure and temperature.
template <typename Apply>
void forEachSourceFluid(Cel &cell, Apply &&apply);

/// Reads and converts the numeric values following the named PVTSim row.
/// A missing row or fewer than two values terminates the run with a file error.
template <typename Convert>
std::vector<double> readPvtSimRow(std::ifstream &file, std::string &word, const char *name,
                                   const std::string &fileName, Convert convert);

/// Tests strict axis growth; unordered comparisons, including NaN, fail.
bool strictlyIncreasing(const std::vector<double> &axis);

/// Replaces the output file with a matrix and an optional quantity header.
void writeTable(const std::string &path, const FullMtx<double> &table, const std::string &header = "");

/// Assigns zero-flow multiple, liquid and gas sources using the inlet fluids
/// and, for the multiple source, its temperature.
void attachEmptySources(Cel &cell, Cel &inlet);

}  // namespace

#endif  // SISPRODSETUPINTERNAL_H_
