#ifndef SISPRODTRENDOUTPUT_H_
#define SISPRODTRENDOUTPUT_H_

/// Trend-file output for the production and service lines.
///
/// A trend is a time series recorded at one fixed cell. Eight writers produce
/// them: two roles (header, rows) times four data sources (production line,
/// service line, and the cross-section temperature trend of each). None of them
/// reads the cells or writes simulation state -- they read the trend buffers and
/// write files, and nothing here feeds back into the numerics.
///
/// The writers take their state through TrendState instead of reading it from
/// SProd: it names, in one place, the state a trend writer is allowed to touch,
/// and it makes the writers callable without an SProd, so they can be exercised
/// against synthetic buffers.

// Declared, not included: TrendState only holds a reference and a pointer to
// these, so the header stays free of the input-deck and globals headers.
class Ler;
struct varGlob1D;

namespace trendoutput {

/// One trend group: the sample buffer and the window still to be written, kept
/// together so that no writer can pair a buffer with another group's counters.
struct TrendSeries {
    /// Samples, indexed [series][row][column] -- MatTrend*.
    double ***samples;
    /// Rows recorded so far, per series -- ntrend*.
    const int *count;
    /// Rows already written, per series -- ntrend*B.
    const int *countBase;
};

/// The state a trend writer reads, and the only state it may read.
///
/// Every field is a reference or a pointer into the owner's state, never a
/// copy: the call site advances the counters between the header call and the
/// row call, so a copy would be read at the wrong moment.
struct TrendState {
    /// Input deck and configuration -- SProd::arq.
    const Ler &input;
    /// Shared 1D globals; only sequenciaAP is read -- SProd::vg1dSP.
    const varGlob1D *globals;
    /// Branch index; negative outside a network -- SProd::indTramo.
    const int &branchIndex;
    /// Count of output passes so far; captions go out on the first --
    /// SProd::kimpT.
    const double &printPassCount;

    /// MatTrendP, ntrend, ntrendB.
    TrendSeries production;
    /// MatTrendG, ntrendg, ntrendgB.
    TrendSeries gasLine;
    /// MatTrendTransP, ntrendtrans, ntrendtransB.
    TrendSeries productionCrossSection;
    /// MatTrendTransG, ntrendtransg, ntrendtransgB.
    TrendSeries gasLineCrossSection;
};

/// Writes the column captions of a production-line trend file, truncating it.
void writeProductionTrendHeader(const TrendState &state, int trendIndex, int networkIndex);

/// Appends the buffered rows of a production-line trend file.
void writeProductionTrendRows(const TrendState &state, int trendIndex, int networkIndex);

/// Writes the column captions of a service-line trend file, truncating it.
void writeGasLineTrendHeader(const TrendState &state, int trendIndex, int networkIndex);

/// Appends the buffered rows of a service-line trend file.
void writeGasLineTrendRows(const TrendState &state, int trendIndex, int networkIndex);

/// Writes the captions of a production cross-section temperature trend file.
void writeProductionCrossSectionTrendHeader(const TrendState &state, int trendIndex);

/// Appends the buffered rows of a production cross-section temperature trend.
void writeProductionCrossSectionTrendRows(const TrendState &state, int trendIndex);

/// Writes the captions of a service cross-section temperature trend file.
void writeGasLineCrossSectionTrendHeader(const TrendState &state, int trendIndex);

/// Appends the buffered rows of a service cross-section temperature trend.
void writeGasLineCrossSectionTrendRows(const TrendState &state, int trendIndex);

} // namespace trendoutput

#endif // SISPRODTRENDOUTPUT_H_
