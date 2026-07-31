/* Can this collect support the measurement you want? Asked before processing.
 *
 * WHY THIS IS SEPARATE FROM THE PIPELINE. Focusing a full-scale collect costs
 * tens of minutes and tens of gigabytes, and every failure this project has
 * recorded looked like a complete, well-formed result until it was taken apart.
 * A configuration blind above 0.3 Hz still prints an observable band of 2.53 Hz;
 * a tracker whose correlation peak wanders the whole window still reports a
 * frequency with a confident prominence. None of that is visible in the output.
 *
 * All of it is decidable in advance, from the collect's geometry and the
 * measurement asked for, by arithmetic that costs milliseconds. This does that
 * arithmetic and says which requirements are met, which are not, and which
 * cannot be settled from the data at all.
 *
 * IT IS NOT A PREDICTION OF SUCCESS. Every check here is necessary and none is
 * sufficient. A collect that passes every one of them may still measure nothing,
 * because whether the ground moves is not a property of the file. See
 * RS_VALIDATE_GROUND_TRUTH, which always reports unknown and says why.
 *
 * The thresholds come from measurement rather than from the literature where the
 * two disagree; runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md
 * records which is which.
 */

#ifndef RESONARSAT_VALIDATE_H
#define RESONARSAT_VALIDATE_H

#include <stddef.h>

#include "resonarsat/resonarsat.h"

/* How much a finding should worry the caller.
 *
 * RS_V_FAIL means the requested measurement cannot work at this configuration,
 * for a reason that is arithmetic rather than probabilistic. RS_V_WARN means it
 * is outside the range anyone has validated but not impossible. RS_V_UNKNOWN
 * means the question cannot be answered from the collect, which is a different
 * thing from a pass and is reported separately so it cannot be mistaken for
 * one. */
typedef enum {
    RS_V_PASS = 0,
    RS_V_WARN,
    RS_V_FAIL,
    RS_V_UNKNOWN
} rs_validate_level_t;

/* The checks, in the order they are reported. */
typedef enum {
    RS_VALIDATE_METADATA = 0,     /* the fields the pipeline needs are present */
    RS_VALIDATE_FREQ_RESOLUTION,  /* dwell against the frequency asked for */
    RS_VALIDATE_BAND,             /* Nyquist of the sub-aperture series */
    RS_VALIDATE_OBSERVATION_RATIO,/* eta: does the sub-look resolve the echoes */
    RS_VALIDATE_AVERAGING_NULL,   /* integer eta, where the response vanishes */
    RS_VALIDATE_APERTURE_FRACTION,/* alpha against the validated range */
    RS_VALIDATE_GRID,             /* n_az >= 2*n_looks for the spectral route */
    RS_VALIDATE_PRF_STABILITY,    /* measured from the file's own pulse times */
    RS_VALIDATE_SENSITIVITY,      /* smallest displacement above the floor */
    RS_VALIDATE_AMBIGUITY,        /* is there any amplitude the tracker can see */
    RS_VALIDATE_COHERENCE_GATE,   /* can the mask reject anything at all */
    RS_VALIDATE_PHASE_FLOOR,      /* the phase estimator's noise floor, Eq 1.18 */
    RS_VALIDATE_MEMORY,           /* what a read of this collect will cost */
    RS_VALIDATE_GROUND_TRUTH,     /* always unknown; see the header */
    RS_VALIDATE_N_CHECKS
} rs_validate_check_t;

/* What the collect provides and what the caller intends to do with it.
 *
 * The geometry fields come from the product. 'prf_min' and 'prf_max' are
 * measured from consecutive pulse times rather than taken from the annotated
 * PRF, because the annotation is a single number and the real thing is not --
 * see rs_subap_stack_t. 'worst_gap_s' is the largest interval between
 * consecutive pulses, which detects dropped vectors.
 *
 * 'target_freq_hz' is the vibration the caller is looking for and is the input
 * that makes most of these checks answerable. 'target_amp_m' may be zero, in
 * which case the sensitivity check reports the smallest detectable amplitude
 * instead of judging a given one. */
typedef struct {
    double dwell_s, prf_hz, lambda_m, slant_range_m, v_platform_ms, incidence_rad;
    size_t n_pulse, n_rbin;
    double prf_min_hz, prf_max_hz, worst_gap_s;

    double target_freq_hz;
    double target_amp_m;
    double alpha;        /* aperture fraction, t_sap / dwell */
    double overlap;
    size_t upsample;
    double cell_m;
    size_t grid_n;       /* grid cells along track, for the spectral route */

    /* The correlation/coherence window edge, in cells, and the coherence mask
     * the run intends to apply. Both default to the pipeline's own defaults;
     * set them to check a specific configuration. */
    size_t win;
    double coherence_min;
} rs_validate_req_t;

/* One check's verdict, with the number behind it. */
typedef struct {
    rs_validate_check_t check;
    rs_validate_level_t level;
    const char *name;
    char detail[512];
} rs_validate_finding_t;

/* Fill a request with defaults: the aperture fraction and overlap the published
 * validation uses, this project's default upsampling, and no target. The caller
 * must still supply the geometry and a target frequency. */
void rs_validate_req_default(rs_validate_req_t *req);

/* Run every check and write the findings.
 *
 * 'out' must have room for RS_VALIDATE_N_CHECKS findings. Writes the count to
 * 'n_out' and returns the worst level encountered, so a caller can branch on the
 * verdict without walking the list. Returns RS_V_FAIL on a NULL argument, since
 * that is not a collect anyone should process. */
rs_validate_level_t rs_validate(const rs_validate_req_t *req,
                                rs_validate_finding_t *out, size_t *n_out);

/* The one-word name of a level, for printing. */
const char *rs_validate_level_str(rs_validate_level_t level);

#endif /* RESONARSAT_VALIDATE_H */
