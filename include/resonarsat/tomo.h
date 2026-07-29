/* Tomographic focusing: turning per-pixel vibration observations into a depth
 * profile. This is the research-grade stage, and the caveats below are part of
 * the interface rather than commentary on it.
 *
 * WHAT SETS THE DEPTH AXIS. Not anything this software measures. The depth
 * scale is fixed by an assumed material wave speed and an assumed
 * "investigation frequency" -- 12.5 kHz in Biondi & Malanga (2022), 22 kHz in
 * patent WO2024008365A1 -- and both are operator inputs with no default. The
 * vibration frequencies the micro-motion stage does measure are of order 1-15
 * Hz and must never be substituted: doing so gives acoustic wavelengths of
 * hundreds of metres and depth resolutions of kilometres.
 *
 * WHAT THE BASELINE IS. In the published method the tomographic baselines are
 * the along-orbit separations of sub-aperture phase centres (the source's
 * figure 8d caption: "the tomographic acquisition geometry along the satellite
 * orbit"). An along-track separation is not an elevation baseline. Model A
 * reproduces the published formulation faithfully and labels it as such; Model
 * C, which uses genuine perpendicular baselines across passes, is the
 * uncontested reference to measure it against. Whether the along-track geometry
 * carries depth information is the method's central open question, and no
 * implementation choice here can settle it. */

#ifndef RESONARSAT_TOMO_H
#define RESONARSAT_TOMO_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/geom.h"
#include "resonarsat/microm.h"
#include "resonarsat/resonarsat.h"

/* Which focusing model to apply. */
typedef enum {
    /* Steering-matrix inversion, Biondi & Malanga (2022) Eqs. 21-24, over the
     * sub-aperture index of a single pass. */
    RS_TOMO_MODEL_A = 0,
    /* Standing-wave back-projection: map each spectral peak to a wavelength and
     * accumulate energy at its harmonics. A simpler heuristic that fails
     * differently from Model A, which is what makes it a useful cross-check. */
    RS_TOMO_MODEL_B = 1,
    /* Classic multi-baseline SAR tomography over genuine perpendicular
     * baselines and the electromagnetic wavelength. Established physics. */
    RS_TOMO_MODEL_C = 2,
    /* Resonance mapping: read depth from the MEASURED vibration frequency
     * rather than from a baseline, via z = v / (2 f).
     *
     * This is a materially different reading of the source, and the reason to
     * implement it is that it explains the published depth scale where Model A
     * does not. Model A's depth axis is linear in the transform's bin index,
     * which is the SAR-tomography convention -- the conjugate of a baseline is a
     * height. It reaches tens of metres. The source claims hundreds of metres to
     * kilometres, and the mass-and-spring diagram in the presentation material
     * describes something else entirely: a layer at depth d ringing at
     * f = v/(2d), so that depth is INVERSE in frequency.
     *
     * With the paper's 6000 m/s over a 20 s dwell and 192 sub-looks, the
     * measurable band of 0.05-4.8 Hz maps to 625 m - 60 km. The shallow end
     * lands within 4 % of the 648 m claimed at Giza, which is difficult to read
     * as coincidence.
     *
     * What follows from it matters, and cuts against several findings recorded
     * in this project. Depth here does not come from baselines, so the objection
     * that an along-track aperture provides no elevation diversity does not
     * apply to it. It needs only enough sub-looks to resolve a frequency rather
     * than to synthesise a deep baseline, so the depth-versus-sensitivity trade
     * does not apply either. And the frequency is MEASURED rather than assumed,
     * leaving the velocity as the only free constant -- and a velocity taken
     * from known lithology is a defensible assumption, not an arbitrary one.
     *
     * Implemented so the question can be settled by measurement. */
    RS_TOMO_MODEL_D = 3
} rs_tomo_model_t;

/* How Model A performs the inversion.
 *
 *   RS_TOMO_SOLVER_DFT     Treat the steering matrix as a Fourier operator, as
 *                          both sources say it is: the paper describes A as
 *                          "the best approximation of a matrix operator
 *                          performing the digital Fourier transform of Y", and
 *                          the patent names the focusing block "the DFT
 *                          mathematical operator". Cost is O(k log k) per pixel
 *                          rather than O(k*F^2).
 *
 *   RS_TOMO_SOLVER_LSTSQ   Build the steering matrix explicitly and solve the
 *                          regularised least-squares problem. Slower, and kept
 *                          as a cross-check: on a uniform depth grid the two
 *                          must agree to numerical precision, and a
 *                          disagreement means the Kz mapping is non-uniform and
 *                          the DFT shortcut is invalid for that geometry.
 */
typedef enum {
    RS_TOMO_SOLVER_DFT = 0,
    RS_TOMO_SOLVER_LSTSQ = 1
} rs_tomo_solver_t;

/* How the data vector Y of Eq. 21 is assembled from the tracking output.
 *
 *   RS_TOMO_Y_SHIFTS   Complex, from both coregistrator shift components:
 *                      Y_k = dx_k + j*dy_k, in metres. This is the paper's own
 *                      construction. Eq. 20 gives the motion as a two-component
 *                      vector r(t) = [a*cos(w0*t), b*sin(w0*t)]*exp(-lambda*t/2)
 *                      "where {a, b} are the instantaneous shifts estimated by
 *                      the coregistrator", and Eq. 21 declares Y to be in C^k.
 *                      The two shifts are the real and imaginary parts.
 *
 *   RS_TOMO_Y_LOS      Real, from the line-of-sight displacement alone.
 *
 * The difference is not cosmetic. A real-valued Y has a Hermitian-symmetric
 * transform, so every depth appears together with its negative and the sign of
 * a recovered depth is undetermined -- half the inversion's output is a mirror
 * of the other half. The complex construction resolves that, which is presumably
 * why the source specifies it. Default to RS_TOMO_Y_SHIFTS. */
typedef enum {
    RS_TOMO_Y_SHIFTS = 0,
    RS_TOMO_Y_LOS = 1
} rs_tomo_y_source_t;

/* Focusing parameters. 'velocity' and 'frequency' have no defaults by design;
 * rs_tomo_params_check() rejects a parameter set in which they are unset. */
typedef struct {
    rs_tomo_model_t  model;
    rs_tomo_solver_t solver;
    rs_tomo_y_source_t y_source;
    rs_wavelength_convention_t convention;

    double velocity;    /* assumed material wave speed, m/s -- REQUIRED */
    double frequency;   /* assumed investigation frequency, Hz -- REQUIRED */

    double depth_max;   /* deepest cell, m */
    double depth_cell;  /* cell size, m */

    double slant_range; /* R, m */
    double incidence;   /* theta, rad */

    /* The baseline extent the model inverts over, in metres.
     *
     * Its meaning differs by model, and the difference is the whole dispute.
     * For Models A and B it is the along-orbit arc flown during the dwell --
     * tens of kilometres. For Model C it is the spread of genuine perpendicular
     * baselines across a repeat-pass stack -- hundreds of metres to a couple of
     * kilometres. Resolution improves with this quantity, so treating an
     * along-track arc as though it were an elevation baseline buys an
     * improvement of one to two orders of magnitude on paper. */
    double aperture;

    /* Radar wavelength in metres, used by Model C. Model C makes no acoustic
     * substitution: it is conventional tomography over real baselines and uses
     * the electromagnetic wavelength throughout. Zero selects an X-band
     * default. */
    double radar_wavelength;

    double regularisation; /* LSTSQ only: Tikhonov factor relative to trace */

    /* Experimental literal scaling of the rendered Eq. 22 exponent.
     *
     * Zero (default) uses the conventional TomoSAR steering term
     * exp(j*Kz*z), where Kz is in rad/m. A positive value supplies a fixed
     * scalar t to the patent's printed exp(j*2*pi*Kz*t*z). Earlier source text
     * defines t as acquisition time but Eq. 22 neither indexes it nor gives a
     * sampling rule. This option is therefore one testable interpretation, not
     * a uniquely literal time model. It is never inferred from acquisition
     * time, sub-aperture dt, or depth spacing.
     *
     * Supported only by Model A with the LSTSQ solver; the DFT shortcut assumes
     * conventional uniform spatial-wavenumber sampling. */
    double eq22_literal_t;

    /* Two conditioning steps this project added, switchable so that a run can
     * be made without anything the source does not specify.
     *
     * 'window' applies a Hann taper along the sub-aperture index before the
     * depth transform. Without it the profile carries the sidelobes of a
     * rectangular k-point transform, which appear as evenly spaced layers under
     * every bright scatterer -- the artefact most easily mistaken for
     * stratigraphy. With it, depth resolution broadens by about half again and
     * genuinely narrow peaks are attenuated. The source specifies neither.
     *
     * 'remove_y_mean' subtracts the mean of the data vector. A constant offset
     * across sub-looks is registration bias rather than motion, and a DC term
     * maps to zero depth by construction, so leaving it deposits a feature at
     * zero depth in every window. The source does not mention it.
     *
     * Both default on. Turning them off is what a faithful-to-the-source run
     * needs, and is expected to be noisier rather than better. */
    int window;          /* nonzero: Hann taper before the depth transform */
    int remove_y_mean;   /* nonzero: subtract the mean of Y */

    /* Exact-source audit fields.
     *
     * These describe stages outside rs_tomo_focus() but still needed before a
     * sidecar may claim a literal patent/paper reproduction. They do not affect
     * focusing arithmetic; the CLI fills them from the measurement chain. */
    int patent_exact;        /* nonzero: caller requested the exact source path */
    int subap_window;        /* Doppler subaperture taper used before tracking */
    double coherence_min;    /* tracking coherence mask used before inversion */

    /* Nonzero: block 7 tracked each look's slave against its own master, the
     * pair held B_shift apart -- WO2024008365A1 [0004] and its Figure 0.3.
     *
     * This belongs in the audit set for the same reason the two fields above
     * do, and it is the most consequential of the three: it decides what the
     * displacement series IS. A common master gives every sample as an offset
     * from one fixed look; the pair gives each as a difference across a fixed
     * lag. Eq. 21's Y is the second, so a run built on the first is not the
     * published method however exactly Eqs. 22-24 are then evaluated. Without
     * this field the sidecar would certify such a run as literal. */
    int pair_reference;

    /* Provenance of the measurement this cube was inverted from.
     *
     * These do not affect the arithmetic. They are recorded because a depth
     * cube is not interpretable without them and the sidecar previously carried
     * none of them: a reader could not tell whether a product came from the
     * paper's sub-aperture decomposition or a substitute, from tracked shifts
     * or pixel phase, or through what coherence gate. Every real-data result in
     * this project's documents was produced before this field existed, and
     * reconstructing which route each used cost more than recording it would
     * have.
     *
     * Set by the caller that ran the micro-motion stage; an empty string means
     * unrecorded, which the sidecar prints as such rather than guessing. */
    char provenance[256];
} rs_tomo_params_t;

/* A focused depth profile per window, plus the parameters that produced it.
 *
 * Every field of 'params' and every derived constant below must accompany the
 * data wherever it goes. A tomogram without its (v, f, lambda_ac, A, dT, model)
 * is not interpretable by anyone, including its author six months later. */
typedef struct {
    double *profile;   /* [n_win][n_depth] focused amplitude */
    double *depth;     /* [n_depth] depth axis, m */
    double *quality;   /* [n_win] carried through from the tracking mask */

    size_t n_win, n_win_az, n_win_rg, n_depth;

    double lambda_ac;      /* acoustic wavelength actually used, m */
    double dT;             /* tomographic resolution, m */
    double z_unambiguous;  /* greatest depth representable without folding, m */
    size_t n_looks;        /* sub-apertures the inversion ran over */

    rs_tomo_params_t params;
} rs_tomo_t;

/* Fill 'params' with everything except the two required assumptions, which are
 * left at zero so that rs_tomo_params_check() will reject the set until the
 * caller supplies them. */
void rs_tomo_params_default(rs_tomo_params_t *params);

/* Validate a parameter set, describing the first problem via rs_set_error().
 *
 * Enforces that 'velocity' and 'frequency' are set and physically plausible,
 * that the depth grid is positive and finite, and -- the check that matters --
 * that the requested cell size is not finer than the resolution the geometry
 * supports. A grid finer than dT does not reveal more structure; it interpolates
 * the same information onto more cells and makes the result look sharper than
 * it is. */
resonarsat_status_t rs_tomo_params_check(const rs_tomo_params_t *params);

/* Return the greatest depth in metres this geometry can represent without
 * ambiguity, given 'n_looks' sub-apertures.
 *
 * The sub-aperture baselines sample the tomographic wavenumber Kz at a spacing
 * of Kz_max/n_looks, and a sampled transform is periodic: depths beyond
 * 2*pi*n_looks/Kz_max fold back onto shallower ones exactly as a frequency
 * above Nyquist folds. This is the depth-domain equivalent of the height of
 * ambiguity in conventional SAR tomography.
 *
 * The number is smaller than intuition suggests. For a 75 km aperture at 500 km
 * range with an acoustic wavelength of 0.12 m, sixteen sub-apertures reach only
 * about 3.7 m. It scales linearly with the look count, so reaching tens of
 * metres takes many tens of looks -- which in turn costs azimuth resolution.
 * That trade is real and is not stated anywhere in the source material.
 *
 * Returns 0.0 if the geometry is degenerate. */
double rs_tomo_max_depth(const rs_tomo_params_t *params, size_t n_looks);

/* Focus vibration observations into depth profiles.
 *
 * 'm' supplies the per-window displacement series that form the data vector Y
 * of Eq. 21, and 'spec' the spectra Model B consumes. Model C additionally
 * requires 'baselines', an array of 'm->n_looks' genuine perpendicular
 * baselines in metres; pass NULL for Models A and B, which derive their
 * along-orbit baselines from the aperture and look count.
 *
 * On success '*out' owns its arrays and must be released with rs_tomo_free().
 *
 * Returns RS_ERR_ARG for an invalid parameter set (see rs_tomo_params_check())
 * or a missing input a model requires. */
resonarsat_status_t rs_tomo_focus(const rs_microm_t *m,
                                  const rs_spectrum_t *spec,
                                  const rs_tomo_params_t *params,
                                  const double *baselines,
                                  rs_tomo_t *out);

/* Release everything a tomogram owns. */
void rs_tomo_free(rs_tomo_t *t);

/* One row of a parameter-sensitivity sweep. */
typedef struct {
    double velocity;    /* assumed wave speed used, m/s */
    double frequency;   /* assumed investigation frequency used, Hz */
    double lambda_ac;   /* resulting acoustic wavelength, m */
    double dT;          /* resulting depth resolution, m */
    double peak_depth;  /* depth of the strongest feature found, m */
    double peak_value;  /* its amplitude */
} rs_tomo_sweep_row_t;

/* Re-run tomographic focusing across a grid of assumed (velocity, frequency)
 * pairs and report where the strongest feature lands in each case.
 *
 * THIS IS THE DECISIVE EXPERIMENT for the depth stage, and the reason is worth
 * stating plainly. The depth axis is scaled by two constants that nothing in
 * this pipeline measures. If the recovered features move in exact proportion to
 * v/f -- as the resolution relation predicts they must, if they are simply that
 * scaling applied to a fixed spectral feature -- then the depths carry no
 * independent information and the tomogram is a relabelled vibration spectrum.
 * If instead features stay put while the assumptions vary, something in the data
 * is pinning them, and that is a genuinely interesting result.
 *
 * The sweep multiplies 'params->velocity' and 'params->frequency' by each of
 * 'n_scale' logarithmically spaced factors between 'scale_min' and 'scale_max',
 * holding everything else fixed. Results are written to 'rows', which the caller
 * allocates with at least n_scale*n_scale entries; '*n_rows' receives the number
 * written.
 *
 * A depth grid finer than a candidate parameter set can support is skipped
 * rather than silently accepted, so some combinations may not appear.
 *
 * Publish the output of this next to the null test. Together they are the two
 * experiments that distinguish this implementation from the one it reproduces. */
resonarsat_status_t rs_tomo_sweep(const rs_microm_t *m,
                                  const rs_spectrum_t *spec,
                                  const rs_tomo_params_t *params,
                                  double scale_min, double scale_max, size_t n_scale,
                                  rs_tomo_sweep_row_t *rows, size_t *n_rows);

/* Merge sweep rows from several scene realisations, averaging the recovered
 * depth at each distinct acoustic wavelength.
 *
 * A single scene gives one peak depth per wavelength, and that number moves
 * substantially with where the scatterers happen to sit -- enough that a fit
 * through one realisation's points says little. Averaging over realisations is
 * what turns the sweep from a suggestive trend into a measurement with a spread
 * attached.
 *
 * 'rows' holds the concatenated output of several rs_tomo_sweep() calls, which
 * need not be in any particular order. Rows sharing an acoustic wavelength are
 * grouped; 'out' receives one row per distinct wavelength with 'peak_depth' and
 * 'peak_value' averaged, and 'sd_out' the sample standard deviation of the
 * depths within each group. Both arrays must hold at least 'n_rows' entries.
 * '*n_out' receives the number of distinct wavelengths found.
 *
 * Feed the merged rows to rs_tomo_sweep_summary() for a fit that reflects the
 * averaged data rather than one scene's noise. */
resonarsat_status_t rs_tomo_sweep_merge(const rs_tomo_sweep_row_t *rows, size_t n_rows,
                                        rs_tomo_sweep_row_t *out, double *sd_out,
                                        size_t *n_out);

/* Summarise how strongly a sweep's recovered depths track the assumed
 * constants, by regressing log(peak_depth) on log(lambda_ac).
 *
 * A slope of 1 with a correlation near 1 means the depths are exactly the
 * assumed scaling and carry no independent depth information -- the outcome the
 * along-track-baseline physics predicts. A slope near 0 means the features are
 * pinned by something in the data. Anything between is ambiguous and should be
 * reported as such rather than rounded to whichever conclusion is preferred.
 *
 * Writes the fitted slope and the Pearson correlation. Returns RS_ERR_ARG if
 * fewer than three usable rows are supplied. */
resonarsat_status_t rs_tomo_sweep_summary(const rs_tomo_sweep_row_t *rows, size_t n_rows,
                                          double *slope, double *correlation);

/* As rs_tomo_sweep_summary(), but restricted to wavelengths whose recovered
 * depth is REPRODUCIBLE across scene realisations.
 *
 * 'sd' holds the per-wavelength standard deviation from rs_tomo_sweep_merge().
 * A wavelength is used only if sd/depth is at most 'max_rel_sd'; '*n_used'
 * receives how many survived.
 *
 * The filter matters, and so does the fact that it is a reproducibility
 * criterion rather than a goodness-of-fit one. On a coarse depth grid the peak
 * occasionally lands a cell away from where it lands in other realisations, and
 * such a point is not a measurement of anything -- it carries the grid's
 * quantisation, not the scene's. Excluding points because they scatter across
 * repeats is defensible; excluding them because they miss the line would be
 * choosing the answer. Decide 'max_rel_sd' before looking at the fit.
 *
 * Returns RS_ERR_ARG if fewer than three wavelengths survive, in which case the
 * sweep was not run over enough realisations or the grid is too coarse to
 * support the question being asked of it. */
resonarsat_status_t rs_tomo_sweep_summary_reliable(const rs_tomo_sweep_row_t *rows,
                                                   const double *sd, size_t n_rows,
                                                   double max_rel_sd,
                                                   double *slope, double *correlation,
                                                   size_t *n_used);

/* Write the parameters and derived constants of a tomogram as text, one field
 * per line, to an open stream. Used for the sidecar file that accompanies every
 * exported product and for the CLI's summary. */
void rs_tomo_write_metadata(const rs_tomo_t *t, void *stream);

#endif /* RESONARSAT_TOMO_H */
