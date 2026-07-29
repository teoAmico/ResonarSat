/* Geometry and resolution relations shared by the whole pipeline.
 *
 * Every quantity in this header is a small closed-form expression, and every
 * one of them has been got wrong at least once somewhere in the literature this
 * project reproduces. They live in one place, with units stated on every
 * argument and a regression test (tests/test_constants.c) that pins them to the
 * published figures, so that no later "cleanup" can quietly reintroduce a
 * factor of two. Nothing here allocates or fails. */

#ifndef RESONARSAT_GEOM_H
#define RESONARSAT_GEOM_H

#define RS_C_LIGHT 299792458.0   /* speed of light in vacuum, m/s */

/* Which of the two published acoustic-wavelength conventions to apply.
 *
 * The two primary sources for this method disagree, and the disagreement scales
 * every depth in the output by a factor of two, so it is exposed as an explicit
 * choice rather than silently resolved:
 *
 *   RS_WAVELEN_PAPER   lambda = v / (2f). Biondi & Malanga (2022) section 4
 *                      writes "a frequency of investigation set by us equal to
 *                      12,500 Hz" but then computes 6000/25000 = 0.24 m. The
 *                      factor of two is consistent with the 4*pi in the
 *                      tomographic wavenumber, and this reading reproduces the
 *                      paper's own published 0.92 m resolution.
 *
 *   RS_WAVELEN_PATENT  lambda = v / f. Patent WO2024008365A1 gives v = 6600 m/s,
 *                      f = 22000 Hz and lambda = 0.30 m, and 6600/22000 = 0.30
 *                      exactly, with no factor of two.
 *
 * The enum values are the divisor itself, so the implementation is a division
 * by (double)convention. Default to RS_WAVELEN_PAPER; whichever is used must be
 * written into the output metadata. */
typedef enum {
    RS_WAVELEN_PAPER  = 2,
    RS_WAVELEN_PATENT = 1
} rs_wavelength_convention_t;

/* Return the acoustic (seismic) wavelength in metres for an assumed material
 * wave speed 'v' in m/s and an assumed investigation frequency 'f' in Hz, under
 * the selected convention.
 *
 * Both inputs are assumptions supplied by the operator, not measurements. This
 * is the single most important fact about the depth axis of any tomogram this
 * software produces: the wavelength returned here scales depth linearly, and no
 * part of the pipeline measures either 'v' or 'f'. The vibration frequencies
 * that the micro-motion stage does measure are of order 1-15 Hz, three orders
 * of magnitude below the kHz "investigation frequency" the method assumes, and
 * must never be substituted here.
 *
 * Returns 0.0 for non-positive 'v' or 'f' rather than dividing by zero; callers
 * that accept these from a user should validate first and report a real error. */
double rs_acoustic_wavelength(double v, double f, rs_wavelength_convention_t convention);

/* Return the synthetic aperture length in metres actually observed, given the
 * platform ground velocity 'v_platform' in m/s and the target illumination time
 * 't_dwell' in seconds.
 *
 * This is the 'A' of the tomographic resolution relation. Biondi & Malanga
 * (2022) section 4 describes it in prose as "about 42,000 m" but substitutes
 * 84,000 m into the arithmetic that yields the published 0.92 m; 84 km is what
 * 7 km/s over a 12 s dwell gives, so the arithmetic is self-consistent and the
 * prose figure is A/2. Compute it here from the platform state rather than
 * inheriting either published number. */
double rs_orbital_aperture(double v_platform, double t_dwell);

/* Return the tomographic (depth) resolution in metres: lambda_ac * R / (2 * A),
 * for acoustic wavelength 'lambda_ac' in metres, slant range 'R' in metres, and
 * observed aperture 'A' in metres.
 *
 * Reproduces Biondi & Malanga (2022) section 4. Because lambda_ac carries the
 * assumed velocity and frequency, so does this: the number returned is the
 * resolution the method claims under those assumptions, not an independently
 * established one. Callers should refuse to build a depth grid finer than this
 * value, and should print it next to the (v, f) that produced it.
 *
 * Returns 0.0 if 'A' is non-positive. */
double rs_tomo_resolution(double lambda_ac, double R, double A);

/* Return the azimuth resolution in metres of an image formed from an aperture
 * of duration 't' seconds: lambda * R / (2 * v_platform * t), with 'lambda' the
 * radar wavelength in metres and 'R' the slant range in metres.
 *
 * Pass the full acquisition time for the full-resolution image, or a
 * sub-aperture duration for one sub-look. The second use is the one that
 * matters: it is the price paid for temporal sampling of vibration, and it is
 * why a wide observable vibration band and a fine image cannot be had together.
 *
 * This is Eq. 1 of Vattulainen et al. (2026), stated there for spotlight mode. */
double rs_azimuth_resolution(double lambda, double R, double v_platform, double t);

/* Return the highest vibration frequency in Hz that a sub-aperture stack can
 * represent without aliasing: n_sub / (2 * t_dwell), the Nyquist limit of a
 * series sampled every t_dwell/n_sub seconds.
 *
 * Note what does NOT appear: the radar's pulse repetition frequency. The PRF
 * bounds how finely the dwell may be subdivided before sub-apertures stop
 * containing distinct pulses, but it is not itself the sampling rate of the
 * vibration measurement. Confusing the two leads to the expectation that a kHz
 * PRF permits kHz vibration measurement, which it does not.
 *
 * With overlapping sub-apertures the effective sampling rate is higher than
 * n_sub/t_dwell; pass the effective count in that case, or use
 * rs_vibration_fmax_overlap(). */
double rs_vibration_fmax(int n_sub, double t_dwell);

/* As rs_vibration_fmax(), but for overlapping sub-aperture windows: given the
 * sub-aperture duration 't_sap' in seconds and the fractional window overlap
 * 'overlap' in [0, 1), the sampling interval is t_sap * (1 - overlap) and the
 * Nyquist limit is half its reciprocal.
 *
 * Overlap is how the published work raises the sampling rate without shortening
 * the sub-aperture and thereby coarsening azimuth resolution; reported working
 * values are 38-43 percent. Returns 0.0 for 'overlap' outside [0, 1). */
double rs_vibration_fmax_overlap(double t_sap, double overlap);

/* Return the observation ratio: sub-aperture duration 't_sap' divided by the
 * period 'period' of the highest vibration frequency of interest, both in
 * seconds and Hz-derived seconds respectively.
 *
 * Keep this below 0.5. Above it the sub-aperture integrates over more than half
 * a cycle of the motion being measured, and the displacement averages toward
 * zero within the window -- the measurement smears itself away. Published
 * successful extractions run 0.39 to 0.57, the highest of those succeeding only
 * because window overlap kept the sampling rate above Nyquist. */
double rs_observation_ratio(double t_sap, double period);

/* Return the tomographic wavenumber Kz for one sub-aperture:
 * 4 * pi * B_perp / (lambda_ac * r * sin(theta)), with 'b_perp' the orthogonal
 * baseline in metres, 'lambda_ac' the acoustic wavelength in metres, 'r' the
 * slant range in metres and 'theta' the incidence angle in radians.
 *
 * This is the conventional spatial-wavenumber part of Eq. 22 of Biondi &
 * Malanga (2022), not the complete rendered exponent. The patent PDF prints
 *
 *     exp(j * 2*pi * Kz_i * t * z)
 *
 * while defining Kz_i itself with 4*pi. The extra 't' is undefined there and
 * dimensionally incompatible if it means time in seconds as it does earlier.
 * Callers therefore use the conventional exp(j*Kz_i*z). This is a documented
 * repair of an inconsistent source, not Eq. 22 "as written".
 *
 * The second caveat is geometric: in the published method the baselines are the
 * along-orbit separations of sub-aperture phase centres (the paper's Figure 8d
 * caption reads "the tomographic acquisition geometry along the satellite
 * orbit"), and an along-track separation is not an elevation baseline. This
 * function reproduces the printed Kz relation; whether the geometry it encodes
 * carries depth information is the open question the null test exists to probe.
 *
 * Returns 0.0 if any denominator factor is zero. */
double rs_tomo_wavenumber(double b_perp, double lambda_ac, double r, double theta);

#endif /* RESONARSAT_GEOM_H */
