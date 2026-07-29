/* Regression test for the resolution constants chain.
 *
 * This test exists because the published sources for this method contain, in
 * three separate documents, arithmetic that does not match their own prose. The
 * quantities pinned here scale the depth axis of every tomogram the software
 * produces, and the acceptance tolerance of the synthetic depth test is itself
 * one of them, so a factor of two loose here invalidates the check that is
 * supposed to catch it.
 *
 * The rule this test enforces: resolution is computed from explicit parameters
 * every time, and no published value is ever inherited as a constant. */

#include "resonarsat/geom.h"
#include "rs_test.h"

int main(void)
{
    /* --------------------------------------------------------------------
     * Biondi & Malanga (2022), section 4 -- the Great Pyramid parameters.
     *
     * The paper's prose says "a frequency of investigation set by us equal to
     * 12,500 Hz" and writes lambda = v/f, but the arithmetic printed alongside
     * is 6000/25,000 = 0.24 m, i.e. it divides by 2f. Taking the arithmetic at
     * face value is what reproduces the paper's own published resolution, so
     * that is the reading RS_WAVELEN_PAPER implements.
     * -------------------------------------------------------------------- */
    RS_CASE("paper convention: acoustic wavelength at Giza");
    const double lambda_giza = rs_acoustic_wavelength(6000.0, 12500.0, RS_WAVELEN_PAPER);
    RS_CHECK_NEAR(lambda_giza, 0.24, 1e-9);

    /* The prose gives the orbital aperture as "about 42,000 m" but substitutes
     * 84,000 m into the resolution arithmetic. 84 km is what the platform
     * actually flies: 7 km/s over a 12 s dwell. The prose figure is A/2.
     *
     * Note a further inconsistency not resolved here: Table 1 of the same paper
     * gives the acquisition duration as 15 s, which would put A at 105 km. We
     * follow the 84 km that the published resolution arithmetic requires. */
    RS_CASE("orbital aperture from platform state");
    const double A_giza = rs_orbital_aperture(7000.0, 12.0);
    RS_CHECK_NEAR(A_giza, 84000.0, 1e-6);

    /* Closing the chain must reproduce the paper's published 0.92 m. */
    RS_CASE("tomographic resolution reproduces the published Giza figure");
    const double dT_giza = rs_tomo_resolution(lambda_giza, 650000.0, A_giza);
    RS_CHECK_NEAR(dT_giza, 0.9286, 0.02);

    /* --------------------------------------------------------------------
     * Patent WO2024008365A1 -- the other convention, from the same author.
     *
     * The patent gives v = 3000 m/s, f = 22 kHz and lambda = 0.136 m. Since
     * 3000/22000 = 0.13636, the patent divides by f alone. The two primary
     * sources therefore disagree by exactly a factor of two on the quantity
     * that sets the depth scale, which is why the convention is an explicit
     * runtime choice rather than a decision buried in the source.
     * -------------------------------------------------------------------- */
    RS_CASE("patent convention: acoustic wavelength");
    const double lambda_patent = rs_acoustic_wavelength(6600.0, 22000.0, RS_WAVELEN_PATENT);
    RS_CHECK_NEAR(lambda_patent, 0.30, 1e-9);

    /* And the resolution the patent computes from those numbers, which is the
     * check that the whole chain reproduces its arithmetic rather than just one
     * division: 0.30 * 650000 / (2 * 75000) = 1.30 m. */
    RS_CHECK_NEAR(rs_tomo_resolution(lambda_patent, 650000.0, 75000.0), 1.30, 1e-9);

    RS_CASE("the two conventions differ by exactly a factor of two");
    const double same_inputs_paper = rs_acoustic_wavelength(6600.0, 22000.0, RS_WAVELEN_PAPER);
    RS_CHECK_NEAR(lambda_patent / same_inputs_paper, 2.0, 1e-12);

    /* --------------------------------------------------------------------
     * The Umbra operating point this project actually targets.
     * -------------------------------------------------------------------- */
    RS_CASE("Umbra X-band operating point");
    const double lambda_umbra = rs_acoustic_wavelength(3000.0, 12500.0, RS_WAVELEN_PAPER);
    const double A_umbra = rs_orbital_aperture(7500.0, 10.0);
    const double dT_umbra = rs_tomo_resolution(lambda_umbra, 500000.0, A_umbra);
    RS_CHECK_NEAR(lambda_umbra, 0.12, 1e-9);
    RS_CHECK_NEAR(A_umbra, 75000.0, 1e-6);
    RS_CHECK_NEAR(dT_umbra, 0.4, 0.01);

    /* --------------------------------------------------------------------
     * The vibration / resolution trade. These two relations together are the
     * project's central feasibility constraint, and the numbers below are the
     * ones quoted in the implementation plan's trade table -- pinned so the
     * documentation and the code cannot drift apart.
     * -------------------------------------------------------------------- */
    RS_CASE("vibration Nyquist limit ignores the PRF");
    RS_CHECK_NEAR(rs_vibration_fmax(16, 10.0), 0.8, 1e-12);
    RS_CHECK_NEAR(rs_vibration_fmax(64, 10.0), 3.2, 1e-12);
    RS_CHECK_NEAR(rs_vibration_fmax(300, 10.0), 15.0, 1e-12);

    RS_CASE("sub-look azimuth resolution degrades as the band widens");
    const double lam_x = 0.031, R_x = 500000.0, vp_x = 7500.0, dwell_x = 10.0;
    const double az16  = rs_azimuth_resolution(lam_x, R_x, vp_x, dwell_x / 16.0);
    const double az64  = rs_azimuth_resolution(lam_x, R_x, vp_x, dwell_x / 64.0);
    const double az300 = rs_azimuth_resolution(lam_x, R_x, vp_x, dwell_x / 300.0);
    RS_CHECK_REL(az16,  1.65, 0.05);
    RS_CHECK_REL(az64,  6.62, 0.05);
    RS_CHECK_REL(az300, 31.0, 0.05);

    /* Resolution degrades as the square of the frequency reach: quadrupling
     * the sub-aperture count quadruples f_max and coarsens resolution by the
     * same factor, so their product scales as N^2. */
    RS_CASE("f_max times resolution scales as N squared");
    const double prod16 = rs_vibration_fmax(16, dwell_x) * az16;
    const double prod64 = rs_vibration_fmax(64, dwell_x) * az64;
    RS_CHECK_REL(prod64 / prod16, 16.0, 1e-6);

    /* --------------------------------------------------------------------
     * Overlap and the observation ratio, from the independently validated
     * processing literature (Vattulainen et al. 2026).
     * -------------------------------------------------------------------- */
    RS_CASE("overlap raises the sampling rate without shortening the window");
    const double fmax_no_overlap = rs_vibration_fmax_overlap(0.30, 0.0);
    const double fmax_40_overlap = rs_vibration_fmax_overlap(0.30, 0.40);
    RS_CHECK_NEAR(fmax_no_overlap, 1.0 / 0.6, 1e-12);
    RS_CHECK(fmax_40_overlap > fmax_no_overlap);

    RS_CASE("observation ratio flags windows that smear the motion away");
    /* A 0.30 s window against a 2 Hz signal (0.5 s period) gives 0.6: above the
     * 0.5 limit, so the window integrates over more than half a cycle. */
    RS_CHECK_NEAR(rs_observation_ratio(0.30, 0.5), 0.6, 1e-12);
    RS_CHECK(rs_observation_ratio(0.30, 0.5) > 0.5);
    /* The same window against a 1 Hz signal is comfortably inside the limit. */
    RS_CHECK(rs_observation_ratio(0.30, 1.0) < 0.5);

    /* --------------------------------------------------------------------
     * Degenerate inputs return zero rather than infinity or NaN, so that a
     * caller passing an unset parameter gets an obviously wrong answer rather
     * than a plausible one.
     * -------------------------------------------------------------------- */
    RS_CASE("degenerate inputs are refused, not propagated");
    RS_CHECK_NEAR(rs_acoustic_wavelength(0.0, 12500.0, RS_WAVELEN_PAPER), 0.0, 0.0);
    RS_CHECK_NEAR(rs_acoustic_wavelength(3000.0, 0.0, RS_WAVELEN_PAPER), 0.0, 0.0);
    RS_CHECK_NEAR(rs_orbital_aperture(-1.0, 10.0), 0.0, 0.0);
    RS_CHECK_NEAR(rs_tomo_resolution(0.24, 650000.0, 0.0), 0.0, 0.0);
    RS_CHECK_NEAR(rs_vibration_fmax(0, 10.0), 0.0, 0.0);
    RS_CHECK_NEAR(rs_vibration_fmax_overlap(0.3, 1.0), 0.0, 0.0);

    /* --------------------------------------------------------------------
     * The tomographic wavenumber, Eq. 22.
     * -------------------------------------------------------------------- */
    RS_CASE("tomographic wavenumber scales as expected");
    const double kz = rs_tomo_wavenumber(100.0, 0.24, 650000.0, 0.5918);
    const double kz_double_baseline = rs_tomo_wavenumber(200.0, 0.24, 650000.0, 0.5918);
    RS_CHECK(kz > 0.0);
    RS_CHECK_NEAR(kz_double_baseline / kz, 2.0, 1e-12);
    RS_CHECK_NEAR(rs_tomo_wavenumber(100.0, 0.24, 650000.0, 0.0), 0.0, 0.0);

    /* --------------------------------------------------------------------
     * Biondi (2022), Remote Sens. 14:3828 -- Vesuvius, the same method applied
     * to a volcano. Pinned here because it settles which reading of the
     * wavelength formula the author intends, and because the contrast with the
     * Giza figures is itself the finding.
     *
     * That paper states a seismic velocity of 3500 km/h ("approximately
     * 972 m/s"), "a frequency of investigation set by us equal to 200 Hz", and
     * computes lambda = v/f = 4.86 m. Unlike the Giza arithmetic above, this
     * one MATCHES the formula both papers print. So v/f is what the author
     * means, and the Giza sub-metre figure depends on a step that contradicts
     * the formula stated beside it.
     *
     * The consequence is checked below: under the author's own stated formula
     * Giza resolves at 1.86 m rather than the published 0.92 m.
     * -------------------------------------------------------------------- */
    RS_CASE("patent convention reproduces the Vesuvius wavelength");
    {
        const double v_vesuvius = 3500.0 * 1000.0 / 3600.0;   /* 3500 km/h */
        const double lambda = rs_acoustic_wavelength(v_vesuvius, 200.0,
                                                     RS_WAVELEN_PATENT);
        RS_CHECK_NEAR(lambda, 4.861, 1e-3);

        /* Resolution: 4.86 * 650,000 / (2 * 42,000), which the paper rounds to
         * "about 36 m". The aperture there is described as half the total orbit
         * length, the opposite of the Giza paper's substitution -- see above. */
        const double dz = rs_tomo_resolution(lambda, 650000.0, 42000.0);
        RS_CHECK(dz > 36.0 && dz < 39.0);
        printf("    Vesuvius: lambda %.3f m, resolution %.1f m\n", lambda, dz);
    }

    RS_CASE("the two papers' constants differ by more than an order of magnitude");
    {
        /* Same method, same author, same year. Only the assumed constants
         * differ, and they set the resolution entirely. This is the project's
         * central claim, stated by the source material rather than by us. */
        const double dz_giza = rs_tomo_resolution(
            rs_acoustic_wavelength(6000.0, 12500.0, RS_WAVELEN_PAPER),
            650000.0, 84000.0);
        const double dz_vesuvius = rs_tomo_resolution(
            rs_acoustic_wavelength(3500.0 * 1000.0 / 3600.0, 200.0, RS_WAVELEN_PATENT),
            650000.0, 42000.0);
        printf("    Giza %.2f m vs Vesuvius %.1f m: a factor of %.0f\n",
               dz_giza, dz_vesuvius, dz_vesuvius / dz_giza);
        RS_CHECK(dz_vesuvius / dz_giza > 30.0);

        /* And under the formula both papers print, Giza does not reach a metre. */
        const double dz_giza_as_written = rs_tomo_resolution(
            rs_acoustic_wavelength(6000.0, 12500.0, RS_WAVELEN_PATENT),
            650000.0, 84000.0);
        printf("    Giza under its own stated formula lambda=v/f: %.2f m, "
               "not the published %.2f m\n", dz_giza_as_written, dz_giza);
        RS_CHECK_NEAR(dz_giza_as_written, 2.0 * dz_giza, 1e-9);
    }

    RS_TEST_END();
}
