/* Geometry and resolution relations. See include/resonarsat/geom.h for the
 * full contract of each function; the comments here cover implementation
 * detail only. */

#include "resonarsat/geom.h"

#include <math.h>

/* Acoustic wavelength from an assumed material wave speed and an assumed
 * investigation frequency. The convention enum doubles as the divisor, so the
 * two published readings differ only in that cast. */
double rs_acoustic_wavelength(double v, double f, rs_wavelength_convention_t convention)
{
    if (v <= 0.0 || f <= 0.0) return 0.0;
    return v / ((double)convention * f);
}

/* Aperture actually flown during the dwell. Trivial, but it exists as a named
 * function so that no call site is tempted to paste in one of the two
 * inconsistent published constants. */
double rs_orbital_aperture(double v_platform, double t_dwell)
{
    if (v_platform <= 0.0 || t_dwell <= 0.0) return 0.0;
    return v_platform * t_dwell;
}

/* Tomographic depth resolution, Biondi & Malanga (2022) section 4. */
double rs_tomo_resolution(double lambda_ac, double R, double A)
{
    if (A <= 0.0) return 0.0;
    return lambda_ac * R / (2.0 * A);
}

/* Azimuth resolution for an aperture of duration t. */
double rs_azimuth_resolution(double lambda, double R, double v_platform, double t)
{
    const double denom = 2.0 * v_platform * t;
    if (denom <= 0.0) return 0.0;
    return lambda * R / denom;
}

/* Vibration Nyquist limit for a stack of n_sub non-overlapping sub-apertures. */
double rs_vibration_fmax(int n_sub, double t_dwell)
{
    if (n_sub <= 0 || t_dwell <= 0.0) return 0.0;
    return (double)n_sub / (2.0 * t_dwell);
}

/* Vibration Nyquist limit for overlapping windows: the step between window
 * centres is t_sap*(1-overlap), which is the actual sampling interval. */
double rs_vibration_fmax_overlap(double t_sap, double overlap)
{
    if (t_sap <= 0.0 || overlap < 0.0 || overlap >= 1.0) return 0.0;
    const double dt = t_sap * (1.0 - overlap);
    if (dt <= 0.0) return 0.0;
    return 1.0 / (2.0 * dt);
}

/* Observation ratio t_sap/T, the smearing criterion. */
double rs_observation_ratio(double t_sap, double period)
{
    if (period <= 0.0) return 0.0;
    return t_sap / period;
}

/* Tomographic wavenumber, Biondi & Malanga (2022) Eq. 22. */
double rs_tomo_wavenumber(double b_perp, double lambda_ac, double r, double theta)
{
    const double s = sin(theta);
    const double denom = lambda_ac * r * s;
    if (denom == 0.0) return 0.0;
    return 4.0 * M_PI * b_perp / denom;
}
