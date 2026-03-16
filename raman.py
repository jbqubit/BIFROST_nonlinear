"""
raman.py  –  Spontaneous Raman Scattering Module for BIFROST
=============================================================
Implements three Raman response-function models selectable via the module-level
``RAMAN_MODEL`` variable or a per-call ``model=`` keyword:

  'bw'        Blow–Wood (1989) two-parameter damped-oscillator model.
              Fast, analytic closed form.  Good near 13 THz; overestimates
              the tail beyond ~20 THz and misses the 15 THz shoulder.

  'tabulated' Lin & Agrawal (2006) 50-point KK-normalised cubic spline.
              Accurate over 0–25 THz; correct tail fall-off.  Requires
              raman_tabulated.py on the Python path.

  'hc'        Hollenbeck & Cantrell (2002) intermediate-broadening model.
              13 vibrational modes (Eq. 9 + Table 1 of HC2002), each a
              Gaussian-damped harmonic oscillator.  Key properties:
                • Peak at 13.2 THz ✓
                • 15.2 THz shoulder ratio ≈ 0.81 (BW: 0.84; tabulated: 0.90)
                • Tail at 25 THz ≈ 0.15 — non-zero because the model
                  includes modes at 32 and 36 THz (physical, unlike the
                  Lin-Agrawal spline which truncates at 25 THz)
                • Response function captures multi-mode interference,
                  ~160 fs phase reversal, and structure up to ~1 ps ✓
              Precomputed at import time via vectorised sine transform;
              subsequent calls use a cubic-spline interpolator.

Usage
-----
    import raman

    # Module-level default (affects all calls that omit model=)
    raman.RAMAN_MODEL = 'hc'

    # Or per-call override
    gR = raman.g_R(Omega, gamma, model='tabulated')

    # Convenience: compute noise in a channel
    N_dot = raman.spRam_noise_in_channel(
        lambda_pump=1451e-9, lambda_channel=1550e-9,
        delta_lambda=0.5e-9, P_pump=1e-3, L=25e3,
        gamma=1.3e-3, T=300.0, model='hc')

References
----------
[1] K. J. Blow & D. Wood, IEEE J. Quantum Electron. 25, 2665 (1989).
[2] G. P. Agrawal, Nonlinear Fiber Optics, 6th ed. Academic Press (2019).
[3] Q. Lin & G. P. Agrawal, Opt. Lett. 31, 3086 (2006).
[4] D. Hollenbeck & C. D. Cantrell, J. Opt. Soc. Am. B 19, 2886 (2002).
"""

from math import pi
import warnings

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson

# ── Physical constants ────────────────────────────────────────────────────────
_hbar = 1.0545718e-34   # J·s
_kB   = 1.380649e-23    # J/K
_c    = 2.99792458e8    # m/s
_c_cm = 2.99792458e10   # cm/s

# ── Module-level model selector ───────────────────────────────────────────────
#: Active Raman model.  Set to 'bw', 'tabulated', or 'hc'.
RAMAN_MODEL: str = 'bw'

VALID_MODELS = ('bw', 'tabulated', 'hc')


def set_model(model: str) -> None:
    """Set the module-level default Raman model ('bw', 'tabulated', or 'hc')."""
    global RAMAN_MODEL
    if model not in VALID_MODELS:
        raise ValueError(f"model must be one of {VALID_MODELS}, got {model!r}")
    RAMAN_MODEL = model


# ── Blow–Wood (BW) constants ──────────────────────────────────────────────────
RAMAN_TAU1_SI = 12.2e-15   # s  phonon oscillation period; sets peak ~13.2 THz
RAMAN_TAU2_SI = 32.0e-15   # s  phonon lifetime; sets linewidth
RAMAN_FR_SI   = 0.18        # fractional Raman contribution (Agrawal Table 2.1)
RAMAN_OMEGA_R = 2 * pi * 13.2e12  # rad/s  Raman peak frequency


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1: Blow–Wood model
# ═══════════════════════════════════════════════════════════════════════════════

def h_R_time(t):
    """Blow–Wood impulse response hR(t).

    Causal (= 0 for t ≤ 0).  Normalised: ∫₀^∞ hR(t) dt = 1.

    Parameters
    ----------
    t : array_like, seconds

    Returns
    -------
    ndarray  (same shape as *t*)
    """
    tau1, tau2 = RAMAN_TAU1_SI, RAMAN_TAU2_SI
    prefac = (tau1**2 + tau2**2) / (tau1 * tau2**2)
    return np.where(t > 0,
                    prefac * np.exp(-t / tau2) * np.sin(t / tau1),
                    0.0)


def h_R_freq(Omega):
    """Analytical Fourier transform of the Blow–Wood hR(t).

    Convention: h̃R(Ω) = ∫₀^∞ hR(t) e^{+iΩt} dt  (Agrawal sign convention).
    Im[h̃R(Ω)] > 0 for Ω > 0 (Stokes gain).

    Parameters
    ----------
    Omega : array_like, rad/s

    Returns
    -------
    ndarray, complex
    """
    tau1, tau2 = RAMAN_TAU1_SI, RAMAN_TAU2_SI
    num = (tau1**2 + tau2**2) / (tau1**2 * tau2**2)
    den = (1 / tau2 - 1j * np.asarray(Omega, dtype=complex))**2 + (1 / tau1)**2
    return num / den


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2: Hollenbeck–Cantrell (HC) model
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Model (Eq. 9 of Hollenbeck & Cantrell 2002):
#
#    hR(t) = Σᵢ Aᵢ exp(−γᵢ t) exp(−Γᵢ² t²/4) sin(ωᵥ,ᵢ t) Θ(t)
#
#  where Θ(t) is the unit step function.
#  Parameters: Table 1 of Hollenbeck & Cantrell (2002).
#  Note: Aᵢ = Aᵢ'/ωᵥ,ᵢ per paper footnote, so amplitudes are simply Aᵢ.
#  Unit conversions from wavenumbers (cm⁻¹):
#    ωᵥ,ᵢ [rad/s] = 2π c_cm × ν_i [cm⁻¹]
#    Γᵢ   [rad/s] = π  c_cm × Gaussian_FWHM [cm⁻¹]   (paper footnote)
#    γᵢ   [rad/s] = π  c_cm × Lorentzian_FWHM [cm⁻¹] (paper footnote)

# Table 1 — component positions (cm⁻¹)
_HC_NU_CM = np.array([
    56.25, 100.00, 231.25, 362.50, 463.00, 497.00, 611.50,
    691.67, 793.67, 835.50, 930.00, 1080.00, 1215.00
])
# Peak intensities Aᵢ (= Aᵢ'/ωᵥ,ᵢ)
_HC_A = np.array([
    1.00, 11.40, 36.67, 67.67, 74.00, 4.50, 6.80,
    4.60, 4.20, 4.50, 2.70, 3.10, 3.00
])
# Gaussian FWHM (cm⁻¹)
_HC_G_FWHM_CM = np.array([
    52.10, 110.42, 175.00, 162.50, 135.33, 24.50, 41.50,
    155.00, 59.50, 64.30, 150.00, 91.00, 160.00
])
# Lorentzian FWHM (cm⁻¹)
_HC_L_FWHM_CM = np.array([
    17.37, 38.81, 58.33, 54.17, 45.11, 8.17, 13.83,
    51.67, 19.83, 21.43, 50.00, 30.33, 53.33
])

# Convert to rad/s
_HC_OMEGA_V = 2.0 * pi * _c_cm * _HC_NU_CM       # angular mode frequencies
_HC_GAMMA   = pi * _c_cm * _HC_G_FWHM_CM          # Gaussian width parameter
_HC_GAMMA_L = pi * _c_cm * _HC_L_FWHM_CM          # Lorentzian damping rate


def _hc_h_R_unnorm(t):
    """Raw (unnormalised) HC impulse response on a time array."""
    t = np.asarray(t, dtype=float)
    h = np.zeros_like(t)
    mask = t > 0
    t_pos = t[mask]
    for A, omega_v, Gamma, gamma_L in zip(_HC_A, _HC_OMEGA_V, _HC_GAMMA, _HC_GAMMA_L):
        h[mask] += (A
                    * np.exp(-gamma_L * t_pos)
                    * np.exp(-Gamma**2 * t_pos**2 / 4.0)
                    * np.sin(omega_v * t_pos))
    return h


# ── Precompute HC spline at import time ───────────────────────────────────────
#   Time grid: 0–3 ps in 0.5 fs steps (6000 points) captures ~1 ps recurrence.
#   Frequency grid: 0–40 THz in 0.05 THz steps (800 points).
#   Im[h̃R(Ω)] = ∫₀^∞ hR(t) sin(Ωt) dt  (Agrawal +iΩt convention).

_HC_T_MAX  = 3.0e-12    # 3 ps
_HC_DT     = 0.5e-15    # 0.5 fs
_t_hc      = np.arange(0.0, _HC_T_MAX, _HC_DT)

_h_hc_raw  = _hc_h_R_unnorm(_t_hc)
_HC_NORM   = np.trapezoid(_h_hc_raw, _t_hc)
_h_hc_norm = _h_hc_raw / _HC_NORM          # normalised: ∫ hR dt = 1

# Sine transform over a positive-frequency grid; use antisymmetry for Ω < 0.
_HC_OMEGA_MAX = 2.0 * pi * 40.0e12          # 40 THz in rad/s
_HC_N_OMEGA   = 801
_omega_hc_grid = np.linspace(0.0, _HC_OMEGA_MAX, _HC_N_OMEGA)

# Vectorised trapezoid: (N_ω, N_t) outer product
_sin_mat_hc = np.sin(np.outer(_omega_hc_grid, _t_hc))   # (801, 6000)
_im_hR_hc_pos = np.trapezoid(_h_hc_norm[np.newaxis, :] * _sin_mat_hc, _t_hc, axis=1)

# Cubic spline on [0, 40 THz]; antisymmetric extension handles Ω < 0.
_HC_SPLINE = CubicSpline(_omega_hc_grid, _im_hR_hc_pos, extrapolate=False)

del _sin_mat_hc   # free ~38 MB


def im_h_R_hc(Omega):
    """Im[h̃R(Ω)] for the Hollenbeck–Cantrell (2002) 13-mode model.

    Normalised so that ∫₀^∞ hR(t) dt = 1.
    Antisymmetric: Im[h̃R(−Ω)] = −Im[h̃R(Ω)].
    Zero for |Ω|/(2π) > 40 THz.

    Parameters
    ----------
    Omega : array_like, rad/s

    Returns
    -------
    ndarray, float
    """
    Omega_arr = np.asarray(Omega, dtype=float)
    abs_Om = np.abs(Omega_arr)
    in_range = abs_Om <= _HC_OMEGA_MAX
    spline_vals = _HC_SPLINE(abs_Om)
    # NaN outside spline range → 0; clip any spline overshoot near zero.
    spline_vals = np.where(np.isfinite(spline_vals), spline_vals, 0.0)
    spline_vals = np.maximum(spline_vals, 0.0)
    val = np.where(in_range, spline_vals, 0.0)
    return np.where(Omega_arr >= 0.0, val, -val)


def h_R_time_hc(t):
    """HC normalised impulse response hR(t) in the time domain.

    Parameters
    ----------
    t : array_like, seconds

    Returns
    -------
    ndarray
    """
    raw = _hc_h_R_unnorm(np.asarray(t, dtype=float))
    return raw / _HC_NORM


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3: Unified g_R dispatcher
# ═══════════════════════════════════════════════════════════════════════════════

def g_R(Omega, gamma, fR=RAMAN_FR_SI, model=None):
    """Raman gain coefficient gR(Ω) [W⁻¹ m⁻¹].

    gR > 0 for Ω > 0 (Stokes gain); gR < 0 for Ω < 0 (anti-Stokes).

    Parameters
    ----------
    Omega : array_like, rad/s
        Angular frequency shift from the pump.
    gamma : float, W⁻¹ m⁻¹
        Fiber nonlinear coefficient.
    fR : float, optional
        Raman fraction (default 0.18).
    model : {'bw', 'tabulated', 'hc'} or None
        Which response function to use.  ``None`` reads ``RAMAN_MODEL``.

    Returns
    -------
    ndarray, W⁻¹ m⁻¹
    """
    m = model if model is not None else RAMAN_MODEL
    if m == 'bw':
        im_hR = np.imag(h_R_freq(Omega))
    elif m == 'tabulated':
        from raman_tabulated import im_h_R_tabulated
        im_hR = im_h_R_tabulated(Omega)
    elif m == 'hc':
        im_hR = im_h_R_hc(Omega)
    else:
        raise ValueError(f"Unknown model {m!r}. Choose from {VALID_MODELS}.")
    return 2.0 * gamma * fR * im_hR


def g_R_from_wavelengths(lambda_pump, lambda_signal, gamma, fR=RAMAN_FR_SI, model=None):
    """Convenience wrapper: evaluate gR at the shift between two wavelengths.

    Parameters
    ----------
    lambda_pump, lambda_signal : float, metres
    gamma : float, W⁻¹ m⁻¹
    fR : float, optional
    model : str or None

    Returns
    -------
    float, W⁻¹ m⁻¹
    """
    Omega = 2.0 * pi * _c * (1.0 / lambda_pump - 1.0 / lambda_signal)
    return float(g_R(Omega, gamma, fR=fR, model=model))


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 4: Thermal phonon occupancy
# ═══════════════════════════════════════════════════════════════════════════════

def thermal_photon_number(Omega, T):
    """Bose–Einstein mean phonon occupancy nth(Ω, T).

    Uses ``np.expm1`` for numerical stability at small arguments.
    Clamps very large arguments to 0 (nth → 0).

    Parameters
    ----------
    Omega : array_like, rad/s  (magnitude used)
    T : float, K

    Returns
    -------
    ndarray
    """
    x = _hbar * np.abs(Omega) / (_kB * T)
    return np.where(x > 500.0, 0.0, 1.0 / np.expm1(x))


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 5: Spontaneous Raman photon-rate density
# ═══════════════════════════════════════════════════════════════════════════════

def _spRam_rate_density_core(gR, P_pump, L, factor, pump_depletion=False):
    """Core rate density formula (shared by Stokes and anti-Stokes).

    Parameters
    ----------
    gR : ndarray, W⁻¹ m⁻¹
    P_pump : float, W
    L : float, m
    factor : ndarray   nth+1 for Stokes; nth for anti-Stokes
    pump_depletion : bool

    Returns
    -------
    ndarray, photons s⁻¹ (rad/s)⁻¹
    """
    if not pump_depletion:
        return P_pump * L * gR * factor / (2.0 * pi)
    else:
        x = gR * P_pump * L
        # L'Hôpital: expm1(x)/x → 1 as gR → 0, recovering P·L
        corrected_PL = np.where(
            gR > 1e-25,
            np.expm1(x) / gR,
            P_pump * L
        )
        return corrected_PL * factor / (2.0 * pi)


def spRam_photon_rate_density(Omega, omega_pump, P_pump, L, gamma, T,
                              sideband='stokes', pump_depletion=False, model=None):
    """Spectral density dṄ/dΩ of spontaneous Raman photons [photons s⁻¹ (rad/s)⁻¹].

    Undepleted-pump CW approximation (pump_depletion=False) or corrected
    exponential formula (pump_depletion=True).

    Sign convention
    ---------------
    Omega = omega_pump - omega_signal.
    • Stokes:      Omega > 0  (pump bluer than signal)  → gR(Omega) > 0
    • Anti-Stokes: Omega < 0  (pump redder than signal) → gR(Omega) < 0,
                  but the physical rate is proportional to |gR|, so we
                  evaluate gR at |Omega| for the anti-Stokes branch.

    Parameters
    ----------
    Omega : array_like, rad/s
        Signed frequency shift from pump (positive = Stokes, negative = anti-Stokes).
    omega_pump : float, rad/s
        Pump angular frequency.
    P_pump : float, W
    L : float, m
        Effective fiber length (Leff_fwd or Leff_back as appropriate).
    gamma : float, W⁻¹ m⁻¹
    T : float, K
    sideband : {'stokes', 'antistokes'}
    pump_depletion : bool
    model : str or None

    Returns
    -------
    ndarray, photons s⁻¹ (rad/s)⁻¹
    """
    # Use |Omega| for gR: gR is symmetric in magnitude, antisymmetric in sign.
    # For anti-Stokes the physical gain magnitude is gR(|Omega|) > 0.
    Omega_for_gR = np.abs(Omega) if sideband == 'antistokes' else Omega
    gR  = g_R(Omega_for_gR, gamma, model=model)
    nth = thermal_photon_number(Omega, T)   # uses |Omega| internally
    factor = (nth + 1.0) if sideband == 'stokes' else nth
    return np.maximum(
        _spRam_rate_density_core(gR, P_pump, L, factor, pump_depletion),
        0.0
    )


def spRam_noise_in_channel(lambda_pump, lambda_channel, delta_lambda,
                           P_pump, L, gamma, T,
                           sideband='stokes', pump_depletion=False,
                           n_points=1001, model=None):
    """Total spontaneous Raman photon rate into a quantum channel [photons/s].

    Integrates the spectral density over the channel bandwidth Δλ centred
    at λ_channel using Simpson's rule.

    Parameters
    ----------
    lambda_pump : float, m
    lambda_channel : float, m
    delta_lambda : float, m  (full optical bandwidth of the channel)
    P_pump : float, W
    L : float, m
    gamma : float, W⁻¹ m⁻¹
    T : float, K
    sideband : {'stokes', 'antistokes'}
    pump_depletion : bool
    n_points : int  (integration points; must be odd for Simpson)
    model : str or None

    Returns
    -------
    float, photons/s
    """
    omega_pump = 2.0 * pi * _c / lambda_pump
    omega_ch   = 2.0 * pi * _c / lambda_channel
    delta_omega = 2.0 * pi * _c * delta_lambda / lambda_channel**2

    omega_lo = omega_ch - delta_omega / 2.0
    omega_hi = omega_ch + delta_omega / 2.0
    omega_s_arr = np.linspace(omega_lo, omega_hi, n_points)
    Omega_arr   = omega_pump - omega_s_arr   # positive = Stokes

    rate_density = spRam_photon_rate_density(
        Omega_arr, omega_pump, P_pump, L, gamma, T,
        sideband=sideband, pump_depletion=pump_depletion, model=model
    )
    return float(simpson(rate_density, x=omega_s_arr))


def check_depletion_validity(P_pump, L, gamma, fR=RAMAN_FR_SI, threshold=0.1):
    """Warn if the pump-depletion correction exceeds *threshold* (default 10 %).

    Returns the correction ratio [exp(x)−1]/x where x = gR_peak · P · L.
    """
    gR_peak = g_R(RAMAN_OMEGA_R, gamma, fR=fR, model='bw')
    x = gR_peak * P_pump * L
    ratio = float(np.expm1(x) / x) if x > 1e-10 else 1.0
    if ratio - 1.0 > threshold:
        warnings.warn(
            f"Pump-depletion correction is {(ratio-1)*100:.1f}% "
            f"(x = gR·P·L = {x:.3f}). Pass pump_depletion=True.",
            UserWarning, stacklevel=2
        )
    return ratio
