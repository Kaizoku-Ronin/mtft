"""estimator_standards.py — A.7 discipline for log-log slope fits (mtft repo).

Guards against the two artifacts identified July 2026 (L1/L2 draft, Appendix A.7):
1. Terminal-bin leverage: a grid endpoint on a bin boundary creates an
   under-populated bin whose phase-lottery RMS has maximal regression leverage
   (observed bias up to 0.035 in slope).
2. Stride-oscillation near-resonance: window strides near an integer multiple
   of the oscillation period in decades (2*pi / (gamma * ln 10)) repeat the
   endpoint phase, converting jitter into apparent systematic drift
   (gamma_3 case: 54.994 cycles per 6-decade stride, 0.01% off integer).

Usage: replace ad-hoc binned fits (e.g. falsify.envelope_slope) with
binned_log_slope(); call stride_resonance_check() before any sliding-window
sweep whose signal contains a known frequency.
"""
import math

def binned_log_slope(ys, vals, bin_width=0.5, min_bin=10):
    """Slope of log10(bin-RMS |vals|) vs log10(y). Drops bins with < min_bin
    samples (terminal-bin guard). Returns (slope, n_bins_used, n_bins_dropped).
    ys, vals: sequences of floats/mpf; vals may be signed."""
    bins = {}
    for y, v in zip(ys, vals):
        b = math.floor(math.log10(float(y)) / bin_width) * bin_width
        bins.setdefault(b, []).append(abs(float(v)))
    used, dropped = [], 0
    for b, vs in bins.items():
        if len(vs) >= min_bin:
            rms = math.sqrt(sum(x * x for x in vs) / len(vs))
            used.append((b, math.log10(rms)))
        else:
            dropped += 1
    if len(used) < 3:
        raise ValueError(f"only {len(used)} usable bins (need >= 3); "
                         f"densify the grid or widen the window")
    used.sort()
    n = len(used)
    sx = sum(p[0] for p in used); sy = sum(p[1] for p in used)
    sxx = sum(p[0] * p[0] for p in used); sxy = sum(p[0] * p[1] for p in used)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return slope, n, dropped

def stride_resonance_check(gamma, stride_decades, tol=0.05):
    """Warn if a sliding-window stride is near-resonant with an oscillation of
    frequency gamma (in ln X). Returns (cycles_per_stride, fractional_distance,
    resonant: bool). Resonant when the fractional distance to the nearest
    integer is < tol."""
    cycles = gamma * stride_decades * math.log(10) / (2 * math.pi)
    frac = abs(cycles - round(cycles))
    return cycles, frac, frac < tol

def recommended_samples_per_decade(gamma, per_period=10):
    """Grid density giving >= per_period samples per oscillation period."""
    periods_per_decade = gamma * math.log(10) / (2 * math.pi)
    return int(math.ceil(per_period * periods_per_decade))

if __name__ == "__main__":
    # regression anchors from the July 2026 sessions
    c, f, r = stride_resonance_check(25.0109, 6.0)
    print(f"gamma_3 stride-6 check: {c:.3f} cycles, frac {f:.4f}, resonant={r}")
    assert r, "gamma_3 / 6-decade near-resonance must be flagged"
    print("samples/decade for gamma_3 at 10/period:",
          recommended_samples_per_decade(25.0109))
