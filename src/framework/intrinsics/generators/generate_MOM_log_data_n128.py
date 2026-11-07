#!/usr/bin/env python3
"""Generate lookup tables for MOM_log_data_n128.F90.

Divide [1, 2) into 128 equal bins. For each bin, choose a nearby representative
c and produce four entries: inv_c, log_c, c_hi, and c_lo.

The default reciprocal-first method seeks a representable reciprocal v = 1/c,
near the reciprocal of the bin midpoint, whose -log(v) is especially close to
a multiple of 2**-43.

To find candidates, start from logarithm-grid points L = q * 2**-43. Compute
exp(-L), then round it to binary64 to obtain a candidate reciprocal v. Before
rounding, -log(exp(-L)) would equal L; after rounding, -log(v) generally misses
the grid slightly. Compare that distance for promising candidates and choose
the smallest. Store the winning v as inv_c, the nearby grid point as log_c,
and split c = 1/v into c_hi and c_lo.

For speed, use NumPy extended precision to identify promising candidates,
then evaluate those with Decimal at --precision digits. The search does not
examine every representable reciprocal in the search range.

Both methods adjust c slightly away from the bin midpoint to make the table
values convenient to store:

  * Log-first chooses log(c) exactly on the logarithm grid, then defines
    c = exp(log(c)) in high precision and rounds 1/c for storage.
  * Reciprocal-first chooses a binary64 reciprocal v and defines c = 1/v in
    high precision. Thus the stored reciprocal defines the center exactly;
    the search seeks v for which -log(v) nearly lands on the logarithm grid.

This is why there is no log_c_lo table. Log-first puts the chosen logarithm
directly on the grid. Reciprocal-first leaves a tiny nonzero difference
between -log(v) and the stored grid value, but searches to make that error
small enough for the solver's accuracy target. This does not guarantee correct
rounding of the complete logarithm.

Use --strategy log-first for the alternative method. The --search-radius and
--objective options apply to that method.

The generated tables are printed to stdout as a complete Fortran module.
"""

from __future__ import annotations

import argparse
import struct
from decimal import Decimal, ROUND_HALF_EVEN, getcontext


# Number of lookup-table subdivisions.
LOG_NDIV = 128
# Bit-pattern step between adjacent binary64 exponents.
EXP_STRIDE = 2**52
# Spacing of the grid used to choose log(c).
LOGC_GRID = Decimal(2) ** -43


def float_to_bits(x: float) -> int:
    """Return the binary64 bit pattern of x as an unsigned integer."""

    return struct.unpack(">Q", struct.pack(">d", x))[0]


def bits_to_float(xb: int) -> float:
    """Return the binary64 value represented by xb."""

    return struct.unpack(">d", struct.pack(">Q", xb))[0]


def nearest_integer(x: Decimal) -> int:
    """Return the nearest integer to x, with ties to even."""

    return int(x.to_integral_value(rounding=ROUND_HALF_EVEN))


def bin_midpoint(j: int) -> Decimal:
    """Return the exact midpoint of table bin j."""

    start_bits = float_to_bits(1.0)
    step_bits = EXP_STRIDE // LOG_NDIV
    center_bits = start_bits + j * step_bits + step_bits // 2
    return Decimal.from_float(bits_to_float(center_bits))


def split_high_low(c: Decimal) -> tuple[float, float]:
    """Split c into a binary64 high part and residual low part."""

    c_hi = float(c)
    c_lo = float(c - Decimal.from_float(c_hi))
    return c_hi, c_lo


def choose_logc(c0: Decimal, search_radius: int, objective: str) -> Decimal:
    """Choose a grid-aligned logarithm near log(c0), for midpoint c0."""

    logc0 = getcontext().ln(c0)
    q0 = nearest_integer(logc0 / LOGC_GRID)

    best_q = q0
    best_score = None
    for q in range(q0 - search_radius, q0 + search_radius + 1):
        # Choose a candidate logarithm first, then derive its center.
        logc = Decimal(q) * LOGC_GRID
        c = getcontext().exp(logc)
        c_hi, c_lo = split_high_low(c)

        if objective == "nearest-log":
            score = abs(logc - logc0)
        elif objective == "nearest-c":
            score = abs(c - c0)
        elif objective == "nearest-c-hi":
            score = abs(Decimal.from_float(c_hi) - c0)
        elif objective == "smallest-c-lo":
            score = abs(Decimal.from_float(c_lo))
        else:
            raise ValueError(f"unknown objective {objective!r}")

        if best_score is None or score < best_score:
            best_score = score
            best_q = q

    return Decimal(best_q) * LOGC_GRID


def choose_reciprocal(c0: Decimal, radius: int):
    """Choose a reciprocal v whose -log(v) nearly falls on the 2**-43 grid.

    c0 is the bin midpoint. Search within radius binary64 bit-pattern steps
    of its rounded reciprocal. Return the winning v, its nearby logarithm-grid
    value, and the corresponding high-precision center c = 1/v.
    """

    import numpy as np

    if np.finfo(np.longdouble).nmant < 63:
        raise ValueError(
            "reciprocal searches need extended-precision NumPy longdouble"
        )
    # Start at the midpoint's rounded reciprocal and bound the search range.
    center_bits = float_to_bits(float(1 / c0))
    first_bits = center_bits - radius
    last_bits = center_bits + radius

    # Find logarithm-grid indices covering that range. The order reverses
    # because -log(v) decreases as the reciprocal v increases.
    qlo = nearest_integer(
        -getcontext().ln(Decimal.from_float(bits_to_float(last_bits)))
        / LOGC_GRID
    )
    qhi = nearest_integer(
        -getcontext().ln(Decimal.from_float(bits_to_float(first_bits)))
        / LOGC_GRID
    )

    # Each L on this grid defines an ideal reciprocal exp(-L). Round that
    # reciprocal to binary64; its negative logarithm will now miss L slightly.
    grid = (
        np.arange(qlo, qhi + 1, dtype=np.int64).astype(np.longdouble)
        * np.longdouble(2) ** -43
    )
    target = np.exp(-grid)
    rounded = target.astype(np.float64)

    # Relative reciprocal-rounding error estimates the logarithm-grid distance.
    # Use this inexpensive estimate to retain 256 promising grid points.
    scores = np.abs(target - rounded.astype(np.longdouble)) / target
    nbest = min(256, len(scores))
    shortlist = np.argpartition(scores, nbest - 1)[:nbest]

    # Include the rounded candidates and their neighbors, since the screening
    # calculation has limited precision. Also check the center and endpoints.
    candidates = {first_bits, center_bits, last_bits}
    for index in shortlist:
        bits = float_to_bits(float(rounded[index]))
        candidates.update(
            b
            for b in (bits - 1, bits, bits + 1)
            if first_bits <= b <= last_bits
        )

    # Recompute -log(v) in high precision and measure its distance to the
    # nearest grid point. Break ties by proximity to the starting reciprocal.
    best = None
    best_key = None
    for bits in sorted(candidates):
        v_float = bits_to_float(bits)
        v = Decimal.from_float(v_float)
        exact_logc = -getcontext().ln(v)
        logc = Decimal(nearest_integer(exact_logc / LOGC_GRID)) * LOGC_GRID
        key = (abs(exact_logc - logc), abs(bits - center_bits))
        if best_key is None or key < best_key:
            best_key = key
            best = (v_float, logc, 1 / v)
    return best


def generate_tables(
    search_radius: int,
    objective: str,
    strategy="reciprocal-first",
    reciprocal_radius=2**29,
):
    """Generate log table arrays."""

    invc = []
    logc = []
    c_hi = []
    c_lo = []

    for j in range(LOG_NDIV):
        c0 = bin_midpoint(j)
        if strategy == "reciprocal-first":
            invcj, logcj, cj = choose_reciprocal(c0, reciprocal_radius)
        else:
            logcj = choose_logc(c0, search_radius, objective)
            cj = getcontext().exp(logcj)
            invcj = float(1 / cj)
        c_hi_j, c_lo_j = split_high_low(cj)

        logc.append(float(logcj))
        invc.append(invcj)
        c_hi.append(c_hi_j)
        c_lo.append(c_lo_j)

    return invc, logc, c_hi, c_lo


def format_array(name: str, values: list[float], comment: str) -> str:
    """Return a Fortran parameter array definition."""

    lines = [comment]
    lines.append(f"real, parameter :: {name}(0:ndiv-1) = [ &")
    for i in range(0, len(values), 4):
        chunk = values[i : i + 4]
        suffix = ", &" if i + 4 < len(values) else " &"
        body = ", ".join(f"{v: .17e}" for v in chunk)
        lines.append(f"    {body}{suffix}")
    lines.append("]")
    return "\n".join(lines)


def print_module(
    invc,
    logc,
    c_hi,
    c_lo,
    objective: str,
    search_radius: int,
    strategy="reciprocal-first",
    reciprocal_radius=2**29,
):
    """Print a complete Fortran table module."""

    print("! This file is part of MOM6, the Modular Ocean Model version 6.")
    print("! See the LICENSE file for licensing information.")
    print("! SPDX-License-Identifier: Apache-2.0")
    print()
    print(
        "!> Lookup tables for the reproducible log() function, "
        "assumes real64 precision."
    )
    print("!! Generated by generators/generate_MOM_log_data_n128.py.")
    if strategy == "reciprocal-first":
        print(
            f"!! Reciprocal-first search; radius: {reciprocal_radius} "
            "binary64 steps on each side."
        )
    else:
        print(
            f"!! log(c) search objective: {objective}; "
            f"search radius: {search_radius} grid steps."
        )
    print("!! Regenerate from src/framework/intrinsics/:")
    print("!!   python3 generators/generate_MOM_log_data_n128.py \\")
    print(f"!!     --precision {getcontext().prec} --strategy {strategy} \\")
    if strategy == "reciprocal-first":
        print(
            f"!!     --reciprocal-radius {reciprocal_radius} "
            "> MOM_log_data_n128.F90"
        )
    else:
        print(
            f"!!     --search-radius {search_radius} --objective {objective} "
            "> MOM_log_data_n128.F90"
        )
    print("module MOM_log_data_n128")
    print()
    print("implicit none ; private")
    print()
    print("public :: ndiv, binbits")
    print("public :: inv_c, log_c, c_hi, c_lo")
    print()
    print("!> Number of subdivisions in the log table interval.")
    print(f"integer, parameter :: ndiv = {LOG_NDIV}")
    print()
    print("!> Number of bits selecting a table subdivision, log2(ndiv).")
    print(f"integer, parameter :: binbits = {LOG_NDIV.bit_length() - 1}")
    print()
    print(
        format_array(
            "inv_c",
            invc,
            "!> Lookup table of 1 / c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "log_c",
            logc,
            "!> Lookup table of log(c) for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "c_hi",
            c_hi,
            "!> Lookup table of high parts of c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "c_lo",
            c_lo,
            "!> Lookup table of low parts of c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print("end module MOM_log_data_n128")


def main() -> None:
    """Parse command-line options and print the generated Fortran module."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--precision", type=int, default=100, help="decimal working precision"
    )
    parser.add_argument(
        "--strategy",
        choices=["log-first", "reciprocal-first"],
        default="reciprocal-first",
        help="choose grid-aligned logarithms first, "
        "or search representable reciprocals",
    )
    parser.add_argument(
        "--reciprocal-radius",
        type=int,
        default=2**29,
        help="binary64 reciprocal steps to search on each side "
        "for reciprocal-first",
    )
    parser.add_argument(
        "--search-radius",
        type=int,
        default=0,
        help="number of 2**-43 log(c) grid points to search on each side",
    )
    parser.add_argument(
        "--objective",
        default="nearest-log",
        choices=["nearest-log", "nearest-c", "nearest-c-hi", "smallest-c-lo"],
        help="criterion for selecting among searched log(c) grid points",
    )
    args = parser.parse_args()

    if args.reciprocal_radius < 0 or args.search_radius < 0:
        parser.error("search radii must be nonnegative")
    if args.strategy == "reciprocal-first" and (
        args.search_radius != 0 or args.objective != "nearest-log"
    ):
        parser.error("--search-radius and --objective apply only to log-first")

    if LOG_NDIV <= 0 or LOG_NDIV & (LOG_NDIV - 1):
        raise SystemExit("LOG_NDIV must be a positive power of two")
    if EXP_STRIDE % LOG_NDIV != 0:
        raise SystemExit(
            "LOG_NDIV must evenly divide one binary64 exponent stride"
        )

    getcontext().prec = args.precision
    invc, logc, c_hi, c_lo = generate_tables(
        args.search_radius,
        args.objective,
        args.strategy,
        args.reciprocal_radius,
    )
    print_module(
        invc,
        logc,
        c_hi,
        c_lo,
        args.objective,
        args.search_radius,
        args.strategy,
        args.reciprocal_radius,
    )


if __name__ == "__main__":
    main()
