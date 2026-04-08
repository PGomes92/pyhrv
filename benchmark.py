#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pyHRV Benchmark: Sequential vs Parallel Computation
----------------------------------------------------

Benchmarks time domain, frequency domain, and full HRV computation
in both sequential and parallel modes using synthetic NNI data of
varying sizes.

Usage:
    python benchmark.py
"""

import time
import sys
import os
import multiprocessing
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for benchmarking

warnings.filterwarnings('ignore')

# Ensure pyhrv is importable from the repo root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pyhrv.utils
import pyhrv.time_domain as td
import pyhrv.frequency_domain as fd
import pyhrv.nonlinear as nl
from pyhrv.hrv import hrv
from pyhrv.parallel import (
    parallel_hrv,
    parallel_time_domain,
    parallel_frequency_domain,
)


def generate_synthetic_nni(n_intervals, mean_nni=800, std_nni=50, seed=42):
    """Generate a synthetic NNI series with realistic characteristics."""
    rng = np.random.RandomState(seed)
    nni = rng.normal(loc=mean_nni, scale=std_nni, size=n_intervals)
    nni = np.clip(nni, 400, 1500)  # Physiologically plausible range
    return nni


def benchmark_function(func, *args, n_runs=3, **kwargs):
    """Run a function multiple times and return the median elapsed time."""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return np.median(times), result


def benchmark_time_domain(nni, n_runs=3):
    """Benchmark sequential vs parallel time domain computation."""
    # Sequential
    seq_time, _ = benchmark_function(
        td.time_domain, nni=nni, show=False, plot=False, n_runs=n_runs
    )

    # Parallel
    par_time, _ = benchmark_function(
        parallel_time_domain, nni=nni, plot=False, n_runs=n_runs
    )

    return seq_time, par_time


def benchmark_frequency_domain(nni, n_runs=3):
    """Benchmark sequential vs parallel frequency domain computation."""
    # Sequential
    seq_time, _ = benchmark_function(
        fd.frequency_domain, nni=nni, show=False, n_runs=n_runs
    )

    # Parallel
    par_time, _ = benchmark_function(
        parallel_frequency_domain, nni=nni, show=False, n_runs=n_runs
    )

    return seq_time, par_time


def benchmark_full_hrv(nni, n_runs=3):
    """Benchmark sequential vs parallel full HRV computation."""
    # Sequential
    seq_time, _ = benchmark_function(
        hrv, nni=nni, show=False, plot_ecg=False, plot_tachogram=False,
        kwargs_time={'plot': False}, parallel=False, n_runs=n_runs
    )

    # Parallel
    par_time, _ = benchmark_function(
        hrv, nni=nni, show=False, plot_ecg=False, plot_tachogram=False,
        kwargs_time={'plot': False}, parallel=True, n_runs=n_runs
    )

    return seq_time, par_time


def print_result(label, seq_time, par_time):
    """Print a formatted benchmark result row."""
    speedup = seq_time / par_time if par_time > 0 else float('inf')
    saved = seq_time - par_time
    print(f"  {label:<30s}  {seq_time:8.3f}s  {par_time:8.3f}s  {speedup:6.2f}x  {saved:+8.3f}s")


def main():
    print("=" * 85)
    print("pyHRV Benchmark: Sequential vs Parallel Computation")
    print("=" * 85)
    print(f"CPU cores available: {multiprocessing.cpu_count()}")
    print(f"Python version: {sys.version.split()[0]}")
    print()

    # Dataset sizes to benchmark
    sizes = [500, 2000, 10000]
    n_runs = 3

    print(f"Each test runs {n_runs} times (median reported)")
    print()

    # --- Time Domain Benchmarks ---
    print("-" * 85)
    print("TIME DOMAIN")
    print("-" * 85)
    print(f"  {'Dataset':<30s}  {'Seq':>8s}  {'Par':>8s}  {'Speed':>6s}  {'Saved':>8s}")
    print(f"  {'-'*30}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}")

    for size in sizes:
        nni = generate_synthetic_nni(size)
        seq_t, par_t = benchmark_time_domain(nni, n_runs=n_runs)
        print_result(f"{size} NNI intervals", seq_t, par_t)

    print()

    # --- Frequency Domain Benchmarks ---
    print("-" * 85)
    print("FREQUENCY DOMAIN")
    print("-" * 85)
    print(f"  {'Dataset':<30s}  {'Seq':>8s}  {'Par':>8s}  {'Speed':>6s}  {'Saved':>8s}")
    print(f"  {'-'*30}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}")

    for size in sizes:
        nni = generate_synthetic_nni(size)
        seq_t, par_t = benchmark_frequency_domain(nni, n_runs=n_runs)
        print_result(f"{size} NNI intervals", seq_t, par_t)

    print()

    # --- Full HRV Benchmarks ---
    print("-" * 85)
    print("FULL HRV (time + frequency + nonlinear)")
    print("-" * 85)
    print(f"  {'Dataset':<30s}  {'Seq':>8s}  {'Par':>8s}  {'Speed':>6s}  {'Saved':>8s}")
    print(f"  {'-'*30}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}")

    for size in sizes:
        nni = generate_synthetic_nni(size)
        seq_t, par_t = benchmark_full_hrv(nni, n_runs=n_runs)
        print_result(f"{size} NNI intervals", seq_t, par_t)

    print()

    # --- Sample NNI from pyHRV files ---
    print("-" * 85)
    print("REAL DATA (pyHRV sample NNI series)")
    print("-" * 85)
    print(f"  {'Dataset':<30s}  {'Seq':>8s}  {'Par':>8s}  {'Speed':>6s}  {'Saved':>8s}")
    print(f"  {'-'*30}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}")

    for series_name in ['short', 'long']:
        nni = pyhrv.utils.load_sample_nni(series=series_name)
        seq_t, par_t = benchmark_full_hrv(nni, n_runs=n_runs)
        print_result(f"Sample '{series_name}' ({len(nni)} NNI)", seq_t, par_t)

    print()
    print("=" * 85)
    print("Benchmark complete.")
    print("=" * 85)


if __name__ == '__main__':
    main()
