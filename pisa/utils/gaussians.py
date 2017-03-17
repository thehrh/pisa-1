#!/usr/bin/env python
#
# Copyright 2017 J.L. Lanfranchi for the IceCube/PINGU collaboration.
"""
Multiple implementations of sum-of-gaussians for compatibility and speed
"""


from __future__ import division

from collections import Iterable, OrderedDict
from math import exp, sqrt
import threading
from time import time

import numpy as np
from scipy import stats

from pisa import (FTYPE, OMP_NUM_THREADS, NUMBA_AVAIL, NUMBA_CUDA_AVAIL,
                  numba_jit)
from pisa.utils.log import logging, set_verbosity, tprofile
from pisa.utils import gaussians_cython
if FTYPE == np.float32:
    from pisa.utils.gaussians_cython import gaussian_s as gaussian
elif FTYPE == np.float64:
    from pisa.utils.gaussians_cython import gaussian_d as gaussian
else:
    raise NotImplementedError('`gaussian` function not implemented for'
                              ' pisa.FTYPE=%s' % FTYPE)


# TODO: if the Numba CUDA functions are defined, then other CUDA (e.g. pycuda)
# code doesn't run. Need to fix this behavior. (E.g. context that gets
# destroyed?)


__all__ = ['GAUS_IMPLEMENTATIONS', 'gaussian', 'gaussians', 'test_gaussians']


GAUS_IMPLEMENTATIONS = ['singlethreaded', 'multithreaded', 'cython']
if NUMBA_CUDA_AVAIL:
    GAUS_IMPLEMENTATIONS.append('cuda')

PI = FTYPE(np.pi)
TWOPI = FTYPE(2*PI)
SQRTPI = FTYPE(sqrt(PI))
SQRT2PI = FTYPE(sqrt(TWOPI))
PISQ = FTYPE(PI**2)


def gaussians(x, mu, sigma, weights=None, implementation=None):
    """Sum of multiple Gaussian curves, normalized to have area of 1.

    Parameters
    ----------
    x : array
        Points at which to evaluate the sum of Gaussians

    mu : arrays
        Means of the Gaussians to accumulate

    sigma : array
        Standard deviations of the Gaussians to accumulate

    weights : array or None
        Weights given to each Gaussian

    implementation : None or string
        One of 'singlethreaded', 'multithreaded', or 'cuda'. Passing None, the
        function will try to determine which of the implementations is best to
        call.

    Returns
    -------
    outbuf : array
        Resulting sum of Gaussians

    Notes
    -----
    This function dynamically calls an appropriate implementation depending
    upon the problem size, the hardware available (multi-core CPU or a GPU),
    and whether the user specifies `implementation`.

    """
    # TODO: figure out which is the roughly the best function to call, if more
    # than one function is available to call (i.e., really small problems are
    # fastest in single core, really large `x` are probably best on GPU, and
    # if/where multithreaded CPU beats GPU is still up in the air
    if implementation is not None:
        implementation = implementation.strip().lower()
        if not implementation in GAUS_IMPLEMENTATIONS:
            raise ValueError('`implementation` must be one of %s'
                             % GAUS_IMPLEMENTATIONS)

    # Convert all inputs to arrays of correct datatype
    if not isinstance(x, Iterable):
        x = [x]
    if not isinstance(mu, Iterable):
        mu = [mu]
    if not isinstance(sigma, Iterable):
        sigma = [sigma]
    if weights is None:
        use_weights = False
        weights = [0]
    else:
        use_weights = True
        if not isinstance(weights, Iterable):
            weights = [weights]
    x = np.asarray(x, dtype=FTYPE)
    mu = np.asarray(mu, dtype=FTYPE)
    inv_sigma = 1/np.asarray(sigma, dtype=FTYPE)
    inv_sigma_sq = -0.5 * inv_sigma * inv_sigma
    sigma = np.asarray(sigma, dtype=FTYPE)
    weights = np.asarray(weights, dtype=FTYPE)

    n_points = len(x)
    n_gaussians = len(mu)

    if use_weights:
        norm = 1/(SQRT2PI * np.sum(weights))
    else:
        norm = 1/(SQRT2PI * n_gaussians)

    # Instantiate an empty output buffer
    outbuf = np.empty(shape=n_points, dtype=FTYPE)

    # Default to CUDA since it's generally fastest
    if implementation == 'cuda' or (implementation is None
                                    and NUMBA_CUDA_AVAIL):
        logging.trace('Using CUDA Gaussians implementation')
        _gaussians_cuda(outbuf, x, mu, sigma, weights, n_gaussians)

    # Use cython version if Numba isn't available
    elif implementation == 'cython' or (implementation is None
                                        and not NUMBA_AVAIL):
        logging.trace('Using cython Gaussians implementation')
        if FTYPE == np.float64:
            gaussians_cython.gaussians_d(
                outbuf, x, mu, inv_sigma, inv_sigma_sq, n_gaussians,
                threads=OMP_NUM_THREADS
            )
        elif FTYPE == np.float32:
            gaussians_cython.gaussians_s(
                outbuf, x, mu, inv_sigma, inv_sigma_sq, n_gaussians,
                threads=OMP_NUM_THREADS
            )

    # Use singlethreaded version if OMP_NUM_THREADS is 1
    elif implementation == 'singlethreaded' or (implementation is None and
                                                OMP_NUM_THREADS == 1):
        logging.trace('Using single-threaded Gaussians implementation')
        _gaussians_singlethreaded(outbuf, x, mu, inv_sigma, inv_sigma_sq,
                                  weights, n_gaussians, start=0, stop=n_points)

    # Use multithreaded version otherwise
    elif implementation == 'multithreaded' or OMP_NUM_THREADS > 1:
        logging.trace('Using multi-threaded Gaussians implementation')
        #_gaussians_multithreaded(outbuf, x, mu, sigma, weights, n_gaussians)
        _gaussians_multithreaded(outbuf, x, mu, inv_sigma, inv_sigma_sq,
                                 weights, n_gaussians)

    else:
        raise ValueError(
            'Unhandled value(s): OMP_NUM_THREADS="%s",'
            ' NUMBA_CUDA_AVAIL="%s", `implementation`="%s"'
            % (OMP_NUM_THREADS, NUMBA_CUDA_AVAIL, implementation)
        )

    return outbuf * norm


#def _gaussians_multithreaded(outbuf, x, mu, sigma, weights, n_gaussians):
#    """Sum of multiple guassians, optimized to be run in multiple threads. This
#    dispatches the single-kernel threaded """
#    n_points = len(x)
#    chunklen = n_points // OMP_NUM_THREADS
#    threads = []
#    start = 0
#    for i in range(OMP_NUM_THREADS):
#        stop = n_points if i == (OMP_NUM_THREADS - 1) else start + chunklen
#        thread = threading.Thread(
#            target=_gaussians_singlethreaded,
#            args=(outbuf, x, mu, sigma, weights, n_gaussians, start, stop)
#        )
#        thread.start()
#        threads.append(thread)
#        start += chunklen
#
#    for thread in threads:
#        thread.join()
#
#
#GAUS_ST_FUNCSIG = (
#    (
#        'void({f:s}[:],{f:s}[:],{f:s}[:],{f:s}[:],{f:s}[:],int64,int64,int64)'
#    ).format(f=FTYPE.__name__)
#)
#@numba_jit(GAUS_ST_FUNCSIG, nopython=True, nogil=True, cache=True)
#def _gaussians_singlethreaded(outbuf, x, mu, sigma, weights, n_gaussians,
#                              start, stop):
#    """Sum of multiple guassians, optimized to be run in a single thread"""
#    if weights[0] != 0:
#        for i in range(start, stop):
#            tot = 0.0
#            for j in range(n_gaussians):
#                xlessmu = x[i] - mu[j]
#                tot += (
#                    exp(-0.5 * (xlessmu*xlessmu) / (sigma[j]*sigma[j]))
#                    * weights[j] / sigma[j]
#                )
#            outbuf[i] = tot
#    else:
#        for i in range(start, stop):
#            tot = 0.0
#            for j in range(n_gaussians):
#                xlessmu = x[i] - mu[j]
#                tot += (
#                    exp(-0.5*(xlessmu*xlessmu)/(sigma[j]*sigma[j]))/sigma[j]
#                )
#            outbuf[i] = tot


def _gaussians_multithreaded(outbuf, x, mu, inv_sigma, inv_sigma_sq, weights,
                             n_gaussians):
    """Sum of multiple guassians, optimized to be run in multiple threads. This
    dispatches the single-kernel threaded """
    n_points = len(x)
    chunklen = n_points // OMP_NUM_THREADS
    threads = []
    start = 0
    for i in range(OMP_NUM_THREADS):
        stop = n_points if i == (OMP_NUM_THREADS - 1) else start + chunklen
        thread = threading.Thread(
            target=_gaussians_singlethreaded,
            args=(outbuf, x, mu, inv_sigma, inv_sigma_sq, weights, n_gaussians,
                  start, stop)
        )
        thread.start()
        threads.append(thread)
        start += chunklen

    for thread in threads:
        thread.join()


GAUS_ST_FUNCSIG = (
    (
        'void({f:s}[:],{f:s}[:],{f:s}[:],{f:s}[:],{f:s}[:],{f:s}[:],int64,int64,int64)'
    ).format(f=FTYPE.__name__)
)
@numba_jit(GAUS_ST_FUNCSIG, nopython=True, nogil=True, cache=True)
def _gaussians_singlethreaded(outbuf, x, mu, inv_sigma, inv_sigma_sq, weights,
                              n_gaussians, start, stop):
    """Sum of multiple guassians, optimized to be run in a single thread"""
    if weights[0] != 0:
        for i in range(start, stop):
            tot = 0.0
            for j in range(n_gaussians):
                xlessmu = x[i] - mu[j]
                tot += (
                    exp(-0.5 * (xlessmu*xlessmu) * inv_sigma_sq[j])
                    * weights[j] * inv_sigma[j]
                )
            outbuf[i] = tot
    else:
        for i in range(start, stop):
            tot = 0.0
            for j in range(n_gaussians):
                xlessmu = x[i] - mu[j]
                tot += (
                    exp((xlessmu*xlessmu)*inv_sigma_sq[j])*inv_sigma[j]
                )
            outbuf[i] = tot


if NUMBA_CUDA_AVAIL:
    from numba import cuda

    def _gaussians_cuda(outbuf, x, mu, sigma, weights, n_gaussians):
        n_points = len(x)

        use_weights = True
        if weights[0] == 0:
            use_weights = False

        threads_per_block = 32
        blocks_per_grid = (
            (n_points + (threads_per_block - 1)) // threads_per_block
        )

        if use_weights:
            func = _gaussians_weighted_cuda_kernel[blocks_per_grid,
                                                   threads_per_block]
        else:
            func = _gaussians_cuda_kernel[blocks_per_grid, threads_per_block]

        # Create empty array on GPU
        d_outbuf = cuda.device_array(shape=n_points, dtype=FTYPE, stream=0)

        # Copy other arguments to GPU
        d_x = cuda.to_device(x)
        d_mu = cuda.to_device(mu)
        d_sigma = cuda.to_device(sigma)
        if use_weights:
            d_weights = cuda.to_device(weights)
            func(d_outbuf, d_x, d_mu, d_sigma, d_weights, sum_weights,
                 n_gaussians)
        else:
            d_weights = None
            func(d_outbuf, d_x, d_mu, d_sigma, n_gaussians)

        # Copy contents of GPU result to host's outbuf
        d_outbuf.copy_to_host(ary=outbuf, stream=0)

        del d_x, d_mu, d_sigma, d_weights, d_outbuf


    GAUS_CUDA_FUNCSIG = (
        (
            'void({f:s}[:], {f:s}[:], {f:s}[:], {f:s}[:], int32)'
        ).format(f=FTYPE.__name__)
    )
    @cuda.jit(GAUS_CUDA_FUNCSIG, inline=True)
    def _gaussians_cuda_kernel(outbuf, x, mu, sigma, n_gaussians):
        pt_idx = cuda.grid(1)
        tot = 0.0
        for g_idx in range(n_gaussians):
            s = sigma[g_idx]
            m = mu[g_idx]
            xlessmu = x[pt_idx] - m
            tot += exp(-0.5 * (xlessmu*xlessmu) / (s*s)) / s
        outbuf[pt_idx] = tot

    GAUS_WTD_CUDA_FUNCSIG = (
        (
            'void({f:s}[:], {f:s}[:], {f:s}[:], {f:s}[:], {f:s}[:], int32)'
        ).format(f=FTYPE.__name__)
    )
    @cuda.jit(GAUS_WTD_CUDA_FUNCSIG, inline=True)
    def _gaussians_weighted_cuda_kernel(outbuf, x, mu, sigma, weights,
                                        n_gaussians):
        pt_idx = cuda.grid(1)
        tot = 0.0
        for g_idx in range(n_gaussians):
            s = sigma[g_idx]
            m = mu[g_idx]
            xlessmu = x[pt_idx] - m
            tot += exp(-0.5 * (xlessmu*xlessmu) / (s*s)) / s
        outbuf[pt_idx] = tot


def test_gaussians():
    """Test `gaussians` function"""
    n_gaus = [1, 10, 100, 1000, 10000]
    n_eval = 1e4

    x = np.linspace(-20, 20, n_eval)
    mu_sigma_sets = [(np.linspace(-50, 50, n), np.linspace(0.5, 100, n))
                     for n in n_gaus]

    timings = OrderedDict()
    for impl in GAUS_IMPLEMENTATIONS:
        timings[impl] = []

    for mus, sigmas in mu_sigma_sets:
        if not isinstance(mus, Iterable):
            mus = [mus]
            sigmas = [sigmas]
        ref = np.sum(
            [stats.norm.pdf(x, loc=m, scale=s) for m, s in zip(mus, sigmas)],
            axis=0
        )/len(mus)
        for impl in GAUS_IMPLEMENTATIONS:
            t0 = time()
            test = gaussians(x, mu=mus, sigma=sigmas, implementation=impl)
            dt = time() - t0
            timings[impl].append(np.round(dt*1000, decimals=3))
            if not np.allclose(test, ref):
                logging.error('BAD RESULT, implementation: %s', impl)
                logging.error('max abs fract diff: %s',
                              np.max(np.abs((test/ref - 1))))

    tprofile.debug('gaussians() timings  (Note: OMP_NUM_THREADS=%d; evaluated'
                   ' at %e points)', OMP_NUM_THREADS, n_eval)
    timings_str = '  '.join([format(t, '10d') for t in n_gaus])
    tprofile.debug('         %15s       %s', 'impl.', timings_str)
    timings_str = '  '.join(['-'*10 for t in n_gaus])
    tprofile.debug('         %15s       %s', '-'*15, timings_str)
    for impl in GAUS_IMPLEMENTATIONS:
        timings_str = '  '.join([format(t, '10.3f') for t in timings[impl]])
        tprofile.debug('Timings, %15s (ms): %s', impl, timings_str)
    logging.info('<< PASS : test_gaussians >>')


if __name__ == "__main__":
    set_verbosity(2)
    test_gaussians()
