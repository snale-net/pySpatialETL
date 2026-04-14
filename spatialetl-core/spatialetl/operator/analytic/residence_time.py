#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from __future__ import division, print_function, absolute_import

import numpy as np


def compute_residence_time(sourceAxis, data):
    """
    """

    ## ORIGINAL ##

    # # --- Dimensions ---
    # Nt, Ny, Nx = data.shape
    # Npix = Nx * Ny  # total number of pixels
    #
    # # --- Reshape to 2D: (Nt,Npix) ---
    # Bios2D = data.reshape(Nt,Npix)
    #
    # # --- Initialize ---
    # TR_vec = np.full((Npix,), np.nan)
    #
    # # --- Loop over pixels ---
    # for p in range(Npix):
    #     C = Bios2D[:,p]
    #
    #     # no data
    #     if np.all(np.isnan(C)):
    #         continue
    #
    #     seuil = C[0] / np.exp(1)  # 1/e threshold
    #
    #     # first index where C <= seuil
    #     indices = np.where(C <= seuil)[0]
    #     idx = indices[0] if len(indices) > 0 else None
    #
    #     if idx is not None and idx > 0:
    #         # Linear interpolation (MATLAB interp1 equivalent)
    #         x = C[idx - 1:idx + 1]
    #         y = sourceAxis[idx - 1:idx + 1]
    #
    #         # np.interp requires increasing x → fix order if needed
    #         if x[0] > x[1]:
    #             TR_vec[p] = np.interp(seuil, x[::-1], y[::-1])
    #         else:
    #             TR_vec[p] = np.interp(seuil, x, y)
    #     else:
    #         TR_vec[p] = sourceAxis[-1]
    #
    # # --- Reshape to (Nx, Ny) ---
    # return TR_vec.reshape(Ny, Nx)

    ## OPTIMIZED ##

    # --- Dimensions ---
    Nt, Ny, Nx = data.shape
    Npix = Nx * Ny

    # --- Reshape to (Nt, Npix) ---
    C = data.reshape(Nt, Npix)

    # --- Initialize ---
    TR_vec = np.full(Npix, np.nan)

    # --- Mask valid pixels ---
    valid = ~np.all(np.isnan(C), axis=0)
    C_valid = C[:, valid]

    if C_valid.size == 0:
        return TR_vec.reshape(Ny, Nx)

    # --- Threshold ---
    seuil = C_valid[0, :] / np.exp(1)

    # --- Find crossings ---
    below = C_valid <= seuil[None, :]  # (Nt, Nvalid)

    # First True index along time axis
    idx = np.argmax(below, axis=0)

    # Detect pixels that never cross
    never_cross = ~np.any(below, axis=0)

    # Valid crossings
    valid_cross = (idx > 0) & (~never_cross)

    # --- Prepare output for valid pixels ---
    TR_valid = np.full(C_valid.shape[1], np.nan)

    # --- Interpolation ---
    i = np.where(valid_cross)[0]
    idx_i = idx[i]

    c0 = C_valid[idx_i - 1, i]
    c1 = C_valid[idx_i, i]
    t0 = sourceAxis[idx_i - 1]
    t1 = sourceAxis[idx_i]

    # Linear interpolation
    TR_valid[i] = t0 + (seuil[i] - c0) * (t1 - t0) / (c1 - c0)

    # --- Edge cases ---
    TR_valid[~valid_cross] = sourceAxis[-1]

    # --- Put back ---
    TR_vec[valid] = TR_valid

    # --- Reshape ---
    return TR_vec.reshape(Ny, Nx)

