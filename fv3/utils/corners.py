#!/usr/bin/env python3


def fill_4corners(q, direction, grid):
    if direction == "x":
        for k in range(q.shape[2]):
            if grid.sw_corner:
                q[grid.is_ - 2, grid.js - 1, k] = q[grid.is_ - 1, grid.js + 1, k]
                q[grid.is_ - 1, grid.js - 1, k] = q[grid.is_ - 1, grid.js, k]
            if grid.se_corner:
                q[grid.ie + 2, grid.js - 1, k] = q[grid.ie + 1, grid.js + 1, k]
                q[grid.ie + 1, grid.js - 1, k] = q[grid.ie + 1, grid.js, k]
            if grid.nw_corner:
                q[grid.is_ - 1, grid.je + 1, k] = q[grid.is_ - 1, grid.je, k]
                q[grid.is_ - 2, grid.je + 1, k] = q[grid.is_ - 1, grid.je - 1, k]
            if grid.ne_corner:
                q[grid.ie + 1, grid.je + 1, k] = q[grid.ie + 1, grid.je, k]
                q[grid.ie + 2, grid.je + 1, k] = q[grid.ie + 1, grid.je - 1, k]
    elif direction == "y":
        for k in range(q.shape[2]):
            if grid.sw_corner:
                q[grid.is_ - 1, grid.js - 1, k] = q[grid.is_, grid.js - 1, k]
                q[grid.is_ - 1, grid.js - 2, k] = q[grid.is_ + 1, grid.js - 1, k]
            if grid.se_corner:
                q[grid.ie + 1, grid.js - 1, k] = q[grid.ie, grid.js - 1, k]
                q[grid.ie + 1, grid.js - 2, k] = q[grid.ie - 1, grid.js - 1, k]
            if grid.nw_corner:
                q[grid.is_ - 1, grid.je + 1, k] = q[grid.is_, grid.je + 1, k]
                q[grid.is_ - 1, grid.je + 2, k] = q[grid.is_ + 1, grid.je + 1, k]
            if grid.ne_corner:
                q[grid.ie + 1, grid.je + 1, k] = q[grid.ie, grid.je + 1, k]
                q[grid.ie + 1, grid.je + 2, k] = q[grid.ie - 1, grid.je + 1, k]
    else:
        raise ValueError("Direction not recognized. Specify either x or y")


def fill2_4corners(q1, q2, direction, grid):
    if direction == "x":
        for k in range(q1.shape[2]):
            if grid.sw_corner:
                q1[grid.is_ - 2, grid.js - 1, k] = q1[grid.is_ - 1, grid.js + 1, k]
                q1[grid.is_ - 1, grid.js - 1, k] = q1[grid.is_ - 1, grid.js, k]
                q2[grid.is_ - 2, grid.js - 1, k] = q2[grid.is_ - 1, grid.js + 1, k]
                q2[grid.is_ - 1, grid.js - 1, k] = q2[grid.is_ - 1, grid.js, k]
            if grid.se_corner:
                q1[grid.ie + 2, grid.js - 1, k] = q1[grid.ie + 1, grid.js + 1, k]
                q1[grid.ie + 1, grid.js - 1, k] = q1[grid.ie + 1, grid.js, k]
                q2[grid.ie + 2, grid.js - 1, k] = q2[grid.ie + 1, grid.js + 1, k]
                q2[grid.ie + 1, grid.js - 1, k] = q2[grid.ie + 1, grid.js, k]
            if grid.nw_corner:
                q1[grid.is_ - 1, grid.je + 1, k] = q1[grid.is_ - 1, grid.je, k]
                q1[grid.is_ - 2, grid.je + 1, k] = q1[grid.is_ - 1, grid.je - 1, k]
                q2[grid.is_ - 1, grid.je + 1, k] = q2[grid.is_ - 1, grid.je, k]
                q2[grid.is_ - 2, grid.je + 1, k] = q2[grid.is_ - 1, grid.je - 1, k]
            if grid.ne_corner:
                q1[grid.ie + 1, grid.je + 1, k] = q1[grid.ie + 1, grid.je, k]
                q1[grid.ie + 2, grid.je + 1, k] = q1[grid.ie + 1, grid.je - 1, k]
                q2[grid.ie + 1, grid.je + 1, k] = q2[grid.ie + 1, grid.je, k]
                q2[grid.ie + 2, grid.je + 1, k] = q2[grid.ie + 1, grid.je - 1, k]
    elif direction == "y":
        for k in range(q1.shape[2]):
            if grid.sw_corner:
                q1[grid.is_ - 1, grid.js - 1, k] = q1[grid.is_, grid.js - 1, k]
                q1[grid.is_ - 1, grid.js - 2, k] = q1[grid.is_ + 1, grid.js - 1, k]
                q2[grid.is_ - 1, grid.js - 1, k] = q2[grid.is_, grid.js - 1, k]
                q2[grid.is_ - 1, grid.js - 2, k] = q2[grid.is_ + 1, grid.js - 1, k]
            if grid.se_corner:
                q1[grid.ie + 1, grid.js - 1, k] = q1[grid.ie, grid.js - 1, k]
                q1[grid.ie + 1, grid.js - 2, k] = q1[grid.ie - 1, grid.js - 1, k]
                q2[grid.ie + 1, grid.js - 1, k] = q2[grid.ie, grid.js - 1, k]
                q2[grid.ie + 1, grid.js - 2, k] = q2[grid.ie - 1, grid.js - 1, k]
            if grid.nw_corner:
                q1[grid.is_ - 1, grid.je + 1, k] = q1[grid.is_, grid.je + 1, k]
                q1[grid.is_ - 1, grid.je + 2, k] = q1[grid.is_ + 1, grid.je + 1, k]
                q2[grid.is_ - 1, grid.je + 1, k] = q2[grid.is_, grid.je + 1, k]
                q2[grid.is_ - 1, grid.je + 2, k] = q2[grid.is_ + 1, grid.je + 1, k]
            if grid.ne_corner:
                q1[grid.ie + 1, grid.je + 1, k] = q1[grid.ie, grid.je + 1, k]
                q1[grid.ie + 1, grid.je + 2, k] = q1[grid.ie - 1, grid.je + 1, k]
                q2[grid.ie + 1, grid.je + 1, k] = q2[grid.ie, grid.je + 1, k]
                q2[grid.ie + 1, grid.je + 2, k] = q2[grid.ie - 1, grid.je + 1, k]
    else:
        raise ValueError("Direction not recognized. Specify either x or y")


def copy_sw_corner(q, direction, grid, kslice):
    for j in range(grid.js - grid.halo, grid.js):
        for i in range(grid.is_ - grid.halo, grid.is_):
            if direction == "x":
                q[i, j, kslice] = q[j, grid.is_ - i + 2, kslice]
            if direction == "y":
                q[i, j, kslice] = q[grid.js - j + 2, i, kslice]


def copy_se_corner(q, direction, grid, kslice):
    for j in range(grid.js - grid.halo, grid.js):
        for i in range(grid.ie + 1, grid.ie + grid.halo + 1):
            if direction == "x":
                q[i, j, kslice] = q[grid.je + 1 - j + 2, i - grid.ie + 2, kslice]
            if direction == "y":
                q[i, j, kslice] = q[grid.je + j - 2, grid.ie + 1 - i + 2, kslice]


def copy_ne_corner(q, direction, grid, kslice):
    for j in range(grid.je + 1, grid.je + grid.halo + 1):
        for i in range(grid.ie + 1, grid.ie + grid.halo + 1):
            if direction == "x":
                q[i, j, kslice] = q[j, 2 * (grid.ie + 1) - 1 - i, kslice]
            if direction == "y":
                q[i, j, kslice] = q[2 * (grid.je + 1) - 1 - j, i, kslice]


def copy_nw_corner(q, direction, grid, kslice):
    for j in range(grid.je + 1, grid.je + grid.halo + 1):
        for i in range(grid.is_ - grid.halo, grid.is_):
            if direction == "x":
                q[i, j, kslice] = q[grid.je + 1 - j + 2, i - 2 + grid.ie, kslice]
            if direction == "y":
                q[i, j, kslice] = q[j + 2 - grid.ie, grid.je + 1 - i + 2, kslice]


# can't actually be a stencil because offsets are variable
def copy_corners(q, direction, grid, kslice=slice(0, None)):
    if grid.sw_corner:
        copy_sw_corner(q, direction, grid, kslice)
    if grid.se_corner:
        copy_se_corner(q, direction, grid, kslice)
    if grid.ne_corner:
        copy_ne_corner(q, direction, grid, kslice)
    if grid.nw_corner:
        copy_nw_corner(q, direction, grid, kslice)


# TODO these can definitely be consolidated/made simpler
def fill_sw_corner_2d_bgrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.is_ - i, grid.js - j, kslice] = q[grid.is_ - j, grid.js + i, kslice]
    if direction == "y":
        q[grid.is_ - j, grid.js - i, kslice] = q[grid.is_ + i, grid.js - j, kslice]


def fill_nw_corner_2d_bgrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.is_ - i, grid.je + 1 + j, kslice] = q[grid.is_ - j, grid.je + 1 - i, kslice]
    if direction == "y":
        q[grid.is_ - j, grid.je + 1 + i, kslice] = q[grid.is_ + i, grid.je + 1 + j, kslice]


def fill_se_corner_2d_bgrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.ie + 1 + i, grid.js - j, kslice] = q[grid.ie + 1 + j, grid.js + i, kslice]
    if direction == "y":
        q[grid.ie + 1 + j, grid.js - i, kslice] = q[grid.ie + 1 - i, grid.js - j, kslice]


def fill_ne_corner_2d_bgrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.ie + 1 + i, grid.je + 1 + j, kslice] = q[grid.ie + 1 + j, grid.je + 1 - i, kslice]
    if direction == "y":
        q[grid.ie + 1 + j, grid.je + 1 + i, kslice] = q[grid.ie + 1 - i, grid.je + 1 + j, kslice]


def fill_sw_corner_2d_agrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.is_ - i, grid.js - j, kslice] = q[grid.is_ - j, i, kslice]
    if direction == "y":
        q[grid.is_ - j, grid.js - i, kslice] = q[i, grid.js - j, kslice]


def fill_nw_corner_2d_agrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.is_ - i, grid.je + j, kslice] = q[grid.is_ - j, grid.je - i + 1, kslice]
    if direction == "y":
        q[grid.is_ - j, grid.je + i, kslice] = q[i, grid.je + j, kslice]


def fill_se_corner_2d_agrid(q, i, j, direction, grid, kslice):
    if direction == "x":
        q[grid.ie + i, grid.js - j, kslice] = q[grid.ie + j, i, kslice]
    if direction == "y":
        q[grid.ie + j, grid.js - i, kslice] = q[grid.ie - i + 1, grid.js - j, kslice]


def fill_ne_corner_2d_agrid(q, i, j, direction, grid, kslice, mysign=1.0):
    if direction == "x":
        q[grid.ie + i, grid.je + j, kslice] = q[grid.ie + j, grid.je - i + 1, kslice]
    if direction == "y":
        q[grid.ie + j, grid.je + i, kslice] = q[grid.ie - i + 1, grid.je + j, kslice]


def fill_corners_2d(q, grid, gridtype, direction="x", kstart=0, nk=None):
    if nk is None:
        nk = grid.npz - kstart
    kslice = slice(kstart, kstart + nk)
    for i in range(1, 1 + grid.halo):
        for j in range(1, 1 + grid.halo):
            if gridtype == "B":
                if grid.sw_corner:
                    fill_sw_corner_2d_bgrid(q, i, j, direction, grid, kslice)
                if grid.nw_corner:
                    fill_nw_corner_2d_bgrid(q, i, j, direction, grid, kslice)
                if grid.se_corner:
                    fill_se_corner_2d_bgrid(q, i, j, direction, grid, kslice)
                if grid.ne_corner:
                    fill_ne_corner_2d_bgrid(q, i, j, direction, grid, kslice)
            if gridtype == "A":
                if grid.sw_corner:
                    fill_sw_corner_2d_agrid(q, i, j, direction, grid, kslice)
                if grid.nw_corner:
                    fill_nw_corner_2d_agrid(q, i, j, direction, grid, kslice)
                if grid.se_corner:
                    fill_se_corner_2d_agrid(q, i, j, direction, grid, kslice)
                if grid.ne_corner:
                    fill_ne_corner_2d_agrid(q, i, j, direction, grid, kslice)


def fill_sw_corner_vector_dgrid(x, y, i, j, grid, mysign, kslice):
    x[grid.is_ - i, grid.js - j, kslice] = mysign * y[grid.is_ - j, i + 2, kslice]
    y[grid.is_ - i, grid.js - j, kslice] = mysign * x[j + 2, grid.js - i, kslice]


def fill_nw_corner_vector_dgrid(x, y, i, j, grid, kslice):
    x[grid.is_ - i, grid.je + 1 + j, kslice] = y[grid.is_ - j, grid.je + 1 - i, kslice]
    y[grid.is_ - i, grid.je + j, kslice] = x[j + 2, grid.je + 1 + i, kslice]


def fill_se_corner_vector_dgrid(x, y, i, j, grid, kslice):
    x[grid.ie + i, grid.js - j, kslice] = y[grid.ie + 1 + j, i + 2, kslice]
    y[grid.ie + 1 + i, grid.js - j, kslice] = x[grid.ie - j + 1, grid.js - i, kslice]


def fill_ne_corner_vector_dgrid(x, y, i, j, grid, mysign, kslice):
    x[grid.ie + i, grid.je + 1 + j, kslice] = mysign * y[grid.ie + 1 + j, grid.je - i + 1, kslice]
    y[grid.ie + 1 + i, grid.je + j, kslice] = mysign * x[grid.ie - j + 1, grid.je + 1 + i, kslice]


def fill_corners_dgrid(x, y, grid, vector, kstart=0, nk=None):
    if nk is None:
        nk = grid.npz - kstart
    kslice = slice(kstart, kstart + nk)
    mysign = 1.0
    if vector:
        mysign = -1.0
    for i in range(1, 1 + grid.halo):
        for j in range(1, 1 + grid.halo):
            if grid.sw_corner:
                fill_sw_corner_vector_dgrid(x, y, i, j, grid, mysign, kslice)
            if grid.nw_corner:
                fill_nw_corner_vector_dgrid(x, y, i, j, grid, kslice)
            if grid.se_corner:
                fill_se_corner_vector_dgrid(x, y, i, j, grid, kslice)
            if grid.ne_corner:
                fill_ne_corner_vector_dgrid(x, y, i, j, grid, mysign, kslice)


def corner_ke(ke, u, v, ut, vt, i, j, dt, offsets, vsign, kslice):
    dt6 = dt / 6.0
    ke[i, j, kslice] = dt6 * (
        (ut[i, j, kslice] + ut[i, j - 1, kslice]) * u[i + offsets["io1"], j, kslice]
        + (vt[i, j, kslice] + vt[i - 1, j, kslice]) * v[i, j + offsets["jo1"], kslice]
        + (ut[i, j + offsets["jo1"], kslice] + vsign * vt[i + offsets["io1"], j, kslice])
        * u[i + offsets["io2"], j, kslice]
    )


def fix_corner_ke(ke, u, v, ut, vt, dt, grid, kstart=0, nk=None):
    if nk is None:
        nk = grid.npz - kstart
    kslice = slice(kstart, kstart + nk)
    if grid.sw_corner:
        offsets = {"io1": 0, "jo1": 0, "io2": -1}
        corner_ke(ke, u, v, ut, vt, grid.is_, grid.js, dt, offsets, 1, kslice)
    if grid.se_corner:
        offsets = {"io1": -1, "jo1": 0, "io2": 0}
        corner_ke(ke, u, v, ut, vt, grid.ie + 1, grid.js, dt, offsets, -1, kslice)
    if grid.ne_corner:
        offsets = {"io1": -1, "jo1": -1, "io2": 0}
        corner_ke(ke, u, v, ut, vt, grid.ie + 1, grid.je + 1, dt, offsets, 1, kslice)
    if grid.nw_corner:
        offsets = {"io1": 0, "jo1": -1, "io2": -1}
        corner_ke(ke, u, v, ut, vt, grid.is_, grid.je + 1, dt, offsets, -1, kslice)
