import gt4py as gt
import gt4py.gtscript as gtscript
import numpy as np

import fv3core._config as spec
import fv3core.utils.gt4py_utils as utils
from fv3core.decorators import gtstencil
from fv3core.stencils.a2b_ord4 import (
    a1,
    a2,
    lagrange_interpolation_x,
    lagrange_interpolation_y,
    lagrange_x_func,
    lagrange_y_func,
)
from fv3core.utils.corners import fill_4corners_x, fill_4corners_y


sd = utils.sd
c1 = -2.0 / 14.0
c2 = 11.0 / 14.0
c3 = 5.0 / 14.0
OFFSET = 2


def grid():
    return spec.grid


# almost the same as a2b_ord4's version
@gtscript.function
def lagrange_y_func_p1(qx):
    return a2 * (qx[0, -1, 0] + qx[0, 2, 0]) + a1 * (qx + qx[0, 1, 0])


@gtscript.function
def lagrange_x_func_p1(qy):
    return a2 * (qy[-1, 0, 0] + qy[2, 0, 0]) + a1 * (qy + qy[1, 0, 0])


@gtscript.function
def avg_x(u):
    return 0.5 * (u + u[0, 1, 0])


@gtscript.function
def avg_y(v):
    return 0.5 * (v + v[1, 0, 0])


@gtstencil()
def avg_box(u: sd, v: sd, utmp: sd, vtmp: sd):
    with computation(PARALLEL), interval(...):
        utmp = avg_x(u)
        vtmp = avg_y(v)


@gtscript.function
def contravariant(u, v, cosa, rsin):
    return (u - v * cosa) * rsin


@gtstencil()
def contravariant_stencil(u: sd, v: sd, cosa: sd, rsin: sd, out: sd):
    with computation(PARALLEL), interval(...):
        out = contravariant(u, v, cosa, rsin)


# @gtstencil()
# def contravariant_components(utmp: sd, vtmp: sd, cosa_s: sd, rsin2: sd, ua: sd, va: sd):
#     with computation(PARALLEL), interval(...):
#         ua = contravariant(utmp, vtmp, cosa_s, rsin2)
#         va = contravariant(vtmp, utmp, cosa_s, rsin2)


@gtstencil()
def ut_main(utmp: sd, uc: sd, v: sd, cosa_u: sd, rsin_u: sd, ut: sd):
    with computation(PARALLEL), interval(...):
        uc = lagrange_x_func(utmp)
        ut = contravariant(uc, v, cosa_u, rsin_u)


@gtstencil()
def vt_main(vtmp: sd, vc: sd, u: sd, cosa_v: sd, rsin_v: sd, vt: sd):
    with computation(PARALLEL), interval(...):
        vc = lagrange_y_func(vtmp)
        vt = contravariant(vc, u, cosa_v, rsin_v)


@gtscript.function
def vol_conserv_cubic_interp_func_x(u):
    return c1 * u[-2, 0, 0] + c2 * u[-1, 0, 0] + c3 * u


@gtscript.function
def vol_conserv_cubic_interp_func_x_rev(u):
    return c1 * u[1, 0, 0] + c2 * u + c3 * u[-1, 0, 0]


@gtscript.function
def vol_conserv_cubic_interp_func_y(v):
    return c1 * v[0, -2, 0] + c2 * v[0, -1, 0] + c3 * v


@gtscript.function
def vol_conserv_cubic_interp_func_y_rev(v):
    return c1 * v[0, 1, 0] + c2 * v + c3 * v[0, -1, 0]


@gtstencil()
def vol_conserv_cubic_interp_x(utmp: sd, uc: sd):
    with computation(PARALLEL), interval(...):
        uc = vol_conserv_cubic_interp_func_x(utmp)


@gtstencil()
def vol_conserv_cubic_interp_x_rev(utmp: sd, uc: sd):
    with computation(PARALLEL), interval(...):
        uc = vol_conserv_cubic_interp_func_x_rev(utmp)


@gtstencil()
def vol_conserv_cubic_interp_y(vtmp: sd, vc: sd):
    with computation(PARALLEL), interval(...):
        vc = vol_conserv_cubic_interp_func_y(vtmp)


@gtstencil()
def vt_edge(vtmp: sd, vc: sd, u: sd, cosa_v: sd, rsin_v: sd, vt: sd, rev: int):
    with computation(PARALLEL), interval(...):
        vc = (
            vol_conserv_cubic_interp_func_y(vtmp)
            if rev == 0
            else vol_conserv_cubic_interp_func_y_rev(vtmp)
        )
        vt = contravariant(vc, u, cosa_v, rsin_v)


@gtstencil()
def uc_x_edge1(ut: sd, sin_sg3: sd, sin_sg1: sd, uc: sd):
    with computation(PARALLEL), interval(...):
        uc = ut * sin_sg3[-1, 0, 0] if ut > 0 else ut * sin_sg1


@gtstencil()
def vc_y_edge1(vt: sd, sin_sg4: sd, sin_sg2: sd, vc: sd):
    with computation(PARALLEL), interval(...):
        vc = vt * sin_sg4[0, -1, 0] if vt > 0 else vt * sin_sg2


# TODO make this a stencil?
def edge_interpolate4_x(ua, dxa):
    t1 = dxa[0, :, :] + dxa[1, :, :]
    t2 = dxa[2, :, :] + dxa[3, :, :]
    n1 = (t1 + dxa[1, :, :]) * ua[1, :, :] - dxa[1, :, :] * ua[0, :, :]
    n2 = (t1 + dxa[2, :, :]) * ua[2, :, :] - dxa[2, :, :] * ua[3, :, :]
    return 0.5 * (n1 / t1 + n2 / t2)


def edge_interpolate4_y(va, dxa):
    t1 = dxa[:, 0, :] + dxa[:, 1, :]
    t2 = dxa[:, 2, :] + dxa[:, 3, :]
    n1 = (t1 + dxa[:, 1, :]) * va[:, 1, :] - dxa[:, 1, :] * va[:, 0, :]
    n2 = (t1 + dxa[:, 2, :]) * va[:, 2, :] - dxa[:, 2, :] * va[:, 3, :]
    return 0.5 * (n1 / t1 + n2 / t2)


@gtstencil()
def lagrange_interpolation(u: sd, v: sd, utmp: sd, vtmp: sd):
    with computation(PARALLEL), interval(...):
        utmp = a2 * (u[0, -1, 0] + u[0, 2, 0]) + a1 * (u + u[0, 1, 0])
        vtmp = a2 * (v[-1, 0, 0] + v[2, 0, 0]) + a1 * (v + v[1, 0, 0])


@gtstencil(externals={"HALO": 3})
def d2a2c_stencil1(
    dord4: int,
    u: sd,
    v: sd,
    ua: sd,
    va: sd,
    utc: sd,
    vtc: sd,
    utmp: sd,
    vtmp: sd,
    cosa_s: sd,
    rsin2: sd,
):
    from __externals__ import HALO, namelist
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):

        assert __INLINED(namelist.grid_type < 3)

        # The order of these blocks matters, so they cannot be merged into a
        # single block since then the order is not guaranteed
        with parallel(region[:, : j_start + HALO]):
            utmp = 0.5 * (u + u[0, 1, 0])
            vtmp = 0.5 * (v + v[1, 0, 0])
        with parallel(region[:, j_end - HALO + 1 :]):
            utmp = 0.5 * (u + u[0, 1, 0])
            vtmp = 0.5 * (v + v[1, 0, 0])
        with parallel(region[: i_start + HALO, :]):
            utmp = 0.5 * (u + u[0, 1, 0])
            vtmp = 0.5 * (v + v[1, 0, 0])
        with parallel(region[i_end - HALO + 1 :, :]):
            utmp = 0.5 * (u + u[0, 1, 0])
            vtmp = 0.5 * (v + v[1, 0, 0])

        ua = contravariant(utmp, vtmp, cosa_s, rsin2)
        va = contravariant(vtmp, utmp, cosa_s, rsin2)

        # SW corner
        with parallel(region[i_start - 3, j_start - 1]):
            utmp = -vtmp[2, 3, 0]
        with parallel(region[i_start - 2, j_start - 1]):
            utmp = -vtmp[1, 2, 0]
        with parallel(region[i_start - 1, j_start - 1]):
            utmp = -vtmp[0, 1, 0]

        # SE corner
        with parallel(region[i_end + 1, j_start - 1]):
            utmp = vtmp[0, 1, 0]
        with parallel(region[i_end + 2, j_start - 1]):
            utmp = vtmp[-1, 2, 0]
        with parallel(region[i_end + 3, j_start - 1]):
            utmp = vtmp[-2, 3, 0]

        # NE corner
        with parallel(region[i_end + 1, j_end + 1]):
            utmp = -vtmp[0, -1, 0]
        with parallel(region[i_end + 2, j_end + 1]):
            utmp = -vtmp[-1, -2, 0]
        with parallel(region[i_end + 3, j_end + 1]):
            utmp = -vtmp[-2, -3, 0]

        # NW corner
        with parallel(region[i_start - 3, j_end + 1]):
            utmp = vtmp[2, -3, 0]
        with parallel(region[i_start - 2, j_end + 1]):
            utmp = vtmp[1, -2, 0]
        with parallel(region[i_start - 1, j_end + 1]):
            utmp = vtmp[0, -1, 0]

        ua = fill_4corners_x(ua, va, sw_mult=-1, se_mult=1, ne_mult=-1, nw_mult=1)


@gtstencil()
def d2a2c_stencil2(
    utmp: sd,
    uc: sd,
    utc: sd,
    ua: sd,
    v: sd,
    cosa_u: sd,
    rsin_u: sd,
    dxa: sd,
    sin_sg1: sd,
    sin_sg3: sd,
):
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):
        uc = lagrange_x_func(utmp)
        utc = contravariant(uc, v, cosa_u, rsin_u)


@gtstencil()
def d2a2c_stencil_west(
    utmp: sd,
    uc: sd,
    utc: sd,
    ua: sd,
    v: sd,
    cosa_u: sd,
    rsin_u: sd,
    dxa: sd,
    sin_sg1: sd,
    sin_sg3: sd,
):
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):
        # West
        with parallel(region[i_start - 1, :]):
            uc = vol_conserv_cubic_interp_func_x(utmp)

        with parallel(region[i_start, :]):
            t1 = dxa[-2, 0, 0] + dxa[-1, 0, 0]
            t2 = dxa[0, 0, 0] + dxa[1, 0, 0]
            n1 = (t1 + dxa[-1, 0, 0]) * ua[-1, 0, 0] - dxa[-1, 0, 0] * ua[-2, 0, 0]
            n2 = (t1 + dxa[0, 0, 0]) * ua[0, 0, 0] - dxa[0, 0, 0] * ua[1, 0, 0]
            utc = 0.5 * (n1 / t1 + n2 / t2)

        with parallel(region[i_start, :]):
            uc = utc * sin_sg3[-1, 0, 0] if utc > 0 else utc * sin_sg1

        with parallel(region[i_start + 1, :]):
            uc = vol_conserv_cubic_interp_func_x_rev(utmp)

        with parallel(region[i_start - 1, :]):
            utc = contravariant(uc, v, cosa_u, rsin_u)

        with parallel(region[i_start + 1, :]):
            utc = contravariant(uc, v, cosa_u, rsin_u)


@gtstencil()
def d2a2c_stencil_east(
    utmp: sd,
    uc: sd,
    utc: sd,
    ua: sd,
    v: sd,
    cosa_u: sd,
    rsin_u: sd,
    dxa: sd,
    sin_sg1: sd,
    sin_sg3: sd,
):
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):
        # East
        with parallel(region[i_end, :]):
            uc = vol_conserv_cubic_interp_func_x(utmp)

        with parallel(region[i_end + 1, j_start - 1 : j_end + 2]):
            t1 = dxa[-2, 0, 0] + dxa[-1, 0, 0]
            t2 = dxa[0, 0, 0] + dxa[1, 0, 0]
            n1 = (t1 + dxa[-1, 0, 0]) * ua[-1, 0, 0] - dxa[-1, 0, 0] * ua[-2, 0, 0]
            n2 = (t1 + dxa[0, 0, 0]) * ua[0, 0, 0] - dxa[0, 0, 0] * ua[1, 0, 0]
            utc = 0.5 * (n1 / t1 + n2 / t2)

        with parallel(region[i_end + 1, :]):
            uc = utc * sin_sg3[-1, 0, 0] if utc > 0 else utc * sin_sg1

        with parallel(region[i_end + 2, :]):
            uc = vol_conserv_cubic_interp_func_x_rev(utmp)

        with parallel(region[i_end, :]):
            utc = contravariant(uc, v, cosa_u, rsin_u)

        with parallel(region[i_end + 2, :]):
            utc = contravariant(uc, v, cosa_u, rsin_u)


@gtstencil(externals={"HALO": 3})
def d2a2c_stencil3(
    dord4: int,
    u: sd,
    v: sd,
    ua: sd,
    va: sd,
    utc: sd,
    vtc: sd,
    utmp: sd,
    vtmp: sd,
    cosa_s: sd,
    rsin2: sd,
):
    from __externals__ import HALO, namelist
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):

        assert __INLINED(namelist.grid_type < 3)

        # if grid.sw_corner:
        #     for j in range(-2, 1):
        #         vtmp[grid.is_ - 1, j + 2, :] = -utmp[grid.is_ - j, grid.js - 1, :]
        #     va[i1, j1 - 1, :] = -ua[i1 + 2, j1, :]
        #     va[i1, j1, :] = -ua[i1 + 1, j1, :]

        # SW corner
        with parallel(region[i_start - 1, j_start - 3]):
            vtmp = -utmp[3, 2, 0]
        with parallel(region[i_start - 1, j_start - 2]):
            vtmp = -utmp[2, 1, 0]
        with parallel(region[i_start - 1, j_start - 1]):
            vtmp = -utmp[1, 0, 0]

        # if grid.se_corner:
        #     for j in range(-2, 1):
        #         vtmp[nx, j + 2, :] = utmp[grid.ie + j, grid.js - 1, :]
        #     va[nx, j1, :] = ua[nx - 1, j1, :]
        #     va[nx, j1 - 1, :] = ua[nx - 2, j1, :]

        # SE corner
        with parallel(region[i_end + 1, j_start - 3]):
            vtmp = utmp[-3, 2, 0]
        with parallel(region[i_end + 1, j_start - 2]):
            vtmp = utmp[-2, 1, 0]
        with parallel(region[i_end + 1, j_start - 1]):
            vtmp = utmp[-1, 0, 0]

        # if grid.ne_corner:
        #     for j in range(0, 3):
        #         vtmp[nx, ny + j, :] = -utmp[grid.ie - j, ny, :]
        #     va[nx, ny, :] = -ua[nx - 1, ny, :]
        #     va[nx, ny + 1, :] = -ua[nx - 2, ny, :]

        # NE corner
        with parallel(region[i_end + 1, j_end + 1]):
            vtmp = -utmp[-1, 0, 0]
        with parallel(region[i_end + 1, j_end + 2]):
            vtmp = -utmp[-2, -1, 0]
        with parallel(region[i_end + 1, j_end + 3]):
            vtmp = -utmp[-3, -2, 0]

        # if grid.nw_corner:
        #     for j in range(0, 3):
        #         vtmp[grid.is_ - 1, ny + j, :] = utmp[j + grid.is_, ny, :]
        #     va[i1, ny, :] = ua[i1 + 1, ny, :]
        #     va[i1, ny + 1, :] = ua[i1 + 2, ny, :]

        # NW corner
        with parallel(region[i_start - 1, j_end + 1]):
            vtmp = utmp[1, 0, 0]
        with parallel(region[i_start - 1, j_end + 2]):
            vtmp = utmp[2, -1, 0]
        with parallel(region[i_start - 1, j_end + 3]):
            vtmp = utmp[3, -2, 0]

        va = fill_4corners_y(va, ua, sw_mult=-1, se_mult=1, ne_mult=-1, nw_mult=1)


@gtstencil()
def d2a2c_stencil_south(
    vtmp: sd,
    vc: sd,
    vtc: sd,
    va: sd,
    u: sd,
    cosa_v: sd,
    rsin_v: sd,
    dya: sd,
    sin_sg2: sd,
    sin_sg4: sd,
):
    from __splitters__ import i_end, i_start, j_end, j_start

    with computation(PARALLEL), interval(...):
        with parallel(region[:, j_start - 1]):
            vc = vol_conserv_cubic_interp_func_y(vtmp)
            vtc = contravariant(vc, u, cosa_v, rsin_v)

        with parallel(region[:, j_start]):
            t1 = dya[0, -2, 0] + dya[0, -1, 0]
            t2 = dya[0, 0, 0] + dya[0, 1, 0]
            n1 = (t1 + dya[0, -1, 0]) * va[0, -1, 0] - dya[0, -1, 0] * va[0, -2, 0]
            n2 = (t1 + dya[0, 0, 0]) * va[0, 0, 0] - dya[0, 0, 0] * va[0, 1, 0]
            vtc = 0.5 * (n1 / t1 + n2 / t2)

        with parallel(region[:, j_start]):
            vc = vtc * sin_sg4[0, -1, 0] if vtc > 0 else vtc * sin_sg2

        with parallel(region[:, j_start + 1]):
            vc = vol_conserv_cubic_interp_func_y_rev(vtmp)
            vtc = contravariant(vc, u, cosa_v, rsin_v)


@gtstencil()
def lagrange_interpolation_y_p1(qx: sd, qout: sd):
    with computation(PARALLEL), interval(...):
        qout = lagrange_y_func_p1(qx)


@gtstencil()
def lagrange_interpolation_x_p1(qy: sd, qout: sd):
    with computation(PARALLEL), interval(...):
        qout = lagrange_x_func_p1(qy)


def compute(dord4, uc, vc, u, v, ua, va, utc, vtc):
    grid = spec.grid
    big_number = 1e30  # 1e8 if 32 bit
    nx = grid.ie + 1  # grid.npx + 2
    ny = grid.je + 1  # grid.npy + 2
    i1 = grid.is_ - 1
    j1 = grid.js - 1
    id_ = 1 if dord4 else 0
    npt = 4 if spec.namelist.grid_type < 3 and not grid.nested else 0
    if npt > grid.nic - 1 or npt > grid.njc - 1:
        npt = 0
    utmp = utils.make_storage_from_shape(ua.shape, grid.default_origin())
    vtmp = utils.make_storage_from_shape(va.shape, grid.default_origin())
    utmp[:] = big_number
    vtmp[:] = big_number
    js1 = npt + OFFSET if grid.south_edge else grid.js - 1
    je1 = ny - npt if grid.north_edge else grid.je + 1
    is1 = npt + OFFSET if grid.west_edge else grid.isd
    ie1 = nx - npt if grid.east_edge else grid.ied

    is1 = npt + OFFSET if grid.west_edge else grid.is_ - 1
    ie1 = nx - npt if grid.east_edge else grid.ie + 1
    js1 = npt + OFFSET if grid.south_edge else grid.jsd
    je1 = ny - npt if grid.north_edge else grid.jed

    lagrange_interpolation(
        u,
        v,
        utmp,
        vtmp,
        origin=(is1, js1, 0),
        domain=(ie1 - is1 + 1, je1 - js1 + 1, grid.npz),
    )

    # # tmp edges
    # if grid.south_edge:
    #     avg_box(
    #         u,
    #         v,
    #         utmp,
    #         vtmp,
    #         origin=(grid.isd, grid.jsd, 0),
    #         domain=(grid.nid, npt + OFFSET - grid.jsd, grid.npz),
    #     )
    # if grid.north_edge:
    #     je2 = ny - npt + 1
    #     avg_box(
    #         u,
    #         v,
    #         utmp,
    #         vtmp,
    #         origin=(grid.isd, je2, 0),
    #         domain=(grid.nid, grid.jed - je2 + 1, grid.npz),
    #     )

    js2 = npt + OFFSET if grid.south_edge else grid.jsd
    je2 = ny - npt if grid.north_edge else grid.jed
    jdiff = je2 - js2 + 1
    # if grid.west_edge:
    #     avg_box(
    #         u,
    #         v,
    #         utmp,
    #         vtmp,
    #         origin=(grid.isd, js2, 0),
    #         domain=(npt + OFFSET - grid.isd, jdiff, grid.npz),
    #     )
    # if grid.east_edge:
    #     avg_box(
    #         u,
    #         v,
    #         utmp,
    #         vtmp,
    #         origin=(nx + 1 - npt, js2, 0),
    #         domain=(grid.ied - nx + npt, jdiff, grid.npz),
    #     )

    # contra-variant components at cell center
    pad = 2 + 2 * id_
    # contravariant_components(
    #     utmp,
    #     vtmp,
    #     grid.cosa_s,
    #     grid.rsin2,
    #     ua,
    #     va,
    #     origin=(grid.is_ - 1 - id_, grid.js - 1 - id_, 0),
    #     domain=(grid.nic + pad, grid.njc + pad, grid.npz),
    # )
    # Fix the edges
    # Xdir:
    # TODO, make stencils? need variable offsets

    # if grid.sw_corner:
    #     for i in range(-2, 1):
    #         utmp[i + 2, grid.js - 1, :] = -vtmp[grid.is_ - 1, grid.js - i, :]
    # #     if spec.namelist.grid_type < 3:
    #         ua[i1 - 1, j1, :] = -va[i1, j1 + 2, :]
    #         ua[i1, j1, :] = -va[i1, j1 + 1, :]
    # if grid.se_corner:
    #     for i in range(0, 3):
    #         utmp[nx + i, grid.js - 1, :] = vtmp[nx, i + grid.js, :]
    #     if spec.namelist.grid_type < 3:
    #         ua[nx, j1, :] = va[nx, j1 + 1, :]
    #         ua[nx + 1, j1, :] = va[nx, j1 + 2, :]
    # if grid.ne_corner:
    #     for i in range(0, 3):
    #         utmp[nx + i, ny, :] = -vtmp[nx, grid.je - i, :]
    #     if spec.namelist.grid_type < 3:
    #         ua[nx, ny, :] = -va[nx, ny - 1, :]
    #         ua[nx + 1, ny, :] = -va[nx, ny - 2, :]
    # if grid.nw_corner:
    #     for i in range(-2, 1):
    #         utmp[i + 2, ny, :] = vtmp[grid.is_ - 1, grid.je + i, :]
    #     if spec.namelist.grid_type < 3:
    #         ua[i1 - 1, ny, :] = va[i1, ny - 2, :]
    #         ua[i1, ny, :] = va[i1, ny - 1, :]

    ifirst = grid.is_ + 2 if grid.west_edge else grid.is_ - 1
    ilast = grid.ie - 1 if grid.east_edge else grid.ie + 2
    idiff = ilast - ifirst + 1

    d2a2c_stencil1(
        dord4,
        u,
        v,
        ua,
        va,
        utc,
        vtc,
        utmp,
        vtmp,
        grid.cosa_s,
        grid.rsin2,
        origin=(grid.is_ - 3, grid.js - 3, 0),
        domain=(grid.nic + 6, grid.njc + 6, grid.npz),
    )

    # ut_main(
    #     utmp,
    #     uc,
    #     v,
    #     grid.cosa_u,
    #     grid.rsin_u,
    #     utc,
    #     origin=(ifirst, j1, 0),
    #     domain=(idiff, grid.njc + 2, grid.npz),
    # )

    d2a2c_stencil2(
        utmp,
        uc,
        utc,
        ua,
        v,
        grid.cosa_u,
        grid.rsin_u,
        grid.dxa,
        grid.sin_sg1,
        grid.sin_sg3,
        origin=(i1, j1, 0),
        domain=(grid.nic + 2, grid.njc + 2, grid.npz),
    )

    d2a2c_stencil_west(
        utmp,
        uc,
        utc,
        ua,
        v,
        grid.cosa_u,
        grid.rsin_u,
        grid.dxa,
        grid.sin_sg1,
        grid.sin_sg3,
        origin=(grid.is_ - 1, grid.js - 1, 0),
        domain=(4, grid.njc + 2, grid.npz),
    )

    d2a2c_stencil_east(
        utmp,
        uc,
        utc,
        ua,
        v,
        grid.cosa_u,
        grid.rsin_u,
        grid.dxa,
        grid.sin_sg1,
        grid.sin_sg3,
        origin=(grid.ie - 1, grid.js - 1, 0),
        domain=(4, grid.njc + 2, grid.npz),
    )

    if spec.namelist.grid_type < 3:
        domain_edge_x = (1, grid.njc + 2, grid.npz)
        jslice = slice(grid.js - 1, grid.je + 2)
        if grid.west_edge and not grid.nested:
            # vol_conserv_cubic_interp_x(
            #     utmp, uc, origin=(i1, j1, 0), domain=domain_edge_x
            # )
            islice = slice(grid.is_ - 2, grid.is_ + 2)
            # utc[grid.is_, jslice, :] = edge_interpolate4_x(
            #     ua[islice, jslice, :], grid.dxa[islice, jslice, :]
            # )
            # uc_x_edge1(
            #     utc,
            #     grid.sin_sg3,
            #     grid.sin_sg1,
            #     uc,
            #     origin=(i1 + 1, j1, 0),
            #     domain=domain_edge_x,
            # )
            # vol_conserv_cubic_interp_x_rev(
            #     utmp, uc, origin=(i1 + 2, j1, 0), domain=domain_edge_x
            # )
            # contravariant_stencil(
            #     uc,
            #     v,
            #     grid.cosa_u,
            #     grid.rsin_u,
            #     utc,
            #     origin=(i1, j1, 0),
            #     domain=domain_edge_x,
            # )
            # contravariant_stencil(
            #     uc,
            #     v,
            #     grid.cosa_u,
            #     grid.rsin_u,
            #     utc,
            #     origin=(i1 + 2, j1, 0),
            #     domain=domain_edge_x,
            # )

        if grid.east_edge and not grid.nested:
            # vol_conserv_cubic_interp_x(
            #     utmp, uc, origin=(nx - 1, j1, 0), domain=domain_edge_x
            # )
            islice = slice(nx - 2, nx + 2)
            # utc[nx, jslice, :] = edge_interpolate4_x(
            #     ua[islice, jslice, :], grid.dxa[islice, jslice, :]
            # )
            # uc_x_edge1(
            #     utc,
            #     grid.sin_sg3,
            #     grid.sin_sg1,
            #     uc,
            #     origin=(nx, j1, 0),
            #     domain=domain_edge_x,
            # )
            # vol_conserv_cubic_interp_x_rev(
            #     utmp, uc, origin=(nx + 1, j1, 0), domain=domain_edge_x
            # )
            # contravariant_stencil(
            #     uc,
            #     v,
            #     grid.cosa_u,
            #     grid.rsin_u,
            #     utc,
            #     origin=(grid.ie, j1, 0),
            #     domain=domain_edge_x,
            # )
            # contravariant_stencil(
            #     uc,
            #     v,
            #     grid.cosa_u,
            #     grid.rsin_u,
            #     utc,
            #     origin=(nx + 1, j1, 0),
            #     domain=domain_edge_x,
            # )
    # # Ydir:
    # if grid.sw_corner:
    #     for j in range(-2, 1):
    #         vtmp[grid.is_ - 1, j + 2, :] = -utmp[grid.is_ - j, grid.js - 1, :]
    #     va[i1, j1 - 1, :] = -ua[i1 + 2, j1, :]
    #     va[i1, j1, :] = -ua[i1 + 1, j1, :]
    # if grid.nw_corner:
    #     for j in range(0, 3):
    #         vtmp[grid.is_ - 1, ny + j, :] = utmp[j + grid.is_, ny, :]
    #     va[i1, ny, :] = ua[i1 + 1, ny, :]
    #     va[i1, ny + 1, :] = ua[i1 + 2, ny, :]
    # if grid.se_corner:
    #     for j in range(-2, 1):
    #         vtmp[nx, j + 2, :] = utmp[grid.ie + j, grid.js - 1, :]
    #     va[nx, j1, :] = ua[nx - 1, j1, :]
    #     va[nx, j1 - 1, :] = ua[nx - 2, j1, :]
    # if grid.ne_corner:
    #     for j in range(0, 3):
    #         vtmp[nx, ny + j, :] = -utmp[grid.ie - j, ny, :]
    #     va[nx, ny, :] = -ua[nx - 1, ny, :]
    #     va[nx, ny + 1, :] = -ua[nx - 2, ny, :]

    d2a2c_stencil3(
        dord4,
        u,
        v,
        ua,
        va,
        utc,
        vtc,
        utmp,
        vtmp,
        grid.cosa_s,
        grid.rsin2,
        origin=(grid.is_ - 3, grid.js - 3, 0),
        domain=(grid.nic + 6, grid.njc + 6, grid.npz),
    )

    d2a2c_stencil_south(
        vtmp,
        vc,
        vtc,
        va,
        u,
        grid.cosa_v,
        grid.rsin_v,
        grid.dya,
        grid.sin_sg2,
        grid.sin_sg4,
        origin=(grid.is_ - 1, grid.js - 1, 0),
        domain=(grid.nic + 2, 4, grid.npz),
    )

    # d2a2c_stencil_north(
    #     utmp,
    #     uc,
    #     utc,
    #     ua,
    #     v,
    #     grid.cosa_u,
    #     grid.rsin_u,
    #     grid.dxa,
    #     grid.sin_sg1,
    #     grid.sin_sg3,
    #     origin=(grid.js - 1, grid.ie - 1, 0),
    #     domain=(grid.nic + 2, 4, grid.npz),
    # )

    if spec.namelist.grid_type < 3:
        domain_edge_y = (grid.nic + 2, 1, grid.npz)
        islice = slice(grid.is_ - 1, grid.ie + 2)
        if grid.south_edge:
            # vt_edge(
            #     vtmp,
            #     vc,
            #     u,
            #     grid.cosa_v,
            #     grid.rsin_v,
            #     vtc,
            #     0,
            #     origin=(i1, j1, 0),
            #     domain=domain_edge_y,
            # )
            jslice = slice(grid.js - 2, grid.js + 2)
            # vtc[islice, grid.js, :] = edge_interpolate4_y(
            #     va[islice, jslice, :], grid.dya[islice, jslice, :]
            # )
            # vc_y_edge1(
            #     vtc,
            #     grid.sin_sg4,
            #     grid.sin_sg2,
            #     vc,
            #     origin=(i1, grid.js, 0),
            #     domain=domain_edge_y,
            # )
            # vt_edge(
            #     vtmp,
            #     vc,
            #     u,
            #     grid.cosa_v,
            #     grid.rsin_v,
            #     vtc,
            #     1,
            #     origin=(i1, grid.js + 1, 0),
            #     domain=domain_edge_y,
            # )
        if grid.north_edge:
            vt_edge(
                vtmp,
                vc,
                u,
                grid.cosa_v,
                grid.rsin_v,
                vtc,
                0,
                origin=(i1, ny - 1, 0),
                domain=domain_edge_y,
            )
            jslice = slice(ny - 2, ny + 2)
            vtc[islice, ny, :] = edge_interpolate4_y(
                va[islice, jslice, :], grid.dya[islice, jslice, :]
            )
            vc_y_edge1(
                vtc,
                grid.sin_sg4,
                grid.sin_sg2,
                vc,
                origin=(grid.is_ - 1, ny, 0),
                domain=domain_edge_y,
            )
            vt_edge(
                vtmp,
                vc,
                u,
                grid.cosa_v,
                grid.rsin_v,
                vtc,
                1,
                origin=(i1, ny + 1, 0),
                domain=domain_edge_y,
            )
        jfirst = grid.js + 2 if grid.south_edge else grid.js - 1
        jlast = grid.je - 1 if grid.north_edge else grid.je + 2
        jdiff = jlast - jfirst + 1

        vt_main(
            vtmp,
            vc,
            u,
            grid.cosa_v,
            grid.rsin_v,
            vtc,
            origin=(i1, jfirst, 0),
            domain=(grid.nic + 2, jdiff, grid.npz),
        )
    else:
        raise Exception("unimplemented grid_type >= 3")
