import fv3.utils.gt4py_utils as utils
import gt4py.gtscript as gtscript
import fv3._config as spec
from gt4py.gtscript import computation, interval, PARALLEL


def grid():
    return spec.grid


sd = utils.sd
stencil_corner = True


@utils.stencil()
def main_ut(uc: sd, vc: sd, cosa_u: sd, rsin_u: sd, ut: sd):
    with computation(PARALLEL), interval(...):
        ut[0, 0, 0] = (
            uc - 0.25 * cosa_u * (vc[-1, 0, 0] + vc + vc[-1, 1, 0] + vc[0, 1, 0])
        ) * rsin_u


@utils.stencil()
def ut_y_edge(uc: sd, sin_sg1: sd, sin_sg3: sd, ut: sd, *, dt: float):
    with computation(PARALLEL), interval(0, -1):
        ut[0, 0, 0] = (uc / sin_sg3[-1, 0, 0]) if (uc * dt > 0) else (uc / sin_sg1)


@utils.stencil()
def ut_x_edge(uc: sd, cosa_u: sd, vt: sd, ut: sd):
    with computation(PARALLEL), interval(0, -1):
        ut[0, 0, 0] = uc - 0.25 * cosa_u * (
            vt[-1, 0, 0] + vt[0, 0, 0] + vt[-1, 1, 0] + vt[0, 1, 0]
        )


@utils.stencil()
def main_vt(uc: sd, vc: sd, cosa_v: sd, rsin_v: sd, vt: sd):
    with computation(PARALLEL), interval(...):
        vt[0, 0, 0] = (
            vc - 0.25 * cosa_v * (uc[0, -1, 0] + uc[1, -1, 0] + uc + uc[1, 0, 0])
        ) * rsin_v


@utils.stencil()
def vt_y_edge(vc: sd, cosa_v: sd, ut: sd, vt: sd):
    with computation(PARALLEL), interval(0, -1):
        vt[0, 0, 0] = vc - 0.25 * cosa_v * (
            ut[0, -1, 0] + ut[1, -1, 0] + ut[0, 0, 0] + ut[1, 0, 0]
        )


@utils.stencil()
def vt_x_edge(vc: sd, sin_sg2: sd, sin_sg4: sd, vt: sd, *, dt: float):
    with computation(PARALLEL), interval(0, -1):
        vt[0, 0, 0] = (vc / sin_sg4[0, -1, 0]) if (vc * dt > 0) else (vc / sin_sg2)


@gtscript.function
def ra_x_func(area, xfx_adv):
    return area + xfx_adv - xfx_adv[1, 0, 0]


@utils.stencil()
def xfx_adv_stencil(
    ut: sd,
    rdxa: sd,
    area: sd,
    dy: sd,
    sin_sg1: sd,
    sin_sg3: sd,
    crx_adv: sd,
    xfx_adv: sd,
    ra_x: sd,
    dt: float,
):
    with computation(PARALLEL), interval(...):
        xfx_adv[0, 0, 0] = dt * ut
        crx_adv[0, 0, 0] = xfx_adv * rdxa[-1, 0, 0] if xfx_adv > 0 else xfx_adv * rdxa
        xfx_adv[0, 0, 0] = (
            dy * xfx_adv * sin_sg3[-1, 0, 0] if xfx_adv > 0 else dy * xfx_adv * sin_sg1
        )
        ra_x = ra_x_func(area, xfx_adv)


@gtscript.function
def ra_y_func(area, yfx_adv):
    return area + yfx_adv - yfx_adv[0, 1, 0]


@utils.stencil()
def yfx_adv_stencil(
    vt: sd,
    rdya: sd,
    area: sd,
    dx: sd,
    sin_sg2: sd,
    sin_sg4: sd,
    cry_adv: sd,
    yfx_adv: sd,
    ra_y: sd,
    dt: float,
):
    with computation(PARALLEL), interval(...):
        yfx_adv[0, 0, 0] = dt * vt
        cry_adv[0, 0, 0] = yfx_adv * rdya[0, -1, 0] if yfx_adv > 0 else yfx_adv * rdya
        yfx_adv[0, 0, 0] = (
            dx * yfx_adv * sin_sg4[0, -1, 0] if yfx_adv > 0 else dx * yfx_adv * sin_sg2
        )
        ra_y = ra_y_func(area, yfx_adv)


def compute(uc_in, vc_in, ut_in, vt_in, xfx_adv, yfx_adv, crx_adv, cry_adv, dt, kstart=0, nk=None):
    if nk is None:
        nk = grid().npz - kstart
    #global kstart_fx
    #global nz_fx
    #global kext
    kstart_fx = kstart
    nz_fx = nk + 1 if nk < grid().npz - 4 else nk
    kext = slice(kstart_fx, kstart_fx + nz_fx)
    ut = compute_ut(uc_in, vc_in, grid().cosa_u, grid().rsin_u, ut_in, kstart_fx, nz_fx)
    vt = compute_vt(
        uc_in,
        vc_in,
        grid().cosa_v,
        grid().rsin_v,
        grid().sin_sg2,
        grid().sin_sg4,
        vt_in, kstart_fx, nz_fx
    )
    update_ut_y_edge(uc_in, grid().sin_sg1, grid().sin_sg3, ut, dt, kstart_fx, nz_fx)
    update_vt_y_edge(vc_in, grid().cosa_v, ut, vt, kstart_fx, nz_fx)
    update_vt_x_edge(vc_in, grid().sin_sg2, grid().sin_sg4, vt, dt, kstart_fx, nz_fx)
    update_ut_x_edge(uc_in, grid().cosa_u, vt, ut, kstart_fx, nz_fx)
    corner_shape = (1, 1, nz_fx)
    if grid().sw_corner:
        sw_corner(uc_in, vc_in, ut, vt, grid().cosa_u, grid().cosa_v, corner_shape, kstart_fx, nz_fx)
    if grid().se_corner:
        se_corner(uc_in, vc_in, ut, vt, grid().cosa_u, grid().cosa_v, corner_shape, kstart_fx, nz_fx)
    if grid().ne_corner:
        ne_corner(uc_in, vc_in, ut, vt, grid().cosa_u, grid().cosa_v, corner_shape, kstart_fx, nz_fx)
    if grid().nw_corner:
        nw_corner(uc_in, vc_in, ut, vt, grid().cosa_u, grid().cosa_v, corner_shape, kstart_fx, nz_fx)
    compute_x_origin = (grid().is_, grid().jsd, kstart_fx)
    compute_y_origin = (grid().isd, grid().js, kstart_fx)
    ra_x = utils.make_storage_from_shape(uc_in.shape, compute_x_origin)
    xfx_adv_stencil(
        ut,
        grid().rdxa,
        grid().area,
        grid().dy,
        grid().sin_sg1,
        grid().sin_sg3,
        crx_adv,
        xfx_adv,
        ra_x,
        dt,
        origin=compute_x_origin,
        domain=(grid().nic + 1, grid().njd, nz_fx),
    )
    ra_y = utils.make_storage_from_shape(vc_in.shape, compute_y_origin)
    yfx_adv_stencil(
        vt,
        grid().rdya,
        grid().area,
        grid().dx,
        grid().sin_sg2,
        grid().sin_sg4,
        cry_adv,
        yfx_adv,
        ra_y,
        dt,
        origin=compute_y_origin,
        domain=(grid().nid, grid().njc + 1, nz_fx),
    )
    # TODO remove the need for a copied extra ut and vt variables, edit in place (resolve issue with data getting zeroed out)
    ut_in[:, :, kext] = ut[:, :, kext]
    vt_in[:, :, kext] = vt[:, :, kext]
    return ra_x, ra_y


def compute_ut(uc_in, vc_in, cosa_u, rsin_u, ut_in, kstart_fx, nz_fx):
    ut_origin = (grid().is_ - 1, grid().jsd, kstart_fx)
    ut = utils.make_storage_from_shape(ut_in.shape, ut_origin)
    kext = slice(kstart_fx, kstart_fx + nz_fx)
    main_ut(
        uc_in,
        vc_in,
        cosa_u,
        rsin_u,
        ut,
        origin=ut_origin,
        domain=(grid().nic + 3, grid().njd, nz_fx),
    )
    ut.data[: grid().is_ - 1, :, kext] = ut_in.data[: grid().is_ - 1, :, kext]
    ut.data[grid().ie + 3 :, :, kext] = ut_in.data[grid().ie + 3 :, :, kext]
    # fill in for j /=2 and j/=3
    if grid().south_edge:
        ut.data[:, grid().js - 1 : grid().js + 1, kext] = ut_in.data[
            :, grid().js - 1 : grid().js + 1, kext
        ]
    # fill in for j/=npy-1 and j /= npy
    if grid().north_edge:
        ut.data[:, grid().je : grid().je + 2, kext] = ut_in.data[
            :, grid().je : grid().je + 2, kext
        ]
    return ut


def update_ut_y_edge(uc, sin_sg1, sin_sg3, ut, dt,kstart_fx, nz_fx):
    edge_shape = (1, ut.shape[1], nz_fx)
    if grid().west_edge:
        ut_y_edge(
            uc,
            sin_sg1,
            sin_sg3,
            ut,
            dt=dt,
            origin=(grid().is_, 0, kstart_fx),
            domain=edge_shape,
        )
    if grid().east_edge:
        ut_y_edge(
            uc,
            sin_sg1,
            sin_sg3,
            ut,
            dt=dt,
            origin=(grid().ie + 1, 0, kstart_fx),
            domain=edge_shape,
        )


def update_ut_x_edge(uc, cosa_u, vt, ut, kstart_fx, nz_fx):
    i1 = grid().is_ + 2 if grid().west_edge else grid().is_
    i2 = grid().ie - 1 if grid().east_edge else grid().ie + 1
    edge_shape = (i2 - i1 + 1, 2, nz_fx)
    if grid().south_edge:
        ut_x_edge(uc, cosa_u, vt, ut, origin=(i1, grid().js - 1, kstart_fx), domain=edge_shape)
    if grid().north_edge:
        ut_x_edge(uc, cosa_u, vt, ut, origin=(i1, grid().je, kstart_fx), domain=edge_shape)


def compute_vt(uc_in, vc_in, cosa_v, rsin_v, sin_sg2, sin_sg4, vt_in, kstart_fx, nz_fx):
    vt_origin = (grid().isd, grid().js - 1, kstart_fx)
    vt = utils.make_storage_from_shape(vt_in.shape, vt_origin)
    kext = slice(kstart_fx, nz_fx)
    main_vt(
        uc_in,
        vc_in,
        cosa_v,
        rsin_v,
        vt,
        origin=vt_origin,
        domain=(grid().nid, grid().njc + 3, nz_fx),
    )  # , origin=(0, 2, kstart_fx), domain=(vt.shape[0]-1, main_j_size, nz_fx))
    # cannot pass vt_in array to stencil without it zeroing out data outside specified domain
    # So... for now copying in so the 'undefined' answers match
    vt.data[:, : grid().js - 1, kext] = vt_in.data[:, : grid().js - 1, kext]
    vt.data[:, grid().je + 3, kext] = vt_in.data[:, grid().je + 3, kext]
    if grid().south_edge:
        vt.data[:, grid().js, kext] = vt_in.data[:, grid().js, kext]
    if grid().north_edge:
        vt.data[:, grid().je + 1, kext] = vt_in.data[:, grid().je + 1, kext]
    return vt


def update_vt_y_edge(vc, cosa_v, ut, vt, kstart_fx, nz_fx):
    if grid().west_edge or grid().east_edge:
        j1 = grid().js + 2 if grid().south_edge else grid().js
        j2 = grid().je if grid().north_edge else grid().je + 2
        edge_shape = (2, j2 - j1, nz_fx)
        if grid().west_edge:
            vt_y_edge(
                vc, cosa_v, ut, vt, origin=(grid().is_ - 1, j1, kstart_fx), domain=edge_shape
            )
        if grid().east_edge:
            vt_y_edge(vc, cosa_v, ut, vt, origin=(grid().ie, j1, kstart_fx), domain=edge_shape)


def update_vt_x_edge(vc, sin_sg2, sin_sg4, vt, dt, kstart_fx, nz_fx):
    if grid().south_edge or grid().north_edge:
        edge_shape = (vt.shape[0], 1, nz_fx)
        if grid().south_edge:
            vt_x_edge(
                vc,
                sin_sg2,
                sin_sg4,
                vt,
                dt=dt,
                origin=(0, grid().js, kstart_fx),
                domain=edge_shape,
            )
        if grid().north_edge:
            vt_x_edge(
                vc,
                sin_sg2,
                sin_sg4,
                vt,
                dt=dt,
                origin=(0, grid().je + 1, kstart_fx),
                domain=edge_shape,
            )


# -------------------- CORNERS-----------------


def corner_ut_stencil(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    from __externals__ import vi, vj, ux, uy, vx, vy

    with computation(PARALLEL), interval(...):
        ut[0, 0, 0] = (
            (
                uc[0, 0, 0]
                - 0.25
                * cosa_u[0, 0, 0]
                * (
                    vt[vi, vy, 0]
                    + vt[vx, vy, 0]
                    + vt[vx, vj, 0]
                    + vc[vi, vj, 0]
                    - 0.25
                    * cosa_v[vi, vj, 0]
                    * (ut[ux, 0, 0] + ut[ux, uy, 0] + ut[0, uy, 0])
                )
            )
            * 1.0
            / (1.0 - 0.0625 * cosa_u[0, 0, 0] * cosa_v[vi, vj, 0])
        )


# for the non-stencil version of filling corners
def get_damp(cosa_u, cosa_v, ui, uj, vi, vj, kext):
    return 1.0 / (1.0 - 0.0625 * cosa_u[ui, uj, kext] * cosa_v[vi, vj, kext])


def index_offset(lower, u, south=True):
    if lower == u:
        offset = 1
    else:
        offset = -1
    if south:
        offset *= -1
    return offset


def corner_ut(
    uc,
    vc,
    ut,
    vt,
    cosa_u,
    cosa_v,
    ui,
    uj,
    vi,
    vj,kstart_fx, nz_fx,
    west,
    lower, 
    south=True,
    vswitch=False
):
    if vswitch:
        lowerfactor = 1 if lower else -1
    else:
        lowerfactor = 1
    vx = vi + index_offset(west, False, south) * lowerfactor
    ux = ui + index_offset(west, True, south) * lowerfactor
    vy = vj + index_offset(lower, False, south) * lowerfactor
    uy = uj + index_offset(lower, True, south) * lowerfactor
    if stencil_corner:
        decorator = gtscript.stencil(
            backend=utils.backend,
            externals={
                "vi": vi - ui,
                "vj": vj - uj,
                "ux": ux - ui,
                "uy": uy - uj,
                "vx": vx - ui,
                "vy": vy - uj,
            },
            rebuild=True,
        )
        corner_stencil = decorator(corner_ut_stencil)
        corner_stencil(
            uc,
            vc,
            ut,
            vt,
            cosa_u,
            cosa_v,
            origin=(ui, uj, kstart_fx),
            domain=(1, 1, nz_fx),
        )
    else:
        kext =  slice(kstart_fx, nz_fx)
        damp = get_damp(cosa_u, cosa_v, ui, uj, vi, vj, kext)
        ut[ui, uj, kext] = (
            uc[ui, uj, kext]
            - 0.25
            * cosa_u[ui, uj, kext]
            * (
                vt[vi, vy, kext]
                + vt[vx, vy, kext]
                + vt[vx, vj, kext]
                + vc[vi, vj, kext]
                - 0.25
                * cosa_v[vi, vj, kext]
                * (ut[ux, uj, kext] + ut[ux, uy, kext] + ut[ui, uy, kext])
            )
        ) * damp


def sw_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape, kstart_fx, nz_fx):
    t = grid().is_ + 1
    n = grid().is_
    z = grid().is_ - 1
    corner_ut(uc, vc, ut, vt, cosa_u, cosa_v, t, z, n, z, kstart_fx, nz_fx, west=True, lower=True)
    corner_ut(
        vc, uc, vt, ut, cosa_v, cosa_u, z, t, z, n, kstart_fx, nz_fx, west=True, lower=True, vswitch=True
    )
    corner_ut(uc, vc, ut, vt, cosa_u, cosa_v, t, n, n, t, kstart_fx, nz_fx, west=True, lower=False)
    corner_ut(
        vc, uc, vt, ut, cosa_v, cosa_u, n, t, t, n, kstart_fx, nz_fx, west=True, lower=False, vswitch=True
    )


def se_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape, kstart_fx, nz_fx):
    t = grid().js + 1
    n = grid().js
    z = grid().js - 1
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        grid().ie,
        z,
        grid().ie,
        z, kstart_fx, nz_fx,
        west=False,
        lower=True,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        grid().ie + 1,
        t,
        grid().ie + 2,
        n, kstart_fx, nz_fx,
        west=False,
        lower=True,
        vswitch=True,
    )
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        grid().ie,
        n,
        grid().ie,
        t, kstart_fx, nz_fx,
        west=False,
        lower=False,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        grid().ie,
        t,
        grid().ie,
        n, kstart_fx, nz_fx,
        west=False,
        lower=False,
        vswitch=True,
    )


def ne_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape, kstart_fx, nz_fx,):
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        grid().ie,
        grid().je + 1,
        grid().ie,
        grid().je + 2, kstart_fx, nz_fx,
        west=False,
        lower=False,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        grid().ie + 1,
        grid().je,
        grid().ie + 2,
        grid().je, kstart_fx, nz_fx,
        west=False,
        lower=False,
        south=False,
        vswitch=True,
    )
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        grid().ie,
        grid().je,
        grid().ie,
        grid().je, kstart_fx, nz_fx,
        west=False,
        lower=True,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        grid().ie,
        grid().je,
        grid().ie,
        grid().je, kstart_fx, nz_fx,
        west=False,
        lower=True,
        south=False,
        vswitch=True,
    )


def nw_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape, kstart_fx, nz_fx):
    t = grid().js + 1
    n = grid().js
    z = grid().js - 1
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        t,
        grid().je + 1,
        n,
        grid().je + 2, kstart_fx, nz_fx,
        west=True,
        lower=False,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        z,
        grid().je,
        z,
        grid().je, kstart_fx, nz_fx,
        west=True,
        lower=False,
        south=False,
        vswitch=True,
    )
    corner_ut(
        uc,
        vc,
        ut,
        vt,
        cosa_u,
        cosa_v,
        t,
        grid().je,
        n,
        grid().je, kstart_fx, nz_fx,
        west=True,
        lower=True,
    )
    corner_ut(
        vc,
        uc,
        vt,
        ut,
        cosa_v,
        cosa_u,
        n,
        grid().je,
        t,
        grid().je, kstart_fx, nz_fx,
        west=True,
        lower=True,
        south=False,
        vswitch=True,
    )


# TODO Probably can delete -- but in case we want to do analysis to show it doesn't matter at all
"""
def sw_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape):
    west_corner_ut_lowest(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ + 1, grid().js - 1, kstart_fx), domain=corner_shape)
    west_corner_ut_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ + 1, grid().js, kstart_fx), domain=corner_shape)
    south_corner_vt_left(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ - 1, grid().js + 1, kstart_fx), domain=corner_shape)
    south_corner_vt_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_, grid().js + 1, kstart_fx), domain=corner_shape)


def se_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape):
    east_corner_ut_lowest(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().js - 1, kstart_fx), domain=corner_shape)
    east_corner_ut_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().js, kstart_fx), domain=corner_shape)
    south_corner_vt_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie + 1, grid().js + 1, kstart_fx), domain=corner_shape)
    south_corner_vt_left(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().js + 1, kstart_fx), domain=corner_shape)

def ne_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape):
    east_corner_ut_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().je + 1, kstart_fx), domain=corner_shape)
    east_corner_ut_lowest(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().je, kstart_fx), domain=corner_shape)
    north_corner_vt_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie + 1, grid().je, kstart_fx), domain=corner_shape)
    north_corner_vt_left(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().ie, grid().je, kstart_fx), domain=corner_shape)

def nw_corner(uc, vc, ut, vt, cosa_u, cosa_v, corner_shape):
    west_corner_ut_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ + 1, grid().je+1, kstart_fx), domain=corner_shape)
    west_corner_ut_lowest(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ + 1, grid().je, kstart_fx), domain=corner_shape)
    north_corner_vt_left(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_ - 1, grid().je, kstart_fx), domain=corner_shape)
    north_corner_vt_adjacent(uc, vc, ut, vt, cosa_u, cosa_v, origin=(grid().is_, grid().je, kstart_fx), domain=corner_shape)


@utils.stencil()
def west_corner_ut_lowest(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_u = 1. / (1.0 - 0.0625 * cosa_u[0, 0, 0] * cosa_v[-1, 0, 0])
        ut[0, 0, 0] = (uc[0, 0, 0]-0.25 * cosa_u[0, 0, 0] * (vt[-1, 1, 0] + vt[0, 1, 0] + vt[0, 0, 0] + vc[-1, 0, 0] -
                                                             0.25 * cosa_v[-1, 0, 0] * (ut[-1, 0, 0] + ut[-1, -1, 0] +
                                                                                        ut[0, -1, 0]))) * damp_u


@utils.stencil()
def west_corner_ut_adjacent(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp = 1. / (1. - 0.0625 * cosa_u[0, 0, 0] * cosa_v[-1, 1, 0])
        ut[0, 0, 0] = (uc[0, 0, 0] - 0.25 * cosa_u[0, 0, 0] * (vt[-1, 0, 0] + vt[0, 0, 0] + vt[0, 1, 0] + vc[-1, 1, 0] -
                                                               0.25 * cosa_v[-1, 1, 0] * (ut[-1, 0, 0] + ut[-1, 1, 0] +
                                                                                          ut[0, 1, 0]))) * damp


@utils.stencil()
def south_corner_vt_left(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_v = 1. / (1.0 - 0.0625 * cosa_u[0, -1, 0] * cosa_v[0, 0, 0])
        vt[0, 0, 0] = (vc[0, 0, 0] - 0.25 * cosa_v[0, 0, 0] * (ut[1, -1, 0] + ut[1, 0, 0] + ut[0, 0, 0] + uc[0, -1, 0] -
                                                               0.25 * cosa_u[0, -1, 0] *
                                                               (vt[0, -1, 0] + vt[-1, -1, 0] + vt[-1, 0, 0]))) * damp_v


@utils.stencil()
def south_corner_vt_adjacent(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_v = 1. / (1.0 - 0.0625 * cosa_u[1, -1, 0] * cosa_v[0, 0, 0])
        vt[0, 0, 0] = (vc[0, 0, 0] - 0.25 * cosa_v[0, 0, 0] * (ut[0, -1, 0] + ut[0, 0, 0] + ut[1, 0, 0] + uc[1, -1, 0] -
                                                               0.25 * cosa_u[1, -1, 0] *
                                                               (vt[0, -1, 0] + vt[1, -1, 0] + vt[1, 0, 0]))) * damp_v


@utils.stencil()
def east_corner_ut_lowest(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_u = 1. / (1.0 - 0.0625 * cosa_u[0, 0, 0] * cosa_v[0, 0, 0])
        ut[0, 0, 0] = (uc[0, 0, 0]-0.25 * cosa_u[0, 0, 0] * (vt[0, 1, 0] + vt[-1, 1, 0] + vt[-1, 0, 0] + vc[0, 0, 0] -
                                                             0.25 * cosa_v[0, 0, 0] * (ut[1, 0, 0] + ut[1, -1, 0] +
                                                                                       ut[0, -1, 0]))) * damp_u


@utils.stencil()
def east_corner_ut_adjacent(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp = 1. / (1. - 0.0625 * cosa_u[0, 0, 0] * cosa_v[0, 1, 0])
        ut[0, 0, 0] = (uc[0, 0, 0] - 0.25 * cosa_u[0, 0, 0] * (vt[0, 0, 0] + vt[-1, 0, 0] + vt[-1, 1, 0] + vc[0, 1, 0] -
                                                               0.25 * cosa_v[0, 1, 0] * (ut[1, 0, 0] + ut[1, 1, 0] +
                                                                                         ut[0, 1, 0]))) * damp


@utils.stencil()
def north_corner_vt_left(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_v = 1. / (1.0 - 0.0625 * cosa_u[0, 0, 0] * cosa_v[0, 0, 0])
        vt[0, 0, 0] = (vc[0, 0, 0] - 0.25 * cosa_v[0, 0, 0] * (ut[1, 0, 0] + ut[1, -1, 0] + ut[0, -1, 0] + uc[0, 0, 0] -
                                                               0.25 * cosa_u[0, 0, 0] *
                                                               (vt[0, 1, 0] + vt[-1, 1, 0] + vt[-1, 0, 0]))) * damp_v

@utils.stencil()
def north_corner_vt_adjacent(uc: sd, vc: sd, ut: sd, vt: sd, cosa_u: sd, cosa_v: sd):
    with computation(PARALLEL), interval(...):
        damp_v = 1. / (1.0 - 0.0625 * cosa_u[1, 0, 0] * cosa_v[0, 0, 0])
        vt[0, 0, 0] = (vc[0, 0, 0] - 0.25 * cosa_v[0, 0, 0] * (ut[0, 0, 0] + ut[0, -1, 0] + ut[1, -1, 0] + uc[1, 0, 0] -
                                                               0.25 * cosa_u[1, 0, 0] *
                                                               (vt[0, 1, 0] + vt[1, 1, 0] + vt[1, 0, 0]))) * damp_v
"""
