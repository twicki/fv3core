#!/usr/bin/env python3
import fv3.utils.gt4py_utils as utils
import gt4py.gtscript as gtscript
import fv3._config as spec
from gt4py.gtscript import computation, interval, PARALLEL

sd = utils.sd


@utils.stencil()
def vorticity(u: sd, dx: sd, vt: sd):
    with computation(PARALLEL), interval(...):
        vt[0, 0, 0] = u * dx


@utils.stencil()
def volume_mean_relative_vorticity(ut: sd, vt: sd, rarea: sd, wk: sd):
    with computation(PARALLEL), interval(...):
        wk[0, 0, 0] = rarea * (vt - vt[0, 1, 0] - ut + ut[1, 0, 0])


def compute(u, v, ut, vt, wk, kstart=0, nk=None):
    grid = spec.grid
    if nk is None:
        nk = grid.npz - kstart
    default_origin = (grid.isd, grid.jsd, kstart)
    vorticity(
        u,
        grid.dx,
        vt,
        origin=default_origin,
        domain=(grid.nid, grid.njd + 1, nk),
    )
    vorticity(
        v,
        grid.dy,
        ut,
        origin=default_origin,
        domain=(grid.nid + 1, grid.njd, nk),
    )
    volume_mean_relative_vorticity(
        ut,
        vt,
        grid.rarea,
        wk,
        origin=default_origin,
        domain=(grid.nid, grid.njd, nk),
    )
