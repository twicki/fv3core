import numpy as np

import gt4py.gtscript as gtscript
from gt4py.gtscript import FORWARD, PARALLEL, computation, interval

import fv3core._config as spec
import fv3core.utils.gt4py_utils as utils
from fv3core.decorators import gtstencil


FloatField = utils.FloatField
FloatFieldIJ = utils.FloatFieldIJ
IntFieldIJ = utils.IntFieldIJ
TracerTuple = utils.TracerTuple


@gtstencil()
def fix_tracer(
    q: TracerTuple,
    dp: FloatField,
    dm: FloatField,
    dm_pos: FloatField,
    zfix: IntFieldIJ,
    upper_fix: FloatField,
    lower_fix: FloatField,
    sum0: FloatFieldIJ,
    sum1: FloatFieldIJ,
    fac: FloatFieldIJ,
):
    # reset fields
    with computation(PARALLEL), interval(...):
        zfix = 0
        lower_fix = 0.0
        upper_fix = 0.0
        sum0 = 0.0
        sum1 = 0.0
    # fix_top:
    with computation(PARALLEL):
        with interval(0, 1):
            if q < 0:
                q = 0
            dm = q * dp
        with interval(1, 2):
            if q[0, 0, -1] < 0.0:
                q = (
                    q + q[0, 0, -1] * dp[0, 0, -1] / dp
                )  # move enough mass up so that the top layer isn't negative
    # fix_interior:
    with computation(FORWARD), interval(1, -1):
        # if a higher layer borrowed from this one, account for that here
        if lower_fix[0, 0, -1] != 0.0:
            q = q - (lower_fix[0, 0, -1] / dp)
        dq = q * dp
        if q < 0.0:
            zfix += 1
            if q[0, 0, -1] > 0.0:
                # Borrow from the layer above
                dq = (
                    q[0, 0, -1] * dp[0, 0, -1]
                    if q[0, 0, -1] * dp[0, 0, -1] < -(q * dp)
                    else -(q * dp)
                )
                q = q + dq / dp
                upper_fix = dq
            if (q < 0.0) and (q[0, 0, 1] > 0.0):
                # borrow from the layer below
                dq = (
                    q[0, 0, 1] * dp[0, 0, 1]
                    if q[0, 0, 1] * dp[0, 0, 1] < -(q * dp)
                    else -(q * dp)
                )
                q = q + dq / dp
                lower_fix = dq
    with computation(PARALLEL), interval(...):
        if upper_fix[0, 0, 1] != 0.0:
            # If a lower layer borrowed from this one, account for that here
            q = q - upper_fix[0, 0, 1] / dp
        dm = q * dp
        dm_pos = dm if dm > 0.0 else 0.0
    # fix_bottom:
    with computation(PARALLEL), interval(-1, None):
        # the 2nd-to-last layer borrowed from this one, account for that here
        if lower_fix[0, 0, -1] != 0.0:
            q = q - (lower_fix[0, 0, -1] / dp)
        qup = q[0, 0, -1] * dp[0, 0, -1]
        qly = -q * dp
        dup = qup if qup < qly else qly
        if (q < 0.0) and (q[0, 0, -1] > 0.0):
            zfix += 1
            q = q + (dup / dp)
            upper_fix = dup
        dm = q * dp
        dm_pos = dm if dm > 0.0 else 0.0
    with computation(PARALLEL), interval(-2, -1):
        # if the bottom layer borrowed from this one, adjust
        if upper_fix[0, 0, 1] != 0.0:
            q = q - (upper_fix[0, 0, 1] / dp)
            dm = q * dp
            dm_pos = dm if dm > 0.0 else 0.0  # now we gotta update these too
    with computation(FORWARD), interval(1, None):
        sum0 += dm
        sum1 += dm_pos
    with computation(PARALLEL), interval(1, None):
        fac = sum0 / sum1 if sum0 > 0.0 else 0.0
        if zfix > 0 and fac > 0.0:
            q = fac * dm / dp if fac * dm / dp > 0.0 else 0.0


def compute(dp2, tracers, im, km, nq, jslice):
    # Same as above, but with multiple tracer fields
    i1 = spec.grid.is_
    js = jslice.start
    jext = jslice.stop - jslice.start

    tracer_list = [tracers[q] for q in utils.tracer_variables[0:nq]]
    shape = tracer_list[0].shape
    shape_ij = shape[0:2]

    dm = utils.make_storage_from_shape(shape, origin=(0, 0, 0))
    dm_pos = utils.make_storage_from_shape(shape, origin=(0, 0, 0))
    upper_fix = utils.make_storage_from_shape(shape, origin=(0, 0, 0))
    lower_fix = utils.make_storage_from_shape(shape, origin=(0, 0, 0))
    zfix = utils.make_storage_from_shape(shape_ij, dtype=np.int, origin=(0, 0))
    fac = utils.make_storage_from_shape(shape_ij, origin=(0, 0))
    sum0 = utils.make_storage_from_shape(shape_ij, origin=(0, 0))
    sum1 = utils.make_storage_from_shape(shape_ij, origin=(0, 0))
    # TODO: implement dev_gfs_physics ifdef when we implement compiler defs

    # for tracer in tracer_list:
    fix_tracer(
        *tracer_list,
        # tracer_list[n],
        dp2,
        dm,
        dm_pos,
        zfix,
        upper_fix,
        lower_fix,
        sum0,
        sum1,
        fac,
        origin=(i1, js, 0),
        domain=(im, jext, km),
    )
    return tracer_list
