from typing import Any, Dict, Optional

import gt4py.gtscript as gtscript
from gt4py.gtscript import PARALLEL, computation, interval

import fv3core._config as spec
import fv3core.utils.gt4py_utils as utils
from fv3core.decorators import gtstencil
from fv3core.utils.typing import FloatField


@gtscript.function
def p_grad_c_u(uc_in, wk, pkc, gz, rdxc, dt2):
    return uc_in + dt2 * rdxc / (wk[-1, 0, 0] + wk) * (
        (gz[-1, 0, 1] - gz) * (pkc[0, 0, 1] - pkc[-1, 0, 0])
        + (gz[-1, 0, 0] - gz[0, 0, 1]) * (pkc[-1, 0, 1] - pkc)
    )


@gtscript.function
def get_wk(pkc, delpc, hydrostatic):
    return pkc[0, 0, 1] - pkc if hydrostatic else delpc


@gtscript.function
def p_grad_c_u_wk(uc_in, delpc, pkc, gz, rdxc, hydrostatic, dt2):
    wk = get_wk(pkc, delpc, hydrostatic)
    return p_grad_c_u(uc_in, wk, pkc, gz, rdxc, dt2)


@gtscript.function
def p_grad_c_v(vc_in, wk, pkc, gz, rdyc, dt2):
    return vc_in + dt2 * rdyc / (wk[0, -1, 0] + wk) * (
        (gz[0, -1, 1] - gz) * (pkc[0, 0, 1] - pkc[0, -1, 0])
        + (gz[0, -1, 0] - gz[0, 0, 1]) * (pkc[0, -1, 1] - pkc)
    )


@gtscript.function
def p_grad_c_v_wk(vc_in, delpc, pkc, gz, rdyc, hydrostatic, dt2):
    wk = get_wk(pkc, delpc, hydrostatic)
    return p_grad_c_v(vc_in, wk, pkc, gz, rdyc, dt2)


@gtscript.function
def p_grad_c_fn(uc_in, vc_in, delpc, pkc, gz, rdxc, rdyc, hydrostatic, dt2):
    wk = get_wk(pkc, delpc, hydrostatic)
    uc_in = p_grad_c_u(uc_in, wk, pkc, gz, rdxc, dt2)
    vc_in = p_grad_c_v(vc_in, wk, pkc, gz, rdyc, dt2)
    return uc_in, vc_in


@gtstencil()
def p_grad_c(
    uc_in: FloatField,
    vc_in: FloatField,
    delpc: FloatField,
    pkc: FloatField,
    gz: FloatField,
    rdxc: FloatField,
    rdyc: FloatField,
    hydrostatic: int,
    dt2: float,
):
    with computation(PARALLEL), interval(0, -1):
        uc_in, vc_in = p_grad_c_fn(
            uc_in,
            vc_in,
            delpc,
            pkc,
            gz,
            rdxc,
            rdyc,
            hydrostatic,
            dt2,  # TODO: add [0, 0, 0] when gt4py bug is fixed
        )


@gtstencil()
def p_grad_c_ustencil(
    uc_in: FloatField,
    delpc: FloatField,
    pkc: FloatField,
    gz: FloatField,
    rdxc: FloatField,
    *,
    hydrostatic: int,
    dt2: float,
):
    with computation(PARALLEL), interval(0, -1):
        uc_in = p_grad_c_u_wk(
            uc_in, delpc, pkc, gz, rdxc, hydrostatic, dt2
        )  # TODO: add [0, 0, 0] when gt4py bug is fixed


@gtstencil()
def p_grad_c_vstencil(
    vc_in: FloatField,
    delpc: FloatField,
    pkc: FloatField,
    gz: FloatField,
    rdyc: FloatField,
    hydrostatic: int,
    dt2: float,
):
    with computation(PARALLEL), interval(0, -1):
        vc_in = p_grad_c_v_wk(
            vc_in, delpc, pkc, gz, rdyc, hydrostatic, dt2
        )  # TODO: add [0, 0, 0] when gt4py bug is fixed


def compute(
    uc: FloatField,
    vc: FloatField,
    delpc: FloatField,
    pkc: FloatField,
    gz: FloatField,
    dt2: float,
    **kwargs: Optional[Any],
):
    if "dfg" not in kwargs:
        kwargs["dfg"] = dict()

    grid = spec.grid
    p_grad_c_ustencil(
        uc,
        delpc,
        pkc,
        gz,
        grid.rdxc,
        hydrostatic=int(spec.namelist.hydrostatic),
        dt2=dt2,
        origin=grid.compute_origin(),
        domain=(grid.nic + 1, grid.njc, grid.npz + 1),
        **kwargs,
    )
    p_grad_c_vstencil(
        vc,
        delpc,
        pkc,
        gz,
        grid.rdyc,
        hydrostatic=int(spec.namelist.hydrostatic),
        dt2=dt2,
        origin=grid.compute_origin(),
        domain=(grid.nic, grid.njc + 1, grid.npz + 1),
        **kwargs,
    )
    dfg: Dict[str, Set[Any]] = kwargs["dfg"]
    stop=1
