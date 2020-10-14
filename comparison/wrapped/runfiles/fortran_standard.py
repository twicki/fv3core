import mpi4py
import yaml
import fv3gfs.wrapper as wrapper
from fv3gfs.util import (
    io,
    QuantityFactory,
    SubtileGridSizer,
)

# May need to run 'ulimit -s unlimited' before running this example
# If you're running in our prepared docker container, you definitely need to do this
# sets the stack size to unlimited

# Run using mpirun -n 6 python3 basic_model.py
# mpirun flags that may be useful:
#     for docker:  --allow-run-as-root
#     for CircleCI: --oversubscribe
#     to silence a certain inconsequential MPI error: --mca btl_vader_single_copy_mechanism none

# All together:
# mpirun -n 6 --allow-run-as-root --oversubscribe --mca btl_vader_single_copy_mechanism none python3 basic_model.py

if __name__ == "__main__":

    # get another namelist for the communicator??
    nml2 = yaml.safe_load(
        open("/fv3core/comparison/wrapped/config/c12_6ranks_baroclinic.yml", "r")
    )["namelist"]

    sizer = SubtileGridSizer.from_namelist(nml2)
    allocator = QuantityFactory.from_backend(sizer, "numpy")

    # MPI stuff
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Set the names of quantities in State. This is everything coming from wrapper.initialize
    initial_names = [
        "cloud_fraction",
        "turbulent_kinetic_energy",
        "specific_humidity",
        "cloud_water_mixing_ratio",
        "rain_mixing_ratio",
        "snow_mixing_ratio",
        "cloud_ice_mixing_ratio",
        "graupel_mixing_ratio",
        "ozone_mixing_ratio",
        "air_temperature",
        "pressure_thickness_of_atmospheric_layer",
        "vertical_thickness_of_atmospheric_layer",
        "logarithm_of_interface_pressure",
        "x_wind",
        "y_wind",
        "vertical_wind",
        "x_wind_on_c_grid",
        "y_wind_on_c_grid",
        "total_condensate_mixing_ratio",
        "interface_pressure",
        "surface_geopotential",
        "interface_pressure_raised_to_power_of_kappa",
        "surface_pressure",
        "vertical_pressure_velocity",
        "atmosphere_hybrid_a_coordinate",
        "atmosphere_hybrid_b_coordinate",
        "accumulated_x_mass_flux",
        "accumulated_y_mass_flux",
        "accumulated_x_courant_number",
        "accumulated_y_courant_number",
        "dissipation_estimate_from_heat_source",
        "eastward_wind",
        "northward_wind",
        "layer_mean_pressure_raised_to_power_of_kappa",
        "time",
    ]

    wrapper.initialize()

    for i in range(wrapper.get_step_count()):
        print("STEP IS ", i)
        if i == 0:
            state = wrapper.get_state(allocator=allocator, names=initial_names)
            io.write_state(state, "instate_{0}.nc".format(rank))
        else:
            state = wrapper.get_state(allocator=allocator, names=initial_names)

        wrapper.step_dynamics()
        wrapper.step_physics()
        wrapper.save_intermediate_restart_if_enabled()

    state = wrapper.get_state(allocator=allocator, names=initial_names)
    io.write_state(state, "outstate_{0}.nc".format(rank))

    wrapper.cleanup()

