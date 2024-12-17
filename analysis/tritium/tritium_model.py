from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
import numpy as np
import json
from libra_toolbox.tritium.lsc_measurements import LIBRARun
from helpers import create_sample, create_gas_streams
from datetime import datetime

lsc_data_folder = "../../data/tritium_detection"
with open("../../data/general.json", "r") as f:
    general_data = json.load(f)


# Make samples
with open(f"{lsc_data_folder}/data.json", "r") as f:
    lsc_data = json.load(f)
list_of_samples = [
    create_sample(sample["label"], f"{lsc_data_folder}/{sample["filename"]}")
    for sample in lsc_data["samples"]
]

# read start time from general.json
all_start_times = []
for generator in general_data["timestamps"]["generators"]:
    start_time = datetime.strptime(generator["on"], "%m/%d/%Y %H:%M")
    all_start_times.append(start_time)
start_time = min(all_start_times)

streams = create_gas_streams(
    samples=list_of_samples,
    start_time=start_time,
    general_data=general_data,
)
IV_stream = streams["IV"]
# OV_stream = streams["OV"]  # TODO add this back

# create run
run = LIBRARun(streams=[IV_stream], start_time=start_time)  # TODO add , OV_stream

# check that background is always substracted
for stream in run.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}"


replacement_times_top = sorted(IV_stream.relative_times_as_pint)
replacement_times_walls = (
    [] * ureg.h
)  # sorted(OV_stream.relative_times_as_pint)# TODO add this back


# tritium model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# read irradiation times from general.json


irradiations = []
for generator in general_data["timestamps"]["generators"]:
    irr_start_time = datetime.strptime(generator["on"], "%m/%d/%Y %H:%M") - start_time
    irr_stop_time = datetime.strptime(generator["off"], "%m/%d/%Y %H:%M") - start_time
    irr_start_time = irr_start_time.total_seconds() * ureg.second
    irr_stop_time = irr_stop_time.total_seconds() * ureg.second
    irradiations.append([irr_start_time, irr_stop_time])

# Neutron rate
# calculated from Kevin's activation foil analysis from run 100 mL #7
# TODO replace for values for this run
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
neutron_rate = (
    P383_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)

# TBR from OpenMC

from pathlib import Path

filename = "../neutron/statepoint.100.h5"
filename = Path(filename)

if not filename.exists():
    raise FileNotFoundError(f"{filename} does not exist, run OpenMC first")

import openmc

sp = openmc.StatePoint(filename)
tally_df = sp.get_tally(name="TBR").get_pandas_dataframe()
calculated_TBR = tally_df["mean"].iloc[0] * ureg.particle * ureg.neutron**-1
calculated_TBR_std_dev = (
    tally_df["std. dev."].iloc[0] * ureg.particle * ureg.neutron**-1
)

# TBR from measurements

total_irradiation_time = sum([irr[1] - irr[0] for irr in irradiations])

T_consumed = neutron_rate * total_irradiation_time
T_produced = (
    IV_stream.get_cumulative_activity("total")[-1]
    # + OV_stream.get_cumulative_activity("total")[-1]  # TODO add this back
)

measured_TBR = (T_produced / quantity_to_activity(T_consumed)).to(
    ureg.particle * ureg.neutron**-1
)

optimised_ratio = 1.7e-2
k_top = 8.9e-8 * ureg.m * ureg.s**-1
k_wall = optimised_ratio * k_top


baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,  # TODO replace by measured_TBR
    neutron_rate=neutron_rate,
    irradiations=irradiations,
    k_top=k_top,
    k_wall=k_wall,
)
