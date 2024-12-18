from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
import numpy as np
import json
from libra_toolbox.tritium.lsc_measurements import (
    LIBRARun,
    LSCFileReader,
    GasStream,
    LSCSample,
    LIBRASample,
)

from datetime import datetime


all_file_readers = []


def create_sample(label: str, filename: str) -> LSCSample:
    # check if a LSCFileReader has been created for this filename
    found = False
    for file_reader in all_file_readers:
        if file_reader.filename == filename:
            found = True
            break

    # if not, create it and add it to the list of LSCFileReaders
    if not found:
        file_reader = LSCFileReader(filename, labels_column="SMPL_ID")

    file_reader.read_file()

    # create the sample
    sample = LSCSample.from_file(file_reader, label)

    # try to find the background sample from the file
    try:
        background_label = "1L-BL-1"
        background_sample = LSCSample.from_file(file_reader, background_label)
    except ValueError:
        background_label = "1L-BL-2"
        background_sample = LSCSample.from_file(file_reader, background_label)

    # substract background
    sample.substract_background(background_sample)

    return sample


lsc_data_folder = "../../data/tritium_detection"
with open("../../data/general.json", "r") as f:
    general_data = json.load(f)

run_nb = general_data["general_data"]["run_nb"]


# read start time from general.json
all_start_times = []
for generator in general_data["generators"]:
    if generator["enabled"] is False:
        continue
    for irradiation_period in generator["periods"]:
        start_time = datetime.strptime(irradiation_period["start"], "%m/%d/%Y %H:%M")
        all_start_times.append(start_time)
start_time = min(all_start_times)


gas_streams = {}
for stream, samples in general_data["tritium_detection"].items():
    stream_samples = []
    for sample_nb, sample_dict in samples.items():
        libra_samples = []
        if sample_dict["actual_sample_time"] is None:
            continue
        for vial_nb, filename in sample_dict["lsc_vials_filenames"].items():
            sample = create_sample(
                label=f"1L-{stream}_{run_nb}-{sample_nb}-{vial_nb}",
                filename=f"{lsc_data_folder}/{filename}",
            )
            libra_samples.append(sample)

        time_sample = datetime.strptime(
            sample_dict["actual_sample_time"], "%m/%d/%Y %H:%M"
        )
        stream_samples.append(LIBRASample(libra_samples, time=time_sample))
    gas_streams[stream] = GasStream(stream_samples, start_time=start_time)


# create run
run = LIBRARun(streams=list(gas_streams.values()), start_time=start_time)

# check that background is always substracted
for stream in run.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}"

IV_stream = gas_streams["IV"]
OV_stream = gas_streams["OV"]

replacement_times_top = sorted(IV_stream.relative_times_as_pint)
replacement_times_walls = sorted(OV_stream.relative_times_as_pint)


# tritium model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# read irradiation times from general.json


irradiations = []
for generator in general_data["generators"]:
    if generator["enabled"] is False:
        continue
    for irradiation_period in generator["periods"]:
        irr_start_time = (
            datetime.strptime(irradiation_period["start"], "%m/%d/%Y %H:%M")
            - start_time
        )
        irr_stop_time = (
            datetime.strptime(irradiation_period["end"], "%m/%d/%Y %H:%M") - start_time
        )
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
T_produced = sum(
    [stream.get_cumulative_activity("total")[-1] for stream in run.streams]
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
