from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
import numpy as np
import json
from lsc_measurements_mod import (
    LIBRARun,
    LSCFileReader,
    GasStream,
    LSCSample,
    LIBRASample,
)

from datetime import datetime


all_file_readers = []
all_quench = []

def create_sample(
        label: str, 
        filename: str, 
        background_filename: str = None, 
        background_label: str = None
        ) -> LSCSample:
    """
    Create a LSCSample from a LSC file with background substracted.

    Args:
        label: the label of the sample in the LSC file
        filename: the filename of the LSC file
        background_filename: the filename of the LSC file containing background count
        background_label: label of the background vial

    Returns:
        the LSCSample object
    """

    if background_filename is None:
        background_filename = filename

    # check if a LSCFileReader has been created for the filename and background_filename
    found = False
    found_background = False
    for file_reader in all_file_readers:
        if file_reader.filename == filename:
            found = True
            file_reader_main = file_reader
        if file_reader.filename == background_filename:
            found_background = True
            file_reader_background = file_reader
        if found and found_background:
            break
    file_reader = file_reader_main

    # if not, create it or them and add to the list of LSCFileReaders
    if not found:
        file_reader = LSCFileReader(filename, labels_column="SMPL_ID")
    if not found_background:
        file_reader_background = LSCFileReader(background_filename, labels_column="SMPL_ID")

    file_reader.read_file()
    file_reader_background.read_file()

    # create the sample
    sample = LSCSample.from_file(file_reader, label)

   # try to find the background sample from the file
    if background_label is None:
        background_labels = ["1L-BL-1", "1L-BL-2", "1L-BL-3"]
        background_sample = None

        for background_label in background_labels:
            try:
                background_sample = LSCSample.from_file(file_reader_background, background_label)
                break
            except ValueError:
                continue

        if background_sample is None:
            raise ValueError(f"Background sample not found in {background_filename}")
    else:
        try:
            background_sample = LSCSample.from_file(file_reader_background, background_label)
        except ValueError as e:
            raise ValueError(f"Background sample '{background_label}' not found in {background_filename}: {str(e)}")



    # substract background
    sample.substract_background(background_sample)

    # read quench set
    all_quench.append(file_reader.quench_set)

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


# create gas streams
gas_streams = {}
gas_streams_repeat = {}
for stream, samples in general_data["tritium_detection"].items():
    stream_samples = []
    stream_samples_repeat = []
    for sample_nb, sample_dict in samples.items():
        libra_samples = []
        libra_samples_repeat = []
        if sample_dict["actual_sample_time"] is None:
            continue
        for vial_nb, filename in sample_dict["lsc_vials_filenames"].items():
            sample = create_sample(
                label=f"1L-{stream}_{run_nb}-{sample_nb}-{vial_nb}",
                filename=f"{lsc_data_folder}/{filename}",
            )
            libra_samples.append(sample)
        for vial_nb, filename_repeat in sample_dict["repeat_count_filenames"].items():
            if filename_repeat is None:
                # Create a zero sample and tag it as background substracted
                empty_sample = LSCSample(
                    activity=0 * ureg.Bq,
                    name=f"1L-{stream}_{run_nb}-{sample_nb}-{vial_nb}"
                )
                empty_sample.background_substracted = True
                libra_samples_repeat.append(sample)
            else:
                sample = create_sample(
                    label=f"1L-{stream}_{run_nb}-{sample_nb}-{vial_nb}",
                    filename=f"{lsc_data_folder}/{filename_repeat}",
                    background_filename=f"{lsc_data_folder}/{sample_dict["repeat_background_filename"]}"
                )
                libra_samples_repeat.append(sample)

        time_sample = datetime.strptime(
            sample_dict["actual_sample_time"], "%m/%d/%Y %H:%M"
        )
        stream_samples.append(LIBRASample(libra_samples, time=time_sample))
        stream_samples_repeat.append(LIBRASample(libra_samples_repeat, time=time_sample))
    gas_streams[stream] = GasStream(stream_samples, start_time=start_time)
    gas_streams_repeat[stream] = GasStream(stream_samples_repeat, start_time=start_time)


# create run
run = LIBRARun(streams=list(gas_streams.values()), start_time=start_time)
run_repeat = LIBRARun(streams=list(gas_streams_repeat.values()), start_time=start_time)

# check that only one quench set is used
assert len(np.unique(all_quench)) == 1

# check that background is always substracted  # TODO this should be done automatically in LIBRARun
for stream in run.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}"
for stream in run_repeat.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}, repeat count"

IV_stream = gas_streams["IV"]
OV_stream = gas_streams["OV"]

IV_stream_repeat = gas_streams_repeat["IV"]
OV_stream_repeat = gas_streams_repeat["OV"]

sampling_times = {
    "IV": sorted(IV_stream.relative_times_as_pint),
    "OV": sorted(OV_stream.relative_times_as_pint),
}

replacement_times_top = sampling_times["IV"]
replacement_times_walls = sampling_times["OV"]


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
neutron_rate_relative_uncertainty = 0.089  # TODO check with Collin what is the uncertainty on this measurement

neutron_rate = 2.611e+08 * ureg.neutron * ureg.s**-1  # TODO from Collin's foil analysis, replace with more robust method

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
T_produced_repeat = sum(
    [stream.get_cumulative_activity("total")[0][-1] for stream in run_repeat.streams]
)
T_produced_repeat_stdev = LSCSample.stdev_addition(
    [stream.get_cumulative_activity("total")[1][-1] for stream in run_repeat.streams]
)

measured_TBR = (T_produced / quantity_to_activity(T_consumed)).to(
    ureg.particle * ureg.neutron**-1
)
measured_TBR_repeat = (T_produced_repeat / quantity_to_activity(T_consumed)).to(
    ureg.particle * ureg.neutron**-1
)
measured_TBR_repeat_stdev = measured_TBR_repeat * np.sqrt(
    (T_produced_repeat_stdev / T_produced_repeat)**2 + neutron_rate_relative_uncertainty**2
).to(
    ureg.particle * ureg.neutron**-1
)

k_top = 0.7*8.9e-8 * ureg.m * ureg.s**-1
k_wall = 0 * ureg.m * ureg.s**-1


baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=measured_TBR,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
    k_top=k_top,
    k_wall=k_wall,
)
baby_model_repeat = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=measured_TBR_repeat,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
    k_top=k_top,
    k_wall=k_wall,
)

# store processed data
processed_data = {
    "modelled_baby_radius": {
        "value": baby_radius.magnitude,
        "unit": str(baby_radius.units),
    },
    "modelled_baby_height": {
        "value": baby_height.magnitude,
        "unit": str(baby_height.units),
    },
    "irradiations": [
        {
            "start_time": {
                "value": irr[0].magnitude,
                "unit": str(irr[0].units),
            },
            "stop_time": {
                "value": irr[1].magnitude,
                "unit": str(irr[1].units),
            },
        }
        for irr in irradiations
    ],
    "neutron_rate_used_in_model": {
        "value": baby_model.neutron_rate.magnitude,
        "unit": str(baby_model.neutron_rate.units),
    },
    "measured_TBR": {
        "value": measured_TBR.magnitude,
        "unit": str(measured_TBR.units),
    },
    "TBR_used_in_model": {
        "value": baby_model.TBR.magnitude,
        "unit": str(baby_model.TBR.units),
    },
    "k_top": {
        "value": baby_model.k_top.magnitude,
        "unit": str(baby_model.k_top.units),
    },
    "k_wall": {
        "value": baby_model.k_wall.magnitude,
        "unit": str(baby_model.k_wall.units),
    },
    "cumulative_tritium_release": {
        label: {
            **{
                form: {
                    "value": gas_stream.get_cumulative_activity(
                        form
                    ).magnitude.tolist(),
                    "unit": str(gas_stream.get_cumulative_activity(form).units),
                }
                for form in ["total", "soluble", "insoluble"]
            },
            "sampling_times": {
                "value": gas_stream.relative_times_as_pint.magnitude.tolist(),
                "unit": str(gas_stream.relative_times_as_pint.units),
            },
        }
        for label, gas_stream in gas_streams.items()
    },
}

# check if the file exists and load it

processed_data_file = "../../data/processed_data.json"

try:
    with open(processed_data_file, "r") as f:
        existing_data = json.load(f)
except FileNotFoundError:
    print(f"Processed data file not found, creating it in {processed_data_file}")
    existing_data = {}

existing_data.update(processed_data)

with open(processed_data_file, "w") as f:
    json.dump(existing_data, f, indent=4)

print(f"Processed data stored in {processed_data_file}")
