from typing import Dict
from datetime import datetime
from typing import List, Dict, Tuple
from libra_toolbox.tritium.lsc_measurements import (
    LSCFileReader,
    LIBRASample,
    LIBRARun,
    LSCSample,
    GasStream,
)

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
    background_label = "1L-BL-1"
    background_sample = LSCSample.from_file(file_reader, background_label)

    # substract background
    sample.substract_background(background_sample)

    return sample


def create_gas_streams(
    samples: List[LSCSample], start_time: str | datetime, general_data: Dict
) -> Tuple[GasStream, GasStream]:
    # group samples by stream IV and OV based on the labels (eg. 1L-IV-..., 1L-OV-...)
    all_IV_samples, all_OV_samples = [], []

    for sample in samples:
        stream = extract_info_from_label(sample.name)["stream"]
        if stream == "IV":
            all_IV_samples.append(sample)
        elif stream == "OV":
            all_OV_samples.append(sample)
        else:
            raise ValueError(
                f"Sample {sample} does not have a valid stream label {sample.name}"
            )

    # group the samples by LIBRA sample number (eg. 1L-IV-1-1-..., 1L-IV-1-2...)
    run_nb = extract_info_from_label(all_IV_samples[0].name)["run_nb"]
    assert all(
        extract_info_from_label(sample.name)["run_nb"] == run_nb
        for sample in all_IV_samples
    ), "All samples should have the same run number"

    IV_samples, OV_samples = {}, {}

    for sample in all_IV_samples:
        # find the LIBRA sample number the format is the following: 1L-IV-[run nb]-[sample nb]-[vial nb]
        # we want to group them by sample nb
        sample_nb = extract_info_from_label(sample.name)["sample_nb"]
        # if the sample number does not exist, create a new LIBRA sample
        if sample_nb not in IV_samples:
            IV_samples[sample_nb] = []

        # add the sample to the LIBRA sample
        IV_samples[sample_nb].append(sample)

    for sample in all_OV_samples:
        # find the LIBRA sample number the format is the following: 1L-OV-[run nb]-[sample nb]-[vial nb]
        # we want to group them by sample nb
        sample_nb = extract_info_from_label(sample.name)["sample_nb"]
        # if the sample number does not exist, create a new LIBRA sample
        if sample_nb not in OV_samples:
            OV_samples[sample_nb] = []

        # add the sample to the LIBRA sample
        OV_samples[sample_nb].append(sample)

    # make LIBRASample objects
    IV_samples_as_libra_samples, OV_samples_as_libra_samples = [], []
    for sample_nb, samples in IV_samples.items():
        time_sample = datetime.strptime(
            general_data["timestamps"]["lsc_sample_times"]["IV"][
                f"{run_nb}-{sample_nb}-x"
            ]["actual"],
            "%m/%d/%Y %H:%M",
        )
        IV_samples_as_libra_samples.append(LIBRASample(samples, time=time_sample))

    for sample_nb, samples in OV_samples.items():
        time_sample = datetime.strptime(
            general_data["timestamps"]["lsc_sample_times"]["OV"][
                f"{run_nb}-{sample_nb}-x"
            ]["actual"],
            "%m/%d/%Y %H:%M",
        )
        OV_samples_as_libra_samples.append(LIBRASample(samples, time=time_sample))

    # return IV and OV gas streams

    IV_stream = GasStream(
        IV_samples_as_libra_samples,
        start_time=start_time,
    )
    OV_stream = GasStream(
        OV_samples_as_libra_samples,
        start_time=start_time,
    )
    return IV_stream, OV_stream


def extract_info_from_label(label: str) -> Dict[str, str]:
    # check that it has the following format
    # [volume]-[stream]_[run nb]-[sample nb]-[vial nb]
    if len(label.split("-")) != 4 or len(label.split("_")) != 2:
        raise ValueError(
            f"Label {label} does not have the correct format [volume]-[stream]_[run nb]-[sample nb]-[vial nb] (eg. 1L-IV_1-1-1)"
        )
    # extract volume, stream, run_nb, sample_nb, vial_nb from the label
    volume_stream, run_sample_vial = label.split("_")
    volume, stream = volume_stream.split("-")
    run_nb, sample_nb, vial_nb = run_sample_vial.split("-")
    return {
        "volume": volume,
        "stream": stream,
        "run_nb": run_nb,
        "sample_nb": sample_nb,
        "vial_nb": vial_nb,
    }
