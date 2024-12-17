from typing import Dict
from datetime import datetime
import re
from typing import List, Dict, Tuple
from libra_toolbox.tritium.lsc_measurements import (
    LSCFileReader,
    LIBRASample,
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
    try:
        background_label = "1L-BL-1"
        background_sample = LSCSample.from_file(file_reader, background_label)
    except ValueError:
        background_label = "1L-BL-2"
        background_sample = LSCSample.from_file(file_reader, background_label)

    # substract background
    sample.substract_background(background_sample)

    return sample


def create_gas_streams(
    samples: List[LSCSample], start_time: str | datetime, general_data: Dict
) -> Dict[str, GasStream]:

    def group_samples_by_stream(samples: List[LSCSample]) -> Dict[str, List[LSCSample]]:
        grouped_samples = {}
        for sample in samples:
            stream = extract_info_from_label(sample.name)["stream"]
            if stream not in grouped_samples:
                grouped_samples[stream] = []
            grouped_samples[stream].append(sample)
        return grouped_samples

    def group_samples_by_number(samples: List[LSCSample]) -> Dict[str, List[LSCSample]]:
        grouped_samples = {}
        for sample in samples:
            sample_nb = extract_info_from_label(sample.name)["sample_nb"]
            if sample_nb not in grouped_samples:
                grouped_samples[sample_nb] = []
            grouped_samples[sample_nb].append(sample)
        return grouped_samples

    def create_libra_samples(
        grouped_samples: Dict[str, List[LSCSample]], stream: str
    ) -> List[LIBRASample]:
        libra_samples = []
        for sample_nb, samples in grouped_samples.items():
            stream_dict = general_data["timestamps"]["lsc_sample_times"][stream]
            time_sample = datetime.strptime(
                stream_dict[f"{run_nb}-{sample_nb}-x"]["actual"],
                "%m/%d/%Y %H:%M",
            )
            libra_samples.append(LIBRASample(samples, time=time_sample))
        return libra_samples

    grouped_by_stream = group_samples_by_stream(samples)

    run_nb = extract_info_from_label(samples[0].name)["run_nb"]
    assert all(
        extract_info_from_label(sample.name)["run_nb"] == run_nb for sample in samples
    ), "All samples should have the same run number"

    gas_streams = {}
    for stream, stream_samples in grouped_by_stream.items():
        grouped_by_number = group_samples_by_number(stream_samples)
        libra_samples = create_libra_samples(grouped_by_number, stream)
        gas_streams[stream] = GasStream(libra_samples, start_time=start_time)

    return gas_streams


def extract_info_from_label(label: str, pattern: str = None) -> Dict[str, str]:
    # Define a regex pattern for the label format
    if not pattern:
        pattern = r"^(?P<volume>\d+[A-Za-z]*)-(?P<stream>[A-Za-z]+)_(?P<run_nb>\d+)-(?P<sample_nb>\d+)-(?P<vial_nb>\d+)$"

    # Match the label against the pattern
    match = re.match(pattern, label)

    if not match:
        raise ValueError(
            f"Label '{label}' does not match the expected format. "
            f"Pattern used: '{pattern}'"
        )

    # Return the matched components as a dictionary
    return match.groupdict()
