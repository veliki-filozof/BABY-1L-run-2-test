import numpy as np
from libra_toolbox.neutron_detection.diamond.process_data import DataProcessor
from pathlib import Path

path_to_data_folder = "../../../data/neutron_detection"


def add_files_from_folder(folder, data_proc: DataProcessor):
    """
    Add all csv files from a folder to the data processor
    """
    for filename in Path(folder).rglob("*.CSV"):
        data_proc.add_file(
            filename, time_column=2, energy_column=3, delimiter=";", skip_header=1
        )


if __name__ == "__main__":
    # day 1 data part 1
    data_day_1_part_1 = DataProcessor()
    add_files_from_folder(
        data_proc=data_day_1_part_1,
        folder=f"{path_to_data_folder}/Single diamond/20241210_LIBRA/UNFILTERED/",
    )
    np.save(
        f"{path_to_data_folder}/binary_data_day_1_part_1.npy",
        np.array([data_day_1_part_1.time_values, data_day_1_part_1.energy_values]),
    )

    # day 1 data part 2
    data_proc_day_1_part2 = DataProcessor()
    add_files_from_folder(
        data_proc=data_proc_day_1_part2,
        folder=f"{path_to_data_folder}/Single diamond/20241210_LIBRA_part2/UNFILTERED/",
    )
    np.save(
        f"{path_to_data_folder}/binary_data_day_1_part_2.npy",
        np.array(
            [data_proc_day_1_part2.time_values, data_proc_day_1_part2.energy_values]
        ),
    )

    # day 2 data
    data_proc_day_2 = DataProcessor()
    add_files_from_folder(
        data_proc=data_proc_day_2,
        folder=f"{path_to_data_folder}/Single diamond/20241211_LIBRA/DAQ/LIBRA_run2/UNFILTERED",
    )
    np.save(
        f"{path_to_data_folder}/binary_data_day_2.npy",
        np.array([data_proc_day_2.time_values, data_proc_day_2.energy_values]),
    )
