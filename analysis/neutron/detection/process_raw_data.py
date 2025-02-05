import numpy as np
from libra_toolbox.neutron_detection.diamond.process_data import DataProcessor
from pathlib import Path


def process_and_save_data(folder, output_file):
    """
    Add all csv files from a folder to the data processor
    """
    data_proc = DataProcessor()
    for filename in Path(folder).rglob("*.CSV"):
        data_proc.add_file(
            filename, time_column=2, energy_column=3, delimiter=";", skip_header=1
        )
    np.save(output_file, np.array([data_proc.time_values, data_proc.energy_values]))


if __name__ == "__main__":
    path_to_data_folder = "../../../data/neutron_detection"

    # day 1 data part 1
    process_and_save_data(
        folder=f"{path_to_data_folder}/20241210_part1/UNFILTERED",
        output_file=f"{path_to_data_folder}/binary_data_day_1_part_1.npy",
    )

    # day 1 data part 2
    process_and_save_data(
        folder=f"{path_to_data_folder}/20241210_part2/UNFILTERED",
        output_file=f"{path_to_data_folder}/binary_data_day_1_part_2.npy",
    )

    # day 2 data
    process_and_save_data(
        folder=f"{path_to_data_folder}/20241211/UNFILTERED",
        output_file=f"{path_to_data_folder}/binary_data_day_2.npy",
    )
