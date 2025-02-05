import numpy as np
import matplotlib.pyplot as plt
from libra_toolbox.neutron_detection.diamond.process_data import DataProcessor
from pathlib import Path


def add_files_from_folder(folder, data_proc):
    for filename in Path(folder).rglob("*.CSV"):
        data_proc.add_file(
            filename, time_column=2, energy_column=3, delimiter=";", skip_header=1
        )


def add_files_from_npy(filename, data_proc):
    data = np.load(filename)
    data_proc.time_values = np.append(data_proc.time_values, data[0])
    data_proc.energy_values = np.append(data_proc.energy_values, data[1])


# day 1 data part 1
data_day_1_part_1 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_1.npy",
    data_day_1_part_1,
)
# add_files_from_folder(
#     data_proc=data_day_1_part_1,
#     folder="../../../data/neutron_detection/Single diamond/20241210_LIBRA/UNFILTERED/",
# )
# np.save(
#     f"binary_data_day_1_part_1.npy",
#     np.array([data_day_1_part_1.time_values, data_day_1_part_1.energy_values]),
# )

# day 1 data part 2
data_proc_day_1_part2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_2.npy",
    data_proc_day_1_part2,
)

# add_files_from_folder(
#     data_proc=data_proc_day_1_part2,
#     folder="../../../data/neutron_detection/Single diamond/20241210_LIBRA_part2/UNFILTERED/",
# )
# np.save(
#     f"binary_data_day_1_part_2.npy",
#     np.array([data_proc_day_1_part2.time_values, data_proc_day_1_part2.energy_values]),
# )

# day 2 data
data_proc_day_2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_2.npy",
    data_proc_day_2,
)
# add_files_from_folder(
#     data_proc=data_proc_day_2,
#     folder="../../../data/neutron_detection/Single diamond/20241211_LIBRA/DAQ/LIBRA_run2/UNFILTERED",
# )
# np.save(
#     f"binary_data_day_2.npy",
#     np.array([data_proc_day_2.time_values, data_proc_day_2.energy_values]),
# )

s_to_h = 1 / 3600

fig, axs = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i, data_proc in enumerate(
    [data_day_1_part_1, data_proc_day_1_part2, data_proc_day_2]
):
    plt.sca(axs[i])
    rates, bins = data_proc.get_count_rate(bin_time=2)
    plt.plot(
        bins[:-1] * s_to_h,
        rates,
        label=data_proc,
    )
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.ylabel("Counts per second")


plt.xlabel("Time (h)")

plt.savefig("out.png")
