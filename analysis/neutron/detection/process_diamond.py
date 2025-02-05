import numpy as np
import matplotlib.pyplot as plt
from libra_toolbox.neutron_detection.diamond.process_data import DataProcessor
from datetime import datetime


def add_files_from_npy(filename, data_proc):
    data = np.load(filename)
    data_proc.time_values = np.append(data_proc.time_values, data[0])
    data_proc.energy_values = np.append(data_proc.energy_values, data[1])


# day 1 data part 1
start_day_1_part_1 = datetime.strptime("12/10: 10:53am", "%m/%d: %I:%M%p")
data_day_1_part_1 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_1.npy",
    data_day_1_part_1,
)

# day 1 data part 2
start_day_1_part_2 = datetime.strptime("12/10: 7:52pm", "%m/%d: %I:%M%p")

# convert start time to seconds after start of part 1
start_day_1_part_2 = (start_day_1_part_2 - start_day_1_part_1).total_seconds()
data_proc_day_1_part2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_2.npy",
    data_proc_day_1_part2,
)
data_proc_day_1_part2.time_values += start_day_1_part_2

# day 2 data
start_day_2 = datetime.strptime("12/11: 10:19am", "%m/%d: %I:%M%p")
start_day_2 = (start_day_2 - start_day_1_part_1).total_seconds()
data_proc_day_2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_2.npy",
    data_proc_day_2,
)
data_proc_day_2.time_values += start_day_2


s_to_h = 1 / 3600

for i, data_proc in enumerate(
    [data_day_1_part_1, data_proc_day_1_part2, data_proc_day_2]
):
    rates, bins = data_proc.get_count_rate(bin_time=2)
    plt.plot(
        bins[:-1] * s_to_h,
        rates,
        label=data_proc,
    )
    plt.fill_between(
        bins[:-1] * s_to_h,
        rates,
        alpha=0.3,
    )
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.ylabel("Counts per second")

plt.xlabel("Time (h)")

plt.savefig("out.png")
