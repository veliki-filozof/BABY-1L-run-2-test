import numpy as np
import matplotlib.pyplot as plt
from libra_toolbox.neutron_detection.diamond.process_data import DataProcessor
from datetime import datetime
import json


def add_files_from_npy(filename, data_proc):
    data = np.load(filename)
    data_proc.time_values = np.append(data_proc.time_values, data[0])
    data_proc.energy_values = np.append(data_proc.energy_values, data[1])


# day 1 data part 1
start_day_1_part_1 = datetime.strptime("2024/12/10: 10:53am", "%Y/%m/%d: %I:%M%p")
data_day_1_part_1 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_1.npy",
    data_day_1_part_1,
)

# day 1 data part 2
start_day_1_part_2 = datetime.strptime("2024/12/10: 7:52pm", "%Y/%m/%d: %I:%M%p")

# convert start time to seconds after start of part 1
start_day_1_part_2 = (start_day_1_part_2 - start_day_1_part_1).total_seconds()
data_proc_day_1_part2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_1_part_2.npy",
    data_proc_day_1_part2,
)
data_proc_day_1_part2.time_values += start_day_1_part_2

# day 2 data
start_day_2 = datetime.strptime("2024/12/11: 10:19am", "%Y/%m/%d: %I:%M%p")
start_day_2 = (start_day_2 - start_day_1_part_1).total_seconds()
data_proc_day_2 = DataProcessor()
add_files_from_npy(
    "binary_data_day_2.npy",
    data_proc_day_2,
)
data_proc_day_2.time_values += start_day_2


# -------- read irradiations from logbook ----------------

with open("../../../data/general.json", "r") as f:
    general_data = json.load(f)

start_time = start_day_1_part_1


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

        lag = (
            30 * 60
        )  # 30 minutes lag due to the time it takes to save one file in COMPASS.
        irr_start_time = irr_start_time.total_seconds() + lag
        irr_stop_time = irr_stop_time.total_seconds() + lag

        irradiations.append([irr_start_time, irr_stop_time])


# -------- plot data ----------------

s_to_h = 1 / 3600
plt.figure(figsize=(13, 5))
for i, data_proc in enumerate(
    [data_day_1_part_1, data_proc_day_1_part2, data_proc_day_2]
):
    rates, bins = data_proc.get_count_rate(bin_time=2)
    plt.plot(bins[:-1] * s_to_h, rates, color="tab:grey", linewidth=0.1)
    plt.fill_between(bins[:-1] * s_to_h, rates, alpha=0.3, color="tab:grey")

for i, irr in enumerate(irradiations):
    plt.fill_betweenx(
        [1600, 1800],
        irr[0] * s_to_h,
        irr[1] * s_to_h,
        color="#EF5B5B",
        alpha=0.3,
    )
plt.annotate(
    "from logbook",
    (7, 1600),
    xytext=(8, 1700),
    ha="center",
    va="center",
    color="#EF5B5B",
)


plt.hlines(y=1000, xmin=5, xmax=9, linestyles="dashed", color="black")
plt.annotate("no aquisition", (7, 1000), xytext=(7, 1100), ha="center", va="center")

plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.ylabel("Counts per second")

plt.xlabel("Time (h)")

plt.ylim(bottom=0)

plt.savefig("out.png")
