import pandas as pd
import numpy as np
from typing import List, Dict, Union
import pint
import pint.facets
from libra_toolbox.tritium import ureg
from datetime import datetime, timedelta
import warnings
from typing import Literal
from collections import defaultdict

# config
DATE_FORMAT = "%m/%d/%Y %I:%M %p"


class LSCFileReader:
    quench_set: str | None
    data: pd.DataFrame | None

    def __init__(
        self,
        file_path: str,
        vial_labels: List[str | None] = None,
        labels_column: str = None,
    ):
        """Reads a LSC file and extracts the Bq:1 values and labels

        Args:
            file_path: Path to the LSC file.
            vial_labels: List of vial labels. Defaults to None.
            labels_column: Column name in the file that contains the vial labels. Defaults to None.
        """
        self.file_path = file_path
        self.vial_labels = vial_labels
        self.labels_column = labels_column
        self.data = None
        self.header_content = None
        self.quench_set = None

    def read_file(self):
        """
        Reads the LSC file and extracts the data in self.data and the vial labels
        (if provided) in self.vial_labels

        Raises:
            ValueError: If both vial_labels and labels_column are provided or if none of them are provided
        """

        # check if vial_labels or labels_column is provided
        if (self.labels_column is None and self.vial_labels is None) or (
            self.labels_column is not None and self.vial_labels is not None
        ):
            raise ValueError("Provide either vial_labels or labels_column")

        # first read the file without dataframe to find the line starting with S#
        header_lines = []
        with open(self.file_path, "r") as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if line.startswith("Quench Set:"):
                    quench_set_line_idx = i + 1
                if line.startswith("S#"):
                    start = i
                    break
                header_lines.append(line)
        self.header_content = "".join(header_lines)
        self.quench_set = lines[quench_set_line_idx].strip()
        # read the file with dataframe starting from the line with S#
        self.data = pd.read_csv(self.file_path, skiprows=start)

        # check if last column is all NaN
        if self.data[self.data.columns[-1]].isnull().all():
            warnings.warn(
                "There seem to be an issue with the last column. Is the format of the file correct?"
            )

        if self.labels_column is not None:
            self.vial_labels = self.data[self.labels_column].tolist()

    def get_bq1_values(self) -> List[float]:
        assert self.data is not None
        return self.data["Bq:1"].tolist()
    
    def get_bq1_values_with_labels(self) -> Dict[str | None, Union[float, List[float]]]:
        """Returns a dictionary with vial labels as keys and Bq:1 values as values.
    
        If all labels have exactly one associated value, returns Dict[str | None, float].
        If at least one label has multiple values, returns Dict[str | None, List[float]].

        Raises:
            ValueError: If vial labels are not provided.

        Returns:
            Dict[str | None, Union[float, List[float]]]: Dictionary with vial labels as keys 
            and either single float values or, for repeat counts, lists of values.
        """
        if self.vial_labels is None:
            raise ValueError("Vial labels must be provided")

        assert len(self.vial_labels) == len(
            self.get_bq1_values()
        ), "Vial labels and Bq:1 values are not equal in length, remember to give None as a label for missing vials"
        
        # Collect values in a defaultdict (list) to handle multiple occurrences
        values = self.get_bq1_values()
        labelled_values = defaultdict(list)
        for label, value in zip(self.vial_labels, values):
            labelled_values[label].append(value)

        # Convert lists to single float values if all labels have exactly one value
        if all(len(v) == 1 for v in labelled_values.values()):
            return {k: v[0] for k, v in labelled_values.items()}  # Dict[str | None, float]

        return dict(labelled_values)  # Dict[str | None, List[float]]

    def get_count_times(self) -> List[float]:
        assert self.data is not None
        return self.data["Count Time"].tolist()

    def get_lum(self) -> List[float]:
        assert self.data is not None
        return self.data["LUM"].tolist()


class LSCSample:
    activity: pint.Quantity
    stdev: pint.Quantity
    origin_file: str | None

    def __init__(self, activity: pint.Quantity | List[pint.Quantity], name: str):
        self.name = name
        self.repeated = isinstance(activity.magnitude, np.ndarray) and activity.magnitude.size > 1
        if self.repeated:
            self.activity = np.mean(activity)
            self.stdev = np.std(activity, ddof=1)
        else:
            self.activity = activity
            self.stdev = 0.0 * ureg.Bq
        # TODO add other attributes available in LSC file
        self.background_substracted = False
        self.origin_file = None
    
    def __str__(self):
        return f"Sample {self.name}"
    
    @staticmethod
    def stdev_addition(stdevs: List[pint.Quantity]) -> pint.Quantity:
        """Calculates standard deviation of summed variables

        Args:
            stdevs (List[pint.Quantity]): List of standard deviations of individual variables
        """
        result = np.sqrt(np.sum(std**2 for std in stdevs))
        return result

    def substract_background(self, background_sample: "LSCSample"):
        """Substracts the background activity from the sample activity

        Args:
            background_sample (LSCSample): Background sample

        Raises:
            ValueError: If background has already been substracted
        """
        if self.background_substracted:
            raise ValueError("Background already substracted")
        self.activity -= background_sample.activity
        if self.activity.magnitude < 0:
            warnings.warn(
                f"Activity of {self.name} is negative after substracting background. Setting to zero."
            )
            self.activity = 0 * ureg.Bq
        self.stdev = self.stdev_addition([self.stdev, background_sample.stdev])
        self.background_substracted = True

    @staticmethod
    def from_file(file_reader: LSCFileReader, vial_name: str) -> "LSCSample":
        """Creates an LSCSample object from a LSCFileReader object

        Args:
            file_reader (LSCFileReader): LSCFileReader object
            vial_name: Name of the vial

        Raises:
            ValueError: If vial_name is not found in the file reader

        Returns:
            the LSCSample object
        """
        values = file_reader.get_bq1_values_with_labels()
        if vial_name not in values:
            raise ValueError(f"Vial {vial_name} not found in the file reader.")
        activity = values[vial_name] * ureg.Bq
        sample = LSCSample(activity, vial_name)
        sample.origin_file = file_reader.file_path
        return sample


class LIBRASample:
    samples: List[LSCSample]

    def __init__(self, samples: List[LSCSample], time: str | datetime):
        self.samples = samples

        if isinstance(time, str):
            self._time = datetime.strptime(time, DATE_FORMAT)
        else:
            self._time = time

    def get_relative_time(self, start_time: str | datetime) -> timedelta:

        if isinstance(start_time, str):
            start_time = datetime.strptime(start_time, DATE_FORMAT)

        sample_time = self._time
        return sample_time - start_time

    def substract_background(self, background_sample: LSCSample):
        for sample in self.samples:
            sample.substract_background(background_sample)

    def get_soluble_activity(self):
        """Returns total activity of soluble samples
        If counts are repeated, returns standard deviation as well
        Backwards compatible with single counts

        Returns:
            tuple containing activity and standard deviation if repeated,
            activity only otherwise
        """
        is_repeated = any(
            sample.repeated
            for sample in self.samples
        )
        act = 0
        stdev = 0
        for sample in self.samples[:2]:
            act += sample.activity
            stdev = LSCSample.stdev_addition([stdev, sample.stdev])

        return (act, stdev) if is_repeated else act

    def get_insoluble_activity(self):
        """Returns total activity of insoluble samples
        If counts are repeated, returns standard deviation as well
        Backwards compatible with single counts
        
        Returns:
            tuple containing activity and standard deviation if repeated,
            activity only otherwise
        """
        is_repeated = any(
            sample.repeated
            for sample in self.samples
        )
        act = 0
        stdev = 0
        for sample in self.samples[2:]:
            act += sample.activity
            stdev = LSCSample.stdev_addition([stdev, sample.stdev])

        return (act, stdev) if is_repeated else act

    def get_total_activity(self):
        """Returns total activity of samples
        If counts are repeated, returns standard deviation as well
        Backwards compatible with single counts

        Returns:
            tuple containing activity and standard deviation if repeated,
            activity only otherwise
        """
        is_repeated = any(
            sample.repeated
            for sample in self.samples
        )
        act_s = 0 * ureg.Bq
        act_i = 0 * ureg.Bq
        stdev_s = 0 * ureg.Bq
        stdev_i = 0 * ureg.Bq
        act = 0 * ureg.Bq
        stdev = 0 * ureg.Bq

        if is_repeated:
            act_s, stdev_s = self.get_soluble_activity()
            act_i, stdev_i = self.get_insoluble_activity()
            act = act_s + act_i
            stdev = LSCSample.stdev_addition([stdev_s, stdev_i])
        else:
            act = self.get_soluble_activity() + self.get_insoluble_activity()
        
        return (act, stdev) if is_repeated else act


class GasStream:
    samples: List[LIBRASample]
    name: str | None

    def __init__(
        self,
        samples: List[LIBRASample],
        start_time: str | datetime,
        name: str | None = None,
    ):
        """Creates a GasStream object from a list of LIBRASample objects

        Args:
            samples: a list of LIBRASample objects
            start_time: the start time at which the gas stream was collected
            name: a name for the GasStream. Defaults to None.
        """
        self.samples = samples
        if isinstance(start_time, str):
            self.start_time = datetime.strptime(start_time, DATE_FORMAT)
        else:
            self.start_time = start_time

        self.name = name


    @staticmethod
    def append_stdev(stdev_list: List, stdev_val: pint.Quantity) -> List:
        """If the list is empty, appends value
        Otherwise uses LSCSample.stdev_addition to calculate the stdev when new value is added, then appends
        """
        if not stdev_list:
            stdev_list.append(stdev_val)
        else:
            stdev_list.append(LSCSample.stdev_addition([stdev_list[-1], stdev_val]))

        return stdev_list

    def get_cumulative_activity(
        self, form: Literal["total", "soluble", "insoluble"] = "total"
    ):
        """Calculates the cumulative activity of the gas stream

        Args:
            form: the form in which the cumulative activity is calculated.

        Raises:
            ValueError: If background has not been substracted

        Returns:
            if counts repeated:
                tuple containing the cumulative activity and standard deviations as pint.Quantity object
            else:
                cumulative activity as a pint.Quantity
        """
        # check if counts are repeated
        is_repeated = any(
            lsc_sample.repeated
            for sample in self.samples
            for lsc_sample in sample.samples
        )
        # check that background has been substracted
        if any(
            not lsc_sample.background_substracted
            for sample in self.samples
            for lsc_sample in sample.samples
        ):
            raise ValueError(
                "Background must be substracted before calculating cumulative activity"
            )
        cumulative_activity = []
        stdevs = []
        if is_repeated:
            for sample in self.samples:
                if form == "total":
                    act, stdev = sample.get_total_activity()
                    cumulative_activity.append(act)
                    self.append_stdev(stdevs, stdev)
                elif form == "soluble":
                    act, stdev = sample.get_soluble_activity()
                    cumulative_activity.append(act)
                    self.append_stdev(stdevs, stdev)
                elif form == "insoluble":
                    act, stdev = sample.get_insoluble_activity()
                    cumulative_activity.append(act)
                    self.append_stdev(stdevs, stdev)
        else:
            for sample in self.samples:
                if form == "total":
                    cumulative_activity.append(sample.get_total_activity())
                elif form == "soluble":
                    cumulative_activity.append(sample.get_soluble_activity())
                elif form == "insoluble":
                    cumulative_activity.append(sample.get_insoluble_activity())
        cumulative_activity = ureg.Quantity.from_list(cumulative_activity)
        if stdevs:
            stdevs = ureg.Quantity.from_list(stdevs)
        cumulative_activity = cumulative_activity.cumsum()
        return (cumulative_activity, stdevs) if is_repeated else cumulative_activity

    @property
    def relative_times(self) -> List[timedelta]:
        """
        The relative times of the samples in the GasStream based on the start_time
        as a list of timedelta objects
        """
        return [sample.get_relative_time(self.start_time) for sample in self.samples]

    @property
    def relative_times_as_pint(self) -> pint.facets.plain.PlainQuantity:
        """
        The relative times of the samples in the GasStream based on the start_time
        as a pint.Quantity object
        """
        times = [t.total_seconds() * ureg.s for t in self.relative_times]
        return ureg.Quantity.from_list(times).to(ureg.day)


class LIBRARun:
    def __init__(self, streams: List[GasStream], start_time: str | datetime):
        self.streams = streams
        if isinstance(start_time, str):
            self.start_time = datetime.strptime(start_time, DATE_FORMAT)
        else:
            self.start_time = start_time


class BABY100mLRun(LIBRARun):
    def __init__(self, inner_vessel_stream, start_time):
        super().__init__([inner_vessel_stream], start_time)


class BABY1LRun(LIBRARun):
    def __init__(self, inner_vessel_stream, outer_vessel_stream, start_time):
        super().__init__([inner_vessel_stream, outer_vessel_stream], start_time)