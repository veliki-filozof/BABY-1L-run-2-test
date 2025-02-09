import pandas as pd
from typing import List, Dict
import pint
import pint.facets
from libra_toolbox.tritium import ureg
from datetime import datetime, timedelta
import warnings
from typing import Literal


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

    def get_bq1_values_with_labels(self) -> Dict[str | None, float]:
        """Returns a dictionary with vial labels as keys and Bq:1 values as values

        Raises:
            ValueError: If vial labels are not provided

        Returns:
            Dictionary with vial labels as keys and Bq:1 values as values
        """
        if self.vial_labels is None:
            raise ValueError("Vial labels must be provided")

        assert len(self.vial_labels) == len(
            self.get_bq1_values()
        ), "Vial labels and Bq:1 values are not equal in length, remember to give None as a label for missing vials"

        values = self.get_bq1_values()
        labelled_values = {label: val for label, val in zip(self.vial_labels, values)}

        return labelled_values

    def get_count_times(self) -> List[float]:
        assert self.data is not None
        return self.data["Count Time"].tolist()

    def get_lum(self) -> List[float]:
        assert self.data is not None
        return self.data["LUM"].tolist()


class LSCSample:
    activity: pint.Quantity
    origin_file: str | None

    def __init__(self, activity: pint.Quantity, name: str):
        self.activity = activity
        self.name = name
        # TODO add other attributes available in LSC file
        self.background_substracted = False
        self.origin_file = None

    def __str__(self):
        return f"Sample {self.name}"

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
        act = 0
        for sample in self.samples[:2]:
            act += sample.activity

        return act

    def get_insoluble_activity(self):
        act = 0
        for sample in self.samples[2:]:
            act += sample.activity

        return act

    def get_total_activity(self):
        return self.get_soluble_activity() + self.get_insoluble_activity()


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

    def get_cumulative_activity(
        self, form: Literal["total", "soluble", "insoluble"] = "total"
    ):
        """Calculates the cumulative activity of the gas stream

        Args:
            form: the form in which the cumulative activity is calculated.

        Raises:
            ValueError: If background has not been substracted

        Returns:
            the cumulative activity as a pint.Quantity object
        """
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
        for sample in self.samples:
            if form == "total":
                cumulative_activity.append(sample.get_total_activity())
            elif form == "soluble":
                cumulative_activity.append(sample.get_soluble_activity())
            elif form == "insoluble":
                cumulative_activity.append(sample.get_insoluble_activity())
        cumulative_activity = ureg.Quantity.from_list(cumulative_activity)
        cumulative_activity = cumulative_activity.cumsum()
        return cumulative_activity

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