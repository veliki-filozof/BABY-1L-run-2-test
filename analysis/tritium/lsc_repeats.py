import pint
from typing import List, Dict, DefaultDict, Literal
import numpy as np
import warnings
from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
from datetime import datetime, timedelta
import json
from collections import defaultdict
from libra_toolbox.tritium.lsc_measurements import (
    LIBRARun,
    LSCFileReader,
    GasStream,
    LSCSample,
    LIBRASample,
)

# config
DATE_FORMAT = "%m/%d/%Y %I:%M %p"

def stdev_addition(stdevs: List[pint.Quantity]) -> pint.Quantity:
    """Calculates stdev of summed variables 
    """
    result = np.sqrt(np.sum(std**2 for std in stdevs))
    return result

class LSCRepeatSample:
    name: str | None
    counted_activities: List[pint.Quantity] | None
    rep_number: int | None
    avg_activity: pint.Quantity | None
    stdev_activity: pint.Quantity | None
    origin_file: str | None

    def __init__(self, name: str, counted_activities: List[pint.Quantity]):
        self.name = name
        self.counted_activities = counted_activities
        self.rep_number = len(self.counted_activities)
        self._perform_statistics()
        self.background_substracted = False
        self.origin_file = None

    # Calculates mean and stdev of repeated counts
    def _perform_statistics(self):
        self.avg_activity = np.mean(self.counted_activities)
        self.stdev_activity = np.std(self.counted_activities, ddof=1) if self.rep_number > 1 else 0.0 

    def __str__(self):
        return f"Sample {self.name}"
    
    def substract_background(self, background_sample: "LSCRepeatSample"):
        """Substracts the background activity from the sample activity, and updates standard deviation.
        Note: Does NOT substract background from individual counts, only from average

        Args:
            background_sample (LSCRepeatSample): Background sample

        Raises:
            ValueError: If background has already been substracted
        """
        if self.background_substracted:
            raise ValueError("Background already substracted")
        self.avg_activity -= background_sample.avg_activity
        if self.avg_activity.magnitude < 0:
            warnings.warn(
                f"Activity of {self.name} is negative after substracting background. Setting to zero."
            )
            self.activity = 0 * ureg.Bq
        self.stdev_activity = stdev_addition((self.stdev_activity, background_sample.stdev_activity))
        self.background_substracted = True

    @staticmethod
    def from_file(file_reader: LSCFileReader, vial_name: str) -> "LSCRepeatSample":
        """Creates an LSCRepeatSample object from a LSCFileReader object

        Args:
            file_reader (LSCFileReader): LSCFileReader object
            vial_name: Name of the vial

        Raises:
            ValueError: If vial_name is not found in the file reader

        Returns:
            the LSCSample object
        """
        values = file_reader.get_bq1_values()
        labels = file_reader.vial_labels
        if vial_name not in labels:
            raise ValueError(f"Vial {vial_name} not found in the file reader.")
        # Group values by labels
        data: DefaultDict[str, List[float]] = defaultdict(list)
        for label, value in zip(labels, values):
            data[label].append(value)
        activities = data[vial_name] * ureg.Bq
        sample = LSCRepeatSample(vial_name, activities)
        sample.origin_file = file_reader.file_path
        return sample
    
class LIBRARepeatSample:
    samples = List[LSCRepeatSample]

    def __init__(self, samples: List[LSCRepeatSample], time: str | datetime):
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

    def substract_background(self, background_sample: LSCRepeatSample):
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
    
class RepeatGasStream:
    samples: List[LIBRARepeatSample]
    name: str | None

    def __init__(
        self,
        samples: List[LIBRARepeatSample],
        start_time: str | datetime,
        name: str | None = None,
    ):
        """Creates a RepeatGasStream object from a list of LIBRARepeatSample objects

        Args:
            samples: a list of LIBRArepeatSample objects
            start_time: the start time at which the gas stream was collected
            name: a name for the RepeatGasStream. Defaults to None.
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

