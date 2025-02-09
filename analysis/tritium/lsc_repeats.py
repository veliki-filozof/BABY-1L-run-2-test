import pint
from typing import List, Dict, DefaultDict
import numpy as np
import warnings
from libra_toolbox.tritium.model import ureg, Model, quantity_to_activity
import json
from collections import defaultdict
from libra_toolbox.tritium.lsc_measurements import (
    LIBRARun,
    LSCFileReader,
    GasStream,
    LSCSample,
    LIBRASample,
)

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
