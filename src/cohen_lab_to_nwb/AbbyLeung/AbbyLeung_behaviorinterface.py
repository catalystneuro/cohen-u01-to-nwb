"""Primary class for converting experiment-specific behavior."""
import os

import numpy as np
from pynwb import TimeSeries
from pynwb.file import NWBFile

from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from scipy.io import loadmat

behav_data = dict(
    wingL="Matrix containing all 3 left wing angles. Row1 = stroke, row2 = deviation, row3 = pitch.",
    wingR="Matrix containing all 3 right wing angles. Row1 = stroke, row2 = deviation, row3 = pitch",
    bodyPitch="body pitch angles",
    bodyRoll="body roll angles",
    bodyYaw="body yaw angles",
)

class AbbyleungBehaviorInterface(BaseDataInterface):
    """Behavior interface for AbbyLeung conversion"""

    keywords = ["behavior"]

    def __init__(self, file_path: str):
        super().__init__(file_path=file_path)

    def get_metadata(self) -> DeepDict:
        # Automatically retrieve as much metadata as possible from the source files available
        metadata = super().get_metadata()   

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict):
        # All the custom code to add the data the nwbfile

        var_name = os.path.splitext(os.path.split(self.source_data.file_path)[1])[0]

        mat_out = loadmat(self.source_data.file_path, simplify_cells=True)
        for i, x in enumerate(mat_out[var_name]):
            for ts_name, desc in behav_data.items():
                nwbfile.add_acquisition(
                    TimeSeries(
                        name=f"{ts_name}{i:03d}",
                        data=x[ts_name],
                        unit="degrees",
                        rate=0.000125,
                        starting_time=np.nan,
                    )
                )
