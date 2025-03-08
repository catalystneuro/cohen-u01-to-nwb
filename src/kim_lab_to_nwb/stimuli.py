from typing import Optional

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FilePathType, DeepDict
from pynwb import NWBFile
from pynwb.image import ImageSeries
from pymatreader import read_mat


class KimLabStimuliInterface(BaseDataInterface):
    """Data interface for Kim Lab visual stimuli data."""

    def __init__(
        self,
        file_path: FilePathType,
        timestamps: np.ndarray,
        verbose: bool = False,
    ):
        """Initialize the stimuli interface.

        Parameters
        ----------
        file_path : FilePathType
            Path to the visual_stimuli.mat file
        timestamps : np.ndarray
            Timestamps for the visual stimuli. This is usedd to synchronize the
            visual stimuli with the rest of the data.
        verbose : bool, default: False
            Whether to print progress information
        """
        super().__init__(file_path=file_path)
        self.file_path = file_path
        self.timestamps = timestamps
        self.verbose = verbose

    def get_metadata(self) -> dict:
        """Get metadata for the visual stimuli.

        Returns
        -------
        dict
            The metadata dictionary
        """
        metadata = super().get_metadata()
        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[DeepDict] = None):
        """Add the visual stimuli data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add the stimuli to
        metadata : Optional[DeepDict], optional
            Metadata dictionary
        """
        # Load the visual stimuli data
        data = read_mat(self.file_path)["visual_stimuli"]

        # Create an image series for the visual stimuli
        # The data shape is (16, 80, 495000) where:
        # - First two dimensions (16, 80) are the stimulus image at each timepoint
        # - Last dimension (495000) is time
        data = np.moveaxis(data, 2, 0)
        image_series = ImageSeries(
            name="VisualStimuli",
            data=data,
            unit="n.a.",
            timestamps=self.timestamps,
            description="Visual stimuli presented during the experiment",
        )

        # Add to the NWB file's stimulus presentation
        nwbfile.add_stimulus(image_series)
