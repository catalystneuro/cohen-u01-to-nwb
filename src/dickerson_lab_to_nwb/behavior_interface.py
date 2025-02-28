from typing import Optional
from pathlib import Path
import h5py
import numpy as np
from pynwb.base import TimeSeries
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FilePathType


class BehaviorInterface(BaseDataInterface):
    """Interface for behavioral data stored in HDF5 format."""

    def __init__(self, file_path: FilePathType, verbose: bool = True):
        """
        Initialize reading of behavioral data.

        Parameters
        ----------
        file_path : FilePathType
            Path to the HDF5 file containing behavioral data
        verbose : bool, default: True
        """
        super().__init__(file_path=file_path, verbose=verbose)
        self.file_path = Path(file_path)

    def get_metadata(self) -> dict:
        """
        Get metadata for behavioral data.

        Returns
        -------
        metadata : dict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        return metadata

    def add_to_nwbfile(self, nwbfile, metadata: Optional[dict] = None):
        """
        Add behavioral data to NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the behavioral data to
        metadata : dict, optional
            Metadata dictionary
        """
        if metadata is None:
            metadata = self.get_metadata()


        with h5py.File(self.file_path, "r") as f:
            
            # Timestamps
            frame_counter = f["CI"]["FrameCounter"]    
            # From personal communication we know that frames are recorded at 4 Hz
            sampling_rate = 4.0
            timestamps = frame_counter[:].squeeze() / sampling_rate  # (frames / (frames / seconds))
            
            
            # WingBeat
            left_wing_beat_amplitude = f["AI"]["LeftWingBeatAmplitude"][:].squeeze()
            left_minus_right_wing_beat = f["AI"]["LeftMinusRightWingBeatAmplitude"][:].squeeze()

            left_wing_beat_amplitude_ts = TimeSeries(
                name="LeftWingBeatAmplitudeTimeSeries",
                data=left_wing_beat_amplitude,
                unit="volt",
                timestamps=timestamps,
                description="Left wing beat amplitude",
            )

            left_minus_right_wing_beat_ts = TimeSeries(
                name="LeftMinusRightWingBeatAmplitudeTimeSeries",
                data=left_minus_right_wing_beat,
                unit="volt",
                description="Left minus right wing beat amplitude",
                timestamps=timestamps,
            )

            nwbfile.add_acquisition(left_wing_beat_amplitude_ts)
            nwbfile.add_acquisition(left_minus_right_wing_beat_ts)

            visual_stimulus_x = f["AI"]["VisualStimulus1"]
            visual_stimulus_y = f["AI"]["VisualStimulus2"]

            from pynwb.behavior import SpatialSeries
            
            description = (
                "Analog output from the controller of the visual arena feeding "
                "information about the pattern/bar movement in the X and Y direction "
                "to Thorsync using the same breakout box."
            )
            
            visual_stimuli_stack = np.stack((visual_stimulus_x[:], visual_stimulus_y[:]), axis=1).squeeze()
            
            spatial_series = SpatialSeries(
                name="VisualStimulusSpatialSeries",
                data=visual_stimuli_stack,
                unit="volt", 
                reference_frame="unknown", # TODO: figure out what this should be
                description=description,
                timestamps=timestamps,
            )

            
            # add as stimulus
            nwbfile.add_stimulus(spatial_series)