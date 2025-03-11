from typing import Optional
from pathlib import Path

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FilePathType, DeepDict
from pynwb import NWBFile
from pynwb.file import TimeSeries
from pynwb.device import Device
from pymatreader import read_mat


class BehaviorInterface(BaseDataInterface):
    """Data interface for Kim Lab behavior data from NIDQ (National Instruments Data Acquisition)."""

    def __init__(
        self,
        file_path: FilePathType,
        verbose: bool = False,
    ):
        """Initialize the behavior interface.

        Parameters
        ----------
        file_path : FilePathType
            Path to the MATLAB data file containing behavior data
        verbose : bool, default: False
            Whether to print progress information
        """
        super().__init__(file_path=file_path)
        self.file_path = Path(file_path)
        self.verbose = verbose
        
        # Validate file exists
        if not self.file_path.is_file():
            raise FileNotFoundError(f"Matlab data file not found at {self.file_path}")
            
        # Load data
        mat_data = read_mat(self.file_path)
        nidq_device_data = mat_data["data"]
        
        # Extract timestamps and other data
        self.timestamps = nidq_device_data[0]
        self.left_wingbeat = nidq_device_data[1]
        self.left_right_wingbeat = nidq_device_data[2]
        self.x_position = nidq_device_data[3]
        self.y_position = nidq_device_data[4]
        self.two_photon_frame_sync = nidq_device_data[5]
        self.behavior_camera_sync = nidq_device_data[6]
        self.stimulus_start = nidq_device_data[7]
        
        # Store protocol information
        self.protocol = mat_data.get("protocol", "unknown")

    def get_metadata(self) -> dict:
        """Get metadata for the behavior data.

        Returns
        -------
        dict
            The metadata dictionary
        """
        metadata = super().get_metadata()
        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[DeepDict] = None):
        """Add the behavior data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add the behavior data to
        metadata : Optional[DeepDict], optional
            Metadata dictionary
        """
        # Add DAQ device
        ni_daq = Device(
            name="NI 782258-01",
            description=(
                "Multifunction DAQ device with USB 2.0 connectivity. "
                "32 single-ended or 16 differential analog inputs (16-bit resolution), "
                "4 analog outputs (±10V, 16-bit resolution), "
                "32 digital I/O, 16 bidirectional channels, "
                "4 counter/timers, 1 MS/s sampling rate."
            ),
            manufacturer="National Instruments (NI)",
        )
        nwbfile.add_device(ni_daq)

        # Add timeseries data
        timeseries_wingbeat = TimeSeries(
            name="LeftWingbeatSeries",
            data=self.left_wingbeat,
            unit="n.a.",
            timestamps=self.timestamps,
            description="wingbeat",
        )
        nwbfile.add_acquisition(timeseries_wingbeat)

        # Note that this pattern indicates that the timestamps are the same as wingbeat
        common_timestamps = timeseries_wingbeat
        timeseries_left_right_wingbeat = TimeSeries(
            name="LeftRightWingbeatSeries",
            data=self.left_right_wingbeat,
            unit="n.a.",
            timestamps=common_timestamps,
            description="left-right wingbeat",
        )
        nwbfile.add_acquisition(timeseries_left_right_wingbeat)

        timeseries_x_position = TimeSeries(
            name="StimulusPositionSeries",
            data=np.c_[self.x_position, self.y_position],
            unit="n.a.",
            timestamps=common_timestamps,
            description="position of the visual pattern",
        )
        nwbfile.add_acquisition(timeseries_x_position)
        
        # Update session description with protocol information
        if "protocol" not in nwbfile.session_description:
            current_description = nwbfile.session_description or ""
            if current_description:
                nwbfile.session_description = f"{current_description}; protocol: {self.protocol}"
            else:
                nwbfile.session_description = f"protocol: {self.protocol}"
                
        if self.verbose:
            print(f"Added behavior data to the NWB file")
