from pathlib import Path
from typing import Optional, List, Dict, Any

import numpy as np
from pymatreader import read_mat
from pynwb.icephys import (
    IntracellularElectrode,
    PatchClampSeries,
    CurrentClampSeries,
)
from neuroconv.basetemporalalignmentinterface import BaseTemporalAlignmentInterface
from neuroconv.utils import FilePathType


class SealTestInterface(BaseTemporalAlignmentInterface):
    """
    Data interface for Suver Lab seal test data stored in MATLAB files.
    """

    def __init__(
        self,
        file_path: FilePathType,
        verbose: bool = True,
    ):
        """
        Initialize the SealTestInterface.

        Parameters
        ----------
        file_path : FilePathType
            Path to the MATLAB file containing seal test data.
        verbose : bool, default: True
            If True, print verbose output.
        """
        super().__init__(file_path=file_path, verbose=verbose)
        self.matlab_data = None
        self.non_array_keys = []
        self.array_keys = []
        self._recording_count = None

    def get_metadata(self) -> dict:
        """
        Extract metadata from the MATLAB file.

        Returns
        -------
        metadata : dict
            Metadata dictionary.
        """
        metadata = super().get_metadata()
        
        # Extract metadata from the MATLAB file
        if self.matlab_data is None:
            self._read_data()
        
        # Extract basic metadata from the first recording
        first_recording = self.matlab_data[0]
        
        # Add metadata to the NWBFile
        file_path = Path(self.source_data["file_path"])
        date_str = file_path.stem.split("_SealTest_")[0]
        
        metadata["NWBFile"].update(
            session_description=f"Seal test recording from {date_str}",
            identifier=f"{date_str}_SealTest",
            session_id=f"{date_str}_SealTest",
            experimenter=["Suver Lab"],
        )
        
        # Add subject metadata
        metadata["Subject"] = {
            "subject_id": f"{first_recording.get('genotype', 'unknown')}",
            "genotype": first_recording.get("genotype", "unknown"),
        }
        
        return metadata

    def _read_data(self):
        """
        Read data from the MATLAB file.
        """
        # Read the MATLAB file
        matlab_data = read_mat(self.source_data["file_path"])["data"]
        
        # Convert to list of dictionaries
        self.matlab_data = [
            {field: record[field] for field in record.dtype.names} 
            for record in matlab_data
        ]
        
        # Identify array and non-array keys
        self.non_array_keys = []
        self.array_keys = []
        
        for key, value in self.matlab_data[0].items():
            if isinstance(value, np.ndarray):
                self.array_keys.append(key)
            else:
                self.non_array_keys.append(key)
        
        # Store the number of recordings
        self._recording_count = len(self.matlab_data)

    def get_original_timestamps(self) -> List[np.ndarray]:
        """
        Get the original timestamps for each recording.

        Returns
        -------
        timestamps : List[np.ndarray]
            List of timestamps arrays, one for each recording.
        """
        if self.matlab_data is None:
            self._read_data()
        
        timestamps = []
        for recording in self.matlab_data:
            # Create timestamps based on sample rate
            sample_rate = float(recording.get("SAMPLERATE", 10000))
            n_samples = len(recording.get("stimTiming", []))
            timestamps.append(np.arange(n_samples) / sample_rate)
        
        return timestamps

    def get_timestamps(self) -> List[np.ndarray]:
        """
        Get the synchronized timestamps for each recording.

        Returns
        -------
        timestamps : List[np.ndarray]
            List of synchronized timestamps arrays, one for each recording.
        """
        if self.matlab_data is None:
            self._read_data()
        
        if self._timestamps is None:
            return self.get_original_timestamps()
        
        return self._timestamps

    def add_to_nwbfile(self, nwbfile, metadata: Optional[dict] = None):
        """
        Add the seal test data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to.
        metadata : dict, optional
            Metadata dictionary.
        """
        if self.matlab_data is None:
            self._read_data()
        
        # Create device if it doesn't exist
        device_name = "PatchClampDevice"
        if device_name not in nwbfile.devices:
            device = nwbfile.create_device(
                name=device_name,
                description="Patch clamp device used for recording",
                manufacturer="Unknown",
            )
        else:
            device = nwbfile.devices[device_name]
        
        # Create electrode if it doesn't exist
        electrode_name = "SealTestElectrode"
        if electrode_name not in nwbfile.intracellular_electrodes:
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description="Seal test electrode",
                device=device,
                location="Unknown",
            )
        else:
            electrode = nwbfile.intracellular_electrodes[electrode_name]
        
        # Get timestamps
        timestamps_list = self.get_timestamps()
        
        # Add each recording to the NWB file
        for i, recording in enumerate(self.matlab_data):
            trial_id = recording.get("trial", i + 1)
            timestamps = timestamps_list[i]
            
            # Add stimTiming as PatchClampSeries
            if "stimTiming" in recording and len(recording["stimTiming"]) > 0:
                stim_series = PatchClampSeries(
                    name=f"stimTiming_sealTest{trial_id:03d}",
                    data=recording["stimTiming"],
                    electrode=electrode,
                    gain=1.0,  # Unknown, using default
                    starting_time=timestamps[0],
                    rate=float(recording.get("SAMPLERATE", 10000)),
                    description="Optogenetic stimulus timing for seal test",
                )
                nwbfile.add_stimulus(stim_series)
            
            # Add micLeft as PatchClampSeries
            if "micLeft" in recording and len(recording["micLeft"]) > 0:
                mic_series = PatchClampSeries(
                    name=f"micLeft_sealTest{trial_id:03d}",
                    data=recording["micLeft"],
                    electrode=electrode,
                    gain=1.0,  # Unknown, using default
                    starting_time=timestamps[0],
                    rate=float(recording.get("SAMPLERATE", 10000)),
                    description="Microphone recording (left) for seal test",
                )
                nwbfile.add_acquisition(mic_series)
            
            # Create a trial for each recording
            pre_time = float(recording.get("PRE_TRIAL_TIME", 0))
            stim_time = float(recording.get("TRIAL_TIME_LIGHT", 0))
            post_time = float(recording.get("POST_TRIAL_TIME", 0))
            
            # Calculate start and stop times
            start_time = timestamps[0]
            stim_start = start_time + pre_time
            stim_stop = stim_start + stim_time
            stop_time = stim_stop + post_time
            
            # Add trial
            nwbfile.add_trial(
                start_time=start_time,
                stop_time=stop_time,
                trial_id=trial_id,
                stim_start=stim_start,
                stim_stop=stim_stop,
                notes=recording.get("notes", ""),
            )
