from pathlib import Path
from typing import Optional, List, Dict, Any

import numpy as np
from pymatreader import read_mat
from pynwb.icephys import (
    IntracellularElectrode,
    PatchClampSeries,
    CurrentClampSeries,
    VoltageClampSeries,
)
from neuroconv.basetemporalalignmentinterface import BaseTemporalAlignmentInterface
from neuroconv.utils import FilePathType


class PatchClampInterface(BaseTemporalAlignmentInterface):
    """
    Data interface for Suver Lab patch clamp data stored in MATLAB files.
    """

    def __init__(
        self,
        file_path: FilePathType,
        verbose: bool = True,
    ):
        """
        Initialize the PatchClampInterface.

        Parameters
        ----------
        file_path : FilePathType
            Path to the MATLAB file containing patch clamp data.
        verbose : bool, default: True
            If True, print verbose output.
        """
        super().__init__(file_path=file_path, verbose=verbose)
        self.matlab_data = None
        self.non_array_keys = []
        self.array_keys = []
        self.metadata_df = None
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
        metadata["NWBFile"].update(
            session_description=f"Patch clamp recording from {first_recording.get('date', 'unknown date')}",
            identifier=f"{first_recording.get('date', 'unknown')}_{first_recording.get('expNumber', 'unknown')}",
            session_id=f"{first_recording.get('date', 'unknown')}_{first_recording.get('expNumber', 'unknown')}",
            experimenter=["Suver Lab"],
        )
        
        # Add subject metadata
        metadata["Subject"] = {
            "subject_id": f"{first_recording.get('genotype', 'unknown')}",
            "genotype": first_recording.get("genotype", "unknown"),
            "age": first_recording.get("age", "unknown"),
        }
        
        # Add device metadata
        metadata["Ecephys"] = {
            "Device": [
                {
                    "name": "PatchClampDevice",
                    "description": "Patch clamp device used for recording",
                    "manufacturer": "Unknown",
                }
            ],
            "ElectrodeGroup": [
                {
                    "name": "PatchClampElectrode",
                    "description": "Patch clamp electrode",
                    "device": "PatchClampDevice",
                    "location": "Unknown",
                }
            ],
            "Electrodes": [
                {
                    "name": "PatchClampElectrode",
                    "description": "Patch clamp electrode",
                }
            ],
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
            sample_rate = float(recording.get("samplerate", 10000))
            n_samples = len(recording.get("Vm", []))
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
        Add the patch clamp data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to.
        metadata : dict, optional
            Metadata dictionary.
        """
        if self.matlab_data is None:
            self._read_data()
        
        # Create device
        device = nwbfile.create_device(
            name="PatchClampDevice",
            description="Patch clamp device used for recording",
            manufacturer="Unknown",
        )
        
        # Create electrode
        electrode = nwbfile.create_icephys_electrode(
            name="PatchClampElectrode",
            description="Patch clamp electrode",
            device=device,
            location="Unknown",
        )
        
        # Get timestamps
        timestamps_list = self.get_timestamps()
        
        # Add each recording to the NWB file
        for i, recording in enumerate(self.matlab_data):
            trial_id = recording.get("trial", i + 1)
            timestamps = timestamps_list[i]
            
            # Add Vm (membrane potential) as CurrentClampSeries
            if "Vm" in recording and len(recording["Vm"]) > 0:
                vm_series = CurrentClampSeries(
                    name=f"Vm_trial{trial_id:03d}",
                    data=recording["Vm"],
                    electrode=electrode,
                    gain=float(recording.get("scaleVoltage", 1.0)),
                    starting_time=timestamps[0],
                    rate=float(recording.get("samplerate", 10000)),
                    description="Membrane potential recording",
                    bias_current=0.0,  # Unknown, using default
                    bridge_balance=0.0,  # Unknown, using default
                    capacitance_compensation=0.0,  # Unknown, using default
                )
                nwbfile.add_acquisition(vm_series)
            
            # Add filteredVm if available
            if "filteredVm" in recording and len(recording["filteredVm"]) > 0:
                filtered_vm_series = CurrentClampSeries(
                    name=f"filteredVm_trial{trial_id:03d}",
                    data=recording["filteredVm"],
                    electrode=electrode,
                    gain=float(recording.get("scaleVoltage", 1.0)),
                    starting_time=timestamps[0],
                    rate=float(recording.get("samplerate", 10000)),
                    description="Filtered membrane potential recording",
                    bias_current=0.0,  # Unknown, using default
                    bridge_balance=0.0,  # Unknown, using default
                    capacitance_compensation=0.0,  # Unknown, using default
                )
                nwbfile.add_acquisition(filtered_vm_series)
            
            # Add I (current) as VoltageClampSeries
            if "I" in recording and len(recording["I"]) > 0:
                i_series = VoltageClampSeries(
                    name=f"I_trial{trial_id:03d}",
                    data=recording["I"],
                    electrode=electrode,
                    gain=float(recording.get("scaleCurrent", 1.0)),
                    starting_time=timestamps[0],
                    rate=float(recording.get("samplerate", 10000)),
                    description="Current recording",
                    capacitance_fast=0.0,  # Unknown, using default
                    capacitance_slow=0.0,  # Unknown, using default
                    resistance_comp_bandwidth=0.0,  # Unknown, using default
                    resistance_comp_correction=0.0,  # Unknown, using default
                    resistance_comp_prediction=0.0,  # Unknown, using default
                    whole_cell_capacitance_comp=0.0,  # Unknown, using default
                    whole_cell_series_resistance_comp=0.0,  # Unknown, using default
                )
                nwbfile.add_acquisition(i_series)
            
            # Add puffer data if available
            if "puffer" in recording and len(recording["puffer"]) > 0:
                puffer_series = PatchClampSeries(
                    name=f"puffer_trial{trial_id:03d}",
                    data=recording["puffer"],
                    electrode=electrode,
                    gain=1.0,  # Unknown, using default
                    starting_time=timestamps[0],
                    rate=float(recording.get("samplerate", 10000)),
                    description="Puffer recording (stimulus delivery)",
                )
                nwbfile.add_stimulus(puffer_series)
            
            # Add tachometer data if available
            if "tachometer" in recording and len(recording["tachometer"]) > 0:
                tachometer_series = PatchClampSeries(
                    name=f"tachometer_trial{trial_id:03d}",
                    data=recording["tachometer"],
                    electrode=electrode,
                    gain=1.0,  # Unknown, using default
                    starting_time=timestamps[0],
                    rate=float(recording.get("samplerate", 10000)),
                    description="Tachometer recording (wing beat frequency)",
                )
                nwbfile.add_acquisition(tachometer_series)
            
            # Create a trial for each recording
            nwbfile.add_trial(
                start_time=timestamps[0],
                stop_time=timestamps[-1],
                trial_id=trial_id,
                condition=recording.get("condition", "unknown"),
            )
            
            # Add intracellular recording entry
            if hasattr(nwbfile, "add_intracellular_recording"):
                # Add to intracellular recordings table if Vm and I are available
                if "Vm" in recording and "I" in recording:
                    nwbfile.add_intracellular_recording(
                        electrode=electrode,
                        response=vm_series,
                        stimulus=i_series,
                        id=trial_id,
                    )
