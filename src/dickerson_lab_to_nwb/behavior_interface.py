from typing import Optional
from pathlib import Path
import h5py
import numpy as np
from pynwb.behavior import BehavioralTimeSeries
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FilePathType

class BehaviorInterface(BaseDataInterface):
    """Interface for behavioral data stored in HDF5 format."""

    def __init__(
        self,
        file_path: FilePathType,
        verbose: bool = True
    ):
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
        
        # Define data structure
        self.data_structure = {
            'AI': {
                'LeftMinusRightWingBeatAmpltude': 'Difference between left and right wing beat amplitudes',
                'LeftWingBeatAmplitude': 'Left wing beat amplitude',
                'VisualStimulus1': 'Visual stimulus channel 1',
                'VisualStimulus2': 'Visual stimulus channel 2'
            },
            'CI': {
                'FrameCounter': 'Frame counter for synchronization'
            },
            'DI': {
                'FrameOut': 'Frame output signal',
                'LineOut': 'Line output signal'
            },
            'Global': {
                'GCtr': 'Global counter'
            }
        }

    def get_metadata(self) -> dict:
        """
        Get metadata for behavioral data.
        
        Returns
        -------
        metadata : dict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()
        
        # Add behavioral metadata
        with h5py.File(self.file_path, 'r') as f:
            metadata.update(
                Behavior=dict(
                    data_structure=self.data_structure,
                    num_samples=f['AI']['LeftWingBeatAmplitude'].shape[0],
                    available_groups=list(f.keys())
                )
            )
        
        return metadata

    def add_to_nwbfile(
        self, 
        nwbfile,
        metadata: Optional[dict] = None
    ):
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
            
        # Create behavioral processing module
        behavior_module = nwbfile.create_processing_module(
            name="behavior",
            description="Processed behavioral data including wing beat measurements and stimuli"
        )
        
        with h5py.File(self.file_path, 'r') as f:
            # Add analog input data
            for channel, description in self.data_structure['AI'].items():
                data = f['AI'][channel][:]
                behavioral_ts = BehavioralTimeSeries(
                    name=f"AI_{channel}",
                    data=data,
                    timestamps=None,  # Timestamps could be derived from sampling rate if known
                    unit="volt",  # Update if different units are used
                    description=description
                )
                behavior_module.add(behavioral_ts)
            
            # Add counter input data
            frame_counter = f['CI']['FrameCounter'][:]
            behavioral_ts = BehavioralTimeSeries(
                name="frame_counter",
                data=frame_counter,
                timestamps=None,
                unit="frames",
                description=self.data_structure['CI']['FrameCounter']
            )
            behavior_module.add(behavioral_ts)
            
            # Add digital input data
            for channel, description in self.data_structure['DI'].items():
                data = f['DI'][channel][:]
                behavioral_ts = BehavioralTimeSeries(
                    name=f"DI_{channel}",
                    data=data,
                    timestamps=None,
                    unit="digital",
                    description=description
                )
                behavior_module.add(behavioral_ts)
            
            # Add global counter
            gctr = f['Global']['GCtr'][:]
            behavioral_ts = BehavioralTimeSeries(
                name="global_counter",
                data=gctr,
                timestamps=None,
                unit="count",
                description=self.data_structure['Global']['GCtr']
            )
            behavior_module.add(behavioral_ts)
