from typing import Optional
from pathlib import Path

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FilePathType, DeepDict
from pynwb import NWBFile
from pymatreader import read_mat


class KimLabTrialsInterface(BaseDataInterface):
    """Data interface for Kim Lab trials data."""

    def __init__(
        self,
        trial_data_file_path: FilePathType,
        condition_data_file_path: FilePathType,
        timestamps: np.ndarray,
        stimuli_blocks: int = 10.0,
        verbose: bool = False,
    ):
        """Initialize the trials interface.

        Parameters
        ----------
        trial_data_file_path : FilePathType
            Path to the file containing trial data
        condition_data_file_path : FilePathType
            Path to the file containing condition data
        timestamps : np.ndarray
            Timestamps for the trial data. This is used to synchronize the
            trial data with the rest of the data.
        stimuli_blocks : int, default: 10.0
            Number of stimuli blocks. This is used to determine the block index
            for each trial.
        verbose : bool, default: False
            Whether to print progress information
        """
        super().__init__(
            trial_data_file_path=trial_data_file_path,
            condition_data_file_path=condition_data_file_path,
        )
        self.trial_data_file_path = Path(trial_data_file_path)
        self.condition_data_file_path = Path(condition_data_file_path)
        
        # Validate files exist
        if not self.trial_data_file_path.is_file():
            raise FileNotFoundError(f"Trial data file not found at {self.trial_data_file_path}")
        if not self.condition_data_file_path.is_file():
            raise FileNotFoundError(f"Condition data file not found at {self.condition_data_file_path}")
        self.timestamps = timestamps
        self.stimuli_blocks = stimuli_blocks
        self.verbose = verbose

    def get_metadata(self) -> dict:
        """Get metadata for the trials.

        Returns
        -------
        dict
            The metadata dictionary
        """
        metadata = super().get_metadata()
        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[DeepDict] = None):
        """Add the trials data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add the trials to
        metadata : Optional[DeepDict], optional
            Metadata dictionary
        """
        # Load the trial and condition data
        trial_data_matlab = read_mat(self.trial_data_file_path)["trial"]
        condition_data_matlab = read_mat(self.condition_data_file_path)["condition"]
        time = self.timestamps

        # Convert to uint16 to avoid any potential floating point issues
        trial_data = trial_data_matlab.astype("uint16")
        condition_data = condition_data_matlab.astype("uint32")

        # Create paired views of the array just once
        current_value = trial_data[:-1]
        next_value = trial_data[1:]

        # Generate indices array only once
        indices = np.arange(len(trial_data) - 1)

        # Trials start when increasing from 0 to a value
        trials_start_mask = (current_value == 0) & (next_value >= 1)
        # Trials end when the value drops to 0 from any non-zero value
        trials_stop_mask = (current_value != 0) & (next_value == 0)

        trials_start_indices = indices[trials_start_mask] + 1  # +1 to get the trial start
        trials_stop_indices = indices[trials_stop_mask]  # No +1 to get the last value before 0

        trial_starts = time[trials_start_indices]
        trial_stops = time[trials_stop_indices]
        stimuli_condition = condition_data[trials_start_indices]
        # The stimuli block index is 1 for trial from 1 to self.stimuli_blocks, 2 for trial from self.stimuli_blocks + 1 to 2*self.stimuli_blocks, etc.
        stimuli_block_index = trial_data[trials_start_indices] // self.stimuli_blocks + 1

        # Add trial columns if they don't exist
        nwbfile.add_trial_column(name="stimuli_condition", description="Stimuli condition")
        nwbfile.add_trial_column(name="stimuli_block", description="Stimuli block")

        # Add trials to the NWB file
        number_of_trials = len(trial_starts)
        for trial_index in range(number_of_trials):
            nwbfile.add_trial(
                start_time=trial_starts[trial_index],
                stop_time=trial_stops[trial_index],
                stimuli_condition=stimuli_condition[trial_index],
                stimuli_block=stimuli_block_index[trial_index],
            )

        if self.verbose:
            print(f"Added {number_of_trials} trials to the NWB file")
