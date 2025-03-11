from datetime import datetime
from pathlib import Path

import numpy as np
from neuroconv.tools.nwb_helpers import get_default_backend_configuration, configure_and_write_nwbfile
from neuroconv.utils import load_dict_from_file
from neuroconv.datainterfaces import VideoInterface, ScanImageImagingInterface
from neuroconv import ConverterPipe
from pynwb import NWBFile, TimeSeries
from pynwb.file import Subject
from pymatreader import read_mat
from pynwb.device import Device

from kim_lab_to_nwb.ophys import KimLabROIInterface
from kim_lab_to_nwb.stimuli import KimLabStimuliInterface
from kim_lab_to_nwb.trials import KimLabTrialsInterface
from kim_lab_to_nwb.behavior import BehaviorInterface
from cohen_u01_nwb_conversion_utils.utils import detect_threshold_crossings
from zoneinfo import ZoneInfo


def convert_session_to_nwb(
    matlab_data_file_path: str | Path,
    experiment_info_file_path: str | Path,
    output_dir: str | Path,
    video_file_path: str | Path = None,
    tiff_folder_path: str | Path = None,
    df_f_file_path: str | Path = None,
    roi_info_file_path: str | Path = None,
    visual_stimuli_file_path: str | Path = None,
    trial_data_file_path: str | Path = None,
    condition_data_file_path: str | Path = None,
    verbose: bool = False,
) -> Path:
    """
    Convert a single experimental session from Kim Lab data to NWB format.

    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB data file containing the following data:
        - data[0,:] : time (timestamps for all measurements)
        - data[1,:] : left wingbeat amplitude
        - data[2,:] : left-right wingbeat (difference between left and right wingbeat)
        - data[3,:] : x-position of the visual pattern
        - data[4,:] : y-position of the visual pattern
        - data[5,:] : 2-photon frame synchronization signal (1 pulse corresponds to 1 frame)
        - data[6,:] : behavior camera signal (1 pulse corresponds to 1 frame)
        - data[7,:] : stimulus start indicator (may be empty in some examples)
    experiment_info_file_path : str or Path
        Path to the experiment info file containing metadata such as:
        - age: age of the fly
        - cross: genotype/strain information
    output_dir : str or Path
        Directory where the NWB file will be saved
    video_file_path : str or Path, optional
        Path to the video file recording fly behavior
    tiff_folder_path : str or Path, optional
        Path to the folder containing TIFF files from ScanImage acquisition.
        These are two-photon microscopy images from a self-made microscope.
    df_f_file_path : str or Path, optional
        Path to the df_f.mat file containing fluorescence traces (ΔF/F) for each ROI.
        Shape is typically (num_rois, num_timepoints).
    roi_info_file_path : str or Path, optional
        Path to the ROI_info.mat file containing ROI coordinates and reference images.
        Contains fields:
        - x_cor: x-coordinates of ROI vertices
        - y_cor: y-coordinates of ROI vertices
        - reference_image: reference image for ROIs
        - total_ROI: total number of ROIs
    visual_stimuli_file_path : str or Path, optional
        Path to the visual_stimuli.mat file containing visual stimuli presented to the fly.
        Shape is typically (16, 80, num_timepoints) where 16x80 is the size of the displayed image.
    trial_data_file_path : str or Path, optional
        Path to the trial.mat file containing trial information.
        Used to detect trial starts and stops based on value transitions.
    condition_data_file_path : str or Path, optional
        Path to the condition.mat file containing condition information for each trial.
        Conditions represent which pattern is being shown to the fly.
    verbose : bool, optional
        Whether to print progress information, by default False

    Returns
    -------
    Path
        Path to the created NWB file

    Raises
    ------
    FileNotFoundError
        If any of the required input files are not found

    Notes
    -----
    The conversion process creates an NWB file with the following components:
    - Subject metadata (species, strain, age, sex)
    - Acquisition data (wingbeat, left-right wingbeat, stimulus position)
    - Imaging data (two-photon microscopy images)
    - ROI data (regions of interest in the imaging data)
    - Visual stimuli data (patterns presented to the fly)
    - Trial data (trial start/stop times and conditions)
    - Video data (behavior camera recordings)

    Optional components are only included if the corresponding files are provided.
    """
    # Convert all paths to Path objects
    matlab_data_file_path = Path(matlab_data_file_path)
    experiment_info_file_path = Path(experiment_info_file_path)
    output_dir = Path(output_dir)
    
    # Convert optional paths to Path objects if provided
    video_file_path = Path(video_file_path) if video_file_path is not None else None
    tiff_folder_path = Path(tiff_folder_path) if tiff_folder_path is not None else None
    df_f_file_path = Path(df_f_file_path) if df_f_file_path is not None else None
    roi_info_file_path = Path(roi_info_file_path) if roi_info_file_path is not None else None
    visual_stimuli_file_path = Path(visual_stimuli_file_path) if visual_stimuli_file_path is not None else None
    trial_data_file_path = Path(trial_data_file_path) if trial_data_file_path is not None else None
    condition_data_file_path = Path(condition_data_file_path) if condition_data_file_path is not None else None

    # File validation is handled by the respective interfaces
    # Validate experiment info file exists (not handled by any interface)
    if not experiment_info_file_path.is_file():
        raise FileNotFoundError(f"Experiment info file not found at {experiment_info_file_path}")

    # Load experiment info
    experiment_info = read_mat(experiment_info_file_path)
    age = experiment_info["age"]
    genotype = experiment_info["cross"]

    # Extract session ID from matlab file path
    session_id = matlab_data_file_path.stem.split("data_")[1]
    if verbose:
        print(f"Processing {session_id=} ({age=}, {genotype=})")

    session_start_time = datetime.strptime(session_id[:8], "%Y%m%d")
    if verbose:
        print(f"Session start time: {session_start_time}")

    # Initialize behavior interface
    behavior_interface = BehaviorInterface(
        file_path=matlab_data_file_path,
        verbose=verbose,
    )
    
    # Initialize data interfaces dictionary
    data_interfaces = {
        "behavior": behavior_interface
    }
    
    # Set up imaging interfaces if tiff folder is provided
    if tiff_folder_path is not None:
        tiff_file_path = tiff_folder_path / f"{session_id}_00001.tif"
        # Set up the imaging interfaces for two channels
        scan_image_interface_channel_1 = ScanImageImagingInterface(
            file_path=tiff_file_path,
            channel_name="Channel 1",
        )

        # scan_image_interface_channel_2 = ScanImageImagingInterface(
        #     file_path=tiff_file_path,
        #     channel_name="Channel 2",
        # )

        num_frames = scan_image_interface_channel_1.imaging_extractor.get_num_frames()
        two_photon_timestamps = behavior_interface.timestamps[detect_threshold_crossings(behavior_interface.two_photon_frame_sync, 0.5)]
        scan_image_interface_channel_1.set_aligned_timestamps(aligned_timestamps=two_photon_timestamps[:num_frames])
#        scan_image_interface_channel_2.set_aligned_timestamps(aligned_timestamps=two_photon_timestamps[:num_frames])
        
        data_interfaces["imaging_channel_1"] = scan_image_interface_channel_1
 #       data_interfaces["imaging_channel_2"] = scan_image_interface_channel_2

    # Set up ROI interface if df_f_file_path and roi_info_file_path are provided
    if df_f_file_path is not None and roi_info_file_path is not None:
        roi_interface = KimLabROIInterface(
            file_path=df_f_file_path,
            roi_info_file_path=roi_info_file_path,
            image_shape=(64, 64),
            timestamps=behavior_interface.timestamps,
            verbose=verbose,
        )
        data_interfaces["roi"] = roi_interface

    # Set up stimuli interface if visual_stimuli_file_path is provided
    if visual_stimuli_file_path is not None:
        stimuli_interface = KimLabStimuliInterface(
            file_path=visual_stimuli_file_path,
            timestamps=behavior_interface.timestamps,
            verbose=verbose,
        )
        data_interfaces["stimuli"] = stimuli_interface
    
    # Set up trials interface if trial_data_file_path and condition_data_file_path are provided
    if trial_data_file_path is not None and condition_data_file_path is not None:
        trials_interface = KimLabTrialsInterface(
            trial_data_file_path=trial_data_file_path,
            condition_data_file_path=condition_data_file_path,
            stimuli_blocks=10.0,
            timestamps=behavior_interface.timestamps,
            verbose=verbose,
        )
        data_interfaces["trials"] = trials_interface

    # Set up video interface if video_file_path is provided
    if video_file_path is not None:
        video_interface = VideoInterface(file_paths=[video_file_path])
        video_timestamps = behavior_interface.timestamps[detect_threshold_crossings(behavior_interface.behavior_camera_sync, 0.5)]
        video_interface.set_aligned_timestamps([video_timestamps])
        data_interfaces["video"] = video_interface

    converter_pipe = ConverterPipe(data_interfaces=data_interfaces)

    california_tz = ZoneInfo("America/Los_Angeles")
    session_start_time = session_start_time.replace(tzinfo=california_tz)

    # Load metadata from YAML file
    metadata_path = Path(__file__).parent / "metadata.yaml"
    metadata = load_dict_from_file(metadata_path)
    
    # Update metadata with session-specific information
    metadata.update(converter_pipe.get_metadata())
    metadata["NWBFile"]["session_start_time"] = session_start_time
    metadata["NWBFile"]["session_description"] = f"{metadata['NWBFile']['session_description']} protocol: {behavior_interface.protocol}"
    metadata["NWBFile"]["session_id"] = session_id
    metadata["Subject"]["strain"] = genotype
    metadata["Subject"]["age"] = f"P{age}D"
    metadata["Subject"]["subject_id"] = session_id

    # Set conversion options for interfaces that are present
    conversion_options = {}
    if "imaging_channel_1" in data_interfaces:
        conversion_options["imaging_channel_1"] = dict(stub_test=False, photon_series_index=0)
    if "imaging_channel_2" in data_interfaces:
        conversion_options["imaging_channel_2"] = dict(stub_test=False, photon_series_index=1)

    # Run the conversion
    nwbfile_path = output_dir / f"{session_id}.nwb"
    converter_pipe.run_conversion(
        nwbfile_path=nwbfile_path,
        metadata=metadata,
        conversion_options=conversion_options,
        overwrite=True,
    )

    if verbose:
        print(f"Created NWB file: {nwbfile_path}")

    return nwbfile_path


if __name__ == "__main__":
    # Example usage with the actual paths
    data_folder_path = Path("/home/heberto/cohen_project/Sample data/Kim Lab/20250301b")
    
    # Define required input paths
    matlab_data_file_path = data_folder_path / "data_20250301b_00007.mat"
    experiment_info_file_path = data_folder_path / "exp_info.mat"
    # Modify this to your desired output directory
    output_dir = Path("/home/heberto/cohen_project/Sample data/Kim Lab/nwb")
    
    # Define optional input paths (set to None if not available)
    video_file_path = data_folder_path / "20250301b_00007.avi"  # Behavior video recording
    tiff_folder_path = data_folder_path  # Contains tiff files like 20250301b_00007_00001.tif
    visual_stimuli_file_path = data_folder_path / "visual_stimuli.mat"  # Visual stimuli presented to the fly
    trial_data_file_path = data_folder_path / "trial.mat"  # Trial information
    condition_data_file_path = data_folder_path / "condition.mat"  # Condition information for each trial
    
    # These files are not available for this dataset
    df_f_file_path = None  # No fluorescence traces available
    roi_info_file_path = None  # No ROI information available

    output_dir.mkdir(exist_ok=True, parents=True)

    nwbfile_path = convert_session_to_nwb(
        # Required parameters
        matlab_data_file_path=matlab_data_file_path,
        experiment_info_file_path=experiment_info_file_path,
        output_dir=output_dir,
        
        # Optional parameters (available for this dataset)
        video_file_path=video_file_path,
        tiff_folder_path=tiff_folder_path,
        visual_stimuli_file_path=visual_stimuli_file_path,
        trial_data_file_path=trial_data_file_path,
        condition_data_file_path=condition_data_file_path,
        
        # Optional parameters (not available for this dataset)
        df_f_file_path=df_f_file_path,  # No fluorescence traces available
        roi_info_file_path=roi_info_file_path,  # No ROI information available
        
        verbose=True,  # Enable verbose output for demonstration
    )
