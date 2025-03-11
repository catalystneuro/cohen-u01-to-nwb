from datetime import datetime
from pathlib import Path
import uuid

import numpy as np
from neuroconv.tools.nwb_helpers import get_default_backend_configuration, configure_and_write_nwbfile
from neuroconv.datainterfaces import VideoInterface, ScanImageImagingInterface
from neuroconv import ConverterPipe
from pynwb import NWBFile, TimeSeries
from pynwb.file import Subject
from pymatreader import read_mat
from pynwb.device import Device

from kim_lab_to_nwb.ophys import KimLabROIInterface
from kim_lab_to_nwb.stimuli import KimLabStimuliInterface
from kim_lab_to_nwb.trials import KimLabTrialsInterface
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
        Path to the MATLAB data file
    experiment_info_file_path : str or Path
        Path to the experiment info file
    output_dir : str or Path
        Directory where the NWB file will be saved
    video_file_path : str or Path, optional
        Path to the video file
    tiff_folder_path : str or Path, optional
        Path to the folder containing TIFF files
    df_f_file_path : str or Path, optional
        Path to the df_f.mat file
    roi_info_file_path : str or Path, optional
        Path to the ROI_info.mat file
    visual_stimuli_file_path : str or Path, optional
        Path to the visual_stimuli.mat file
    trial_data_file_path : str or Path, optional
        Path to the trial_data.mat file
    condition_data_file_path : str or Path, optional
        Path to the condition_data.mat file
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
    """
    # Convert all paths to Path objects
    matlab_data_file_path = Path(matlab_data_file_path)
    experiment_info_file_path = Path(experiment_info_file_path)
    output_dir = Path(output_dir)
    
    # Convert optional paths to Path objects if provided
    if video_file_path is not None:
        video_file_path = Path(video_file_path)
    if tiff_folder_path is not None:
        tiff_folder_path = Path(tiff_folder_path)
    if df_f_file_path is not None:
        df_f_file_path = Path(df_f_file_path)
    if roi_info_file_path is not None:
        roi_info_file_path = Path(roi_info_file_path)
    if visual_stimuli_file_path is not None:
        visual_stimuli_file_path = Path(visual_stimuli_file_path)
    if trial_data_file_path is not None:
        trial_data_file_path = Path(trial_data_file_path)
    if condition_data_file_path is not None:
        condition_data_file_path = Path(condition_data_file_path)

    # Validate required input files exist
    if not matlab_data_file_path.is_file():
        raise FileNotFoundError(f"Matlab data file not found at {matlab_data_file_path}")
    if not experiment_info_file_path.is_file():
        raise FileNotFoundError(f"Experiment info file not found at {experiment_info_file_path}")
    
    # Validate optional input files exist if provided
    if video_file_path is not None and not video_file_path.is_file():
        raise FileNotFoundError(f"Video file not found at {video_file_path}")
    if tiff_folder_path is not None and not tiff_folder_path.exists():
        raise FileNotFoundError(f"Tiff folder not found at {tiff_folder_path}")
    if df_f_file_path is not None and not df_f_file_path.is_file():
        raise FileNotFoundError(f"df_f.mat file not found at {df_f_file_path}")
    if roi_info_file_path is not None and not roi_info_file_path.is_file():
        raise FileNotFoundError(f"ROI_info.mat file not found at {roi_info_file_path}")
    if visual_stimuli_file_path is not None and not visual_stimuli_file_path.is_file():
        raise FileNotFoundError(f"visual_stimuli.mat file not found at {visual_stimuli_file_path}")
    if trial_data_file_path is not None and not trial_data_file_path.is_file():
        raise FileNotFoundError(f"trial_data.mat file not found at {trial_data_file_path}")
    if condition_data_file_path is not None and not condition_data_file_path.is_file():
        raise FileNotFoundError(f"condition_data.mat file not found at {condition_data_file_path}")

    # Load data
    mat_data = read_mat(matlab_data_file_path)
    data = mat_data["data"]



    protocol = mat_data["protocol"]
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

    # Unpack data
    time = data[0]
    left_wingbeat = data[1]
    left_right_wingbeat = data[2]
    x_position = data[3]
    y_position = data[4]
    two_photon_frame_sync = data[5]
    behavior_camera_sync = data[6]
    stimulus_start = data[7]  # Not used in this example, empty for the example data

    # Initialize data interfaces dictionary
    data_interfaces = {}
    
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
        two_photon_timestamps = time[detect_threshold_crossings(two_photon_frame_sync, 0.5)]
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
            timestamps=time,
            verbose=verbose,
        )
        data_interfaces["roi"] = roi_interface

    # Set up stimuli interface if visual_stimuli_file_path is provided
    if visual_stimuli_file_path is not None:
        stimuli_interface = KimLabStimuliInterface(
            file_path=visual_stimuli_file_path,
            timestamps=time,
            verbose=verbose,
        )
        data_interfaces["stimuli"] = stimuli_interface
    
    # Set up trials interface if trial_data_file_path and condition_data_file_path are provided
    if trial_data_file_path is not None and condition_data_file_path is not None:
        trials_interface = KimLabTrialsInterface(
            trial_data_file_path=trial_data_file_path,
            condition_data_file_path=condition_data_file_path,
            stimuli_blocks=10.0,
            timestamps=time,
            verbose=verbose,
        )
        data_interfaces["trials"] = trials_interface

    # Set up video interface if video_file_path is provided
    if video_file_path is not None:
        video_interface = VideoInterface(file_paths=[video_file_path])
        video_timestamps = time[detect_threshold_crossings(behavior_camera_sync, 0.5)]
        video_interface.set_aligned_timestamps([video_timestamps])
        data_interfaces["video"] = video_interface

    converter_pipe = ConverterPipe(data_interfaces=data_interfaces)

    california_tz = ZoneInfo("America/Los_Angeles")
    session_start_time = session_start_time.replace(tzinfo=california_tz)

    metadata = converter_pipe.get_metadata()
    metadata["NWBFile"]["session_start_time"] = session_start_time
    metadata["NWBFile"]["session_description"] = f"protocol: {protocol}"
    metadata["NWBFile"]["session_id"] = session_id
    metadata["Subject"]["species"] = "Drosophila melanogaster"
    metadata["Subject"]["strain"] = genotype
    metadata["Subject"]["age"] = f"P{age}D"
    metadata["Subject"]["subject_id"] = session_id
    metadata["Subject"]["sex"] = "F"

    # Set conversion options for interfaces that are present
    conversion_options = {}
    if "imaging_channel_1" in data_interfaces:
        conversion_options["imaging_channel_1"] = dict(stub_test=False, photon_series_index=0)
    if "imaging_channel_2" in data_interfaces:
        conversion_options["imaging_channel_2"] = dict(stub_test=False, photon_series_index=1)

    nwbfile = converter_pipe.create_nwbfile(
        metadata=metadata,
        conversion_options=conversion_options,
    )

    # Add timeseries data
    timeseries_wingbeat = TimeSeries(
        name="wingbeat",
        data=left_wingbeat,
        unit="n.a.",
        timestamps=time,
        description="wingbeat",
    )
    nwbfile.add_acquisition(timeseries_wingbeat)

    # Note that this pattern indicates that the timestamps are the same as wingbeat
    common_timestamps = timeseries_wingbeat
    timeseries_left_right_wingbeat = TimeSeries(
        name="left_right_wingbeat",
        data=left_right_wingbeat,
        unit="n.a.",
        timestamps=common_timestamps,
        description="left-right wingbeat",
    )
    nwbfile.add_acquisition(timeseries_left_right_wingbeat)

    timeseries_x_position = TimeSeries(
        name="stimulus_position",
        data=np.c_[x_position, y_position],
        unit="n.a.",
        timestamps=common_timestamps,
        description="position of the visual pattern",
    )
    nwbfile.add_acquisition(timeseries_x_position)

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

    # Configure and save the NWB file
    nwbfile_path = output_dir / f"{session_id}.nwb"
    configure_and_write_nwbfile(nwbfile=nwbfile, output_filepath=nwbfile_path, backend="hdf5")

    if verbose:
        print(f"Created NWB file: {nwbfile_path}")

    return nwbfile_path


if __name__ == "__main__":
    # Example usage with the actual paths
    data_folder_path = Path("/home/heberto/cohen_project/Sample data/Kim Lab/20250301b")
    # Define input paths
    matlab_data_file_path = data_folder_path / "data_20250301b_00007.mat"
    video_file_path = data_folder_path / "20250301b_00007.avi"
    experiment_info_file_path = data_folder_path / "exp_info.mat"
    tiff_folder_path = data_folder_path  # Contains tiff files like 20250301b_00007_00001.tif
    visual_stimuli_file_path = data_folder_path / "visual_stimuli.mat"
    trial_data_file_path = data_folder_path / "trial.mat"
    condition_data_file_path = data_folder_path / "condition.mat"
    output_dir = Path("/home/heberto/cohen_project/Sample data/Kim Lab/nwb")

    output_dir.mkdir(exist_ok=True, parents=True)

    nwbfile_path = convert_session_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        experiment_info_file_path=experiment_info_file_path,
        video_file_path=video_file_path,
        tiff_folder_path=tiff_folder_path,
        visual_stimuli_file_path=visual_stimuli_file_path,
        trial_data_file_path=trial_data_file_path,
        condition_data_file_path=condition_data_file_path,
        output_dir=output_dir,
        verbose=True,  # Enable verbose output for demonstration
    )
