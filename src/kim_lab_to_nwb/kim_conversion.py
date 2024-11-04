from datetime import datetime
from pathlib import Path
import uuid

import numpy as np
from neuroconv.tools.nwb_helpers import get_default_backend_configuration, configure_backend
from neuroconv.datainterfaces import VideoInterface
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from pynwb.file import Subject
from pymatreader import read_mat

from kim_lab_to_nwb.ophys import MultiTiffMultiPageTiffImagingInterface
from cohen_u01_nwb_conversion_utils.utils import detect_threshold_crossings


def convert_session_to_nwb(
    matlab_data_file_path: str | Path,
    video_file_path: str | Path,
    experiment_info_file_path: str | Path,
    tiff_folder_path: str | Path,
    output_dir: str | Path,
    verbose: bool = False,
) -> Path:
    """
    Convert a single experimental session from Kim Lab data to NWB format.

    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB data file
    video_file_path : str or Path
        Path to the video file
    experiment_info_file_path : str or Path
        Path to the experiment info file
    tiff_folder_path : str or Path
        Path to the folder containing TIFF files
    output_dir : str or Path
        Directory where the NWB file will be saved
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
    video_file_path = Path(video_file_path)
    experiment_info_file_path = Path(experiment_info_file_path)
    tiff_folder_path = Path(tiff_folder_path)
    output_dir = Path(output_dir)

    # Validate input files exist
    if not matlab_data_file_path.is_file():
        raise FileNotFoundError(f"Matlab data file not found at {matlab_data_file_path}")
    if not video_file_path.is_file():
        raise FileNotFoundError(f"Video file not found at {video_file_path}")
    if not experiment_info_file_path.is_file():
        raise FileNotFoundError(f"Experiment info file not found at {experiment_info_file_path}")
    if not tiff_folder_path.exists():
        raise FileNotFoundError(f"Tiff folder not found at {tiff_folder_path}")

    # Load data
    mat_data = read_mat(matlab_data_file_path)
    data = mat_data["data"]
    protocol = mat_data["protocol"]

    experiment_info = read_mat(experiment_info_file_path)
    age = experiment_info["age"]
    genotype = experiment_info["cross"]

    # Extract session ID from matlab file path
    session_id = matlab_data_file_path.stem.split('data_')[1]
    if verbose:
        print(f"Processing {session_id=} ({age=}, {genotype=})")

    session_start_time = datetime.strptime(session_id[:8], '%Y%m%d')
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
    stimulus_start = data[7]

    # Create NWB file
    nwbfile = NWBFile(
        session_description=f"protocol: {protocol}",
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time,
        session_id=session_id,
    )

    # Set up imaging interface
    ophys_interface = MultiTiffMultiPageTiffImagingInterface(
        tiff_folder_path,
        pattern=session_id + "_{frame:05d}.tif",
        sampling_frequency=30.0,
        verbose=verbose
    )

    aligned_timestamps = time[detect_threshold_crossings(two_photon_frame_sync, 0.5)]
    aligned_timestamps = aligned_timestamps[:ophys_interface.imaging_extractor.get_num_frames()]
    ophys_interface.set_aligned_timestamps(aligned_timestamps=aligned_timestamps)
    ophys_interface.add_to_nwbfile(nwbfile, metadata=dict())

    # Set up video interface
    video_interface = VideoInterface(
        file_paths=[video_file_path],
    )

    video_timestamps = time[detect_threshold_crossings(behavior_camera_sync, 0.5)]
    video_interface.set_aligned_timestamps([video_timestamps])
    video_interface.add_to_nwbfile(nwbfile, metadata=dict())

    # Add subject information
    nwbfile.subject = Subject(
        subject_id=session_id,
        genotype=genotype,
        age=f"P{age}D",
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

    timeseries_left_right_wingbeat = TimeSeries(
        name="left_right_wingbeat",
        data=left_right_wingbeat,
        unit="n.a.",
        timestamps=timeseries_wingbeat,
        description="left-right wingbeat",
    )
    nwbfile.add_acquisition(timeseries_left_right_wingbeat)

    timeseries_x_position = TimeSeries(
        name="stimulus_position",
        data=np.c_[x_position, y_position],
        unit="n.a.",
        timestamps=timeseries_wingbeat,
        description="position of the visual pattern",
    )
    nwbfile.add_acquisition(timeseries_x_position)

    # Configure and save the NWB file
    backend_configuration = get_default_backend_configuration(nwbfile, backend="hdf5")
    configure_backend(nwbfile=nwbfile, backend_configuration=backend_configuration)

    nwbfile_path = output_dir / f"{session_id}.nwb"
    with NWBHDF5IO(nwbfile_path, mode="w") as io:
        io.write(nwbfile)
    
    if verbose:
        print(f"Created NWB file: {nwbfile_path}")
    
    return nwbfile_path


if __name__ == "__main__":
    # Example usage with the original paths
    data_folder_path = Path("/Users/heberto/project_data/Sample data-selected/Kim Lab")
    
    # Define input paths
    matlab_data_file_path = data_folder_path / "raw data" / "data_20240108b_00003.mat"
    video_file_path = data_folder_path / "raw data" / "20240108b_00003.avi"
    experiment_info_file_path = data_folder_path / "raw data" / "exp_info.mat"
    tiff_folder_path = data_folder_path / "raw data"
    output_dir = data_folder_path / "nwb"

    output_dir.mkdir(exist_ok=True, parents=True)
    
    nwbfile_path = convert_session_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        video_file_path=video_file_path,
        experiment_info_file_path=experiment_info_file_path,
        tiff_folder_path=tiff_folder_path,
        output_dir=output_dir,
        verbose=True  # Enable verbose output for demonstration
    )
