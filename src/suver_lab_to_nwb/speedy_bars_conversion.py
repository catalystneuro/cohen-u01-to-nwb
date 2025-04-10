from pathlib import Path
import time
from datetime import datetime
import uuid
import zoneinfo
from typing import Optional, Union

import numpy as np
from neuroconv.datainterfaces import DeepLabCutInterface, ExternalVideoInterface
from pymatreader import read_mat
from pynwb import NWBFile
from pynwb.file import Subject
from pynwb.base import TimeSeries
from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile


def convert_speedy_bars_to_nwb(
    matlab_data_file_path: Union[str, Path],
    output_dir: Union[str, Path],
    dlc_file_path: Optional[Union[str, Path]] = None,
    video_file_path: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    stub_test: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Convert Suver Lab SpeedyBars behavioral data to NWB format.

    This function converts behavioral data from the SpeedyBars experiment to the NWB format.
    In this experiment, each fly experienced 60 optic flow trials ranging from 0–50 cm/s in 
    5 cm/s increments along with 100 cm/s (5 trials at each optic flow speed, pseudo-randomized 
    order). Each trial consisted of an 8 second period of constant optic flow.

    Hardware Configuration:
    ----------------------
    - National Instruments DAQ system (Dev2)
    - Input channels:
      - Tachometer signal: analog input channel 0 (ai0)
      - Photodiode signal: analog input channel 2 (ai2) with SingleEnded terminal configuration
      - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
    - Output channels:
      - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
      - Solenoid valve control: analog output channel 2 (ao2)
      - Camera trigger: analog output channel 3 (ao3)

    Visual Stimulation System:
    ------------------------
    - Uses a TCP/IP server (port 5000) to communicate with a separate visual stimulus program
    - Optic flow speeds: 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, and 100 cm/s
    - 5 trials per stimulus speed, randomized order (60 trials total)
    - Trial timing:
      - Adjustment time: 2 seconds
      - Recording time: 6 seconds
      - Total trial time: 8 seconds

    Synchronization:
    --------------
    - Uses a timer object for precise synchronization between visual stimuli and data acquisition
    - 6-second delay between trial setup and start
    - Sends start time to visual stimulus program via TCP/IP
    - Waits for confirmation from visual stimulus program before proceeding

    Video Recording:
    --------------
    - Camera: Dorsal view (labeled as "Coronal (x) Camera")
    - Resolution: 640x480 pixels, Y8 format
    - Frame rate: 60 fps
    - Manual shutter mode with shutter value of 800
    - Videos saved as Motion JPEG AVI files
    - Camera is hardware-triggered via the DAQ system

    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB file containing SpeedyBars data.
    output_dir : str or Path
        Path to the directory where the NWB file will be saved.
    dlc_file_path : str or Path, optional
        Path to the DeepLabCut file containing pose estimation data.
    video_file_path : str or Path, optional
        Path to the video file to be included as an external file reference.
    overwrite : bool, default: False
        Whether to overwrite existing NWB files.
    stub_test : bool, default: False
        If True, only the first five trials will be processed and saved to a stub NWB file.
    verbose : bool, default: True
        Whether to print verbose output.

    Returns
    -------
    nwbfile_path : Path
        Path to the generated NWB file.
    """
    # Start timing
    start_time = time.time()

    # Convert paths to Path objects
    matlab_data_file_path = Path(matlab_data_file_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if dlc_file_path is not None:
        dlc_file_path = Path(dlc_file_path)

    # Read the MATLAB data
    session_data = read_mat(matlab_data_file_path)["data"]

    # Get the first trial data to extract common information
    first_trial_data = {key: value[0] for key, value in session_data.items()}
    
    # Extract metadata
    session_date = first_trial_data.get("date")
    experiment_number = first_trial_data.get("expnumber")
    condition = first_trial_data.get("condition")
    min_age = first_trial_data.get("min_age")
    max_age = first_trial_data.get("max_age")
    sample_rate = float(first_trial_data.get("samplerate"))
    adjust_time = float(first_trial_data.get("adjust_time"))
    record_time = float(first_trial_data.get("record_time"))
    fps = float(first_trial_data.get("fps"))
    
    # Format age range
    age_range = f"P{min_age}D/P{max_age}D"
    
    # Transform session date (format may vary, assuming YYYY_MM_DD if available)
    if session_date:
        try:
            session_start_time = datetime.strptime(session_date, "%Y_%m_%d")
        except (ValueError, TypeError):
            # If date is not available or in unexpected format, use current date
            session_start_time = datetime.now()
    else:
        session_start_time = datetime.now()
        
    # Set timezone to Nashville (Vanderbilt University)
    vanderbilt_time_zone = zoneinfo.ZoneInfo("America/Chicago")
    session_start_time = session_start_time.replace(tzinfo=vanderbilt_time_zone)

    # Create session ID and description
    session_id = f"SpeedyBars_{condition}_{min_age}-{max_age}days"
    num_trials = len(session_data["trial"])

    # Create a detailed session description
    session_description = (
        f"SpeedyBars experiment with Drosophila melanogaster. "
        f"Each fly experienced 60 optic flow trials ranging from 0–50 cm/s in 5 cm/s "
        f"increments along with 100 cm/s (5 trials at each optic flow speed, pseudo-randomized "
        f"order). Each trial consisted of an 8 second period of constant optic flow. "
        f"For each optic flow speed, averages were taken over the last 6 seconds of each trial. "
        f"Condition: {condition}, Age range: {min_age}-{max_age} days."
    )

    # Create subject information
    subject = Subject(
        species="Drosophila melanogaster",
        description=f"Drosophila melanogaster, condition: {condition}",
        age=age_range,
        sex="U",  # Unknown sex
    )

    # Create NWB file
    identifier = str(uuid.uuid4())
    nwbfile = NWBFile(
        session_start_time=session_start_time,
        session_description=session_description,
        identifier=identifier,
        session_id=session_id,
        institution="Vanderbilt University",
        lab="Suver Lab",
        experimenter=["Suver, Marie"],
        subject=subject,
        keywords=["Drosophila", "behavior", "optic flow", "SpeedyBars"],
    )

    # Create file path for the NWB file
    if session_date and experiment_number:
        nwbfile_path = output_dir / f"{session_date}_E{experiment_number}_SpeedyBars.nwb"
    else:
        nwbfile_path = output_dir / f"{session_id}.nwb"
        
    # Check if file exists and handle overwrite
    if nwbfile_path.exists() and not overwrite:
        raise FileExistsError(f"File {nwbfile_path} already exists. Set overwrite=True to overwrite.")

    # Process trials
    num_trials_to_process = min(5, num_trials) if stub_test else num_trials
    
    # Create a trials table with columns for stimulus (optic flow speed) and timing information
    nwbfile.add_trial_column(name="stimulus", description="Optic flow speed in cm/s")
    nwbfile.add_trial_column(name="adjust_time", description="Adjustment time at each optic speed (seconds)")
    nwbfile.add_trial_column(name="record_time", description="Recording time used for averages (seconds)")
    nwbfile.add_trial_column(name="tcp_command", description="Command sent to the TCP/IP server for visual stimulation")
    
    # Define optic flow speeds
    optic_flow_speeds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]  # cm/s
    
    # Calculate number of digits needed for zero-padding based on total trials
    num_digits = len(str(num_trials))
    
    for trial_index in range(num_trials_to_process):
        # Extract trial data
        trial_data = {key: value[trial_index] for key, value in session_data.items()}
        trial_num = trial_data["trial"]
        stimulus = trial_data["stimulus"]  # Optic flow speed in cm/s
        
        # Format trial number with leading zeros based on total number of trials
        formatted_trial = f"{trial_num:0{num_digits}d}"
        
        # Calculate trial start time (assuming trials are sequential)
        trial_start_time = float(trial_index * 10.0)  # Approximate start time in seconds
        
        # Add trial to the trials table with additional columns
        nwbfile.add_trial(
            start_time=trial_start_time,
            stop_time=trial_start_time + (len(trial_data["pufferSignal"]) / sample_rate),
            stimulus=stimulus,
            adjust_time=adjust_time,
            record_time=record_time,
            tcp_command=stimulus  # In the original MATLAB code, the TCP command is the same as the speed
        )
        
        # Add puffer signal (air puff sensory stimulus)
        puffer_stimulus = TimeSeries(
            name=f"TimeSeriesPufferTrial{formatted_trial}",
            description="Air puff sensory stimulus delivery. Increases in volts indicate an air puff (sensory input).",
            data=trial_data["pufferSignal"],
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_stimulus(puffer_stimulus)
        
        # Add tachometer signal (wingbeat data)
        tachometer_series = TimeSeries(
            name=f"TimeSeriesTachometerTrial{formatted_trial}",
            description="Tachometer recording wingbeat data from fly. Recorded in volts.",
            data=trial_data["tachometerSignal"],
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_acquisition(tachometer_series)
        
        # Add photodiode signal (visual stimuli)
        photodiode_series = TimeSeries(
            name=f"TimeSeriesPhotodiodeTrial{formatted_trial}",
            description="Voltage reading from an amplified photodiode pointed at the visual stimuli. Signal increases when white bar is present relative to when black bar is present.",
            data=trial_data["photodiodeSignal"],
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_acquisition(photodiode_series)
    
    # Add DeepLabCut data if available
    if dlc_file_path is not None and dlc_file_path.exists():
        # Calculate timestamps for DLC data based on video frame rate
        # This is an approximation and would need to be adjusted based on actual synchronization
        dlc_timestamps = np.arange(0, num_trials_to_process * 10, 1/fps)
        
        dlc_interface = DeepLabCutInterface(file_path=dlc_file_path)
        dlc_interface.set_aligned_timestamps(dlc_timestamps)
        dlc_interface.add_to_nwbfile(nwbfile=nwbfile, container_name="PoseEstimation")
    
    # Add video data if available
    if video_file_path is not None and Path(video_file_path).exists():
        # Calculate timestamps for video data based on video frame rate
        video_timestamps = np.arange(0, num_trials_to_process * 10, 1/fps)
        
        video_interface = ExternalVideoInterface(file_paths=[video_file_path])
        video_interface.set_aligned_timestamps(video_timestamps)
        video_interface.add_to_nwbfile(nwbfile=nwbfile, container_name="BehavioralVideo")

    # Save the NWB file
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path)

    # Print conversion time
    if verbose:
        end_time = time.time()
        duration = end_time - start_time
        if duration < 60:
            print(f"Conversion completed in {duration:.2f} seconds")
        elif duration < 3600:
            print(f"Conversion completed in {duration/60:.2f} minutes")
        else:
            print(f"Conversion completed in {duration/3600:.2f} hours")

        print(f"NWB file saved to: {nwbfile_path}")

    return nwbfile_path


if __name__ == "__main__":
    # Example usage
    data_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/behavioralDataExamples/")
    output_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/nwb")
    
    # Path to the SpeedyBars data file
    matlab_data_file_path = data_dir / "2024_05_29_E3.mat"
    
    # Path to DeepLabCut data (commented out as it's not available yet)
    # dlc_file_path = data_dir / "DeepLabCut" / "speedy_bars_dlc.h5"
    
    # Path to video file (commented out as it's not available yet)
    # video_file_path = data_dir / "videos" / "2024_05_29_E3_Video_Dorsal.avi"
    
    # Run conversion
    convert_speedy_bars_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        output_dir=output_dir,
        # dlc_file_path=dlc_file_path,  # Uncomment when DLC data becomes available
        # video_file_path=video_file_path,  # Uncomment when video data becomes available
        overwrite=True,
        verbose=True,
    )
