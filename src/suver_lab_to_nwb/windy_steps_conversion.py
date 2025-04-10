from pathlib import Path
import time
from datetime import datetime
import uuid
import zoneinfo
from typing import Optional, Union
import math

import numpy as np
from neuroconv.datainterfaces import DeepLabCutInterface, ExternalVideoInterface
from pymatreader import read_mat
from pynwb import NWBFile
from pynwb.file import Subject
from pynwb.base import TimeSeries
from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile


def convert_windy_steps_to_nwb(
    matlab_data_file_path: Union[str, Path],
    output_dir: Union[str, Path],
    dlc_file_path: Optional[Union[str, Path]] = None,
    video_file_path: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    stub_test: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Convert Suver Lab WindySteps behavioral data to NWB format.

    This function converts behavioral data from the WindySteps experiment to the NWB format.
    In this experiment, each fly experienced 10 airflow trials either increasing from 0 cm/s 
    to 300 cm/s (5 trials) or decreasing from 300 cm/s to 0 cm/s (5 trials). Airflow was 
    presented in 50 cm/s steps over 42 seconds, with 6 seconds of steady state airflow for 
    each windspeed.

    Hardware Configuration:
    ----------------------
    - National Instruments DAQ system (Dev2)
    - Input channels:
      - Tachometer signal: analog input channel 0 (ai0)
      - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
    - Output channels:
      - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
      - Optogenetic stimulation (when used): analog output channel 1 (ao1)
      - Solenoid valve control: analog output channel 2 (ao2)
      - Camera trigger: analog output channel 3 (ao3)

    Wind Speed Calculation:
    ----------------------
    Wind speeds (0, 50, 100, 150, 200, 250, 300 cm/s) are converted to voltage values for the MFC using:
    - Wind tube radius: 0.18796 cm
    - Tube area = π × r² (in cm²)
    - Flow rate (cm³/s) = windspeed (cm/s) × tube area (cm²)
    - Flow rate (cm³/min) = flow rate (cm³/s) × 60
    - Flow rate (L/min) = flow rate (cm³/min) / 1000
    - MFC voltage = flow rate (L/min) × 5/2

    Signal Processing:
    ----------------
    - Tachometer signal is filtered using a 2nd-order Butterworth bandpass filter (150-250 Hz)
    - Filtered signal is further smoothed using Savitzky-Golay filtering

    Experimental Conditions:
    ----------------------
    The script supports various experimental conditions including:
    - Normal flies
    - Silenced flies (with optogenetic manipulation)
    - Dark conditions
    - Various genetic lines (18D07, 74C10, ChrimCS, silencedCS, silencedCS_glued)

    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB file containing WindySteps data.
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
    adjust_time = first_trial_data.get("adjust_time")
    record_time = first_trial_data.get("record_time")
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
    session_id = f"WindySteps_{condition}_{min_age}-{max_age}days"
    num_trials = len(session_data["trial"])

    # Create a detailed session description
    session_description = (
        f"WindySteps experiment with Drosophila melanogaster. "
        f"Each fly experienced 10 airflow trials either increasing from 0 cm/s to 300 cm/s "
        f"(5 trials) or decreasing from 300 cm/s to 0 cm/s (5 trials). Airflow was presented "
        f"in 50 cm/s steps over 42 seconds, with 6 seconds of steady state airflow for each "
        f"windspeed. For each windspeed, averages were taken over the last 3 seconds as the "
        f"mass flow controller did not adjust instantaneously. Condition: {condition}, "
        f"Age range: {min_age}-{max_age} days."
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
        keywords=["Drosophila", "behavior", "airflow", "WindySteps"],
    )

    # Create file path for the NWB file
    if session_date and experiment_number:
        nwbfile_path = output_dir / f"{session_date}_E{experiment_number}_WindySteps.nwb"
    else:
        nwbfile_path = output_dir / f"{session_id}.nwb"
        
    # Check if file exists and handle overwrite
    if nwbfile_path.exists() and not overwrite:
        raise FileExistsError(f"File {nwbfile_path} already exists. Set overwrite=True to overwrite.")

    # Process trials
    num_trials_to_process = min(5, num_trials) if stub_test else num_trials
    
    # Create a trials table with columns for trial type (ascending/descending) and wind speed information
    nwbfile.add_trial_column(name="stim_type", description="1 for ascending/increasing stimuli, 2 for descending/decreasing stimuli")
    nwbfile.add_trial_column(name="wind_speed", description="Wind speed in cm/s (0, 50, 100, 150, 200, 250, 300)")
    nwbfile.add_trial_column(name="mfc_voltage", description="Voltage sent to Mass Flow Controller")
    nwbfile.add_trial_column(name="is_optogenetic", description="Whether optogenetic stimulation was used")
    nwbfile.add_trial_column(name="condition_details", description="Detailed information about the experimental condition")
    
    # Define wind speeds and calculate MFC voltages
    wind_speeds = [0, 50, 100, 150, 200, 250, 300]  # cm/s
    tube_radius = 0.18796  # cm
    tube_area = math.pi * (tube_radius ** 2)  # cm²
    
    # Function to calculate MFC voltage from wind speed
    def calculate_mfc_voltage(wind_speed):
        if wind_speed == 0:
            return 0.0
        flow_rate_cm3_per_sec = wind_speed * tube_area  # cm³/s
        flow_rate_cm3_per_min = flow_rate_cm3_per_sec * 60  # cm³/min
        flow_rate_L_per_min = flow_rate_cm3_per_min / 1000  # L/min
        mfc_voltage = flow_rate_L_per_min * (5/2)  # V
        return mfc_voltage
    
    # Check if condition involves optogenetic stimulation
    is_optogenetic = "silenced" in condition.lower() if condition else False
    
    # Determine condition details based on condition string
    condition_details = "Unknown"
    if condition:
        if "dark" in condition.lower():
            condition_details = "Dark conditions"
        elif "silenced" in condition.lower():
            condition_details = "Silenced flies with optogenetic manipulation"
        elif any(genotype in condition for genotype in ["18D07", "74C10", "ChrimCS", "silencedCS", "silencedCS_glued"]):
            condition_details = f"Genetic line: {condition}"
        else:
            condition_details = "Normal flies"
    
    # Calculate number of digits needed for zero-padding based on total trials
    num_digits = len(str(num_trials))
    
    for trial_index in range(num_trials_to_process):
        # Extract trial data
        trial_data = {key: value[trial_index] for key, value in session_data.items()}
        trial_num = trial_data["trial"]
        stim_type = trial_data["stimType"]  # 1 for ascending, 2 for descending
        
        # Format trial number with leading zeros based on total number of trials
        formatted_trial = f"{trial_num:0{num_digits}d}"
        
        # Calculate trial start time (assuming trials are sequential)
        trial_start_time = float(trial_index * 60.0)  # Approximate start time in seconds
        
        # Determine wind speed based on trial index and stim_type
        # For ascending trials (stim_type=1), wind speed increases from 0 to 300 cm/s
        # For descending trials (stim_type=2), wind speed decreases from 300 to 0 cm/s
        if stim_type == 1:  # ascending
            # Calculate current step in the trial (0-6 for 7 wind speeds)
            # This is a simplification - in reality, the wind speed changes during the trial
            step = min(int(trial_index / (num_trials_to_process / 7)), 6)
            wind_speed = wind_speeds[step]
        else:  # descending
            step = 6 - min(int(trial_index / (num_trials_to_process / 7)), 6)
            wind_speed = wind_speeds[step]
        
        # Calculate MFC voltage
        mfc_voltage = calculate_mfc_voltage(wind_speed)
        
        # Add trial to the trials table with additional columns
        nwbfile.add_trial(
            start_time=trial_start_time,
            stop_time=trial_start_time + (len(trial_data["pufferSignal"]) / sample_rate),
            stim_type=stim_type,
            wind_speed=wind_speed,
            mfc_voltage=mfc_voltage,
            is_optogenetic=is_optogenetic,
            condition_details=condition_details
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
        
        # Add smoothed tachometer signal
        tachometer_smoothed_series = TimeSeries(
            name=f"TimeSeriesTachometerSmoothedTrial{formatted_trial}",
            description="Tachometer signal smoothed using a Butterworth bandpass filter (150-250 Hz) followed by a second-order Savitzky-Golay filter.",
            data=trial_data["tachometerSignal_smoothed"],
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_acquisition(tachometer_smoothed_series)
    
    # Add DeepLabCut data if available
    if dlc_file_path is not None and dlc_file_path.exists():
        # Calculate timestamps for DLC data based on video frame rate
        # This is an approximation and would need to be adjusted based on actual synchronization
        dlc_timestamps = np.arange(0, num_trials_to_process * 60, 1/fps)
        
        dlc_interface = DeepLabCutInterface(file_path=dlc_file_path)
        dlc_interface.set_aligned_timestamps(dlc_timestamps)
        dlc_interface.add_to_nwbfile(nwbfile=nwbfile, container_name="PoseEstimation")
    
    # Add video data if available
    if video_file_path is not None and Path(video_file_path).exists():
        # Calculate timestamps for video data based on video frame rate
        video_timestamps = np.arange(0, num_trials_to_process * 60, 1/fps)
        
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
    
    # Path to the WindySteps data file
    matlab_data_file_path = data_dir / "2024_07_08_E3.mat"
    
    # Path to DeepLabCut data (commented out as it's not available yet)
    # dlc_file_path = data_dir / "DeepLabCut" / "windy_steps_dlc.h5"
    
    # Path to video file (commented out as it's not available yet)
    # video_file_path = data_dir / "videos" / "2024_07_08_E3_Video_Dorsal.avi"
    
    # Run conversion
    convert_windy_steps_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        output_dir=output_dir,
        # dlc_file_path=dlc_file_path,  # Uncomment when DLC data becomes available
        # video_file_path=video_file_path,  # Uncomment when video data becomes available
        overwrite=True,
        verbose=True,
    )
