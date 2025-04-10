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


def convert_coco_to_nwb(
    matlab_data_file_path: Union[str, Path],
    output_dir: Union[str, Path],
    dlc_file_path: Optional[Union[str, Path]] = None,
    video_file_path: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    stub_test: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Convert Suver Lab Coco behavioral data to NWB format.

    This function converts behavioral data from the Coco experiment to the NWB format.
    In this experiment, each fly experienced 36 trials at 3 oscillatory frequencies 
    (0.3, 1.3, 2.3 Hz, 12 trials at each frequency) with windspeeds ranging from 
    100–200 cm/s and optic flow speeds from 5–35 cm/s.

    Hardware Configuration:
    ----------------------
    - National Instruments DAQ system (Dev2)
    - Input Channels:
      - Tachometer signal: analog input channel 0 (ai0)
      - Puffer signal: analog input channel 5 (ai5) with SingleEnded terminal configuration
      - Photodiode signal: analog input channel 6 (ai6) with SingleEnded terminal configuration
    - Output Channels:
      - Mass Flow Controller (MFC) signal: analog output channel 0 (ao0)
      - Solenoid valve control: analog output channel 2 (ao2)
      - Camera trigger: analog output channel 3 (ao3)
    - Uses external digital trigger (PFI0) with rising edge condition

    Experimental Conditions:
    ----------------------
    - Four possible conditions: 'wind', 'visual', 'both', or 'none'
    - Condition determines whether wind oscillation, visual oscillation, both, or neither is applied

    Oscillation Parameters:
    ---------------------
    - Three oscillation frequencies: 0.3, 1.3, and 2.3 Hz
    - Corresponding TCP frequencies for visual system: 3, 13, and 23
    - Amplitude voltages: 0.835, 1.335, and 2.09 V
    - Base voltage (150 cm/s wind speed): 2.4972571 V
    - Wind oscillation formula: ((amplitude * sin(2*pi*t*frequency))+baseV)

    Trial Structure:
    --------------
    - Total trial time: 19 seconds
    - Video recording time: 16 seconds
    - Baseline time: 6 seconds (constant wind speed)
    - Oscillation time: 12 seconds
    - End time: 1 second

    Experimental Design:
    -----------------
    - 12 blocks with 3 trials per block (36 trials total)
    - Each block contains all three frequencies in random order
    - Block and trial-within-block structure is explicitly tracked

    Visual Stimulation System:
    ------------------------
    - Uses a TCP/IP server (port 5000) to communicate with a separate visual stimulus program
    - Visual stimulus frequencies are mapped to TCP frequencies (3, 13, 23)
    - For 'wind' or 'none' conditions, a special value (100) is sent to indicate no visual stimulus

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
        Path to the MATLAB file containing Coco data.
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
    fps = float(first_trial_data.get("fps"))
    trial_length = float(first_trial_data.get("trialLength"))
    
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
    session_id = f"Coco_{condition}_{min_age}-{max_age}days"
    num_trials = len(session_data["trialnum"])

    # Create a detailed session description
    session_description = (
        f"Coco experiment with Drosophila melanogaster. "
        f"Each fly experienced 36 trials at 3 oscillatory frequencies "
        f"(0.3, 1.3, 2.3 Hz, 12 trials at each frequency) with windspeeds ranging from "
        f"100–200 cm/s and optic flow speeds from 5–35 cm/s. "
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
        keywords=["Drosophila", "behavior", "oscillatory", "Coco"],
    )

    # Create file path for the NWB file
    if session_date and experiment_number:
        nwbfile_path = output_dir / f"{session_date}_E{experiment_number}_Coco.nwb"
    else:
        nwbfile_path = output_dir / f"{session_id}.nwb"
        
    # Check if file exists and handle overwrite
    if nwbfile_path.exists() and not overwrite:
        raise FileExistsError(f"File {nwbfile_path} already exists. Set overwrite=True to overwrite.")

    # Process trials
    num_trials_to_process = min(5, num_trials) if stub_test else num_trials
    
    # Create a trials table with columns for block, block_trial, stimulus, and oscillation parameters
    nwbfile.add_trial_column(name="block", description="Current block for current trial")
    nwbfile.add_trial_column(name="block_trial", description="Trial within the current block (1, 2, or 3)")
    nwbfile.add_trial_column(name="stimulus", description="Current oscillation frequency being presented (Hz)")
    nwbfile.add_trial_column(name="oscillation_type", description="Type of oscillation ('wind', 'visual', 'both', 'none')")
    nwbfile.add_trial_column(name="tcp_frequency", description="Corresponding TCP frequency for visual system (3, 13, 23)")
    nwbfile.add_trial_column(name="amplitude_voltage", description="Amplitude voltage used (0.835, 1.335, 2.09 V)")
    nwbfile.add_trial_column(name="base_voltage", description="Base voltage (2.4972571 V)")
    nwbfile.add_trial_column(name="baseline_time", description="Baseline time with constant wind speed (seconds)")
    nwbfile.add_trial_column(name="oscillation_time", description="Oscillation time (seconds)")
    
    # Define oscillation parameters
    oscillation_frequencies = [0.3, 1.3, 2.3]  # Hz
    tcp_frequencies = [3, 13, 23]  # Corresponding TCP frequencies
    amplitude_voltages = [0.835, 1.335, 2.09]  # V
    base_voltage = 2.4972571  # V (150 cm/s wind speed)
    baseline_time = 6.0  # seconds
    oscillation_time = 12.0  # seconds
    
    # Determine oscillation type based on condition
    oscillation_type = condition.lower() if condition else "unknown"
    if oscillation_type not in ["wind", "visual", "both", "none"]:
        oscillation_type = "unknown"
    
    # Calculate number of digits needed for zero-padding based on total trials
    num_digits = len(str(num_trials))
    
    for trial_index in range(num_trials_to_process):
        # Extract trial data
        trial_data = {key: value[trial_index] for key, value in session_data.items()}
        trial_num = trial_data["trialnum"]
        block = trial_data["block"]
        block_trial = trial_data["block_trial"]
        stimulus = trial_data["stimulus"]  # Oscillation frequency in Hz
        
        # Format trial number with leading zeros based on total number of trials
        formatted_trial = f"{trial_num:0{num_digits}d}"
        
        # Calculate trial start time (assuming trials are sequential)
        trial_start_time = float(trial_index * 20.0)  # Approximate start time in seconds
        
        # Find the index of the current oscillation frequency in the list
        try:
            freq_index = oscillation_frequencies.index(stimulus)
            tcp_frequency = tcp_frequencies[freq_index]
            amplitude_voltage = amplitude_voltages[freq_index]
        except ValueError:
            # If the frequency is not in the list, use default values
            tcp_frequency = 0
            amplitude_voltage = 0
        
        # Add trial to the trials table with additional columns
        nwbfile.add_trial(
            start_time=trial_start_time,
            stop_time=trial_start_time + (len(trial_data["pufferSignal"]) / sample_rate),
            block=block,
            block_trial=block_trial,
            stimulus=stimulus,
            oscillation_type=oscillation_type,
            tcp_frequency=tcp_frequency,
            amplitude_voltage=amplitude_voltage,
            base_voltage=base_voltage,
            baseline_time=baseline_time,
            oscillation_time=oscillation_time
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
        dlc_timestamps = np.arange(0, num_trials_to_process * 20, 1/fps)
        
        dlc_interface = DeepLabCutInterface(file_path=dlc_file_path)
        dlc_interface.set_aligned_timestamps(dlc_timestamps)
        dlc_interface.add_to_nwbfile(nwbfile=nwbfile, container_name="PoseEstimation")
    
    # Add video data if available
    if video_file_path is not None and Path(video_file_path).exists():
        # Calculate timestamps for video data based on video frame rate
        video_timestamps = np.arange(0, num_trials_to_process * 20, 1/fps)
        
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
    
    # Path to the Coco data file
    matlab_data_file_path = data_dir / "2024_10_28_E4.mat"
    
    # Path to DeepLabCut data (commented out as it's not available yet)
    # dlc_file_path = data_dir / "DeepLabCut" / "coco_dlc.h5"
    
    # Path to video file (commented out as it's not available yet)
    # video_file_path = data_dir / "videos" / "2024_10_28_E4_Video_Dorsal.avi"
    
    # Run conversion
    convert_coco_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        output_dir=output_dir,
        # dlc_file_path=dlc_file_path,  # Uncomment when DLC data becomes available
        # video_file_path=video_file_path,  # Uncomment when video data becomes available
        overwrite=True,
        verbose=True,
    )
