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
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries
from pynwb.base import TimeSeries
from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile


def calculate_seal_test_metrics(current, voltage):
    """
    Calculate metrics for a seal test based on current and voltage traces.
    
    This is a simplified version of the ComputeCellStats.m algorithm.
    
    Parameters
    ----------
    current : numpy.ndarray
        Current trace data
    voltage : numpy.ndarray
        Voltage trace data
    
    Returns
    -------
    dict
        Dictionary containing calculated metrics
    """
    # Constants from ComputeCellStats.m
    THRESH = 3  # Threshold for detecting pulse edges
    SCALE = 10  # Scale factor for resistance calculations
    
    # Find pulse onset/offset
    on_indices = np.where(np.diff(voltage) > THRESH)[0]
    off_indices = np.where(np.diff(voltage) < -THRESH)[0]
    
    # Skip if we can't find clear pulses
    if len(on_indices) == 0 or len(off_indices) == 0:
        return {
            "ra": np.nan,
            "rinput": np.nan,
            "i0": np.nan,
            "im_pulse": np.nan,
            "vm_pulse": np.nan
        }
    
    # Ensure off_indices come after on_indices
    while off_indices[0] < on_indices[0]:
        off_indices = off_indices[1:]
        if len(off_indices) == 0:
            return {
                "ra": np.nan,
                "rinput": np.nan,
                "i0": np.nan,
                "im_pulse": np.nan,
                "vm_pulse": np.nan
            }
    
    # Get baseline values (before pulse)
    pre_samples = 50  # Number of samples before pulse to use for baseline
    if on_indices[0] > pre_samples:
        base_vm = np.mean(voltage[:on_indices[0]-5])
        base_im = np.mean(current[:on_indices[0]-5])
    else:
        base_vm = voltage[0]
        base_im = current[0]
    
    # Get steady-state values during pulse
    pulse_length = off_indices[0] - on_indices[0]
    center_pulse = pulse_length // 2
    fraction = 0.25  # Use middle 50% of pulse
    
    section_begin = on_indices[0] + center_pulse - int(fraction * pulse_length)
    section_end = on_indices[0] + center_pulse + int(fraction * pulse_length)
    
    # Ensure indices are within bounds
    section_begin = max(0, section_begin)
    section_end = min(len(voltage) - 1, section_end)
    
    vm_pulse = np.mean(voltage[section_begin:section_end]) - base_vm
    im_pulse = np.mean(current[section_begin:section_end]) - base_im
    
    # Find peak current
    i0_idx = np.argmax(current[on_indices[0]:off_indices[0]]) + on_indices[0]
    i0 = current[i0_idx] - base_im
    
    # Calculate resistances
    if im_pulse != 0:
        rinput = vm_pulse / im_pulse * SCALE  # Input resistance in MOhm
    else:
        rinput = np.nan
    
    if i0 != 0:
        ga = i0 / vm_pulse  # Conductance in nS
        ra = 1 / ga * SCALE if ga != 0 else np.nan  # Access resistance in MOhm
    else:
        ra = np.nan
    
    return {
        "ra": ra,
        "rinput": rinput,
        "i0": i0,
        "im_pulse": im_pulse,
        "vm_pulse": vm_pulse
    }


def convert_session_to_nwb(
    matlab_data_file_path: Union[str, Path],
    output_dir: Union[str, Path],
    data_dir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    stub_test: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Convert Suver Lab data to NWB format.

    This function converts patch clamp data, videos, and DeepLabCut pose estimation
    from the Suver Lab to the NWB format. All trials from the same experiment are
    included in a single NWB file, as they were acquired in the same run from the
    same fly and together comprise a single experiment.

    The function uses specific file naming conventions:
    - Video files: {date}_E{experiment_number}_Video_lateral_flyLeft{trial}.avi
                {date}_E{experiment_number}_Video_lateral_flyRight_{trial}.avi
                {date}_E{experiment_number}_Video_lateral_ventral_{trial}.avi
    - DeepLabCut files: {date}_E{experiment_number}_Video_lateral_flyLeft{trial}DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5

    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB file containing patch clamp data.
    output_dir : str or Path
        Path to the directory where the NWB file will be saved.
    data_dir : str or Path, optional
        Path to the directory containing the data files. If None, it will be inferred
        from the matlab_data_file_path.
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

    # Infer data_dir if not provided
    if data_dir is None:
        data_dir = matlab_data_file_path.parent
    else:
        data_dir = Path(data_dir)

    # Read the MATLAB data
    session_data = read_mat(matlab_data_file_path)["data"]

    # Get the first trial data to extract common information
    first_trial_data = {key: value[0] for key, value in session_data.items()}
    session_date = first_trial_data["date"]
    experiment_number = first_trial_data["expNumber"]
    condition = first_trial_data["condition"]
    age = first_trial_data["age"]
    genotype = first_trial_data["genotype"]

    # Transform session date (format: 2024_10_24) to datetime object
    session_start_time = datetime.strptime(session_date, "%Y_%m_%d")
    vanderbilt_time_zone = zoneinfo.ZoneInfo("America/Chicago")
    session_start_time = session_start_time.replace(tzinfo=vanderbilt_time_zone)

    # Create session ID and description
    session_id = f"{condition}_{age}_{genotype}"
    num_trials = len(session_data["trial"])

    # Create a more detailed session description based on the email information
    session_description = (
        f"In-vivo whole-cell patch clamp recordings in Drosophila melanogaster "
        f"(UAS-10xGFP;GMR-24C06-GAL4), while simultaneously recording videos of the "
        f"animal's antennae with cameras, wingbeat data with the tachometer, and "
        f"sensory input with the puffer. All trials were acquired in the same run "
        f"from the same fly on {session_date}, experiment {experiment_number}."
    )

    # Create subject information
    subject = Subject(
        species="Drosophila melanogaster",
        description="UAS-10xGFP;GMR-24C06-GAL4 raised on standard food and on a 12h light/dark cycle",
        age=age,
        genotype=genotype,
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
    )

    # Create device with detailed information from the email
    device = nwbfile.create_device(
        name="AMMC Patch Clamp Amplifier 2400",
        description="AM-Systems 2400 patch clamp amplifier with 10xVm port for filtered and amplified recordings",
    )

    # Define the DeepLabCut suffix based on the provided file patterns
    dlc_suffix = "DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5"

    # First, process seal test 15 to get the seal resistance value
    seal_tests_dir = data_dir / "seal tests"
    seal_test_files = list(seal_tests_dir.glob(f"{session_date}_SealTest_*.mat"))
    
    # Default seal value
    seal_value = "TBD"
    
    # Look for seal test 15 to get the seal resistance
    for seal_test_file in seal_test_files:
        seal_test_num = seal_test_file.stem.split('_')[-1]
        if seal_test_num == "15":
            # Read the seal test data
            seal_test_trial_data = read_mat(seal_test_file)["data"]
            
            # Extract the first (and only) trial data            
            # Get the current and voltage data
            seal_test_current = seal_test_trial_data["I"]
            seal_test_voltage = seal_test_trial_data["Vm"]
            seal_test_sample_rate = float(seal_test_trial_data["SAMPLERATE"])
            
            # Calculate metrics
            metrics = calculate_seal_test_metrics(seal_test_current, seal_test_voltage)
            
            # Set seal value if rinput is available
            if not np.isnan(metrics["rinput"]):
                seal_value = f"{metrics['rinput']:.2f} MOhm"
            
            break
    
    # Create the electrode with the correct seal resistance value
    intracellular_electrode = nwbfile.create_icephys_electrode(
        name=f"IcephysElectrode",
        description=f"Patch clamp electrode for trial",
        device=device,
        slice="Drosophila brain",  # Information about the slice preparation
        seal=seal_value,  # Seal resistance from seal test 15
    )

    # Now process all seal tests to add them to the NWB file
    
    for seal_test_file in seal_test_files:
        # Extract seal test number from filename
        seal_test_num = seal_test_file.stem.split('_')[-1]
        
        # Read the seal test data
        seal_test_trial_data = read_mat(seal_test_file)["data"]
                
        # Get the current and voltage data
        seal_test_current = seal_test_trial_data["I"]
        seal_test_voltage = seal_test_trial_data["Vm"]
        seal_test_sample_rate = float(seal_test_trial_data["SAMPLERATE"])
        
        # Get detailed description based on seal test number
        if seal_test_num == "14":
            description = "Electrode resistance measurement (Rinput) performed before the experiment to assess electrode quality"
            timing = "pre-experiment"
        elif seal_test_num == "15":
            description = "Seal quality measurement (Rinput) performed before the experiment to assess the quality of the seal between the cell and the electrode"
            timing = "pre-experiment"
        elif seal_test_num == "16":
            description = "Intracellular access quality measurement performed before the experiment by comparing access resistance (Raccess) to input resistance (Rinput)"
            timing = "pre-experiment"
        elif seal_test_num == "17":
            description = "Intracellular access quality measurement performed after the experiment to assess cell health during recording"
            timing = "post-experiment"
        else:
            description = f"Seal test {seal_test_num}"
            timing = "unknown"
        
        # Calculate metrics if possible
        metrics = calculate_seal_test_metrics(seal_test_current, seal_test_voltage)
        
        # Add metrics to description if available
        metrics_text = ""
        if not np.isnan(metrics["rinput"]):
            metrics_text += f" Input resistance: {metrics['rinput']:.2f} MOhm."
        if not np.isnan(metrics["ra"]):
            metrics_text += f" Access resistance: {metrics['ra']:.2f} MOhm."
        
        # Create current clamp stimulus series for the seal test
        seal_test_stimulus = CurrentClampStimulusSeries(
            name=f"CurrentClampStimulusSeriesSealTest{seal_test_num}",
            description=f"{description} ({timing}).{metrics_text} Current clamp stimulus with 10mV peak-peak pulse at 100 Hz.",
            data=seal_test_current,
            starting_time=np.nan,  # The time is unclear
            rate=seal_test_sample_rate,
            electrode=intracellular_electrode,
        )
        
        # Add the seal test stimulus to the NWB file
        
        # Also add the seal test response
        seal_test_response = CurrentClampSeries(
            name=f"CurrentClampSeriesSealTest{seal_test_num}",
            description=f"{description} ({timing}).{metrics_text} Membrane potential recording in response to 10mV peak-peak pulse at 100 Hz.",
            data=seal_test_voltage,
            starting_time=np.nan,  # Using np.nan because this information is unavailable
            rate=seal_test_sample_rate,
            electrode=intracellular_electrode,
        )
        
        # Add the seal test response to the NWB file
        seal_test_id = int(seal_test_num)  + num_trials
        nwbfile.add_intracellular_recording(
            electrode=intracellular_electrode,
            stimulus=seal_test_stimulus,
            response=seal_test_response,
            id=seal_test_id,
        )
    
    # Process each trial
    # Calculate number of digits needed for zero-padding based on total trials
    num_digits = len(str(num_trials))

    num_trials = 5 if stub_test else num_trials  # Limit to first 5 trials if stub_test is True
    for trial_index in range(num_trials):
        # Extract trial data
        trial_data = {key: value[trial_index] for key, value in session_data.items()}
        trial = trial_data["trial"]
        # Format trial number with leading zeros based on total number of trials
        sample_rate = float(trial_data["samplerate"])
        scale_current = trial_data["scaleCurrent"]
        scale_voltage = trial_data["scaleVoltage"]

        # Apply scaling factors to convert raw values to appropriate units
        # The filteredVm is recorded from the 10xVm port on the AM-Systems 2400 amplifier
        # which is already filtered and amplified by a factor of 10
        # The scale_voltage factor (100) is then applied to obtain accurate Vm values
        voltage = trial_data["Vm"] * scale_voltage
        current = trial_data["I"] * scale_current
        filtered_voltage = trial_data["filteredVm"] * scale_voltage
        puffer_data = trial_data["puffer"]
        tachometer_data = trial_data["tachometer"]

        # Calculate trial start time (10 seconds apart)
        trial_start_time = float(trial_index * 10.0)

        # Create timestamps for the video and pose estimation
        nframes = trial_data["nframes"]  # Number of video frames
        fps = trial_data["fps"]  # Video frame rate
        video_timestamps = np.arange(0, nframes) / fps  # Timestamps in seconds
        video_timestamps += trial_start_time  # Offset by trial start time

        formatted_trial = f"{trial:0{num_digits}d}"

        # Create electrode for this trial
        # Create stimulus and response series
        stimulus = CurrentClampStimulusSeries(
            name=f"CurrentClampStimulusSeriesTrial{formatted_trial}",
            description="Current injection stimulus",
            data=current,
            starting_time=trial_start_time,
            rate=sample_rate,
            electrode=intracellular_electrode,
        )

        response = CurrentClampSeries(
            name=f"CurrentClampSeriesTrial{formatted_trial}",
            description="Membrane potential recording",
            data=voltage,
            starting_time=trial_start_time,
            rate=sample_rate,
            electrode=intracellular_electrode,
        )

        # Add the current clamp series to the NWB file
        rowindex = nwbfile.add_intracellular_recording(
            electrode=intracellular_electrode,
            stimulus=stimulus,
            response=response,
            id=trial,
        )

        # Add filtered voltage as a separate series
        filtered_response = CurrentClampSeries(
            name=f"CurrentClampFilteredResponseTrial{formatted_trial}",
            description="Filtered membrane potential recording from 10xVm port",
            data=filtered_voltage,
            starting_time=trial_start_time,
            rate=sample_rate,
            electrode=intracellular_electrode,
        )
        nwbfile.add_acquisition(filtered_response)

        # Add puffer data as a stimulus
        puffer_stimulus = TimeSeries(
            name=f"TimeSeriesPufferStimulusTrial{formatted_trial}",
            description="Air puff sensory stimulus delivery. Increases in volts indicate an air puff (sensory input).",
            data=puffer_data,
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_stimulus(puffer_stimulus)

        # Add tachometer data
        time_series_tachometer = TimeSeries(
            name=f"TimeSeriesTachometer{formatted_trial}",
            description="Light diode sensor that detects changes in wing position, used for detecting flight",
            data=tachometer_data,
            unit="V",
            starting_time=trial_start_time,
            rate=sample_rate,
        )
        nwbfile.add_acquisition(time_series_tachometer)

        # Add video data for each view
        videos_dir = data_dir / "videos"
        dlc_dir = data_dir / "DeepLabCut"

        # 1. Video_lateral_flyLeft
        angle_name = "LateralFlyLeft"
        fly_left_video_path = (
            videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_flyLeft{formatted_trial}.avi"
        )
        if fly_left_video_path.exists():
            video_interface = ExternalVideoInterface(
                file_paths=[fly_left_video_path],
                video_name=f"Video{angle_name}_trial{formatted_trial}",
            )
            video_interface.set_aligned_timestamps([video_timestamps])
            video_interface.add_to_nwbfile(nwbfile=nwbfile)

            # Add DeepLabCut data for the flyLeft view
            dlc_file_path = dlc_dir / f"{fly_left_video_path.stem}{dlc_suffix}"
            if dlc_file_path.exists():
                dlc_interface = DeepLabCutInterface(file_path=dlc_file_path)
                container_name = f"PoseEstimation{angle_name}_trial{formatted_trial}"
                dlc_interface.set_aligned_timestamps(video_timestamps)
                dlc_interface.add_to_nwbfile(nwbfile=nwbfile, container_name=container_name)

        # 2. Video_lateral_flyRight_
        angle_name = "LateralFlyRight"
        fly_right_video_path = (
            videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_flyRight_{formatted_trial}.avi"
        )
        if fly_right_video_path.exists():

            video_interface = ExternalVideoInterface(
                file_paths=[fly_right_video_path],
                video_name=f"Video{angle_name}_trial{formatted_trial}",
            )
            video_interface.set_aligned_timestamps([video_timestamps])
            video_interface.add_to_nwbfile(nwbfile=nwbfile)

        # 3. Video_lateral_ventral_
        angle_name = "LateralVentral"
        ventral_video_path = (
            videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_ventral_{formatted_trial}.avi"
        )
        if ventral_video_path.exists():

            video_interface = ExternalVideoInterface(
                file_paths=[ventral_video_path],
                video_name=f"Video{angle_name}_trial{formatted_trial}",
            )
            video_interface.set_aligned_timestamps([video_timestamps])
            video_interface.add_to_nwbfile(nwbfile=nwbfile)

    # Save the NWB file
    nwbfile_path = output_dir / f"{session_id}_experiment{experiment_number}.nwb"

    # Check if file exists and handle overwrite
    if nwbfile_path.exists() and not overwrite:
        raise FileExistsError(f"File {nwbfile_path} already exists. Set overwrite=True to overwrite.")

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
    # Set base paths
    data_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/E.phys/data")
    output_dir = Path("/home/heberto/cohen_project/Sample data/Suver Lab/nwb")

    # Set specific file paths
    session_date = "2024_10_24"
    experiment_number = "2"

    # Determine file paths
    matlab_data_file_path = data_dir / f"{session_date}_E{experiment_number}.mat"

    # Run conversion
    convert_session_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        output_dir=output_dir,
        data_dir=data_dir,
        overwrite=True,
        stub_test=False,
        verbose=True,
    )
