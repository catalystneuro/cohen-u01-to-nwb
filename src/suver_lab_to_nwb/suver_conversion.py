from pathlib import Path
import os
import time
from datetime import datetime
import re
import uuid
import zoneinfo
from typing import Dict, Any, Optional, Union, List

import numpy as np
from neuroconv.datainterfaces import DeepLabCutInterface, ExternalVideoInterface
from pymatreader import read_mat
from pynwb import NWBFile, NWBHDF5IO
from pynwb.file import Subject
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries
from pynwb.base import TimeSeries


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
        genotype=genotype
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
        subject=subject
    )
    
    # Create device with detailed information from the email
    device = nwbfile.create_device(
        name="AMMC Patch Clamp Amplifier 2400",
        description="AM-Systems 2400 patch clamp amplifier with 10xVm port for filtered and amplified recordings"
    )
    
    # Define the DeepLabCut suffix based on the provided file patterns
    dlc_suffix = "DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5"

    intracellular_electrode = nwbfile.create_icephys_electrode(
        name=f"IcephysElectrode",
        description=f"Patch clamp electrode for trial",
        device=device,
        slice="Drosophila brain",  # Information about the slice preparation
        seal="TBD",  # Seal resistance
    )

    # Process each trial
    num_trials = len(session_data["trial"])
    # Calculate number of digits needed for zero-padding based on total trials
    num_digits = len(str(num_trials))
    
    
    num_trials = 5 if stub_test else num_trials  # Limit to first 5 trials if stub_test is True
    for trial_index in range(num_trials):
        # Extract trial data
        trial_data = {key: value[trial_index] for key, value in session_data.items()}
        trial = trial_data["trial"]
        # Format trial number with leading zeros based on total number of trials
        sample_rate = float(trial_data["samplerate"])
        fps = trial_data["fps"]  # Video frame rate
        scale_current = trial_data["scaleCurrent"]
        scale_voltage = trial_data["scaleVoltage"]
        nframes = trial_data["nframes"]  # Number of video frames
        
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
            name=f"CurrentClampResponseTrial{formatted_trial}",
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
        
        # Explicitly handle each video view
        
        # 1. Video_lateral_flyLeft
        angle_name = "LateralFlyLeft"
        fly_left_video_path = videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_flyLeft{formatted_trial}.avi"
        if fly_left_video_path.exists():
            video_interface = ExternalVideoInterface(
                file_paths=[fly_left_video_path],
            )
            video_interface.add_to_nwbfile(nwbfile=nwbfile)
            
            # Add DeepLabCut data for the flyLeft view
            dlc_file_path = dlc_dir / f"{fly_left_video_path.stem}{dlc_suffix}"
            if dlc_file_path.exists():                
                dlc_interface = DeepLabCutInterface(file_path=dlc_file_path)
                container_name = f"PoseEstimation{angle_name}_trial{formatted_trial}"
                dlc_interface.add_to_nwbfile(nwbfile=nwbfile, container_name=container_name)
        
        # 2. Video_lateral_flyRight_
        fly_right_video_path = videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_flyRight_{formatted_trial}.avi"
        if fly_right_video_path.exists():

            
            video_interface = ExternalVideoInterface(
                file_paths=[fly_right_video_path],
            )
            video_interface.add_to_nwbfile(nwbfile=nwbfile)
        
        # 3. Video_lateral_ventral_
        ventral_video_path = videos_dir / f"{session_date}_E{experiment_number}_Video_lateral_ventral_{formatted_trial}.avi"
        if ventral_video_path.exists():

            video_interface = ExternalVideoInterface(
                file_paths=[ventral_video_path],
            )
            video_interface.add_to_nwbfile(nwbfile=nwbfile)
    
    # Save the NWB file
    nwbfile_path = output_dir / f"{session_id}_experiment{experiment_number}.nwb"
    
    # Check if file exists and handle overwrite
    if nwbfile_path.exists() and not overwrite:
        raise FileExistsError(f"File {nwbfile_path} already exists. Set overwrite=True to overwrite.")
    
    with NWBHDF5IO(nwbfile_path, mode="w") as io:
        io.write(nwbfile)
    
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
        stub_test=True,
        verbose=True
    )
