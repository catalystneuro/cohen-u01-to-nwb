from pathlib import Path
import os
import time
from datetime import datetime
import re
import uuid
import zoneinfo
from typing import Dict, Any, Optional, Union, List

import numpy as np
from neuroconv.datainterfaces import DeepLabCutInterface, VideoInterface
from pymatreader import read_mat
from pynwb import NWBFile, NWBHDF5IO
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries
from pynwb.base import TimeSeries


def convert_session_to_nwb(
    matlab_data_file_path: Union[str, Path],
    video_file_path: Union[str, Path],
    deeplab_cut_file_path: Union[str, Path],
    output_dir: Union[str, Path],
    trial_index: int = 0,
    overwrite: bool = False,
    verbose: bool = True
) -> Path:
    """
    Convert Suver Lab data to NWB format.
    
    This function converts patch clamp data, videos, and DeepLabCut pose estimation
    from the Suver Lab to the NWB format.
    
    Parameters
    ----------
    matlab_data_file_path : str or Path
        Path to the MATLAB file containing patch clamp data.
    video_file_path : str or Path
        Path to the video file.
    deeplab_cut_file_path : str or Path
        Path to the DeepLabCut file.
    output_dir : str or Path
        Path to the directory where the NWB file will be saved.
    trial_index : int, default: 0
        Index of the trial to convert.
    overwrite : bool, default: False
        Whether to overwrite existing NWB files.
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
    video_file_path = Path(video_file_path)
    deeplab_cut_file_path = Path(deeplab_cut_file_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read the MATLAB data
    session_data = read_mat(matlab_data_file_path)["data"]
    # Extract trial data
    trial_data_index = {key: value[trial_index] for key, value in session_data.items()}
    session_date = trial_data_index["date"]
    trial = trial_data_index["trial"]
    experiment_number = trial_data_index["expNumber"]
    condition = trial_data_index["condition"]
    age = trial_data_index["age"]
    genotype = trial_data_index["genotype"]
    sample_rate = float(trial_data_index["samplerate"])
    fps = trial_data_index["fps"]  # Video frame rate
    scale_current = trial_data_index["scaleCurrent"]
    scale_voltage = trial_data_index["scaleVoltage"]
    nframes = trial_data_index["nframes"]  # Number of video frames
    voltage = trial_data_index["Vm"]
    current = trial_data_index["I"]
    filtered_voltage = trial_data_index["filteredVm"]
    puffer_data = trial_data_index["puffer"]
    tachometer_data = trial_data_index["tachometer"]
    # Transform session date (format: 2024_10_24) to datetime object
    session_start_time = datetime.strptime(session_date, "%Y_%m_%d")
    vanderbilt_time_zone = zoneinfo.ZoneInfo("America/Chicago")
    session_start_time = session_start_time.replace(tzinfo=vanderbilt_time_zone)
    
    # Create session ID and description
    session_id = f"{condition}_{trial}_{age}_{genotype}"
    session_description = f"Patch clamp recording from {session_date}, trial {trial}"
    
    # Create NWB file
    identifier = str(uuid.uuid4())
    nwbfile = NWBFile(
        session_start_time=session_start_time,
        session_description=session_description,
        identifier=identifier,
        session_id=session_id
    )

    # Create device and electrode
    device = nwbfile.create_device(name="Patch clamp amplifier")
    
    intracellular_electrode = nwbfile.create_icephys_electrode(
        name=f"Electrode_trial{trial}",
        description=f"Patch clamp electrode for trial {trial}",
        device=device,
        slice="Drosophila brain",  # Information about the slice preparation
        seal="TBD",  # Seal resistance
    )
    
    # Create stimulus and response series
    stimulus = CurrentClampStimulusSeries(
        name=f"stimulus_trial{trial}",
        description="Current injection stimulus",
        data=current,
        starting_time=0.0,
        rate=sample_rate,
        electrode=intracellular_electrode,
    )
    
    response = CurrentClampSeries(
        name=f"response_trial{trial}",
        description="Membrane potential recording",
        data=voltage,
        starting_time=0.0,
        rate=sample_rate,
        electrode=intracellular_electrode,
    )
    
    # Add the current clamp series to the NWB file
    rowindex = nwbfile.add_intracellular_recording(
        electrode=intracellular_electrode, 
        stimulus=stimulus, 
        response=response, 
        id=trial
    )
    
    # Add filtered voltage as a separate series
    filtered_response = CurrentClampSeries(
        name=f"filtered_response_trial{trial}",
        description="Filtered membrane potential recording",
        data=filtered_voltage,
        starting_time=0.0,
        rate=sample_rate,
        electrode=intracellular_electrode,
    )
    nwbfile.add_acquisition(filtered_response)
    
    # Add puffer and tachometer data
    time_series_puffer = TimeSeries(
        name=f"puffer_trial{trial}",
        description="Puffer stimulus delivery",
        data=puffer_data,
        unit="V",
        starting_time=0.0,
        rate=sample_rate,
    )
    
    time_series_tachometer = TimeSeries(
        name=f"tachometer_trial{trial}",
        description="Wing beat frequency measurement",
        data=tachometer_data,
        unit="V",
        starting_time=0.0,
        rate=sample_rate,
    )
    
    nwbfile.add_acquisition(time_series_puffer)
    nwbfile.add_acquisition(time_series_tachometer)
    
    # Add video data
    video_interface = VideoInterface(
        file_path=[video_file_path],
    )
    video_interface.add_to_nwbfile(nwbfile=nwbfile)
    
    # Add DeepLabCut data
    dlc_interface = DeepLabCutInterface(file_path=deeplab_cut_file_path)
    dlc_interface.add_to_nwbfile(nwbfile=nwbfile)
    
    # Save the NWB file
    nwbfile_path = output_dir / f"{session_id}_trial{trial_index}.nwb"
    
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
    trial_index = 0  # Use the first trial by default
    
    # Determine file paths
    matlab_data_file_path = data_dir / f"{session_date}_E{experiment_number}.mat"
    
    # Read the MATLAB data to get trial information
    session_data = read_mat(matlab_data_file_path)["data"]
    trial = session_data["trial"][trial_index]
    
    # Determine video and DeepLabCut file paths
    view = "Video_lateral_flyLeft"
    video_file_path = data_dir / "videos" / f"{session_date}_E{experiment_number}_{view}{trial}.avi"
    
    deeplab_cut_prefix = "DLC_resnet50_2x_JONsEffect_noGlueOct29shuffle1_1030000.h5"
    deeplab_cut_file_path = data_dir / "DeepLabCut" / f"{video_file_path.stem}{deeplab_cut_prefix}"
    
    # Run conversion
    convert_session_to_nwb(
        matlab_data_file_path=matlab_data_file_path,
        video_file_path=video_file_path,
        deeplab_cut_file_path=deeplab_cut_file_path,
        output_dir=output_dir,
        trial_index=trial_index,
        overwrite=True,
        verbose=True
    )
