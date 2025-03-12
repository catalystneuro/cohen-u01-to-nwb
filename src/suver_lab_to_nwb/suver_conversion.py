from pathlib import Path
import os
import time
from datetime import datetime
import re

import numpy as np
from neuroconv import ConverterPipe
from neuroconv.datainterfaces import DeepLabCutInterface, VideoInterface
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBHDF5IO, NWBFile
from tqdm import tqdm

from suver_lab_to_nwb.patch_clamp_interface import PatchClampInterface
from suver_lab_to_nwb.seal_test_interface import SealTestInterface


def convert_suver_lab_data(
    data_dir: str | Path,
    output_dir: str | Path,
    session_id: str = None,
    verbose: bool = True
) -> Path:
    """
    Convert Suver Lab data to NWB format.
    
    This function converts patch clamp data, videos, DeepLabCut pose estimation,
    and seal test data from the Suver Lab to the NWB format.
    
    Parameters
    ----------
    data_dir : str or Path
        Path to the directory containing the data.
    output_dir : str or Path
        Path to the directory where the NWB file will be saved.
    session_id : str, optional
        Session ID to use for the conversion. If None, it will be extracted from
        the directory name.
    verbose : bool, default: True
        Whether to print verbose output.
        
    Returns
    -------
    nwbfile_path : Path
        Path to the generated NWB file.
        
    Notes
    -----
    The function expects the following data structure:
    - A main patch clamp data file: {session_id}.mat
    - Video files in a 'videos' subdirectory or the main directory
    - DeepLabCut files in a 'DeepLabCut' subdirectory or the main directory
    - Seal test files in a 'seal tests' subdirectory or the main directory
    
    The function will create a single NWB file containing all the data.
    """
    # Start timing
    start_time = time.time()
    
    # Convert paths to Path objects
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Get session ID from directory name if not provided
    if session_id is None:
        # Extract date and experiment number from directory name
        match = re.search(r"(\d{4}_\d{2}_\d{2})_E(\d+)", data_dir.name)
        if match:
            date_str, exp_num = match.groups()
            session_id = f"{date_str}_E{exp_num}"
        else:
            session_id = data_dir.name
    
    # Find the main patch clamp data file
    patch_clamp_file = data_dir / f"{session_id}.mat"
    if not patch_clamp_file.exists():
        raise FileNotFoundError(f"Patch clamp file not found: {patch_clamp_file}")
    
    # Find video files
    video_dir = data_dir / "videos"
    if not video_dir.exists():
        video_dir = data_dir  # Videos might be in the main directory
    
    video_files = {}
    for view in ["lateral_flyLeft", "lateral_flyRight", "lateral_ventral"]:
        pattern = f"{session_id}_Video_{view}_*.avi"
        matching_files = sorted(list(video_dir.glob(pattern)))
        if matching_files:
            video_files[view] = matching_files
    
    # Find DeepLabCut files
    dlc_dir = data_dir / "DeepLabCut"
    if not dlc_dir.exists():
        dlc_dir = data_dir  # DLC files might be in the main directory
    
    pattern = f"{session_id}_Video_*DLC_*.h5"
    dlc_files = sorted(list(dlc_dir.glob(pattern)))
    
    # Find seal test files
    seal_test_dir = data_dir / "seal tests"
    if not seal_test_dir.exists():
        seal_test_dir = data_dir  # Seal test files might be in the main directory
    
    date_str = session_id.split("_E")[0]
    pattern = f"{date_str}_SealTest_*.mat"
    seal_test_files = sorted(list(seal_test_dir.glob(pattern)))

    # Initialize interfaces
    data_interfaces = {}
    
    # Add patch clamp interface
    patch_clamp_interface = PatchClampInterface(file_path=patch_clamp_file, verbose=verbose)
    data_interfaces["PatchClamp"] = patch_clamp_interface
    
    # Add video interfaces
    for view, files in video_files.items():
        video_interface = VideoInterface(file_paths=files, verbose=verbose)
        data_interfaces[f"Video_{view}"] = video_interface
    
    # Add DeepLabCut interfaces
    for i, dlc_file in enumerate(dlc_files):
        dlc_interface = DeepLabCutInterface(file_path=dlc_file, verbose=verbose)
        # Extract view from filename
        view_match = re.search(r"Video_(lateral_\w+)", dlc_file.name)
        view = view_match.group(1) if view_match else f"view{i+1}"
        data_interfaces[f"DeepLabCut_{view}"] = dlc_interface
    
    # Add seal test interfaces
    for i, seal_test_file in enumerate(seal_test_files):
        seal_test_interface = SealTestInterface(file_path=seal_test_file, verbose=verbose)
        data_interfaces[f"SealTest_{i+1}"] = seal_test_interface
    
    # Create converter
    converter = ConverterPipe(data_interfaces=data_interfaces, verbose=verbose)
    
    # Get metadata
    metadata = converter.get_metadata()
    
    # Load and merge custom metadata
    default_metadata_file = Path(__file__).parent / "metadata.yaml"
    if default_metadata_file.exists():
        custom_metadata = load_dict_from_file(default_metadata_file)
        metadata = dict_deep_update(metadata, custom_metadata)
    
    # Set up conversion options
    conversion_options = {}
    
    # Add conversion options for DeepLabCut interfaces
    for interface_name in data_interfaces:
        if interface_name.startswith("DeepLabCut_"):
            view = interface_name.replace("DeepLabCut_", "")
            conversion_options[interface_name] = {"container_name": f"PoseEstimation_{view}"}
    
    # Run conversion
    # Create output file path
    nwbfile_path = output_dir / f"{session_id}.nwb"
    converter.run_conversion(
        nwbfile_path=nwbfile_path,
        metadata=metadata,
        conversion_options=conversion_options,
        overwrite=True,
    )
    
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
    # Set paths
    data_dir = "/home/heberto/cohen_project/Sample data/Suver Lab/E.phys/data"
    output_dir = "/home/heberto/cohen_project/Sample data/Suver Lab/nwb"
    
    # Run conversion
    convert_suver_lab_data(
        data_dir=data_dir,
        output_dir=output_dir,
        session_id="2024_10_24_E2",
        overwrite=True,
        verbose=True
    )
