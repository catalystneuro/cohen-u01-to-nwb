from pathlib import Path
from datetime import datetime
from typing import Optional, Union
from neuroconv.utils import load_dict_from_file, dict_deep_update
from pynwb import NWBHDF5IO, NWBFile

from .behavior_interface import BehaviorInterface
from .thor_interface import ThorTiffImagingInterface

def convert_session(
    session_path: Union[str, Path],
    output_folder_path: Union[str, Path],
    metadata: Optional[dict] = None,
    stub_test: bool = False,
) -> None:
    """
    Convert a session to NWB format.

    Parameters
    ----------
    session_path : Union[str, Path]
        Path to session folder
    nwbfile_path : Union[str, Path]
        Path to save NWB file
    metadata : dict, optional
        Metadata dictionary
    stub_test : bool, default: False
        If True, only write metadata to test conversion
    """
    session_path = Path(session_path)
    output_folder_path = Path(output_folder_path)
    
    # Extract session info from folder name
    # Example: Tshx18D07_240124_115923_f3_r1
    session_id = session_path.name
    
    # Initialize interfaces
    behavior_interface = BehaviorInterface(
        file_path=session_path / f"{session_id}.h5",  # HDF5 file with behavioral data
        verbose=True
    )
    
    # Find first OME TIFF file
    first_tiff = session_path / "ChanA_001_001_001_001.tif"
    if not first_tiff.exists():
        raise ValueError(f"Could not find first OME TIFF file: {first_tiff}")
        
    thor_interface = ThorTiffImagingInterface(
        file_path=first_tiff,
        verbose=True
    )
    
    # Get metadata from interfaces
    metadata = {}
    metadata.update(behavior_interface.get_metadata())
    metadata.update(thor_interface.get_metadata())
    
    # Create NWB file
    nwbfile = NWBFile(
        session_description=f"Payel session {session_id}",
        identifier=session_id,
        session_start_time=metadata.get("NWBFile", {}).get("session_start_time", datetime.now()),
        experimenter="Payel",
        lab="Dickerson Lab",
        institution="Case Western Reserve University",
        experiment_description="Optogenetic stimulation during behavior",
    )
    
    # Add data from interfaces
    behavior_interface.add_to_nwbfile(nwbfile, metadata=metadata)
    thor_interface.add_to_nwbfile(nwbfile, metadata=metadata)
    
    # Write NWB file
    nwbfile_path = output_folder_path / f"{session_id}.nwb"
    nwbfile_path.parent.mkdir(parents=True, exist_ok=True)
    with NWBHDF5IO(nwbfile_path, mode="w") as io:
        io.write(nwbfile)

if __name__ == "__main__":
    # Test the conversion with a sample session
    base_folder_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/data_google_drive/")
    session_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/data_google_drive/Sample_trial/Sample_1")
    output_folder_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/nwb_files")
    
    h5_file_behavior_path = base_folder_path / "Sample_trial/Sample_1/Animal_1_Trial_6000/" / "Episode001.h5.h5"
    assert h5_file_behavior_path.is_file()
    
    convert_session(
        session_path=session_path,
        output_folder_path=output_folder_path,
        stub_test=True  # Set to False for full conversion
    )
