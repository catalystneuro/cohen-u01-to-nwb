from pathlib import Path
from datetime import datetime
from typing import Optional, Union
from neuroconv import ConverterPipe

from dickerson_lab_to_nwb.behavior_interface import BehaviorInterface
from dickerson_lab_to_nwb.thor_interface import ThorImagingInterface

from neuroconv.utils import dict_deep_update, load_dict_from_file

def convert_session(
    thor_first_tiff_file_path: Union[str, Path],
    behavior_hdf5_file_path: Union[str, Path],
    output_folder_path: Union[str, Path],
    stub_test: bool = False,
    verbose: bool = False,
) -> None:
    """
    Convert a session to NWB format.

    Parameters
    ----------
    session_path : Union[str, Path]
        Path to session folder
    nwbfile_path : Union[str, Path]
        Path to save NWB file
    stub_test : bool, default: False
        If True, only write metadata to test conversion
    """
    thor_first_tiff_file_path = Path(thor_first_tiff_file_path)
    output_folder_path = Path(output_folder_path)

    session_id = thor_first_tiff_file_path.parent.name

    # Initialize interfaces
    behavior_interface = BehaviorInterface(file_path=behavior_hdf5_file_path, verbose=True)

    thor_interface_channel_A = ThorImagingInterface(
        file_path=thor_first_tiff_file_path,
        channel_name="ChanA",
        verbose=verbose,
    )

    thor_interface_channel_B = ThorImagingInterface(
        file_path=thor_first_tiff_file_path,
        channel_name="ChanB",
        verbose=verbose,
    )

    converter = ConverterPipe(
        data_interfaces={
            "ThorChanA": thor_interface_channel_A,
            "ThorChanB": thor_interface_channel_B,
            "Behavior": behavior_interface,
        },
        verbose=verbose,
    )

    conversion_options = {
        "ThorChanA": dict(stub_test=stub_test, photon_series_index=0),
        "ThorChanB": dict(stub_test=stub_test, photon_series_index=1),
    }

    converter_metadata = converter.get_metadata()

    editable_metadata_path = Path(__file__).parent / "metadata.yaml"
    editable_metadata = load_dict_from_file(editable_metadata_path)
    metadata = dict_deep_update(converter_metadata, editable_metadata)


    nwbfile_path = output_folder_path / f"{session_id}.nwb"
    nwbfile_path.parent.mkdir(parents=True, exist_ok=True)

    converter.run_conversion(
        nwbfile_path=nwbfile_path,
        metadata=metadata,
        conversion_options=conversion_options,
        overwrite=True,
    )
    
    if verbose:
        print(f"Session {session_id} converted successfully and saved to {nwbfile_path}")
    


if __name__ == "__main__":
    # Test the conversion with a sample session
    base_folder_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/data_google_drive/")
    session_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/data_google_drive/Sample_trial/Sample_1")
    output_folder_path = Path("/home/heberto/cohen_project/Sample data/Dickerson Lab/nwb_files")

    thor_first_tiff_file_path = base_folder_path / "Sample_trial/Sample_1/sample/ChanA_001_001_001_001.tif"
    assert thor_first_tiff_file_path.is_file()

    behavior_hdf5_file_path = base_folder_path / "Sample_trial/Sample_1/Animal_1_Trial_6000/" / "Episode001.h5"
    assert behavior_hdf5_file_path.is_file()

    convert_session(
        thor_first_tiff_file_path=thor_first_tiff_file_path,
        behavior_hdf5_file_path=behavior_hdf5_file_path,
        output_folder_path=output_folder_path,
        verbose=True,
    )
