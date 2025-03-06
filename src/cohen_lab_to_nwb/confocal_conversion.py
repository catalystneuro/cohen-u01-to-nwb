from pathlib import Path
import datetime
import uuid
import zoneinfo

from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile
from pynwb import NWBFile
from pynwb.file import Subject

from cohen_lab_to_nwb.zeiss_confocal_interface import ZeissConfocalInterface


def convert_confocal_to_nwb(
    file_path: str | Path,
    output_dir: str | Path,
    verbose: bool = False,
):
    """
    Convert Zeiss confocal microscopy data to NWB format.

    Parameters
    ----------
    file_path : str | Path
        Path to the .czi file
    output_dir : str | Path
        Directory where the NWB file will be saved
    verbose : bool, default: False
        If True, print detailed information during conversion
    """
    file_path = Path(file_path)
    assert file_path.is_file(), f"CZI file not found at {file_path}"

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Create interfaces for each channel
    interfaces = []
    channel_names = ["Alexa Fluor 488", "Alexa Fluor 633"]
    for channel_name in channel_names:
        interface = ZeissConfocalInterface(
            file_path=file_path,
            channel_name=channel_name,
            verbose=verbose
        )
        interfaces.append(interface)

    # Get metadata from the first interface
    metadata = interfaces[0].get_metadata()
    
    # Set timezone to Eastern Time (Cornell)
    cornell_timezone = zoneinfo.ZoneInfo("America/New_York")
    session_start_time = metadata["NWBFile"]["session_start_time"].replace(tzinfo=cornell_timezone)

    # Create a detailed session description
    session_description = (
        f"Confocal microscopy imaging using {', '.join(channel_names)}. "
        f"Data acquired as Z-stack with {interfaces[0].imaging_extractor._num_z} slices. "
        f"Image dimensions: {interfaces[0].imaging_extractor._num_rows}x{interfaces[0].imaging_extractor._num_columns} pixels."
    )

    # Create NWB file
    nwbfile = NWBFile(
        session_start_time=session_start_time,
        identifier=str(uuid.uuid4()),
        session_id=file_path.stem,
        session_description=session_description,
        experimenter=["Cohen Lab"],
        keywords=["Confocal microscopy", "Fluorescence imaging", "Z-stack"],
    )

    # Add subject information
    nwbfile.subject = Subject(
        subject_id=file_path.stem,
        species="Drosophila melanogaster",
        description="Adult Drosophila melanogaster",
    )

    # Add each channel's data to the NWB file
    for interface in interfaces:
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)

    # Save the NWB file
    nwbfile_path = output_dir / f"{file_path.stem}.nwb"
    configure_and_write_nwbfile(nwbfile=nwbfile, output_filepath=nwbfile_path)

    if verbose:
        print(f"\nCreated NWB file: {nwbfile_path.name}")
        print(f"Session ID: {nwbfile.session_id}")
        print(f"Subject: {nwbfile.subject.species}")
        print(f"Channels: {', '.join(channel_names)}")
        print(f"Z-stack slices: {interfaces[0].imaging_extractor._num_z}")
        print(f"Image dimensions: {interfaces[0].imaging_extractor._num_rows}x{interfaces[0].imaging_extractor._num_columns}")

    return nwbfile_path


if __name__ == "__main__":
    # Example usage
    folder_path = Path("/home/heberto/cohen_project/Sample data/Cohen Lab/Confocal images")
    assert folder_path.is_dir()

    file_path = folder_path / "tp1-alt_VNC2.czi"
    assert file_path.is_file()

    output_dir = folder_path.parent / "nwb_files"

    nwbfile_path = convert_confocal_to_nwb(
        file_path=file_path,
        output_dir=output_dir,
        verbose=True,
    )
