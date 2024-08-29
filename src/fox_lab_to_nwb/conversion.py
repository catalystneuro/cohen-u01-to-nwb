from pathlib import Path
import time
from datetime import datetime 
from typing import Optional
from zoneinfo import ZoneInfo
from neuroconv.utils import dict_deep_update, load_dict_from_file

from neuroconv import ConverterPipe
from fox_lab_to_nwb.behavior import BehaviorInterface


def run_trial_conversion(trial_data_folder: Path, output_dir_path: Optional[Path] = None, verbose: bool = True):

    if verbose:
        start_time = time.time()

    if output_dir_path is None:
        output_dir_path = Path.home() / "conversion_nwb"

    # session_id = f"{subject}_{session_date}_{session_time}"
    session_id = "session"
    nwbfile_path = output_dir_path / f"{session_id}.nwb"

    file_path = trial_data_folder / "Tshx18D07_240124_115923_f3_r1.fly2"
    assert file_path.exists(), f"File {file_path} does not exist"
    interface = BehaviorInterface(file_path=file_path)

    converter = ConverterPipe(data_interfaces={"behavior": interface})


    # Add datetime to conversion
    metadata = converter.get_metadata()
    session_start_time = datetime.now().astimezone(ZoneInfo("America/New_York"))
    metadata["NWBFile"]["session_start_time"] = session_start_time

    # Update default metadata with the editable in the corresponding yaml file
    editable_metadata_path = Path(__file__).parent / "metadata.yaml"
    editable_metadata = load_dict_from_file(editable_metadata_path)
    metadata = dict_deep_update(metadata, editable_metadata)

    subject_metadata = metadata["Subject"]
    subject = "subject"
    subject_metadata["subject_id"] = f"{subject}"

    # Run conversion, this adds the basic data to the NWBFile
    conversion_options = {}
    converter.run_conversion(
        metadata=metadata,
        nwbfile_path=nwbfile_path,
        conversion_options=conversion_options,
        overwrite=True,
    )

    if verbose:
        stop_time = time.time()
        conversion_time_seconds = stop_time - start_time
        if conversion_time_seconds <= 60 * 3:
            print(f"Conversion took {conversion_time_seconds:.2f} seconds")
        elif conversion_time_seconds <= 60 * 60:
            print(f"Conversion took {conversion_time_seconds / 60:.2f} minutes")
        else:
            print(f"Conversion took {conversion_time_seconds / 60 / 60:.2f} hours")


if __name__ == "__main__":

    verbose = True
    data_path = Path("/home/heberto/cohen_project/Sample data/Fox Lab")
    assert data_path.exists(), f"Folder {data_path} does not exist"
    trial_data_folder = data_path / "Tshx18D07_240124_115923_f3_r1"
    assert trial_data_folder.exists(), f"Folder {trial_data_folder} does not exist"

    run_trial_conversion(trial_data_folder=trial_data_folder, verbose=verbose)
