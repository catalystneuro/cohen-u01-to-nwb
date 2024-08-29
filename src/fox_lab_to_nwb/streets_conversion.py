from pathlib import Path
import time
from datetime import datetime 
from typing import Optional
from zoneinfo import ZoneInfo

import parse 
from neuroconv.utils import dict_deep_update, load_dict_from_file
from neuroconv import ConverterPipe

from fox_lab_to_nwb.behavior import BehaviorInterface


def run_trial_conversion(trial_data_folder: Path, output_dir_path: Optional[Path] = None, verbose: bool = True):

    if verbose:
        start_time = time.time()

    if output_dir_path is None:
        output_dir_path = Path.home() / "conversion_nwb"



    # Parse metadta from the trial_data_folder name
    # Define the correct parsing pattern
    pattern = r"{cross}_{date}_{time}_f{fly_number}_r{trial_repeat_number}"
    parsed_metadata = parse.parse(pattern, trial_data_folder.name)
    
    trial_repeat_number = parsed_metadata["trial_repeat_number"]
    cross = parsed_metadata["cross"] # TODO, where in subject metadata should this be? Example Tshx18D07
    fly_number = parsed_metadata["fly_number"]
    session_date = parsed_metadata["date"]
    session_time = parsed_metadata["time"]
    
    subject = f"Fly{fly_number}"    

    session_id = f"{subject}_{session_date}_{session_time}"
    nwbfile_path = output_dir_path / f"{session_id}.nwb"

    # Behavior interface
    format = "fly2"  # Authors said in email this might change
    daq_file_name = f"{cross}_{session_date}_{session_time}_f{fly_number}_r{trial_repeat_number}.{format}"
    file_path = trial_data_folder / daq_file_name
    interface = BehaviorInterface(file_path=file_path)

    converter = ConverterPipe(data_interfaces={"behavior": interface})


    # Add datetime to conversion
    metadata = converter.get_metadata()
    session_datetime = datetime.strptime(f"{session_date}_{session_time}", "%y%m%d_%H%M%S")
    session_start_time = session_datetime.astimezone(ZoneInfo("America/New_York"))
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
