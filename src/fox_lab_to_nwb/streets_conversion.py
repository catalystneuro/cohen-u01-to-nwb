from pathlib import Path
import time
from datetime import datetime
from typing import Optional
from zoneinfo import ZoneInfo

import parse
from neuroconv.utils import dict_deep_update, load_dict_from_file
from neuroconv import ConverterPipe
from neuroconv.datainterfaces import DeepLabCutInterface, VideoInterface

from fox_lab_to_nwb.behavior import BehaviorInterface


def run_trial_conversion(trial_folder_path: Path, output_dir_path: Optional[Path] = None, verbose: bool = True):

    if verbose:
        start_time = time.time()

    if output_dir_path is None:
        output_dir_path = Path.home() / "conversion_nwb"

    # Parse metadta from the trial_folder_path name
    # Define the correct parsing pattern
    pattern = r"{cross}_{date}_{time}_f{fly_number}_r{trial_repeat_number}"
    parsed_metadata = parse.parse(pattern, trial_folder_path.name)

    trial_repeat_number = parsed_metadata["trial_repeat_number"]
    cross = parsed_metadata["cross"]  # TODO, where in subject metadata should this be? Example Tshx18D07
    fly_number = parsed_metadata["fly_number"]
    session_date = parsed_metadata["date"]
    session_time = parsed_metadata["time"]

    subject = f"Fly{fly_number}"

    session_id = f"{subject}_{session_date}_{session_time}"
    nwbfile_path = output_dir_path / f"{session_id}.nwb"

    # Behavior interface
    format = "fly2"  # Authors said in email this might change
    daq_file_name = f"{cross}_{session_date}_{session_time}_f{fly_number}_r{trial_repeat_number}.{format}"
    file_path = trial_folder_path / daq_file_name
    behavior_interface = BehaviorInterface(file_path=file_path)

    # Video Interface
    # TODO: are the names of the video files always the same per trial? we need
    # More trials to find out

    side_cam_name = "SideCam_000000.avi"
    file_path = trial_folder_path / side_cam_name
    side_cam_interface = VideoInterface(file_paths=[file_path], metadata_key_name="SideCam")

    top_camera_name = "TOPCAM_000000.avi"
    file_path = trial_folder_path / top_camera_name

    top_cam_interface = VideoInterface(file_paths=[file_path], metadata_key_name="TopCam")

    haltere_camera_name = "XZ_1_186.mp4"

    file_path = trial_folder_path / haltere_camera_name
    haltere_cam_interface = VideoInterface(file_paths=[file_path], metadata_key_name="BackCam")

    # DLC interface
    top_cam_dlc_file_name = "TOPCAM_000000DLC_resnet50_antennatrackingMar11shuffle1_100000.h5"
    top_cam_file_path = trial_folder_path / top_cam_dlc_file_name
    top_cam_dlc_interface = DeepLabCutInterface(file_path=top_cam_file_path)

    haltere_cam_dlc_file_name = "XZ_1_186DLC_resnet50_haltereMar13shuffle1_100000.h5"
    halter_cam_file_path = trial_folder_path / haltere_cam_dlc_file_name
    haltere_cam_dlc_interface = DeepLabCutInterface(file_path=halter_cam_file_path)

    data_interface = {
        "Behavior": behavior_interface,
        "SideCam": side_cam_interface,
        "TopCam": top_cam_interface,
        "BackCam": haltere_cam_interface,
        "DeepLabCutTopCam": top_cam_dlc_interface,
        "DeepLabCutHalterCam": haltere_cam_dlc_interface,
    }

    converter = ConverterPipe(data_interfaces=data_interface)

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
    conversion_options = {
        "DeepLabCutTopCam": {"container_name": "PoseEstimationTopCam"},
        "DeepLabCutHalterCam": {"container_name": "PoseEstimationHaltereCam"},
    }

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
    trial_folder_path = data_path / "Tshx18D07_240124_115923_f3_r1"
    assert trial_folder_path.exists(), f"Folder {trial_folder_path} does not exist"

    run_trial_conversion(trial_folder_path=trial_folder_path, verbose=verbose)
