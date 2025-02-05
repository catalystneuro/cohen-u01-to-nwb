from pathlib import Path
from pymatreader import read_mat
from neuroconv.datainterfaces import VideoInterface
from pynwb import NWBFile, TimeSeries
from pynwb.file import Subject
import datetime
import uuid

import zoneinfo
from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile

def convert_stimuli_experiment_to_nwb(
    matlab_struct_file_path: str | Path,
    video_folder_path: str | Path,
    output_dir: str | Path,
    verbose: bool = False,
):

    matlab_struct_file_path = Path(matlab_struct_file_path)
    video_folder_path = Path(video_folder_path)
    
    assert matlab_struct_file_path.is_file(), f"Matlab struct file not found at {matlab_struct_file_path}"
    assert video_folder_path.exists(), f"Video folder not found at {video_folder_path}"

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    

    data_matlab_struct = read_mat(matlab_struct_file_path)["datastruct"]

    # Simple fields
    effector = data_matlab_struct["effector"]
    driver = data_matlab_struct["driver"]
    experiment_date_in_ddmmyy_format = data_matlab_struct["ExprDate"]
    movie_number = data_matlab_struct["MovNum"]
    opto_stimuli_range_start_stop_ms = data_matlab_struct["OptoStimRange"]
    opto_intensity_in_ampers = data_matlab_struct["optoIntensity"]  # TODO confirm units
    
    # Fields with vectors
    timestamps = data_matlab_struct["t"]
    body_pitch_angle = data_matlab_struct["bodyPitch"]
    body_roll_angle = data_matlab_struct["bodyRoll"]
    body_yaw_angle = data_matlab_struct["bodyYaw"]
    
    # Fields with matrices
    left_wing_angles = data_matlab_struct["wingL"]
    right_ring_angles = data_matlab_struct["wingR"]

    number_of_experiments = len(timestamps)
    nwbfile_path_list = []
    for experiment_index in range(number_of_experiments):
        experiment_timestamp = timestamps[experiment_index]
        experiment_body_pitch_angle = body_pitch_angle[experiment_index]
        experiment_body_roll_angle = body_roll_angle[experiment_index]
        experiment_body_yaw_angle = body_yaw_angle[experiment_index]
        experiment_left_wing_angles = left_wing_angles[experiment_index]
        experiment_right_wing_angles = right_ring_angles[experiment_index]
        
        session_date = experiment_date_in_ddmmyy_format[experiment_index]

        date_in_datetime_format = datetime.datetime.strptime(session_date, "%d%m%Y")
        # Define the timezone for Cornell
        cornell_timezone = zoneinfo.ZoneInfo("America/New_York")
        # Add time zone
        session_start_time = date_in_datetime_format.replace(tzinfo=cornell_timezone)

        session_description = f"TBD"
        experiment_video = movie_number[experiment_index]
        video_number_xxx_format = f"{experiment_video:03}"
        experiment_opto_intensity_in_ampers = opto_intensity_in_ampers[experiment_index]
        
        # TODO Check if this terminology is correct
        genotype = driver[experiment_index]
        experiment_effector = effector[experiment_index]
                
        session_id = f"{session_date}_{genotype}_{experiment_effector}_{experiment_opto_intensity_in_ampers}_{video_number_xxx_format}"
        nwbfile = NWBFile(
            session_start_time = session_start_time,
            identifier=str(uuid.uuid4()),
            session_id=session_id,
            session_description=session_description,
        )
        
        genotype = driver[experiment_index]
        nwbfile.subject = Subject(
            subject_id="generic_fly",
            genotype=genotype,
            age="TBD",
            strain="TBD",
        )

        # Add the video data
        video_file_path = video_folder_path / f"Expr_43_movie_{video_number_xxx_format}.mp4"
        assert video_file_path.is_file(), f"Video file not found at {video_file_path}"
        video_interface = VideoInterface(file_paths=[video_file_path])
        
        # Add the video interface to the NWB file
        video_interface.add_to_nwbfile(nwbfile=nwbfile)
        
        # Add the TimeSeries for the body angles
        body_pitch_angle_timeseries = TimeSeries(
            name="body_pitch_angle",
            data=experiment_body_pitch_angle,
            unit="degrees",
            timestamps=experiment_timestamp,
        )
        
        nwbfile.add_acquisition(body_pitch_angle_timeseries)
        
        body_roll_angle_timeseries = TimeSeries(
            name="body_roll_angle",
            data=experiment_body_roll_angle,
            unit="degrees",
            timestamps=experiment_timestamp,
        )
        
        nwbfile.add_acquisition(body_roll_angle_timeseries)
        
        body_yaw_angle_timeseries = TimeSeries(
            name="body_yaw_angle",
            data=experiment_body_yaw_angle,
            unit="degrees",
            timestamps=experiment_timestamp,
        )
        
        nwbfile.add_acquisition(body_yaw_angle_timeseries)
    
        
        # TODO, extract the datat form the matrices and write each component to the NWBFile
        
        
        nwbfile_path = output_dir / f"{session_id}.nwb"
        configure_and_write_nwbfile(nwbfile=nwbfile, output_filepath=nwbfile_path)        
    
        nwbfile_path_list.append(nwbfile_path)
        
    return nwbfile_path_list

if __name__ == "__main__":
    
    driver = "SS40851"
    effector = "UAS-CsChrimson"
    intensity = "0.67A"
    experiment_type = "Opto Activation"


    folder_path = Path("/home/heberto/cohen_project/Sample data/Cohen Lab/Free Flight Optogenetics Data")
    experiment_folder = folder_path / f"{driver}  {effector}"/ f"{intensity} {experiment_type}"
    assert experiment_folder.is_dir()

    matlab_struct_file_path = experiment_folder / "SS40851_0_67A_activation.mat"
    video_folder_path = experiment_folder / "mp4"
    output_dir = folder_path.parent / "nwb_files" / f"{driver}  {effector}" / f"{intensity} {experiment_type}"
    
    nwbfile_path_list = convert_stimuli_experiment_to_nwb(
        matlab_struct_file_path=matlab_struct_file_path,
        video_folder_path=video_folder_path,
        output_dir=output_dir,
        verbose=True,
    )
    
    print(nwbfile_path_list)

