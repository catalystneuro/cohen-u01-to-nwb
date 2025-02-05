import time
from pathlib import Path
from pymatreader import read_mat
from neuroconv.datainterfaces import VideoInterface
from tqdm import tqdm
from pynwb import NWBFile
from pynwb.behavior import SpatialSeries
from pynwb.file import Subject
from pynwb.behavior import CompassDirection
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

    start_time = time.time()
    if verbose:
        print(f"\nStarting conversion of {matlab_struct_file_path.name}")
        print(f"Found {len(timestamps)} experiments to convert")
    
    number_of_experiments = len(timestamps)
    nwbfile_path_list = []
    
    # Create iterator with optional progress bar
    experiment_iterator = tqdm(range(number_of_experiments), desc="Converting experiments") if verbose else range(number_of_experiments)
    
    for experiment_index in experiment_iterator:
        experiment_timestamp = timestamps[experiment_index]
        experiment_body_pitch_angle = body_pitch_angle[experiment_index]
        experiment_body_roll_angle = body_roll_angle[experiment_index]
        experiment_body_yaw_angle = body_yaw_angle[experiment_index]
        experiment_left_wing_angles = left_wing_angles[experiment_index]
        experiment_right_wing_angles = right_ring_angles[experiment_index]

        # Build session start time as a python datetime object with timezone information
        session_date = experiment_date_in_ddmmyy_format[experiment_index]
        date_in_datetime_format = datetime.datetime.strptime(session_date, "%d%m%Y")
        cornell_timezone = zoneinfo.ZoneInfo("America/New_York")
        session_start_time = date_in_datetime_format.replace(tzinfo=cornell_timezone)

        # Create detailed session description
        experiment_opto_stim_range = opto_stimuli_range_start_stop_ms[experiment_index]
        stim_start, stim_end = experiment_opto_stim_range
        stim_duration = stim_end - stim_start
        

        experiment_video = movie_number[experiment_index]
        video_number_xxx_format = f"{experiment_video:03}"
        experiment_opto_intensity_in_ampers = opto_intensity_in_ampers[experiment_index]
        
        genotype = driver[experiment_index]
        experiment_effector = effector[experiment_index]
                
        session_id = f"{session_date}_{genotype}_{experiment_effector}_{experiment_opto_intensity_in_ampers}_{video_number_xxx_format}"

        session_description = (
            f"Free flight optogenetics experiment using {experiment_effector} with {driver[experiment_index]} driver. "
            f"Optogenetic stimulation applied at {experiment_opto_intensity_in_ampers}A intensity "
            f"from {stim_start}ms to {stim_end}ms (duration: {stim_duration}ms). "
            f"High-speed video recording at 8000fps capturing wing and body kinematics "
            f"during optogenetic manipulation."
        )
        nwbfile = NWBFile(
            session_start_time = session_start_time,
            identifier=str(uuid.uuid4()),
            session_id=session_id,
            session_description=session_description,
        )
        
        # Create detailed subject information
        full_genotype = f"{driver[experiment_index]}/{experiment_effector}"
        nwbfile.subject = Subject(
            subject_id=f"fly_{session_date}_{video_number_xxx_format}",
            sex="F",  # All flies are female
            age="P3-5",  # 3-5 days post-eclosion
            description="Female fly, 3-5 days old post-eclosion",
            genotype=full_genotype,
            strain=f"{driver[experiment_index]}>UAS-{experiment_effector}",  # Standard notation for fly crosses
            species="Drosophila melanogaster"  # Common fruit fly
        )

        # Add the video data
        video_file_path = video_folder_path / f"Expr_43_movie_{video_number_xxx_format}.mp4"
        assert video_file_path.is_file(), f"Video file not found at {video_file_path}"
        video_interface = VideoInterface(file_paths=[video_file_path])
        
        # Add the video interface to the NWB file
        video_interface.add_to_nwbfile(nwbfile=nwbfile)
        
        # Create a CompassDirection container for body orientation angles
        body_orientation = CompassDirection(name="BodyOrientation")
        nwbfile.add_acquisition(body_orientation)
        
        # Store body angles as SpatialSeries
        body_angles = {
            'Pitch': experiment_body_pitch_angle,
            'Roll': experiment_body_roll_angle,
            'Yaw': experiment_body_yaw_angle
        }
        
        for angle_name, angle_data in body_angles.items():
            body_angle_series = SpatialSeries(
                name=f"SpatialSeriesBody{angle_name}",
                description=f"Body {angle_name} angle in the lab reference frame",
                data=angle_data,
                timestamps=experiment_timestamp,
                reference_frame="lab-centered spherical coordinates",
                unit="degrees"
            )
            body_orientation.add_spatial_series(body_angle_series)
    
        
        # Create a CompassDirection container for wing angles
        wing_angles = CompassDirection(name="WingAngles")
        nwbfile.add_acquisition(wing_angles)
        
        # Extract and store left wing angles
        left_wing_components = ['Stroke', 'Deviation', 'Pitch']
        for i, component in enumerate(left_wing_components):
            left_wing_series = SpatialSeries(
                name=f"SpatialSeriesLeftWing{component}",
                description=f"Left wing {component.lower()} angle",
                data=experiment_left_wing_angles[i],
                timestamps=experiment_timestamp,
                reference_frame="body-centered spherical coordinates",
                unit="degrees"
            )
            wing_angles.add_spatial_series(left_wing_series)
        
        # Extract and store right wing angles
        right_wing_components = ['Stroke', 'Deviation', 'Pitch']
        for i, component in enumerate(right_wing_components):
            right_wing_series = SpatialSeries(
                name=f"SpatialSeriesRightWing{component}",
                description=f"Right wing {component.lower()} angle",
                data=experiment_right_wing_angles[i],
                timestamps=experiment_timestamp,
                reference_frame="body-centered spherical coordinates",
                unit="degrees"
            )
            wing_angles.add_spatial_series(right_wing_series)
        
        nwbfile_path = output_dir / f"{session_id}.nwb"
        configure_and_write_nwbfile(nwbfile=nwbfile, output_filepath=nwbfile_path)        
    
        if verbose:
            print(f"\nCreated NWB file: {nwbfile_path.name}")
            print(f"Session ID: {session_id}")
            print(f"Subject: {nwbfile.subject.species}, {nwbfile.subject.genotype}")
            print(f"Experiment date: {session_date}")
            print(f"Stimulation: {stim_duration}ms at {experiment_opto_intensity_in_ampers}A")
        
        nwbfile_path_list.append(nwbfile_path)
    
    end_time = time.time()
    if verbose:
        conversion_time = end_time - start_time
        print(f"\nConversion completed in {conversion_time:.2f} seconds")
        print(f"Created {len(nwbfile_path_list)} NWB files in {output_dir}")
        
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
