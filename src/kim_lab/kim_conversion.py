import os
from datetime import datetime

import h5py
import numpy as np
from neuroconv.tools.nwb_helpers import get_default_backend_configuration, configure_backend
from neuroconv.datainterfaces import VideoInterface
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from pynwb.file import Subject
from scipy.io import loadmat
from tqdm import tqdm
import uuid

from src.kim_lab.ophys import MultiTiffMultiPageTiffImagingInterface
from src.utils import detect_threshold_crossings

data_dir = "/Users/bendichter/Downloads/Kim Lab/input"
output_dir = "/Users/bendichter/Downloads/Kim Lab/nwb"

# /Users/bendichter/Downloads/Kim Lab/20240108b_00003/raw data/data_20240108b_00003.mat

# data(1,:) is time
# data(2,:) is left wingbeat
# data(3,:) is left-right wingbeat
# data(4,:) is x-position of the visual pattern
# data(5,:) is y-position of the visual pattern
# data(6,:) is 2-photon frame synchronization signal (1 pulse corresponds to 1 frame)
# data(7,:) is behavior camera signal (1 pulse corresponds to 1 frame)
# data(8,:) indicates the start of a stimulus (it is empty in this example)


for session_path in tqdm(os.listdir(data_dir)):
    if os.path.isdir(os.path.join(data_dir, session_path)):
        session_id = session_path
        session_dir = os.path.join(data_dir, session_path)
        session_data_dir = os.path.join(session_dir, "raw data")
        session_data_fpath = os.path.join(session_data_dir, f"data_{session_id}.mat")
        if not os.path.exists(session_data_fpath):
            continue

        # Load data
        mat_data = loadmat(session_data_fpath)
        data = mat_data["data"]
        protocol = mat_data["protocol"]

        with h5py.File(os.path.join(session_data_dir, "exp_info.mat"), "r") as f:
            age = f["age"][0,0]
            genotype = ''.join([chr(x) for x in f["cross"][:].ravel()])

        print(f"Processing {session_id=} ({age=}, {genotype=})")

        session_start_time = datetime.strptime(session_id[:8], '%Y%m%d').date()
        print(session_start_time)

        time = data[0]
        left_wingbeat = data[1]
        left_right_wingbeat = data[2]
        x_position = data[3]
        y_position = data[4]
        two_photon_frame_sync = data[5]
        behavior_camera_sync = data[6]
        stimulus_start = data[7]

        # Create NWB file
        nwbfile = NWBFile(
            session_description=f"protocol: {protocol}",
            identifier=str(uuid.uuid4()),
            session_start_time=session_start_time,
            session_id=session_id,
        )

        ophys_interface = MultiTiffMultiPageTiffImagingInterface(
            session_data_dir,
            pattern=session_id + "_{frame:05d}.tif",
            sampling_frequency=30.0,
            verbose=True
        )

        aligned_timestamps = time[detect_threshold_crossings(two_photon_frame_sync, 0.5)]
        aligned_timestamps = aligned_timestamps[:ophys_interface.imaging_extractor.get_num_frames()]
        ophys_interface.set_aligned_timestamps(aligned_timestamps=aligned_timestamps)

        ophys_interface.add_to_nwbfile(nwbfile, metadata=dict())

        video_interface = VideoInterface(
            file_paths=["/Users/bendichter/Downloads/Kim Lab/input/20240108b_00003/raw data/20240108b_00003.avi"],
        )

        video_timestamps = time[detect_threshold_crossings(behavior_camera_sync, 0.5)]
        video_timestamps = video_timestamps[:video_interface.get_num_frames()[0]]
        video_interface.set_aligned_timestamps([video_timestamps])
        video_interface.add_to_nwbfile(nwbfile, metadata=dict())

        nwbfile.subject = Subject(
            subject_id=session_id,
            genotype=genotype,
            age=f"P{age}D",
        )

        # Add data
        timeseries_wingbeat = TimeSeries(
            name="wingbeat",
            data=left_wingbeat,
            unit="n.a.",
            timestamps=time,
            description="wingbeat",
        )
        nwbfile.add_acquisition(timeseries_wingbeat)

        timeseries_left_right_wingbeat = TimeSeries(
            name="left_right_wingbeat",
            data=left_right_wingbeat,
            unit="n.a.",
            timestamps=timeseries_wingbeat,
            description="left-right wingbeat",
        )
        nwbfile.add_acquisition(timeseries_left_right_wingbeat)

        timeseries_x_position = TimeSeries(
            name="stimulus_position",
            data=np.c_[x_position, y_position],
            unit="n.a.",
            timestamps=timeseries_wingbeat,
            description="position of the visual pattern",
        )

        backend_configuration = get_default_backend_configuration(nwbfile, backend="hdf5")
        configure_backend(nwbfile=nwbfile, backend_configuration=backend_configuration)

        with NWBHDF5IO(os.path.join(output_dir, session_id + ".nwb"), mode="w") as io:
            io.write(nwbfile)
