import os
import uuid
from datetime import datetime

import numpy as np
from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile
from pynwb import NWBFile, TimeSeries
from pynwb.image import ImageSeries
from pynwb.ogen import OptogeneticSeries
from scipy.io import loadmat


mat_file = "/Users/bendichter/data/CohenU01_data/Cohen Lab/SS40851_1A_Activation.mat"
video_dir = "/Users/bendichter/data/CohenU01_data/Cohen Lab/videos"

def convert_session(mat_file, video_dir, metadata):
    nwbfile = NWBFile(
        identifier=str(uuid.uuid4()),
        **metadata["NWBFile"]
    )

    # add videos
    for video_file in os.listdir(video_dir):

        image_series = ImageSeries(
            name=os.path.splitext(video_file)[0],
            external_file=[video_file],
            starting_frame=[0],
            format="external",
            starting_time=np.nan,
            rate=8000,
        )

        nwbfile.add_acquisition(image_series)

    # add behavior
    var_name = os.path.splitext(os.path.split(mat_file)[1])[0]
    mat_out = loadmat(mat_file, simplify_cells=True)
    for i, x in enumerate(mat_out[var_name]):
        for ts_name, info in metadata["Behavior"].items():
            nwbfile.add_acquisition(
                TimeSeries(
                    name=f"{ts_name}{i:03d}",
                    description=info["description"],
                    data=x[ts_name],
                    unit=info["unit"],
                    rate=8000,
                    starting_time=np.nan,
                )
            )

    #optogenetic stimulation
    device = nwbfile.create_device(**metadata["Ogen"]["Device"])

    ogen_site = nwbfile.create_ogen_site(
        device=device,
        **metadata["Ogen"]["OgenSite"],
    )

    ogen_series = OptogeneticSeries(
        name="OptogeneticSeries",
        data=,  # watts
        site=ogen_site,
        rate=30.0,  # Hz
    )

    configure_and_write_nwbfile(nwbfile=nwbfile, backend="hdf5", output_filepath="test2.nwb")


if __name__ == "__main__":
    convert_session(mat_file, video_dir, session_start_time=datetime.now())

