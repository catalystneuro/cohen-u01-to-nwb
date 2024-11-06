import os
from datetime import datetime
import yaml

import numpy as np
from neuroconv.tools.nwb_helpers import get_default_backend_configuration, configure_backend
from pynwb import TimeSeries, NWBHDF5IO
from pynwb.file import Subject, NWBFile
from pynwb.image import ImageSeries
from pynwb.ogen import OptogeneticSeries
from scipy.io import loadmat
from tqdm import tqdm


def convert_session(session_dir: str, output_dir: str, metadata: dict):

    _, session_id = os.path.split(session_dir)

    mat_fpath = os.path.join(session_dir, f"{session_id}.mat")
    data = loadmat(mat_fpath, squeeze_me=True)["data"]
    records_list = [{field: record[field] for field in record.dtype.names} for record in data]

    session_start_time = datetime.strptime(records_list[0]["date"], '%Y_%m_%d').date()

    nwbfile = NWBFile(
        session_description=metadata["NWBFile"]["session_description"] or "unknown",
        identifier=f"session_{session_id}",
        session_start_time=session_start_time,
        session_id=session_id,
    )

    ogen_device = nwbfile.create_device(**metadata["ogen_device"])
    ogen_stim_site = nwbfile.create_ogen_site(device=ogen_device, **metadata["ogen_site"])

    nwbfile.subject = Subject(
        subject_id=records_list[0]["notes"],
        genotype=records_list[0]["genotype"],
        **metadata["Subject"],
    )

    for records in records_list:
        mic_series = TimeSeries(
            name=f"micLeft{records['trial']:03d}",
            data=records["micLeft"],
            unit="n.a.",
            rate=float(records["SAMPLERATE"]),
            starting_time=np.nan,
            description="used to monitor the flapping of the fly's wing",
        )
        nwbfile.add_acquisition(mic_series)

        ogen_series = OptogeneticSeries(
            name=f"opto{records['trial']:03d}",
            site=ogen_stim_site,
            data=records["stimTiming"],
            rate=float(records["SAMPLERATE"]),
        )
        nwbfile.add_stimulus(ogen_series)

        for angle in ("frontal", "lateral"):
            img_series = ImageSeries(
                name=f"Video_{angle}_{records['trial']:03d}",
                external_file=[
                    os.path.join(session_dir, f"{session_id}_Video_{angle}_{records['trial']:02d}.avi"),
                ],
                starting_time=np.nan,
                description="unknown",
                format="external",
                rate=float(records["fps"]),
            )
            nwbfile.add_acquisition(img_series)

    backend_configuration = get_default_backend_configuration(nwbfile, backend="hdf5")
    # you can mess with the backend configuration here. I am proceeding with the defaults.
    configure_backend(nwbfile=nwbfile, backend_configuration=backend_configuration)

    with NWBHDF5IO(os.path.join(output_dir, session_id + ".nwb"), mode="w") as io:
        io.write(nwbfile)


#### USAGE ####
data_dir = "/Users/bendichter/Downloads/Suver Lab/input"

metadata_fpath = os.path.join(data_dir, "metadata.yaml")
for session_path in tqdm(os.listdir(data_dir)):
    if os.path.isdir(os.path.join(data_dir, session_path)):
        with open(metadata_fpath, "r") as f:
            metadata = yaml.safe_load(f)
        session_fpath = os.path.join(data_dir, session_path)
        output_dir = "/Users/bendichter/Downloads/Suver Lab/nwb"
        convert_session(session_fpath, output_dir, metadata)
