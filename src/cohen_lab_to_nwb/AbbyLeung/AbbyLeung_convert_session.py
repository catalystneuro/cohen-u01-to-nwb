"""Primary script to run to convert an entire session for of data using the NWBConverter."""
from pathlib import Path
from typing import Union
import datetime

from neuroconv.utils import load_dict_from_file, dict_deep_update

from AbbyLeung.AbbyLeung_nwbconverter import AbbyleungNWBConverter


def session_to_nwb(
    video_filepath: Union[str, Path],
    behav_filepath: Union[str, Path],
    output_filepath: Union[str, Path],
    session_start_time: datetime.datetime,
    session_id,
    subject_id,
    stub_test: bool = False,
):

    video_filepath = Path(video_filepath)
    behav_filepath = Path(behav_filepath)

    nwbfile_path = output_dir_path / f"{session_id}.nwb"

    source_data = dict()
    conversion_options = dict()

    # Add Behavior
    source_data.update(dict(Behavior=dict(file_path=behav_filepath)))
    conversion_options.update(dict(Behavior=dict(stub_test=stub_test)))

    # Add Video
    source_data.update(dict(Video=dict(file_path=video_filepath)))
    conversion_options.update(dict(Video=dict(stub_test=stub_test)))

    converter = AbbyleungNWBConverter(source_data=source_data)

    # Add datetime to conversion
    metadata = converter.get_metadata()
    metadata["NWBFile"]["session_start_time"] = session_start_time

    # Update default metadata with the editable in the corresponding yaml file
    editable_metadata_path = Path(__file__).parent / "AbbyLeung_metadata.yaml"
    editable_metadata = load_dict_from_file(editable_metadata_path)
    metadata = dict_deep_update(metadata, editable_metadata)

    # Run conversion
    converter.run_conversion(metadata=metadata, nwbfile_path=nwbfile_path, conversion_options=conversion_options)


# def convert_all_sessions():
#     if stub_test:
#         output_dir_path = output_dir_path / "nwb_stub"
#     output_dir_path.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":

    # Parameters for conversion
    behav_filepath = Path("/Users/bendichter/data/CohenU01_data/Cohen Lab/SS40851_1A_Activation.mat")
    video_filepath = Path("/Users/bendichter/data/CohenU01_data/Cohen Lab/Expr_81_movie_011.mp4")
    output_dir_path = Path("~/conversion_nwb/")
    stub_test = False

    session_to_nwb(
        video_filepath=video_filepath,
        behav_filepath=behav_filepath,
        output_dir_path=output_dir_path,
        stub_test=stub_test,
        session_start_time=datetime.datetime(2021, 1, 1),
    )
