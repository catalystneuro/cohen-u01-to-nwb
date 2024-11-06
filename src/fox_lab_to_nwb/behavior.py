from pathlib import Path
from typing import Optional

from neuroconv import BaseDataInterface
from pydantic import FilePath
from pynwb import NWBFile
from pynwb import TimeSeries
from pymatreader import read_mat


class BehaviorInterface(BaseDataInterface):
    def __init__(self, file_path: FilePath):
        self.file_path = Path(file_path)

    def add_to_nwbfile(
        self,
        nwbfile: NWBFile,
        metadata: Optional[dict] = None,
    ) -> None:

        mat = read_mat(self.file_path)
        recording_structure = mat["rec"]
        daq_struct = recording_structure["daq"]

        timestamps = daq_struct["tstamps"]

        # Extracted from the matlab manually
        # neither scipy or pymatreader can read string characters from the .mat file
        channel_names = [
            "CamSync",
            "CamTrigger",
            "OptoTrigger",
            "LWingBeatAmp",
            "RWingBeatAmp",
            "WingBeatFreq",
            "LHutchen",
            "RHutchen",
            "PTrigger",
        ]

        # TODO: figure how how to store this
        # Synchronization signals
        cam_sync = daq_struct["data"][:, 0]
        cam_trigger = daq_struct["data"][:, 1]
        opto_trigger = daq_struct["data"][:, 2]
        ptrigger = daq_struct["data"][:, 8]

        # Behavior signals
        left_wing_beat_amplitude = daq_struct["data"][:, 3]
        right_wing_beat_amplitude = daq_struct["data"][:, 4]
        wing_beat_frequency = daq_struct["data"][:, 5]

        unit = "tbd"
        description = "tbd"
        left_wing_beat_amplitude_time_series = TimeSeries(
            name="LeftWingBeatAmplitudeTimeSeries",
            data=left_wing_beat_amplitude,
            unit=unit,
            timestamps=timestamps,
            description=description,
        )

        description = "tbd"
        right_wing_beat_amplitude_time_series = TimeSeries(
            name="RightWingBeatAmplitudeTimeSeries",
            data=right_wing_beat_amplitude,
            unit=unit,
            timestamps=timestamps,
            description=description,
        )

        description = "tbd"
        unit = "Hz"  # TODO: Figure this out, the values in the plot are around 3 but should be higher for flies
        wing_beat_frequency_time_series = TimeSeries(
            name="WingBeatFrequencyTimeSeries",
            data=wing_beat_frequency,
            unit=unit,
            timestamps=timestamps,
            description=description,
        )

        nwbfile.add_acquisition(left_wing_beat_amplitude_time_series)
        nwbfile.add_acquisition(right_wing_beat_amplitude_time_series)
        nwbfile.add_acquisition(wing_beat_frequency_time_series)

        # TODO: Ask Ben if this nesting makes sense? probably not
        # time_series = [left_wing_beat_amplitude_time_series, right_wing_beat_amplitude_time_series, wing_beat_frequency_time_series]
        # behavioral_time_series_container = BehavioralTimeSeries(name="BehavioralTimeSeries", time_series=time_series)
        # nwbfile.add_acquisition(behavioral_time_series_container)

        # Not clear what are those signals, haltere?
        lhutchen = daq_struct["data"][:, 6]
        rhutchen = daq_struct["data"][:, 7]

        unit = "tbd"
        description = "tbd"
        lhutchen_time_series = TimeSeries(
            name="LHutchenTimeSeries", data=lhutchen, unit=unit, timestamps=timestamps, description=description
        )

        description = "tbd"
        rhutchen_time_series = TimeSeries(
            name="RHutchenTimeSeries", data=rhutchen, unit=unit, timestamps=timestamps, description=description
        )

        nwbfile.add_acquisition(lhutchen_time_series)
        nwbfile.add_acquisition(rhutchen_time_series)
