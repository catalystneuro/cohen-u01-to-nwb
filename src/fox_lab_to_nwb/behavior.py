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
            "LHutchen",  # hütchens means little hat and is the shape of the signal on the sensor
            "RHutchen",
            "PTrigger",
        ]

        # Behavior signals
        left_wing_beat_amplitude = daq_struct["data"][:, 3]
        right_wing_beat_amplitude = daq_struct["data"][:, 4]
        wing_beat_frequency = daq_struct["data"][:, 5]

        unit = "tbd"
        description = "Recorded by the WingBeat Analyzer, an LED directly overhead a photodiode that detects changes in light (i.e., changes in how much the fly’s wing is occluding the LED)"
        left_wing_beat_amplitude_time_series = TimeSeries(
            name="LeftWingBeatAmplitudeTimeSeries",
            data=left_wing_beat_amplitude,
            unit=unit,
            timestamps=timestamps,
            description=description,
        )

        description = "Recorded by the WingBeat Analyzer, an LED directly overhead a photodiode that detects changes in light (i.e., changes in how much the fly’s wing is occluding the LED)"
        right_wing_beat_amplitude_time_series = TimeSeries(
            name="RightWingBeatAmplitudeTimeSeries",
            data=right_wing_beat_amplitude,
            unit=unit,
            timestamps=timestamps,
            description=description,
        )

        description = "Recorded by the WingBeat Analyzer "
        unit = "Hz"  # a normal Drosophila will flap between 180 and 220 Hz
        wing_beat_frequency_time_series = TimeSeries(
            name="WingBeatFrequencyTimeSeries",
            data=wing_beat_frequency,
            unit=unit,
            timestamps=timestamps,
            conversion=0.1,
            description=description,
        )

        nwbfile.add_acquisition(left_wing_beat_amplitude_time_series)
        nwbfile.add_acquisition(right_wing_beat_amplitude_time_series)
        nwbfile.add_acquisition(wing_beat_frequency_time_series)

        # Hutchen is the shape of the signal on the sensor, means little hats
        lhutchen = daq_struct["data"][:, 6]
        rhutchen = daq_struct["data"][:, 7]

        unit = "arbitrary units or z-score"
        lhutchen_description = (
            "Voltage recording of the shadow of the fly’s left wingbeat. "
            "Low voltages indicate the wing is fully extended forward (minimal light occlusion), "
            "and high voltages indicate the wing is fully back (maximal light occlusion). "
            "Voltage values depend on the fly’s positioning over the photodiode."
        )

        rhutchen_description = (
            "Voltage recording of the shadow of the fly’s right wingbeat. "
            "Low voltages indicate the wing is fully extended forward (minimal light occlusion), "
            "and high voltages indicate the wing is fully back (maximal light occlusion). "
            "Voltage values depend on the fly’s positioning over the photodiode."
        )

        # Updated object names and definitions
        left_wing_voltage_series = TimeSeries(
            name="LeftWingVoltageTimeSeries",
            data=lhutchen,
            unit=unit,
            timestamps=timestamps,
            description=lhutchen_description
        )

        right_wing_voltage_series = TimeSeries(
            name="RightWingVoltageTimeSeries",
            data=rhutchen,
            unit=unit,
            timestamps=timestamps,
            description=rhutchen_description
        )
        nwbfile.add_acquisition(left_wing_voltage_series)
        nwbfile.add_acquisition(right_wing_voltage_series)

    def extract_synchronization_signals_info(self):

        mat = read_mat(self.file_path)
        recording_structure = mat["rec"]
        daq_struct = recording_structure["daq"]

        daq_sampling_rate = daq_struct["fs"]

        cam_sync = daq_struct["data"][:, 0]  # This is a defunct channel and not used according to authors
        cam_trigger = daq_struct["data"][:, 1]
        opto_trigger = daq_struct["data"][:, 2]
        ptrigger = daq_struct["data"][:, 8]

        return_dict = {
            "daq_sampling_rate": daq_sampling_rate,
            "cam_sync": cam_sync,
            "cam_trigger": cam_trigger,
            "opto_trigger": opto_trigger,
            "ptrigger": ptrigger,
        }

        return return_dict
