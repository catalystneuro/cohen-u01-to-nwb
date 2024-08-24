"""Primary NWBConverter class for this dataset."""
from neuroconv import NWBConverter
from neuroconv.datainterfaces import VideoInterface

from .AbbyLeung_behaviorinterface import AbbyleungBehaviorInterface


class AbbyleungNWBConverter(NWBConverter):
    """Primary conversion class for my extracellular electrophysiology dataset."""

    data_interface_classes = dict(
        Behavior=AbbyleungBehaviorInterface,
        Video=VideoInterface,
    )
