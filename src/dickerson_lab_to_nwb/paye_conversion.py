import os
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Tuple, Optional

import numpy as np
from PIL import Image
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from neuroconv.utils import FilePathType, DeepDict
from roiextractors.extraction_tools import PathType
from roiextractors.imagingextractor import ImagingExtractor

from ..cohen_u01_utils.utils import match_paths


def extract_experiment_details(xml_file_path: str):
    """
    Extract the frameRate from the LSM element and the start time from the Date element.

    Parameters
    ----------
    xml_file_path : str
        Path to the XML file containing the experiment details.

    Returns
    -------
    dict
        A dictionary containing the frameRate and startTime if available.
    """
    # Dictionary to hold the extracted values
    details = {}

    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Extract frameRate from the LSM element
    lsm_element = root.find(".//LSM")
    if lsm_element is not None and "frameRate" in lsm_element.attrib:
        details["frameRate"] = float(lsm_element.attrib["frameRate"])

    # Extract startTime from the Date element
    date_element = root.find(".//Date")
    if date_element is not None and "date" in date_element.attrib:
        date_str = date_element.attrib["date"]
        details["startTime"] = datetime.strptime(date_str, "%m/%d/%Y %H:%M:%S")

    return details


class ThorTiffImagingExtractor(ImagingExtractor):
    """A ImagingExtractor for multiple TIFF files."""

    extractor_name = "ThorTiffImaging"
    is_writable = False

    def __init__(self, folder_path: PathType, pattern="{channel}_001_001_001_{frame:d}.tif"):
        """Create a ThorTiffImagingExtractor instance from a TIFF file.

        Parameters
        ----------
        folder_path : str
            List of path to each TIFF file.
        """

        super().__init__()
        self.folder_path = folder_path

        paths = match_paths(folder_path, pattern)

        channels = list(set(x["channel"] for x in paths.values()))

        self._video = {}
        for channel in channels:
            data = []
            for fpath in paths:
                img = Image.open(fpath)
                data.append(np.array(img))
            self._video[channel] = np.array(data)

        shape = self._video[channels[0]].shape
        self._num_frames, self._num_rows, self._num_columns = shape
        self._num_channels = len(channels)
        self._channel_names = channels

        extracted_metadata = extract_experiment_details(os.path.join(folder_path, "Experiment.xml"))
        self._sampling_frequency = extracted_metadata.get("frameRate", None)
        self.start_time = extracted_metadata.get("startTime", None)

        self._kwargs = {"folder_path": folder_path}

    def get_frames(self, frame_idxs, channel: int = 0):
        return self._video[channel][frame_idxs, ...]

    def get_video(self, start_frame=None, end_frame=None, channel: Optional[int] = 0) -> np.ndarray:
        return self._video[channel][start_frame:end_frame, ...]

    def get_image_size(self) -> Tuple[int, int]:
        return self._num_rows, self._num_columns

    def get_num_frames(self):
        return self._num_frames

    def get_sampling_frequency(self):
        return self._sampling_frequency

    def get_num_channels(self):
        return self._num_channels

    def get_channel_names(self):
        return self._channel_names


class ThorTiffImagingInterface(BaseImagingExtractorInterface):
    """Interface for multi-page TIFF files."""

    display_name = "ThorLabs TIFF Imaging"
    Extractor = ThorTiffImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        source_schema = super().get_source_schema()
        source_schema["properties"]["folder_path"]["description"] = (
            "Directory that contains the TIFF files and " "Experiment.xml file. "
        )
        return source_schema

    def __init__(self, folder_path: FilePathType, verbose: bool = False):
        """
        Initialize reading of TIFF file.

        Parameters
        ----------
        folder_path : FilePathType
        verbose : bool, default: False
        """
        super().__init__(folder_path=folder_path, verbose=verbose)

    def get_metadata(self, photon_series_type="TwoPhotonSeries") -> DeepDict:
        metadata = super().get_metadata(photon_series_type=photon_series_type)
        metadata["NWBFile"]["session_start_time"] = self.extractor.start_time
        return metadata


folder_path = "/Users/bendichter/Downloads/Dickerson Lab/Sample_trial-20240508T162829Z-001/Sample_trial/sample/"

ThorTiffImagingInterface(folder_path=folder_path)
