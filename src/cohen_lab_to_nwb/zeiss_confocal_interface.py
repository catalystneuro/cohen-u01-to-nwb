from collections import defaultdict, namedtuple
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List, Union

import numpy as np
import czifile
import warnings
import xml.etree.ElementTree as ET

from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType


class ZeissConfocalImagingExtractor(ImagingExtractor):
    """
    An ImagingExtractor for Zeiss .czi confocal microscopy files.
    
    This extractor handles multi-dimensional image data from Zeiss confocal microscopes,
    supporting multiple channels, Z-stacks, and time series data.
    
    The CZI format uses the following dimension order:
    B (Block) - Used for tiled/mosaic images
    V (View) - For multi-view imaging
    C (Channel) - Number of channels
    T (Time) - Time points
    Z (Z-stack) - Focal planes
    Y (Height) - Y dimension
    X (Width) - X dimension
    0 (Sample) - For multi-sample imaging
    """

    extractor_name = "ZeissConfocalImaging"
    is_writable = False

    def __init__(self, file_path: Union[str, Path], channel_name: Optional[str] = None):
        """
        Create a ZeissConfocalImagingExtractor instance from a CZI file.

        Parameters
        ----------
        file_path : Union[str, Path]
            Path to the .czi file
        channel_name : Optional[str]
            Name of the channel to extract (e.g., "Ch1" or "Ch2")
        """
        super().__init__()
        self.file_path = Path(file_path)
        self.channel_name = channel_name
        self._data = None

        # Open and read the CZI file
        with czifile.CziFile(self.file_path) as czi_reader:
            self._metadata = czi_reader.metadata()
            metadata_root = ET.fromstring(self._metadata)
            
            # Extract image dimensions from metadata
            image_info = metadata_root.find(".//Image/Information/Image")
            self._num_channels = int(image_info.find("SizeC").text) if image_info.find("SizeC") is not None else 1
            self._num_frames = int(image_info.find("SizeT").text) if image_info.find("SizeT") is not None else 1
            self._num_z = int(image_info.find("SizeZ").text) if image_info.find("SizeZ") is not None else 1
            self._num_rows = int(image_info.find("SizeY").text)
            self._num_columns = int(image_info.find("SizeX").text)
            
            # Get channel information
            channels = metadata_root.findall(".//Dimensions/Channels/Channel")
            self._channel_names = []
            for channel in channels:
                fluor = channel.find("Fluor")
                if fluor is not None and fluor.text:
                    self._channel_names.append(fluor.text)
                else:
                    self._channel_names.append(f"Ch{len(self._channel_names) + 1}")

            # Filter by channel if specified
            if channel_name is not None:
                if channel_name not in self._channel_names:
                    raise ValueError(f"Channel '{channel_name}' not found. Available channels: {self._channel_names}")
                self._channel_index = self._channel_names.index(channel_name)
                self._channel_names = [channel_name]
                self._num_channels = 1
            else:
                self._channel_index = None

            # Get sampling frequency if available
            laser_scan_info = metadata_root.find(".//LaserScanInfo")
            if laser_scan_info is not None:
                frame_time = float(laser_scan_info.find("FrameTime").text)
                self._sampling_frequency = 1.0 / frame_time if frame_time > 0 else None
            else:
                self._sampling_frequency = None

            # Load the complete data array
            self._data = czi_reader.asarray()
            # Remove singleton dimensions (B, V, 0)
            self._data = np.squeeze(self._data)

            # Store the data type
            self._dtype = self._data.dtype

        self._kwargs = {"file_path": str(file_path)}

    def get_frames(self, frame_idxs: List[int]) -> np.ndarray:
        """
        Get specific frames by their indices.

        Parameters
        ----------
        frame_idxs : List[int]
            List of frame indices to retrieve

        Returns
        -------
        np.ndarray
            Array of shape (n_frames, height, width) if no depth,
            or (n_frames, height, width, n_z) if Z dimension exists
        """
        if self._channel_index is not None:
            if self._num_z > 1:
                return self._data[frame_idxs, self._channel_index, :, :, :]
            else:
                return self._data[frame_idxs, self._channel_index, :, :]
        else:
            if self._num_z > 1:
                return self._data[frame_idxs, :, :, :, :]
            else:
                return self._data[frame_idxs, :, :, :]

    def get_video(self, start_frame: Optional[int] = None, end_frame: Optional[int] = None) -> np.ndarray:
        """
        Get a range of frames.

        Parameters
        ----------
        start_frame : Optional[int]
            Start frame index
        end_frame : Optional[int]
            End frame index

        Returns
        -------
        np.ndarray
            Video data array
        """
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = self._num_frames
        frame_idxs = list(range(start_frame, end_frame))
        return self.get_frames(frame_idxs)

    def get_image_size(self) -> Tuple[int, int]:
        """Return the image dimensions (height, width)."""
        return self._num_rows, self._num_columns

    def get_num_frames(self) -> int:
        """Return the number of frames."""
        return self._num_frames

    def get_sampling_frequency(self) -> Optional[float]:
        """Return the sampling frequency if available."""
        return self._sampling_frequency

    def get_num_channels(self) -> int:
        """Return the number of channels."""
        return self._num_channels

    def get_channel_names(self) -> List[str]:
        """Return the channel names."""
        return self._channel_names

    def get_dtype(self):
        """Return the data type."""
        return self._dtype


class ZeissConfocalInterface(BaseImagingExtractorInterface):
    """
    Interface for Zeiss confocal microscopy CZI files.
    """
    display_name = "Zeiss Confocal Imaging"
    Extractor = ZeissConfocalImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        """Get the source schema for the interface."""
        source_schema = super().get_source_schema()
        source_schema["properties"]["file_path"] = {
            "type": "string",
            "description": "Path to the Zeiss .czi file"
        }
        source_schema["properties"]["channel_name"] = {
            "type": "string",
            "description": "Name of the channel to extract (e.g., 'Alexa Fluor 488')",
            "required": False
        }
        return source_schema

    def __init__(self, file_path: FilePathType, channel_name: Optional[str] = None, verbose: bool = False):
        """Initialize the interface."""
        super().__init__(file_path=file_path, channel_name=channel_name, verbose=verbose)
        self.channel_name = channel_name

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the CZI file.

        Returns
        -------
        DeepDict
            Metadata dictionary containing imaging device and optical channel information
        """
        metadata = super().get_metadata()

        # Parse CZI metadata
        with czifile.CziFile(self.source_data["file_path"]) as czi_reader:
            root = ET.fromstring(czi_reader.metadata())

            # Get creation date
            document = root.find(".//Document")
            if document is not None:
                creation_date = document.find("CreationDate")
                if creation_date is not None:
                    self.session_start_time = datetime.strptime(creation_date.text, "%Y-%m-%dT%H:%M:%S")
                    metadata["NWBFile"]["session_start_time"] = self.session_start_time

            # Get microscope information
            microscope = root.find(".//Microscopes/Microscope")
            if microscope is not None:
                system = microscope.find("System")
                device_description = f"Zeiss {system.text}" if system is not None else "Zeiss Confocal Microscope"
            else:
                device_description = "Zeiss Confocal Microscope"

            metadata.setdefault("Ophys", {})["Device"] = [{
                "name": "ZeissMicroscope",
                "description": device_description
            }]

            # Get objective information
            objective = root.find(".//Objectives/Objective")
            if objective is not None:
                lens_na = float(objective.find("LensNA").text) if objective.find("LensNA") is not None else None
                magnification = float(objective.find("NominalMagnification").text) if objective.find("NominalMagnification") is not None else None
                immersion = objective.find("Immersion").text if objective.find("Immersion") is not None else None

            # Get scaling information
            scaling = root.find(".//Scaling/Items/Distance")
            if scaling is not None:
                pixel_size_x = float(scaling.find("./[@Id='X']/Value").text) if scaling.find("./[@Id='X']/Value") is not None else None
                pixel_size_y = float(scaling.find("./[@Id='Y']/Value").text) if scaling.find("./[@Id='Y']/Value") is not None else None

            # Process channel information
            channels = root.findall(".//Dimensions/Channels/Channel")
            optical_channels = []
            
            for channel in channels:
                channel_name = channel.find("Name").text if channel.find("Name") is not None else None
                if channel_name is None or (self.channel_name is not None and channel_name != self.channel_name):
                    continue

                fluor = channel.find("Fluor").text if channel.find("Fluor") is not None else None
                excitation = float(channel.find("ExcitationWavelength").text) if channel.find("ExcitationWavelength") is not None else None
                emission = float(channel.find("EmissionWavelength").text) if channel.find("EmissionWavelength") is not None else None

                optical_channels.append({
                    "name": channel_name,
                    "description": f"{fluor} channel" if fluor else f"Channel {channel_name}",
                    "emission_lambda": emission
                })

            # Create imaging plane metadata
            channel_name_formatted = self.channel_name.replace(" ", "") if self.channel_name else "Default"
            imaging_plane_name = f"ImagingPlane{channel_name_formatted}"
            
            imaging_plane_metadata = {
                "name": imaging_plane_name,
                "optical_channel": optical_channels,
                "description": "Confocal Imaging Plane",
                "device": "ZeissMicroscope",
                "excitation_lambda": excitation if 'excitation' in locals() else None,
                "indicator": fluor if 'fluor' in locals() else "unknown",
                "location": "unknown",
                "grid_spacing": [pixel_size_x, pixel_size_y] if 'pixel_size_x' in locals() and 'pixel_size_y' in locals() else None,
                "grid_spacing_unit": "meters",
                "imaging_rate": self.imaging_extractor.get_sampling_frequency()
            }
            metadata["Ophys"]["ImagingPlane"] = [imaging_plane_metadata]

            # Create two-photon series metadata
            two_photon_series_name = f"ConfocalSeries{channel_name_formatted}"
            two_photon_series_metadata = {
                "name": two_photon_series_name,
                "imaging_plane": imaging_plane_name,
                "unit": "n.a."
            }
            metadata["Ophys"]["TwoPhotonSeries"] = [two_photon_series_metadata]

        return metadata
