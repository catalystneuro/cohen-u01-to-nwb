from collections import defaultdict, namedtuple
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List, Union

import numpy as np
import czifile
import warnings
import xml.etree.ElementTree as ET

from neuroconv.basedatainterface import BaseDataInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType
from pynwb.ophys import OpticalChannel, ImagingPlane
from pynwb.image import GrayscaleImage
from pynwb.base import Images


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
            
            if image_info is None:
                # Try alternative paths for image dimensions
                self._num_channels = 1
                self._num_frames = 1
                self._num_z = 1
                
                # Try to get dimensions from the shape of the data array
                shape = czi_reader.shape
                axes = czi_reader.axes
                
                if 'C' in axes:
                    self._num_channels = shape[axes.index('C')]
                if 'T' in axes:
                    self._num_frames = shape[axes.index('T')]
                if 'Z' in axes:
                    self._num_z = shape[axes.index('Z')]
                if 'Y' in axes:
                    self._num_rows = shape[axes.index('Y')]
                else:
                    raise ValueError("Could not determine Y dimension from CZI file")
                if 'X' in axes:
                    self._num_columns = shape[axes.index('X')]
                else:
                    raise ValueError("Could not determine X dimension from CZI file")
            else:
                # Extract dimensions from the image_info element
                self._num_channels = int(image_info.find("SizeC").text) if image_info.find("SizeC") is not None else 1
                self._num_frames = int(image_info.find("SizeT").text) if image_info.find("SizeT") is not None else 1
                self._num_z = int(image_info.find("SizeZ").text) if image_info.find("SizeZ") is not None else 1
                self._num_rows = int(image_info.find("SizeY").text)
                self._num_columns = int(image_info.find("SizeX").text)
            
            # Check if there's more than one time point
            if self._num_frames > 1:
                raise ValueError("This interface only supports CZI files with a single time point. "
                                 f"Found {self._num_frames} time points.")
            
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


class ZeissConfocalInterface(BaseDataInterface):
    """
    Interface for Zeiss confocal microscopy CZI files.
    
    This interface inherits from BaseDataInterface and handles the extraction
    of confocal microscopy data from Zeiss CZI files, saving the ImagingPlane
    and images from both channels as an Images container.
    """
    display_name = "Zeiss Confocal Imaging"

    def __init__(self, file_path: FilePathType, channel_names: Optional[List[str]] = None, verbose: bool = False):
        """
        Initialize the interface.
        
        Parameters
        ----------
        file_path : FilePathType
            Path to the Zeiss .czi file
        channel_names : Optional[List[str]], default: None
            List of channel names to extract (e.g., ['Alexa Fluor 488', 'Alexa Fluor 633'])
            If None, all available channels will be extracted
        verbose : bool, default: False
            If True, print verbose output
        """
        super().__init__(file_path=file_path, verbose=verbose)
        self.file_path = Path(file_path)
        self.channel_names = channel_names
        self.extractors = {}
        
        # Create extractors for each channel
        if channel_names:
            for channel_name in channel_names:
                self.extractors[channel_name] = ZeissConfocalImagingExtractor(
                    file_path=file_path, 
                    channel_name=channel_name
                )
        else:
            # If no channels specified, create a single extractor without channel filtering
            extractor = ZeissConfocalImagingExtractor(file_path=file_path)
            # Get all available channels
            available_channels = extractor.get_channel_names()
            for channel_name in available_channels:
                self.extractors[channel_name] = ZeissConfocalImagingExtractor(
                    file_path=file_path, 
                    channel_name=channel_name
                )
            self.channel_names = available_channels

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
        with czifile.CziFile(self.file_path) as czi_reader:
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
            optical_channels_metadata = []
            
            for channel in channels:
                channel_name = channel.find("Name").text if channel.find("Name") is not None else None
                if channel_name is None:
                    continue

                fluor = channel.find("Fluor").text if channel.find("Fluor") is not None else None
                excitation = float(channel.find("ExcitationWavelength").text) if channel.find("ExcitationWavelength") is not None else None
                emission = float(channel.find("EmissionWavelength").text) if channel.find("EmissionWavelength") is not None else None

                optical_channels_metadata.append({
                    "name": channel_name,
                    "description": f"{fluor} channel" if fluor else f"Channel {channel_name}",
                    "emission_lambda": emission
                })

            # Create imaging plane metadata
            imaging_plane_metadata = {
                "name": "ImagingPlane",
                "optical_channel": optical_channels_metadata,
                "description": "Confocal Imaging Plane",
                "device": "ZeissMicroscope",
                "excitation_lambda": excitation if 'excitation' in locals() else None,
                "indicator": "multiple" if len(optical_channels_metadata) > 1 else (fluor if 'fluor' in locals() else "unknown"),
                "location": "unknown",
                "grid_spacing": [pixel_size_x, pixel_size_y] if 'pixel_size_x' in locals() and 'pixel_size_y' in locals() else None,
                "grid_spacing_unit": "meters",
                "imaging_rate": None  # Confocal images typically don't have a frame rate
            }
            metadata["Ophys"]["ImagingPlane"] = [imaging_plane_metadata]

            # Create images metadata for each channel
            images_metadata = []
            for channel_metadata in optical_channels_metadata:
                channel_name = channel_metadata["name"]
                if self.channel_names and channel_name not in self.channel_names:
                    continue
                
                channel_name_formatted = channel_name.replace(" ", "")
                images_metadata.append({
                    "name": f"ConfocalImages{channel_name_formatted}",
                    "description": f"Confocal images for {channel_name}"
                })
            
            metadata["Ophys"]["Images"] = images_metadata

        return metadata
        
    def add_to_nwbfile(self, nwbfile, metadata: Optional[dict] = None):
        """
        Add the confocal data to an NWB file.
        
        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the confocal data to
        metadata : Optional[dict], default: None
            Metadata dictionary
        """
        if metadata is None:
            metadata = self.get_metadata()
            
        # Get device
        device_metadata = metadata["Ophys"]["Device"][0]
        device = nwbfile.create_device(
            name=device_metadata["name"],
            description=device_metadata["description"]
        )
        
        # Create at least one optical channel (required by NWB)
        # Even if we don't have channel metadata, we need to create a default optical channel
        optical_channels = []
        
        if len(metadata["Ophys"]["ImagingPlane"][0]["optical_channel"]) > 0:
            # Use the channel metadata if available
            for channel_metadata in metadata["Ophys"]["ImagingPlane"][0]["optical_channel"]:
                if self.channel_names and channel_metadata["name"] not in self.channel_names:
                    continue
                    
                # Ensure emission_lambda is not None (required by NWB)
                emission_lambda = channel_metadata.get("emission_lambda")
                if emission_lambda is None:
                    # Use a default value if not available
                    emission_lambda = 520.0  # Common emission wavelength for fluorescence microscopy
                    
                optical_channel = OpticalChannel(
                    name=channel_metadata["name"],
                    description=channel_metadata["description"],
                    emission_lambda=emission_lambda
                )
                optical_channels.append(optical_channel)
        
        # If no optical channels were created, create a default one
        if len(optical_channels) == 0:
            optical_channel = OpticalChannel(
                name="OpticalChannel",
                description="Default optical channel",
                emission_lambda=520.0
            )
            optical_channels.append(optical_channel)
            
        # Create imaging plane
        imaging_plane_metadata = metadata["Ophys"]["ImagingPlane"][0]
        grid_spacing = imaging_plane_metadata.get("grid_spacing")
        grid_spacing_unit = imaging_plane_metadata.get("grid_spacing_unit")
        
        # Ensure excitation_lambda is not None (required by NWB)
        excitation_lambda = np.nan
        
        imaging_plane = nwbfile.create_imaging_plane(
            name=imaging_plane_metadata["name"],
            optical_channel=optical_channels,
            description=imaging_plane_metadata["description"],
            device=device,
            excitation_lambda=excitation_lambda,
            indicator=imaging_plane_metadata.get("indicator", "unknown"),
            location=imaging_plane_metadata.get("location", "unknown"),
            grid_spacing=grid_spacing,
            grid_spacing_unit=grid_spacing_unit
        )
        
        # Add image data for each channel as a single Images container
        for channel_name, extractor in self.extractors.items():
            channel_name_formatted = channel_name.replace(" ", "")
            images_name = f"ConfocalImages{channel_name_formatted}"
            
            # Get raw data directly from the extractor
            raw_data = extractor._data
            
            # Get the channel index (it will always be available since we create an extractor per channel)
            channel_index = extractor._channel_index
            
            # Create a list to hold all images for this channel
            grayscale_images = []
            
            # For Z-stack data
            if extractor._num_z > 1:
                # For each Z slice
                for z in range(extractor._num_z):
                    # Extract the image data for this channel and Z slice
                    # The shape of raw_data might vary, so we need to handle different cases
                    if len(raw_data.shape) == 5:  # (T, C, Z, Y, X)
                        image_data = raw_data[0, channel_index, z, :, :]
                    elif len(raw_data.shape) == 4:  # (C, Z, Y, X)
                        image_data = raw_data[channel_index, z, :, :]
                    else:
                        raise ValueError(f"Unexpected data shape: {raw_data.shape}")
                    
                    # Create a GrayscaleImage for this slice
                    image_name = f"{images_name}_z{z}"
                    grayscale_image = GrayscaleImage(
                        name=image_name,
                        data=image_data,
                        description=f"Confocal image for {channel_name}, depth index: {z}"
                    )
                    grayscale_images.append(grayscale_image)
            else:
                # For 2D data (no Z stack)
                if len(raw_data.shape) == 4:  # (T, C, Y, X)
                    image_data = raw_data[0, channel_index, :, :]
                elif len(raw_data.shape) == 3:  # (C, Y, X)
                    image_data = raw_data[channel_index, :, :]
                else:
                    raise ValueError(f"Unexpected data shape: {raw_data.shape}")
                
                # Create a GrayscaleImage
                grayscale_image = GrayscaleImage(
                    name=images_name,
                    data=image_data,
                    description=f"Confocal image for {channel_name}"
                )
                grayscale_images.append(grayscale_image)
            
            # Create a single Images container with all the GrayscaleImage objects for this channel
            images_container = Images(
                name=images_name,
                images=grayscale_images,
                description=f"Z-stack of confocal images for {channel_name} with {len(grayscale_images)} slices"
            )
            
            # Add the Images container to the NWB file
            nwbfile.add_acquisition(images_container)
