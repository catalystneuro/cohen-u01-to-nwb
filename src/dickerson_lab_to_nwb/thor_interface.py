from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path
import numpy as np
import tifffile
from lxml import etree as ET
from datetime import datetime
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType


class ThorTiffImagingExtractor(ImagingExtractor):
    """A ImagingExtractor for multiple TIFF files using OME metadata."""

    extractor_name = "ThorTiffImaging"
    is_writable = False

    def __init__(self, file_path: PathType, channel_name: Optional[str] = None):
        """Create a ThorTiffImagingExtractor instance from a TIFF file.

        Parameters
        ----------
        file_path : str
            Path to first OME TIFF file (e.g., ChanA_001_001_001_001.tif)
        """
        super().__init__()
        self.file_path = Path(file_path)
        self.folder_path = self.file_path.parent
        self.ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2010-06"}
        self._data = None  # Cache for loaded data
        
        # Load experiment metadata first to validate channel name
        xml_path = self.folder_path / "Experiment.xml"
        if xml_path.exists():
            root = ET.parse(str(xml_path)).getroot()
            laser_scanning_microscope = root.xpath("//LSM")[0]

            pixel_x = laser_scanning_microscope.attrib["pixelX"]
            pixel_y = laser_scanning_microscope.attrib["pixelY"]
            pixel_size_um = laser_scanning_microscope.attrib["pixelSizeUM"]
            frame_rate = laser_scanning_microscope.attrib["frameRate"]

            self._sampling_frequency = frame_rate

            wavelenghts = root.xpath("//Wavelength")
            channel_names = [element.attrib["name"] for element in wavelenghts]
            if channel_name is not None and channel_name not in channel_names:
                raise ValueError(
                    f"Channel {channel_name} not available. Channels that are available are {channel_names}"
                )

            self.channel_name = channel_name

        # Read metadata from first file
        with tifffile.TiffFile(self.file_path) as tif:
            self._ome_metadata = tif.ome_metadata

            # Parse OME XML to get file structure
            root = ET.fromstring(self._ome_metadata)

            # Extract dimensions from Pixels element
            pixels = root.find(".//{*}Pixels")
            self._num_channels = int(pixels.get("SizeC"))
            self._num_frames = int(pixels.get("SizeT"))
            self._num_rows = int(pixels.get("SizeY"))
            self._num_columns = int(pixels.get("SizeX"))

            # Get channel names from OME metadata
            self._channel_names = []
            for i in range(self._num_channels):
                channel = pixels.find(f".//ome:Channel[@ID='Channel:0:{i}']", namespaces=self.ns)
                name = f"Channel_{i}" if channel is None else channel.get("Name", f"Channel_{i}")
                self._channel_names.append(name)

            # Update number of channels if filtering
            if self.channel_name is not None:
                # Find channel index from filename pattern
                self._channel_idx = None
                for i, name in enumerate(self._channel_names):
                    if self.channel_name in name:
                        self._channel_idx = i
                        break
                
                if self._channel_idx is not None:
                    self._channel_names = [self._channel_names[self._channel_idx]]
                    self._num_channels = 1

        self._kwargs = {"file_path": str(file_path)}

    def get_frames(self, frame_idxs, channel: int = 0) -> np.ndarray:
        """Get specific frames from the specified channel."""
        # Load data if not already loaded
        if self._data is None:
            with tifffile.TiffFile(self.file_path) as tif:
                self._data = tif.asarray()  # Shape: (frames, channels, height, width)
                if self.channel_name is not None and self._channel_idx is not None:
                    self._data = self._data[:, self._channel_idx:self._channel_idx+1, :, :]

        return self._data[frame_idxs, channel, :, :]

    def get_video(self, start_frame=None, end_frame=None, channel: Optional[int] = 0) -> np.ndarray:
        """Get a range of frames from the specified channel."""
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = self._num_frames

        frame_idxs = range(start_frame, end_frame)
        return self.get_frames(frame_idxs, channel)

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


class ThorImagingInterface(BaseImagingExtractorInterface):
    """Interface for Thor TIFF files with OME metadata."""

    display_name = "ThorLabs TIFF Imaging"
    Extractor = ThorTiffImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        source_schema = super().get_source_schema()
        source_schema["properties"]["file_path"] = dict(
            type="string", description="Path to first OME TIFF file (e.g., ChanA_001_001_001_001.tif)"
        )
        source_schema["properties"]["channel_name"] = dict(
            type="string",
            description="Name of the channel to extract (must match name in Experiment.xml)",
            required=False
        )
        return source_schema

    def __init__(self, file_path: FilePathType, channel_name: Optional[str] = None, verbose: bool = False):
        """
        Initialize reading of TIFF file.

        Parameters
        ----------
        file_path : FilePathType
            Path to first OME TIFF file
        channel_name : str, optional
            Name of the channel to extract (must match name in Experiment.xml)
        verbose : bool, default: False
        """
        super().__init__(file_path=file_path, channel_name=channel_name, verbose=verbose)

    def get_metadata(self) -> DeepDict:
        metadata = super().get_metadata()
        
        
        # Parse Experiment.xml for additional metadata
        xml_path = self.source_data["file_path"].parent / "Experiment.xml"
        if xml_path.exists():
            root = ET.parse(str(xml_path)).getroot()
            
            # Device metadata
            software = root.find(".//Software")
            software_version = software.get("version") if software is not None else None
            device_description = f"ThorLabs 2P Microscope running ThorImageLS {software_version}" if software_version else "ThorLabs 2P Microscope"
            
            metadata["Ophys"] = {
                "Device": [{
                    "name": "ThorMicroscope",
                    "description": device_description
                }]
            }
            
            # LSM metadata
            lsm = root.find(".//LSM")
            if lsm is not None:
                pixel_size = float(lsm.get("pixelSizeUM", 0))  # in micrometers
                frame_rate = float(lsm.get("frameRate", 0))
                width_um = float(lsm.get("widthUM", 0))
                height_um = float(lsm.get("heightUM", 0))
                
                # PMT metadata
                pmt = root.find(".//PMT")
                pmt_gains = {}
                if pmt is not None:
                    pmt_gains["ChanA"] = float(pmt.get("gainA", 0))
                    pmt_gains["ChanB"] = float(pmt.get("gainB", 0))
                
                # Get wavelength info and map channels to indicators
                wavelengths = root.find(".//Wavelengths")
                if wavelengths is not None:
                    channel_elements = wavelengths.findall(".//Wavelength")
                    optical_channels = []
                    channel_indicators = {"ChanA": "tdTomato", "ChanB": "GCaMP"}
                    
                    for channel in channel_elements:
                        name = channel.get("name", "")
                        indicator = channel_indicators.get(name, "unknown")
                        optical_channels.append({
                            "name": name,
                            "description": f"{indicator} channel",
                            "emission_lambda": 520.0  # Placeholder
                        })
                
                # ImagingPlane metadata
                metadata["Ophys"]["ImagingPlane"] = [{
                    "name": "ImagingPlane",
                    "optical_channel": optical_channels,
                    "description": "2P Imaging Plane",
                    "device": "ThorMicroscope",
                    "excitation_lambda": 920.0,  # Placeholder
                    "indicator": "tdTomato and GCaMP",
                    "location": "unknown",
                    "grid_spacing": [pixel_size * 1e-6, pixel_size * 1e-6],  # Convert to meters
                    "grid_spacing_unit": "meters",
                    "imaging_rate": frame_rate
                }]
                
                # TwoPhotonSeries metadata
                selected_channel = self.source_data.get("channel_name")
                pmt_gain = pmt_gains.get(selected_channel) if selected_channel else None
                
                metadata["Ophys"]["TwoPhotonSeries"] = [{
                    "name": "TwoPhotonSeries",
                    "imaging_plane": "ImagingPlane",
                    "field_of_view": [width_um * 1e-6, height_um * 1e-6],  # Convert to meters
                    "pmt_gain": pmt_gain,
                    "scan_line_rate": frame_rate * float(lsm.get("pixelY", 0)),
                    "unit": "n.a."
                }]

        return metadata
