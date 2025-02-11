
from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path
import numpy as np
import tifffile
import xml.etree.ElementTree as ET
from datetime import datetime
from neuroconv.datainterfaces.ophys.baseimagingextractor import BaseImagingExtractorInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType

class ThorTiffImagingExtractor(ImagingExtractor):
    """A ImagingExtractor for multiple TIFF files using OME metadata."""

    extractor_name = "ThorTiffImaging"
    is_writable = False

    def __init__(self, file_path: PathType):
        """Create a ThorTiffImagingExtractor instance from a TIFF file.

        Parameters
        ----------
        file_path : str
            Path to first OME TIFF file (e.g., ChanA_001_001_001_001.tif)
        """
        super().__init__()
        self.file_path = Path(file_path)
        self.folder_path = self.file_path.parent
        
        # Read first file and its OME metadata
        with tifffile.TiffFile(self.file_path) as tif:
            self._ome_metadata = tif.ome_metadata
            
            # Parse OME XML to get file structure
            root = ET.fromstring(self._ome_metadata)
            
            # Extract dimensions from Pixels element
            pixels = root.find(".//{*}Pixels")
            self._num_channels = int(pixels.get("SizeC", 1))
            self._num_frames = int(pixels.get("SizeT", 1))
            self._num_rows = int(pixels.get("SizeY"))
            self._num_columns = int(pixels.get("SizeX"))
            
            # Get channel names
            ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2010-06'}
            self._channel_names = []
            for i in range(self._num_channels):
                channel = pixels.find(f".//ome:Channel[@ID='Channel:0:{i}']", namespaces=ns)
                name = f"Channel_{i}" if channel is None else channel.get("Name", f"Channel_{i}")
                self._channel_names.append(name)
            
            # Map files to their positions
            self._file_map = {}
            for tiff_data in pixels.findall(".//ome:TiffData", namespaces=ns):
                t = int(tiff_data.get("FirstT", 0))
                c = int(tiff_data.get("FirstC", 0))
                uuid_elem = tiff_data.find(".//ome:UUID", namespaces=ns)
                if uuid_elem is not None:
                    filename = uuid_elem.get("FileName")
                    self._file_map[(t, c)] = self.folder_path / filename
        
        # Load experiment metadata
        xml_path = self.folder_path / "Experiment.xml"
        if xml_path.exists():
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Extract frame rate
            if root.find(".//frameRate") is not None:
                self._sampling_frequency = float(root.find(".//frameRate").text)
            else:
                self._sampling_frequency = None
                
            # Extract start time
            start_time_elem = root.find(".//startTime")
            if start_time_elem is not None:
                try:
                    self.start_time = datetime.strptime(
                        start_time_elem.text,
                        "%Y-%m-%d %H:%M:%S"
                    )
                except ValueError:
                    self.start_time = None
            else:
                self.start_time = None
        else:
            self._sampling_frequency = None
            self.start_time = None
            
        self._kwargs = {"file_path": str(file_path)}

    def get_frames(self, frame_idxs, channel: int = 0) -> np.ndarray:
        """Get specific frames from the specified channel."""
        frames = []
        for frame_idx in frame_idxs:
            file_path = self._file_map.get((frame_idx, channel))
            if file_path and file_path.exists():
                with tifffile.TiffFile(file_path) as tif:
                    frames.append(tif.asarray())
        return np.array(frames)

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


class ThorTiffImagingInterface(BaseImagingExtractorInterface):
    """Interface for Thor TIFF files with OME metadata."""

    display_name = "ThorLabs TIFF Imaging"
    Extractor = ThorTiffImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        source_schema = super().get_source_schema()
        source_schema["properties"]["file_path"] = dict(
            type="string",
            description="Path to first OME TIFF file (e.g., ChanA_001_001_001_001.tif)"
        )
        return source_schema

    def __init__(self, file_path: FilePathType, verbose: bool = False):
        """
        Initialize reading of TIFF file.

        Parameters
        ----------
        file_path : FilePathType
            Path to first OME TIFF file
        verbose : bool, default: False
        """
        super().__init__(file_path=file_path, verbose=verbose)

    def get_metadata(self, photon_series_type="TwoPhotonSeries") -> DeepDict:
        metadata = super().get_metadata(photon_series_type=photon_series_type)

        # Add session start time if available
        if self.extractor.start_time:
            metadata["NWBFile"]["session_start_time"] = self.extractor.start_time
            
        # Add frame rate if available
        if self.extractor._sampling_frequency:
            metadata["Ophys"]["frame_rate"] = self.extractor._sampling_frequency
            
        # Add OME metadata
        if hasattr(self.extractor, '_ome_metadata'):
            metadata["OME"] = self.extractor._ome_metadata
            
        return metadata
