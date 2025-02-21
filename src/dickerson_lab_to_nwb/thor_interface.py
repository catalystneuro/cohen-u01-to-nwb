from collections import defaultdict, namedtuple
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List, Union

import numpy as np
import tifffile
import warnings
from lxml import etree as ET

from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType


class ThorTiffImagingExtractor(ImagingExtractor):
    """
    An ImagingExtractor for multiple TIFF files using OME metadata.
    
    
    This extractor builds a mapping between the T (time) dimension and the corresponding
    pages/IFD or the tiff files using a named tuple structure:
    
    For each time frame (T), we record a list of PageMapping objects that corresponds to the 
    pages of the tiff file that contain the image data for that frame.
    
    Each PageMapping object contains:
      - page_index: The index of the page in the TIFF file (which holds the complete X and Y image data),
      - channel_index: The coordinate along the channel (C) axis (or None if absent),
      - depth_index: The coordinate along the depth (Z) axis (or None if absent).

    When get_frames() is called, the mapping is used to load only the pages for the requested
    frames into a preallocated NumPy array. 
    
    Note: According to the OME specification (see
    https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html#Pixels_DimensionOrder),
    the spatial dimensions (X and Y) are always stored on a single page.
    """

    extractor_name = "ThorTiffImaging"
    is_writable = False

    # Named tuple to hold page mapping details.
    PageMapping = namedtuple("PageMapping", ["page_index", "channel_index", "depth_index"])

    @staticmethod
    def _parse_ome_metadata(metadata_string: str) -> ET._Element:
        """
        Parse an OME metadata string using lxml.etree.
        Removes XML comments if present and attempts to parse as bytes first.
        """
        if metadata_string.lstrip().startswith("<!--"):
            metadata_string = metadata_string.replace("<!--", "").replace("-->", "")
        try:
            return ET.fromstring(metadata_string.encode("utf-8"))
        except ValueError:
            return ET.fromstring(metadata_string)

    def __init__(self, file_path: Union[str, Path], channel_name: Optional[str] = None):
        """
        Create a ThorTiffImagingExtractor instance from a TIFF file.
        """
        super().__init__()
        self.file_path = Path(file_path)
        self.folder_path = self.file_path.parent
        self.channel_name = channel_name
        self._data = None  # Mapping pages on demand

        # Load Experiment.xml metadata if available.
        self._parse_experiment_xml()

        # Open the TIFF file to extract OME metadata and series information.
        with tifffile.TiffFile(self.file_path) as tiff_reader:
            self._ome_metadata = tiff_reader.ome_metadata
            ome_root = self._parse_ome_metadata(self._ome_metadata)
            pixels_element = ome_root.find(".//{*}Pixels")
            if pixels_element is None:
                raise ValueError("Could not find 'Pixels' element in OME metadata.")

            self._num_channels = int(pixels_element.get("SizeC", "1"))
            self._num_frames = int(pixels_element.get("SizeT", "1"))
            self._num_rows = int(pixels_element.get("SizeY"))
            self._num_columns = int(pixels_element.get("SizeX"))
            self._num_z = int(pixels_element.get("SizeZ", "1"))
            self._dimension_order = pixels_element.get("DimensionOrder")

            series = tiff_reader.series[0]
            self._dtype = series.dtype
            number_of_pages = len(series)
            series_axes = series.axes  # e.g., "XYZTC" or "XYCZT"
            series_shape = series.shape

        # Determine non-spatial axes (remove X and Y).
        non_spatial_axes = [axis for axis in series_axes if axis not in ("X", "Y")]
        non_spatial_shape = [dim for axis, dim in zip(series_axes, series_shape) if axis not in ("X", "Y")]

        if "T" not in non_spatial_axes:
            raise ValueError("The TIFF file must have a T (time) dimension. Static images are not supported.")

        total_expected_pages = np.prod(non_spatial_shape)
        if total_expected_pages != number_of_pages:
            warnings.warn(f"Expected {total_expected_pages} pages but found {number_of_pages} pages in the series.")

        # Identify axis indices.
        self._time_axis_index = non_spatial_axes.index("T")
        self._z_axis_index = non_spatial_axes.index("Z") if "Z" in non_spatial_axes else None
        self._channel_axis_index = non_spatial_axes.index("C") if "C" in non_spatial_axes else None

        self._non_spatial_axes = non_spatial_axes
        self._non_spatial_shape = non_spatial_shape

        # Build the mapping from each time frame (T) to its corresponding pages.
        self._frame_page_mapping: Dict[int, List[ThorTiffImagingExtractor.PageMapping]] = defaultdict(list)
        for page_index in range(number_of_pages):
            page_multi_index = np.unravel_index(page_index, non_spatial_shape, order="C")
            time_index = page_multi_index[self._time_axis_index]
            channel_index = (
                page_multi_index[self._channel_axis_index] if self._channel_axis_index is not None else None
            )
            depth_index = (
                page_multi_index[self._z_axis_index] if self._z_axis_index is not None else None
            )
            mapping_entry = ThorTiffImagingExtractor.PageMapping(
                page_index=page_index, channel_index=channel_index, depth_index=depth_index
            )
            self._frame_page_mapping[time_index].append(mapping_entry)

        self._kwargs = {"file_path": str(file_path)}

    def _parse_experiment_xml(self) -> None:
        """
        Helper function to parse Experiment.xml and extract metadata such as frame rate
        and channel names.
        """
        experiment_xml_path = self.folder_path / "Experiment.xml"
        if experiment_xml_path.exists():
            experiment_tree = ET.parse(str(experiment_xml_path))
            experiment_root = experiment_tree.getroot()
            lsm_element = experiment_root.find(".//LSM")
            if lsm_element is not None:
                frame_rate = lsm_element.attrib.get("frameRate")
                self._sampling_frequency = float(frame_rate) if frame_rate is not None else None
            else:
                raise ValueError("Could not find 'LSM' element in Experiment.xml which contains the frame rate.")

            wavelength_elements = experiment_root.findall(".//Wavelength")
            self._channel_names = [
                wavelength.attrib.get("name") for wavelength in wavelength_elements if "name" in wavelength.attrib
            ]
            if self.channel_name is not None and self.channel_name not in self._channel_names:
                raise ValueError(
                    f"Channel '{self.channel_name}' not available. Available channels: {self._channel_names}"
                )
            # Set channel filter if a channel name is provided.
            if self.channel_name is not None:
                self._channel_index_for_filter = None
                for index, name in enumerate(self._channel_names):
                    if self.channel_name in name:
                        self._channel_index_for_filter = index
                        break
                if self._channel_index_for_filter is None:
                    raise ValueError(f"Channel '{self.channel_name}' not found in Experiment.xml.")
                # Update channel names and count based on the filter.
                self._channel_names = [self._channel_names[self._channel_index_for_filter]]
                self._num_channels = 1
        else:
            # If Experiment.xml is missing, set defaults.
            self._sampling_frequency = None
            self._channel_names = []
    
    def get_frames(self, frame_idxs: List[int]) -> np.ndarray:
        """
        Get specific frames by their time indices.

        Parameters
        ----------
        frame_idxs : List[int]
            List of time/frame indices to retrieve.

        Returns
        -------
        np.ndarray
            Array of shape (n_frames, height, width) if no depth, or
            (n_frames, height, width, n_z) if a Z dimension exists.
        """
        with tifffile.TiffFile(self.file_path) as tiff_reader:
            series = tiff_reader.series[0]
            data_type = series.dtype
            image_height = self._num_rows
            image_width = self._num_columns

            has_z_dimension = self._z_axis_index is not None and self._num_z > 1
            number_of_z_planes = self._num_z if has_z_dimension else 1

            n_frames = len(frame_idxs)
            output_shape = (n_frames, image_height, image_width, number_of_z_planes) if has_z_dimension \
                else (n_frames, image_height, image_width)
            output_array = np.empty(output_shape, dtype=data_type)

            for frame_counter, frame_idx in enumerate(frame_idxs):
                if frame_idx not in self._frame_page_mapping:
                    raise ValueError(f"No pages found for frame {frame_idx}.")
                page_mappings = self._frame_page_mapping[frame_idx]

                # Filter by channel if a channel name was provided.
                if self._channel_axis_index is not None and self.channel_name is not None:
                    filter_index = self._channel_index_for_filter
                    page_mappings = [m for m in page_mappings if m.channel_index == filter_index]

                if has_z_dimension:
                    page_mappings.sort(key=lambda entry: entry.depth_index)
                    if len(page_mappings) != number_of_z_planes:
                        raise ValueError(
                            f"Expected {number_of_z_planes} pages for frame {frame_idx} but got {len(page_mappings)}."
                        )
                    for depth_counter, mapping_entry in enumerate(page_mappings):
                        page_data = series.pages[mapping_entry.page_index].asarray()
                        output_array[frame_counter, :, :, depth_counter] = page_data
                else:
                    if len(page_mappings) != 1:
                        raise ValueError(
                            f"Expected 1 page for frame {frame_idx} but got {len(page_mappings)}."
                        )
                    single_page_index = page_mappings[0].page_index
                    page_data = series.pages[single_page_index].asarray()
                    output_array[frame_counter, :, :] = page_data

        return output_array

    def get_video(self, start_frame: Optional[int] = None, end_frame: Optional[int] = None) -> np.ndarray:
        """
        Get a range of frames.
        """
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = self._num_frames
        frame_indices = list(range(start_frame, end_frame))
        return self.get_frames(frame_indices)

    def get_image_size(self) -> Tuple[int, int]:
        """Return the image dimensions (height, width)."""
        return self._num_rows, self._num_columns

    def get_num_frames(self) -> int:
        """Return the number of frames (time points)."""
        return self._num_frames

    def get_sampling_frequency(self) -> Optional[float]:
        """Return the sampling frequency, if available."""
        return self._sampling_frequency

    def get_num_channels(self) -> int:
        """Return the number of channels."""
        return self._num_channels

    def get_channel_names(self) -> List[str]:
        """Return the channel names."""
        return self._channel_names
    
        """Return the channel names."""
        return self._channel_names

    def get_dtype(self):
        return self._dtype


class ThorImagingInterface(BaseImagingExtractorInterface):
    """
    Interface for Thor TIFF files with OME metadata.
    """
    display_name = "ThorLabs TIFF Imaging"
    Extractor = ThorTiffImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        source_schema = super().get_source_schema()
        source_schema["properties"]["file_path"] = {
            "type": "string",
            "description": "Path to first OME TIFF file (e.g., ChanA_001_001_001_001.tif)"
        }
        source_schema["properties"]["channel_name"] = {
            "type": "string",
            "description": "Name of the channel to extract (must match name in Experiment.xml)",
            "required": False
        }
        return source_schema

    def __init__(self, file_path: FilePathType, channel_name: Optional[str] = None, verbose: bool = False):
        """
        Initialize reading of a TIFF file.
        """
        super().__init__(file_path=file_path, channel_name=channel_name, verbose=verbose)
        self.channel_name = channel_name

    def get_metadata(self) -> DeepDict:
        metadata = super().get_metadata()

        xml_path = Path(self.source_data["file_path"]).parent / "Experiment.xml"
        if xml_path.exists():
            root = ET.parse(str(xml_path)).getroot()

            # Device metadata
            software = root.find(".//Software")
            software_version = software.get("version") if software is not None else None
            device_description = (
                f"ThorLabs 2P Microscope running ThorImageLS {software_version}"
                if software_version
                else "ThorLabs 2P Microscope"
            )

            date_element = root.find(".//Date")
            date_attribute = date_element.attrib["date"]
            self.session_start_time = datetime.strptime(date_attribute, "%m/%d/%Y %H:%M:%S")
            metadata["NWBFile"]["session_start_time"] = self.session_start_time

            metadata.setdefault("Ophys", {})["Device"] = [{
                "name": "ThorMicroscope",
                "description": device_description
            }]

            # LSM metadata
            lsm = root.find(".//LSM")
            if lsm is not None:
                pixel_size = float(lsm.get("pixelSizeUM", 0))
                frame_rate = float(lsm.get("frameRate", 0))
                width_um = float(lsm.get("widthUM", 0))
                height_um = float(lsm.get("heightUM", 0))

                # PMT metadata
                pmt = root.find(".//PMT")
                pmt_gains = {}
                if pmt is not None:
                    pmt_gains["ChanA"] = float(pmt.get("gainA", 0))
                    pmt_gains["ChanB"] = float(pmt.get("gainB", 0))

                # Define channel indicators outside of wavelength processing.
                channel_indicators = {"ChanA": "tdTomato", "ChanB": "GCaMP"}

                wavelengths = root.find(".//Wavelengths")
                optical_channels = []
                if wavelengths is not None:
                    channel_elements = wavelengths.findall(".//Wavelength")
                    for channel in channel_elements:
                        name = channel.get("name", "")
                        indicator = channel_indicators.get(name, "unknown")
                        optical_channels.append({
                            "name": name,
                            "description": f"{indicator} channel",
                            "emission_lambda": 520.0  # Placeholder
                        })

                # Helper to convert to CamelCase
                def to_camel_case(s: str) -> str:
                    parts = s.split('_')
                    return parts[0] + ''.join(word.title() for word in parts[1:])

                ChannelName = to_camel_case(self.channel_name) if self.channel_name else "Default"

                imaging_plane_name = f"ImagingPlane{ChannelName}"
                channel_imaging_plane_metadata = {
                    "name": imaging_plane_name,
                    "optical_channel": optical_channels,
                    "description": "2P Imaging Plane",
                    "device": "ThorMicroscope",
                    "excitation_lambda": 920.0,  # Placeholder
                    "indicator": channel_indicators.get(self.channel_name, "unknown"),
                    "location": "unknown",  # TODO: Extract if possible.
                    "grid_spacing": [pixel_size * 1e-6, pixel_size * 1e-6],  # Convert um to meters
                    "grid_spacing_unit": "meters",
                    "imaging_rate": frame_rate
                }
                metadata["Ophys"]["ImagingPlane"] = [channel_imaging_plane_metadata]

                selected_channel = self.source_data.get("channel_name")
                pmt_gain = pmt_gains.get(selected_channel) if selected_channel else None

                two_photon_series_name = f"TwoPhotonSeries{ChannelName}"
                two_photon_series_metadata = {
                    "name": two_photon_series_name,
                    "imaging_plane": imaging_plane_name,
                    "field_of_view": [width_um * 1e-6, height_um * 1e-6],  # Convert um to meters
                    "pmt_gain": pmt_gain,
                    "scan_line_rate": frame_rate * float(lsm.get("pixelY", 0)),
                    "unit": "n.a."
                }
                metadata["Ophys"]["TwoPhotonSeries"] = [two_photon_series_metadata]

        return metadata
