from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path
import numpy as np
import tifffile
from lxml import etree as ET
from datetime import datetime
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from roiextractors import ImagingExtractor
from neuroconv.utils import FilePathType, DeepDict, PathType


import numpy as np
import tifffile
import warnings
import xml.etree.ElementTree as ET  # For Experiment.xml parsing
from pathlib import Path
from typing import Optional, Tuple, Union
from collections import defaultdict, namedtuple

# Assume ImagingExtractor and PathType are defined elsewhere.
# For example:
# class ImagingExtractor:
#     pass
# PathType = Union[str, Path]

class ThorTiffImagingExtractor(ImagingExtractor):
    """An ImagingExtractor for multiple TIFF files using OME metadata.

    This extractor builds a mapping between the T (time) dimension and the corresponding
    pages using a named tuple structure. For each time frame (T), we record a list of
    PageMapping objects. Each PageMapping object contains:
      - page_index: The index of the page in the TIFF file (which holds the complete X and Y image data),
      - channel_index: The coordinate along the channel (C) axis (or None if absent),
      - depth_index: The coordinate along the depth (Z) axis (or None if absent).

    When get_frames() is called, the mapping is used to load only the pages for the requested
    frames into a preallocated NumPy array. If a Z dimension exists (SizeZ > 1), the output array
    is of shape (n_frames, height, width, n_z) with the depth (Z) as the last dimension.
    Otherwise, it is of shape (n_frames, height, width).

    Note: According to the OME specification (see
    https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html#Pixels_DimensionOrder),
    the spatial dimensions (X and Y) are stored on a single page.
    """

    extractor_name = "ThorTiffImaging"
    is_writable = False

    # According to the OME XML schema for Pixels_DimensionOrder
    # (https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html#Pixels_DimensionOrder),
    # XY and are always present in that order and should be totally contained within a tiff page/IFD.
    # ------------------------


    # Define a named tuple to hold page mapping details.
    # Each PageMapping object contains:
    #   - page_index: The index of the page in the TIFF file (the page contains the full X and Y image data),
    #   - channel_index: The coordinate along the channel (C) axis (or None if no channel dimension exists),
    #   - depth_index: The coordinate along the depth (Z) axis (or None if no depth dimension exists).
    PageMapping = namedtuple("PageMapping", ["page_index", "channel_index", "depth_index"])

    @staticmethod
    def _parse_ome_metadata(metadata_string: str):
        """
        Parse an OME metadata string using lxml.etree.

        This function removes XML comments if present and attempts to parse
        the metadata string. It first attempts to parse the metadata as bytes.
        If that fails, it attempts to parse it as an XML fragment. Older ome metadata was written as a comment
        instead of proper XML.
        """
        import lxml.etree
        if metadata_string.lstrip().startswith("<!--"):
            # Remove XML comments
            metadata_string = metadata_string.replace("<!--", "").replace("-->", "")
        try:
            return lxml.etree.fromstring(metadata_string.encode("utf-8"))
        except ValueError:
            return lxml.etree.fromstring(metadata_string)

    def __init__(self, file_path: Union[str, Path], channel_name: Optional[str] = None):
        """
        Create a ThorTiffImagingExtractor instance from a TIFF file.

        Parameters
        ----------
        file_path : str or Path
            Path to the first OME TIFF file (e.g., ChanA_001_001_001_001.tif)
        channel_name : str, optional
            If provided, filter to the channel with a matching name.
        """
        super().__init__()
        self.file_path = Path(file_path)
        self.folder_path = self.file_path.parent
        self.ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2010-06"}
        self._data = None  # Not used since we map pages on demand
        self.channel_name = channel_name

        # --- Load Experiment Metadata if available ---
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
            experiment_channel_names = [
                wavelength.attrib.get("name") for wavelength in wavelength_elements if "name" in wavelength.attrib
            ]
            if channel_name is not None and channel_name not in experiment_channel_names:
                raise ValueError(
                    f"Channel '{channel_name}' not available. Available channels: {experiment_channel_names}"
                )
            
            self._channel_names = experiment_channel_names
            
            # If filtering by channel name, update the channel index and count.
            if channel_name is not None:
                self._channel_index_for_filter = None
                for index, name in enumerate(self._channel_names):
                    if channel_name in name:
                        self._channel_index_for_filter = index
                        break
                if self._channel_index_for_filter is not None:
                    self._channel_names = [self._channel_names[self._channel_index_for_filter]]
                    self._num_channels = 1
                else:
                    raise ValueError(f"Channel '{channel_name}' not found in Experiment.xml.")
                
        # --- Read OME Metadata from the TIFF file ---
        with tifffile.TiffFile(self.file_path) as tiff_reader:
            self._ome_metadata = tiff_reader.ome_metadata
            ome_root = self._parse_ome_metadata(self._ome_metadata)
            pixels_element = ome_root.find(".//{*}Pixels")
            if pixels_element is None:
                raise ValueError("Could not find 'Pixels' element in OME metadata.")
            
            # Extract dimensions from the Pixels element (defaulting missing values to 1)
            self._num_channels = int(pixels_element.get("SizeC", "1"))
            self._num_frames = int(pixels_element.get("SizeT", "1"))
            self._num_rows = int(pixels_element.get("SizeY"))
            self._num_columns = int(pixels_element.get("SizeX"))
            self._num_z = int(pixels_element.get("SizeZ", "1"))
            self._dimension_order = pixels_element.get("DimensionOrder")
            
            # --- Build the mapping from each time frame (T) to its corresponding pages ---
            # For each page in the TIFF file (accessed via tiff_reader.pages), we compute its multi-index
            # along the non-spatial dimensions (T, Z, and C).
            # For each page we create a PageMapping named tuple that records:
            #   - page_index: the index of the page in the TIFF file,
            #   - channel_index: the coordinate along the C axis (or None if absent),
            #   - depth_index: the coordinate along the Z axis (or None if absent).
            # We then group these PageMapping objects by the time coordinate (T) so that we have a
            # mapping from each time frame to a list of PageMapping entries.
            series = tiff_reader.series[0]
            number_of_pages = len(series)
            series_axes = series.axes  # e.g., "XYZTC" or "XYCZT"
            series_shape = series.shape

        # Remove spatial dimensions "X" and "Y" so that the remaining axes determine the pages.
        non_spatial_axes = [axis for axis in series_axes if axis not in ("X", "Y")]
        non_spatial_shape = [dim for axis, dim in zip(series_axes, series_shape) if axis not in ("X", "Y")]

        # Throw an error if the T (time) dimension is not available.
        if "T" not in non_spatial_axes:
            raise ValueError("The TIFF file must have a T (time) dimension. Static images are not supported.")

        total_expected_pages = np.prod(non_spatial_shape)
        if total_expected_pages != len(series):
            warnings.warn(f"Expected {total_expected_pages} pages but found {len(series)} pages in the series.")

        # Determine the indices for T, Z, and C within the non-spatial axes.
        self._time_axis_index = non_spatial_axes.index("T")
        self._z_axis_index = non_spatial_axes.index("Z") if "Z" in non_spatial_axes else None
        self._channel_axis_index = non_spatial_axes.index("C") if "C" in non_spatial_axes else None

        self._non_spatial_axes = non_spatial_axes
        self._non_spatial_shape = non_spatial_shape

        self._frame_page_mapping = defaultdict(list)
        for page_index in range(number_of_pages):
            page_multi_index = np.unravel_index(page_index, non_spatial_shape, order="C")
            time_index = page_multi_index[self._time_axis_index]
            channel_index = (
                page_multi_index[self._channel_axis_index]
                if self._channel_axis_index is not None
                else None
            )
            depth_index = (
                page_multi_index[self._z_axis_index]
                if self._z_axis_index is not None
                else None
            )
            mapping_entry = ThorTiffImagingExtractor.PageMapping(
                page_index=page_index,
                channel_index=channel_index,
                depth_index=depth_index,
            )
            self._frame_page_mapping[time_index].append(mapping_entry)
                
        self._kwargs = {"file_path": str(file_path)}

    def get_frames(self, frame_idxs, channel: int = 0) -> np.ndarray:
        """
        Get specific frames for the specified channel without loading all data.

        This method uses the precomputed mapping between T (time) frames and the corresponding
        pages (computed in __init__). For each requested frame, it preallocates an output array
        and fills it with the page data. If a Z dimension exists (SizeZ > 1), the output array is
        of shape (n_frames, height, width, n_z) with the depth (Z) as the last dimension.
        Otherwise, the output array is of shape (n_frames, height, width).

        Parameters
        ----------
        frame_idxs : sequence of int
            The time/frame indices to retrieve.
        channel : int, optional
            The channel index to retrieve. Pages are filtered based on their channel_index.
            (If the TIFF file has no channel dimension, this parameter is ignored.)

        Returns
        -------
        np.ndarray
            If no Z dimension is present: an array of shape (n_frames, height, width).
            If a Z dimension exists: an array of shape (n_frames, height, width, n_z).
        """
        with tifffile.TiffFile(self.file_path) as tiff_reader:
            series = tiff_reader.series[0]
            data_type = series.dtype
            image_height = self._num_rows
            image_width = self._num_columns

            # Determine whether a Z dimension exists.
            if self._z_axis_index is not None and self._num_z > 1:
                has_z_dimension = True
                number_of_z_planes = self._num_z
            else:
                has_z_dimension = False
                number_of_z_planes = 1

            # Preallocate the output array.
            number_of_frames = len(frame_idxs)
            if has_z_dimension:
                output_shape = (number_of_frames, image_height, image_width, number_of_z_planes)
            else:
                output_shape = (number_of_frames, image_height, image_width)

            output_array = np.empty(output_shape, dtype=data_type)

            # Loop over the requested frames.
            for frame_counter, time_frame in enumerate(frame_idxs):
                if time_frame not in self._frame_page_mapping:
                    raise ValueError(f"No pages found for frame {time_frame}.")
                # Get the list of PageMapping entries for this frame.
                page_mappings = self._frame_page_mapping[time_frame]
                # If the file contains a channel dimension, filter for the requested channel.
                if self._channel_axis_index is not None:
                    
                    filter_index = self._channel_index_for_filter if self.channel_name is not None else channel
                    page_mappings = [mapping for mapping in page_mappings if mapping.channel_index == filter_index]

                if has_z_dimension:
                    # Sort the pages by the depth (Z) coordinate.
                    page_mappings.sort(key=lambda entry: entry.depth_index)
                    if len(page_mappings) != number_of_z_planes:
                        raise ValueError(
                            f"Expected {number_of_z_planes} pages for frame {time_frame} but got {len(page_mappings)}."
                        )
                    # For each depth plane, load the corresponding page.
                    for depth_counter, mapping_entry in enumerate(page_mappings):
                        page_data = tiff_reader.pages[mapping_entry.page_index].asarray()
                        output_array[frame_counter, :, :, depth_counter] = page_data
                else:
                    # In the case without depth or channel, there is only one TIFF page per frame.
                    # Therefore, we load the single page using the first (and only) mapping entry.
                    if len(page_mappings) != 1:
                        raise ValueError(
                            f"Expected 1 page for frame {time_frame} but got {len(page_mappings)}."
                        )
                    single_page_index = page_mappings[0].page_index
                    page_data = series.pages[single_page_index].asarray()
                    output_array[frame_counter, :, :] = page_data
            return output_array

    def get_video(self, start_frame: Optional[int] = None, end_frame: Optional[int] = None, channel: int = 0) -> np.ndarray:
        """Get a range of frames for the specified channel."""
        if start_frame is None:
            start_frame = 0
        if end_frame is None:
            end_frame = self._num_frames
        frame_indices = list(range(start_frame, end_frame))
        return self.get_frames(frame_indices, channel)

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

    def get_channel_names(self):
        """Return the channel names."""
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
        self.channel_name  = channel_name
        
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
        
            date_element = root.find(".//Date")
            # date_attribute like this date="05/06/2024 13:55:12"
            date_attribute = date_element.attrib["date"]
            self.session_start_time = datetime.strptime(date_attribute, "%m/%d/%Y %H:%M:%S")
        
            metadata["NWBFile"]["session_start_time"] = self.session_start_time
        
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
                
                to_camel_case = lambda s: s.split('_')[0] + ''.join(x.title() for x in s.split('_')[1:])
                ChannelName = to_camel_case(self.channel_name)

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
                indicator = channel_indicators.get(self.channel_name, "unknown")
                
                imaging_plane_name = f"ImagingPlane{ChannelName}"
                channel_imaging_plane_metadata = {
                    "name": imaging_plane_name,
                    "optical_channel": optical_channels,
                    "description": "2P Imaging Plane",
                    "device": "ThorMicroscope",
                    "excitation_lambda": 920.0,  # Placeholder
                    "indicator": indicator,
                    "location": "unknown",  # TODO: can it be extracted from the ome xml?
                    "grid_spacing": [pixel_size * 1e-6, pixel_size * 1e-6],  # Convert to meters
                    "grid_spacing_unit": "meters",
                    "imaging_rate": frame_rate
                }
                metadata["Ophys"]["ImagingPlane"] = [channel_imaging_plane_metadata]
                
                # TwoPhotonSeries metadata
                selected_channel = self.source_data.get("channel_name")
                pmt_gain = pmt_gains.get(selected_channel) if selected_channel else None

                two_photon_series_name = f"TwoPhotonSeries{ChannelName}"
                two_photon_series_metadata = {
                    "name": two_photon_series_name,
                    "imaging_plane": imaging_plane_name,
                    "field_of_view": [width_um * 1e-6, height_um * 1e-6],  # Convert to meters
                    "pmt_gain": pmt_gain,
                    "scan_line_rate": frame_rate * float(lsm.get("pixelY", 0)),
                    "unit": "n.a."
                }
                

                metadata["Ophys"]["TwoPhotonSeries"] = [two_photon_series_metadata]

        return metadata
