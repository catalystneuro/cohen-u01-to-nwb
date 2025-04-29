from typing import Optional, Tuple, Dict, List, Union
from pathlib import Path
from copy import deepcopy

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import FolderPathType, DeepDict, FilePathType, get_json_schema_from_method_signature
from neuroconv.nwbconverter import ConverterPipe
from pydantic import FilePath, validate_call
from pynwb import NWBFile
from pynwb.ophys import PlaneSegmentation, ImageSegmentation, RoiResponseSeries, ImagingPlane, OpticalChannel
from pynwb.device import Device
from pymatreader import read_mat
from pynwb.ophys import DfOverF
from pathlib import Path
from typing import Optional, Tuple
import warnings
from warnings import warn
import numpy as np

from roiextractors.extraction_tools import PathType, DtypeType, get_package
from roiextractors import ImagingExtractor
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface


class ScanImageConverter(ConverterPipe):
    """
    Primary conversion class for handling multiple ScanImage channels.

    This converter automatically detects and adds interfaces for all available channels
    in a ScanImage file.
    """

    display_name = "ScanImage Converter"
    associated_suffixes = (".tiff", ".tif", ".TIFF", ".TIF")
    info = "Converter for multi-channel ScanImage imaging data."

    @classmethod
    def get_source_schema(cls) -> dict:
        """
        Get the schema for the source arguments.

        Returns
        -------
        dict
            The schema dictionary containing input parameters and descriptions
            for initializing the ScanImage converter.
        """
        source_schema = get_json_schema_from_method_signature(method=cls.__init__, exclude=["channels"])
        source_schema["properties"]["file_path"]["description"] = "Path to the ScanImage TIFF file."
        return source_schema

    @classmethod
    def get_available_channel_names(cls, file_path: FilePath) -> list[str]:
        """
        Get the channel names available in the ScanImage file.

        Parameters
        ----------
        file_path : FilePath
            Path to the ScanImage TIFF file.

        Returns
        -------
        list of str
            The names of all available channels in the file.
        """
        return ScanImageImagingExtractor.get_available_channel_names(file_path)

    @validate_call
    def __init__(
        self,
        file_path: FilePath,
        file_paths: Optional[list[FilePath]] = None,
        photon_series_type: str = "TwoPhotonSeries",
        verbose: bool = False,
    ):
        """
        Initialize the ScanImage converter with interfaces for all available channels.

        Parameters
        ----------
        file_path : FilePath
            Path to the ScanImage TIFF file.
        file_paths : list of FilePath, optional
            List of paths to ScanImage TIFF files. If provided, this overrides the automatic file detection.
        photon_series_type : str, default: "TwoPhotonSeries"
            Type of photon series to use. Either "TwoPhotonSeries" or "OnePhotonSeries".
        verbose : bool, default: False
            Whether to print verbose output.
        """
        file_path = Path(file_path)

        # Get available channels
        channel_names = self.get_available_channel_names(file_path=file_path)
        self.photon_series_type = photon_series_type
        if not channel_names:
            raise ValueError(f"No channels found in ScanImage file: {file_path}")

        # Create interfaces for each channel
        data_interfaces = {}
        self.interface_channels = {}
        for channel_name in channel_names:
            interface_name = f"ScanImage{channel_name.replace(' ', '').capitalize()}ImagingInterface"
            interface = ScanImageImagingInterface(
                file_path=file_path,
                channel_name=channel_name,
                file_paths=file_paths,
            )

            self.interface_channels[channel_name] = interface
            # Set photon series type
            interface.photon_series_type = photon_series_type
            data_interfaces[interface_name] = interface

        # Initialize the ConverterPipe with the data interfaces
        super().__init__(data_interfaces=data_interfaces, verbose=verbose)

    def get_conversion_options_schema(self) -> dict:
        """
        Get the schema for the conversion options.

        Returns
        -------
        dict
            The schema dictionary containing conversion options for all interfaces.
        """
        conversion_options_schema = super().get_conversion_options_schema()

        return conversion_options_schema

    def add_to_nwbfile(self, nwbfile, metadata=None, conversion_options=None):

        # All of this is a trick because photon_series index is not availalbe at initilization
        # So we need to find out where the photon series is on the metadata list and then 
        # extract its index to pass it as a conversion option 
        metadata = metadata or self.get_metadata()
        photon_series_metadata: list = metadata["Ophys"][self.photon_series_type]
        interface_name_to_photon_series_index = {}
        for interface_name, interface in self.data_interface_objects.items():
            channel_name = interface.channel_name
            for photon_series_index, series_metadata in enumerate(photon_series_metadata):
                if channel_name in series_metadata["description"]:
                    # If the channel name is found in the description, use that photon_series_index
                    break
            interface_name_to_photon_series_index[interface_name] = photon_series_index

        conversion_options = conversion_options or {key: {} for key in self.data_interface_objects.keys()}
        for interface_name, interface in self.data_interface_objects.items():
            photon_series_index = interface_name_to_photon_series_index[interface_name]
            conversion_options[interface_name]["photon_series_index"] = photon_series_index

        return super().add_to_nwbfile(nwbfile, metadata, conversion_options)


class ScanImageImagingExtractor(ImagingExtractor):
    """
    Specialized extractor for reading TIFF files produced via ScanImage software.

    This extractor is designed to handle the structure of ScanImage TIFF files, which can contain
    multi channel and both planar and volumetric data. It also supports both single-file and multi-file datasets generated
    by ScanImage in various acquisition modes (grab, focus, loop).

    The extractor creates a mapping between each frame in the dataset and its corresponding physical file
    and IFD (Image File Directory) location. This mapping enables efficient retrieval of specific frames
    without loading the entire dataset into memory, making it suitable for large datasets.

    For datasets with multiple frames per slice, a slice_sample parameter must be provided.


    Key features:
    - Handles multi-channel data with channel selection
    - Supports volumetric (multi-plane) imaging data
    - Automatically detects and loads multi-file datasets based on ScanImage naming conventions
    - Extracts and provides access to ScanImage metadata
    - Efficiently retrieves frames using lazy loading
    - Handles flyback frames in volumetric data by ignoring them in the mapping

    """

    extractor_name = "ScanImageImagingExtractor"

    def __init__(
        self,
        file_path: Optional[PathType] = None,
        channel_name: Optional[str] = None,
        file_paths: Optional[list[PathType]] = None,
        slice_sample: Optional[int] = None,
    ):
        """
        Initialize the extractor.

        Parameters
        ----------
        file_path : PathType
            Path to the TIFF file. If this is part of a multi-file series, this should be the first file.
        channel_name : str, optional
            Name of the channel to extract. If None and multiple channels are available, the first channel will be used.
            Check available channels with `get_available_channel_names`.
        file_paths : list[PathType], optional
            list of file paths to use. If provided, this overrides the automatic
            file detection heuristics. Use this if automatic detection does not work correctly and you know
            exactly which files should be included.  The file paths should be provided in an order that
            reflects the temporal order of the frames in the dataset.
        slice_sample : int, optional
            When frames_per_slice > 1 (multiple frames per slice), this parameter specifies which frame to use
            for each slice. Must be between 0 and frames_per_slice-1. If None and frames_per_slice > 1,
            a ValueError will be raised.
        """
        super().__init__()
        self.file_path = file_paths[0] if file_paths is not None else file_path
        assert self.file_path is not None, "file_path or file_paths must be provided"

        # Validate file suffix
        valid_suffixes = [".tiff", ".tif", ".TIFF", ".TIF"]
        if self.file_path.suffix not in valid_suffixes:
            suffix_string = ", ".join(valid_suffixes[:-1]) + f", or {valid_suffixes[-1]}"
            warn(
                f"Suffix ({self.file_path.suffix}) is not of type {suffix_string}! "
                f"The {self.extractor_name} Extractor may not be appropriate for the file."
            )

        # Open the
        tifffile = get_package(package_name="tifffile")
        tiff_reader = tifffile.TiffReader(self.file_path)

        self._general_metadata = tiff_reader.scanimage_metadata
        self._metadata = self._general_metadata["FrameData"]

        self._num_rows, self._num_columns = tiff_reader.pages[0].shape
        self._dtype = tiff_reader.pages[0].dtype

        # This criteria was confirmed by Lawrence Niu, a developer of ScanImage
        self.is_volumetric = self._metadata["SI.hStackManager.enable"]
        if self.is_volumetric:
            self._sampling_frequency = self._metadata["SI.hRoiManager.scanVolumeRate"]
            self._num_planes = self._metadata["SI.hStackManager.numSlices"]

            self._frames_per_slice = self._metadata["SI.hStackManager.framesPerSlice"]
            if self._frames_per_slice > 1:
                if slice_sample is None:
                    error_msg = (
                        f"Multiple frames per slice detected (frames_per_slice = {self._frames_per_slice}). "
                        f"Please specify a slice_sample between 0 and {self._frames_per_slice - 1} to select which frame to use for each slice."
                    )
                    raise ValueError(error_msg)

                if not (0 <= slice_sample < self._frames_per_slice):
                    error_msg = f"slice_sample must be between 0 and {self._frames_per_slice - 1} (frames_per_slice - 1), but got {slice_sample}."
                    raise ValueError(error_msg)

                self._slice_sample = slice_sample
            else:
                self._slice_sample = None

            self._frames_per_volume_per_channel = self._metadata["SI.hStackManager.numFramesPerVolume"]
            self._frames_per_volume_with_flyback = self._metadata["SI.hStackManager.numFramesPerVolumeWithFlyback"]

            self.num_flyback_frames = self._frames_per_volume_with_flyback - self._frames_per_volume_per_channel
        else:
            self._sampling_frequency = self._metadata["SI.hRoiManager.scanFrameRate"]
            self._num_planes = 1
            self._frames_per_slice = 1
            self.num_flyback_frames = 0

        # This piece of the metadata is the indication that the channel is saved on the data
        channels_available = self._metadata["SI.hChannels.channelSave"]
        channels_available = [channels_available] if isinstance(channels_available, int) else channels_available
        self._num_channels = len(channels_available)

        # Determine their name and use matlab 1-indexing
        all_channel_names = self._metadata["SI.hChannels.channelName"]
        self.channel_names = [all_channel_names[channel_index - 1] for channel_index in channels_available]

        # Channel selection checks
        self._is_multi_channel_data = len(self.channel_names) > 1
        if self._is_multi_channel_data and channel_name is None:

            error_msg = (
                f"Multiple channels available in the data {self.channel_names}"
                "Please specify a channel name to extract data from."
            )
            raise ValueError(error_msg)
        elif self._is_multi_channel_data and channel_name is not None:
            if channel_name not in self.channel_names:
                error_msg = (
                    f"Channel name ({channel_name}) not found in available channels ({self.channel_names}). "
                    "Please specify a valid channel name."
                )
                raise ValueError(error_msg)

            self.channel_name = channel_name
            self._channel_index = self.channel_names.index(channel_name)
        else:  # single channel data

            self.channel_name = self.channel_names[0]
            self._channel_index = 0

        # Check if this is a multi-file dataset
        if file_paths is None:
            self.file_paths = self._find_data_files()
        else:
            self.file_paths = file_paths

        # Open all TIFF files and store only file readers for lazy loading
        total_ifds = 0
        self._tiff_readers = []
        for file_path in self.file_paths:
            try:
                tiff_reader = tifffile.TiffFile(file_path)
                self._tiff_readers.append(tiff_reader)
                total_ifds += len(tiff_reader.pages)
            except Exception as e:
                for tiff_reader in self._tiff_readers:
                    tiff_reader.close()
                raise RuntimeError(f"Error opening TIFF file {file_path}: {e}")

        # Calculate total IFDs and samples
        self._ifds_per_file = [len(tiff_reader.pages) for tiff_reader in self._tiff_readers]

        # Note that this includes all the frames for all the channels including flyback frames
        self._num_frames_in_dataset = sum(self._ifds_per_file)

        image_frames_per_cycle = self._num_planes * self._num_channels * self._frames_per_slice
        total_frames_per_cycle = image_frames_per_cycle + self.num_flyback_frames * self._num_channels

        # Note that the acquisition might end without completing the last cycle and we discard those frames
        num_acquisition_cycles = self._num_frames_in_dataset // (total_frames_per_cycle)

        #  Every cycle is a full channel sample either volume or planar
        self._num_samples = num_acquisition_cycles

        # Map IFDs and files to frames, channel, depth, and acquisition cycle
        full_frames_to_ifds_table = self._create_frame_to_ifd_table(
            num_channels=self._num_channels,
            num_planes=self._num_planes,
            num_acquisition_cycles=num_acquisition_cycles,
            num_frames_per_slice=self._frames_per_slice,
            num_flyback_frames=self.num_flyback_frames,
            ifds_per_file=self._ifds_per_file,
        )

        # Filter mapping for the specified channel
        channel_mask = full_frames_to_ifds_table["channel_index"] == self._channel_index
        channel_frames_to_ifd_table = full_frames_to_ifds_table[channel_mask]

        self._frames_to_ifd_table = channel_frames_to_ifd_table

        # Filter mapping for the specified slice_sample
        if self.is_volumetric and self._slice_sample is not None:
            slice_sample_mask = channel_frames_to_ifd_table["slice_sample_index"] == self._slice_sample
            self._frames_to_ifd_table = channel_frames_to_ifd_table[slice_sample_mask]

    @staticmethod
    def _create_frame_to_ifd_table(
        num_channels: int,
        num_planes: int,
        num_acquisition_cycles: int,
        ifds_per_file: list[int],
        num_frames_per_slice: int = 1,
        num_flyback_frames: int = 0,
    ) -> np.ndarray:
        """
        Create a table that describes the data layout of the dataset.

        Every row in the table corresponds to a frame in the dataset and contains:
        - file_index: The index of the file in the series
        - IFD_index: The index of the IFD in the file
        - channel_index: The index of the channel
        - depth_index: The index of the depth
        - acquisition_cycle_index: The index of the time

        The table is represented as a structured numpy array that maps each combination of time,
        channel, and depth to its corresponding physical location in the TIFF files.

        Parameters
        ----------
        num_channels : int
            Number of channels.
        num_planes: int
            The number of planes which corresponds to the depth index or the number of frames per volume
            per channel.
        num_acquisition_cycles : int
            Number of acquisition cycles. For ScanImage, this is the number of samples.
        ifds_per_file : list[int]
            Number of IFDs in each file.
        num_frames_per_slice : int
            Number of frames per slice. This is used to determine the slice_sample index.
        num_flyback_frames : int
            Number of flyback frames.

        Returns
        -------
        np.ndarray
            A structured array mapping all combinations of time, channel, and depth to file
            and IFD indices.
        """
        # Create structured dtype for the table
        mapping_dtype = np.dtype(
            [
                ("file_index", np.uint16),
                ("IFD_index", np.uint16),
                ("channel_index", np.uint8),
                ("depth_index", np.uint8),
                ("slice_sample_index", np.uint8),
                ("acquisition_cycle_index", np.uint16),
            ]
        )

        # Calculate total number of entries
        image_frames_per_cycle = num_planes * num_frames_per_slice * num_channels
        total_frames_per_cycle = image_frames_per_cycle + num_flyback_frames * num_channels

        # Generate global ifd indices for complete cycles only
        # This ensures we only include frames from complete acquisition cycles
        num_frames_in_complete_cycles = num_acquisition_cycles * total_frames_per_cycle
        global_ifd_indices = np.arange(num_frames_in_complete_cycles, dtype=np.uint32)

        # We need to filter out the flyback frames, we create an index within each acquisition cycle
        # And then filter out the non-image frames (flyback frames)
        index_in_acquisition_cycle = global_ifd_indices % total_frames_per_cycle
        is_imaging_frame = index_in_acquisition_cycle < image_frames_per_cycle

        global_ifd_indices = global_ifd_indices[is_imaging_frame]
        index_in_acquisition_cycle = index_in_acquisition_cycle[is_imaging_frame]

        # To find their file index we need file boundaries
        file_boundaries = np.zeros(len(ifds_per_file) + 1, dtype=np.uint32)
        file_boundaries[1:] = np.cumsum(ifds_per_file)

        # Find which file each global index belongs to
        file_indices = np.searchsorted(file_boundaries, global_ifd_indices, side="right") - 1

        # Now, we offset the global IFD indices by the starting position of the file
        # to get local IFD indices that start at 0 for each file
        ifd_indices = global_ifd_indices - file_boundaries[file_indices]

        # Calculate indices for each dimension based on the frame position within the cycle
        # For ScanImage, the order is always CZT which means that the channel index comes first,
        # followed by the frames per slice, then depth and finally the acquisition cycle
        channel_indices = index_in_acquisition_cycle % num_channels
        slice_sample_indices = (index_in_acquisition_cycle // num_channels) % num_frames_per_slice
        depth_indices = (global_ifd_indices // (num_channels * num_frames_per_slice)) % num_planes
        acquisition_cycle_indices = global_ifd_indices // total_frames_per_cycle

        # Create the structured array with the correct size (number of imaging frames after filtering)
        mapping = np.zeros(len(global_ifd_indices), dtype=mapping_dtype)
        mapping["file_index"] = file_indices
        mapping["IFD_index"] = ifd_indices
        mapping["channel_index"] = channel_indices
        mapping["slice_sample_index"] = slice_sample_indices
        mapping["depth_index"] = depth_indices
        mapping["acquisition_cycle_index"] = acquisition_cycle_indices

        return mapping

    def _find_data_files(self) -> list[PathType]:
        """Find additional files in the series based on the file naming pattern.

        This method determines which files to include in the dataset using one of these approaches:

        1. If file_paths is provided: Uses the provided list of file paths directly
        2. If file_pattern is provided: Uses the provided pattern to glob for files
        3. Otherwise, analyzes the file name and ScanImage metadata to determine if the current file
            is part of a multi-file dataset. It uses different strategies based on the acquisition mode:
            - For 'grab' mode with finite frames per file: Uses base_name_acquisition_* pattern
            - For 'loop' mode: Uses base_name_* pattern
            - For 'slow' stack mode with volumetric data: Uses base_name_* pattern
            - Otherwise: Returns only the current file

        This information about ScanImage file naming was shared in a private conversation with
        Lawrence Niu, who is a developer of ScanImage.

        Returns
        -------
        list[PathType]
            list of paths to all files in the series, sorted naturally (e.g., file_1, file_2, file_10)
        """
        # Parse the file name to extract base name, acquisition number, and file index
        file_stem = self.file_path.stem

        # Can be grab, focus or loop, see
        # https://docs.scanimage.org/Basic+Features/Acquisitions.html
        acquisition_state = self._metadata["SI.acqState"]
        frames_per_file = self._metadata["SI.hScan2D.logFramesPerFile"]
        stack_mode = self._metadata["SI.hStackManager.stackMode"]

        # This is the happy path that is well specified in the documentation
        if acquisition_state == "grab" and frames_per_file != float("inf"):
            name_parts = file_stem.split("_")
            base_name, acquisition, file_index = "_".join(name_parts[:-2]), name_parts[-2], name_parts[-1]
            pattern = f"{base_name}_{acquisition}_*{self.file_path.suffix}"
        # Looped acquisitions also divides the files according to Lawrence Niu in private conversation
        elif acquisition_state == "loop":  # This also separates the files
            base_name = "_".join(file_stem.split("_")[:-1])  # Everything before the last _
            pattern = f"{base_name}_*{self.file_path.suffix}"
        # This also divided the files according to Lawrence Niu in private conversation
        elif stack_mode == "slow" and self.is_volumetric:
            base_name = "_".join(file_stem.split("_")[:-1])  # Everything before the last _
            pattern = f"{base_name}_*{self.file_path.suffix}"
        else:
            files_found = [self.file_path]
            return files_found

        from natsort import natsorted

        files_found = natsorted(self.file_path.parent.glob(pattern))
        return files_found

    def get_series(self, start_sample: Optional[int] = None, end_sample: Optional[int] = None) -> np.ndarray:
        """
        Get data as a time series from start_sample to end_sample.

        This method retrieves frames at the specified range from the ScanImage TIFF file(s).
        It uses the mapping created during initialization to efficiently locate and load only
        the requested frames, without loading the entire dataset into memory.

        For volumetric data (multiple planes), the returned array will have an additional dimension
        for the planes. For planar data (single plane), the plane dimension is squeezed out.

        Parameters
        ----------
        start_sample : int
        end_sample : int

        Returns
        -------
        numpy.ndarray
            Array of data with shape (num_samples, height, width) if num_planes is 1,
            or (num_samples, height, width, num_planes) if num_planes > 1.

            For example, for a non-volumetric dataset with 512x512 frames, requesting 3 samples
            would return an array with shape (3, 512, 512).

            For a volumetric dataset with 5 planes and 512x512 frames, requesting 3 samples
            would return an array with shape (3, 512, 512, 5).
        """
        start_sample = int(start_sample) if start_sample is not None else 0
        end_sample = int(end_sample) if end_sample is not None else self.get_num_samples()

        samples_in_series = end_sample - start_sample

        # Preallocate output array as volumetric and squeeze if not volumetric before returning
        num_rows, num_columns, num_planes = self.get_volume_shape()
        dtype = self.get_dtype()
        samples = np.empty((samples_in_series, num_rows, num_columns, num_planes), dtype=dtype)

        for return_index, sample_index in enumerate(range(start_sample, end_sample)):
            for depth_position in range(num_planes):

                # Calculate the index in the mapping table array
                frame_index = sample_index * num_planes + depth_position
                table_row = self._frames_to_ifd_table[frame_index]
                file_index = table_row["file_index"]
                ifd_index = table_row["IFD_index"]

                tiff_reader = self._tiff_readers[file_index]
                image_file_directory = tiff_reader.pages[ifd_index]
                samples[return_index, :, :, depth_position] = image_file_directory.asarray()

        # Squeeze the depth dimension if not volumetric
        if not self.is_volumetric:
            samples = samples.squeeze(axis=3)

        return samples

    def get_image_shape(self) -> Tuple[int, int]:
        """Get the shape of the video frame (num_rows, num_columns).

        Returns
        -------
        tuple
            Shape of the video frame (num_rows, num_columns).
        """
        return (self._num_rows, self._num_columns)

    def get_frame_shape(self) -> Tuple[int, int]:
        """Get the shape of a single frame (num_rows, num_columns).

        Returns
        -------
        tuple
            Shape of a single frame (num_rows, num_columns).
        """
        return (self._num_rows, self._num_columns)

    def get_sample_shape(self):
        """
        Get the shape of a sample.

        Returns
        -------
        tuple of int
            Shape of a single sample. If the data is volumetric, the shape is hape of a single sample (num_rows, num_columns).
            (num_rows, num_columns, num_planes). Otherwise, the shape is
            (num_rows, num_columns).
        """
        if self.is_volumetric:
            return (self._num_rows, self._num_columns, self._num_planes)
        else:
            return (self._num_rows, self._num_columns)

    def get_volume_shape(self) -> Tuple[int, int, int]:
        """Get the shape of a single volume (num_rows, num_columns, num_planes).

        Returns
        -------
        tuple
            Shape of a single volume (num_rows, num_columns, num_planes).
        """
        return (self._num_rows, self._num_columns, self._num_planes)

    def get_num_samples(self) -> int:
        """Get the number of samples in the video.

        Returns
        -------
        int
            Number of samples in the video.
        """
        return self._num_samples

    def get_sampling_frequency(self) -> float:
        """Get the sampling frequency in Hz.

        Returns
        -------
        float
            Sampling frequency in Hz.
        """
        return self._sampling_frequency

    def get_channel_names(self):
        return self.channel_names

    @staticmethod
    def get_available_channel_names(file_path: PathType) -> list:
        """Get the channel names available in a ScanImage TIFF file.

        This static method extracts the channel names from a ScanImage TIFF file
        without needing to create an extractor instance. This is useful for
        determining which channels are available before creating an extractor.

        Parameters
        ----------
        file_path : PathType
            Path to the ScanImage TIFF file.

        Returns
        -------
        list
            list of channel names available in the file.

        Examples
        --------
        >>> channel_names = ScanImageImagingExtractor.get_channel_names('path/to/file.tif')
        >>> print(f"Available channels: {channel_names}")
        """
        from tifffile import read_scanimage_metadata

        with open(file_path, "rb") as fh:
            all_metadata = read_scanimage_metadata(fh)
            non_varying_frame_metadata = all_metadata[0]

        # `channelSave` indicates whether the channel is saved
        # We check `channelSave` first but keep the `channelsActive` check for backward compatibility
        channel_availability_keys = ["SI.hChannels.channelSave", "SI.hChannels.channelsActive"]
        for channel_availability in channel_availability_keys:
            if channel_availability in non_varying_frame_metadata.keys():
                break

        available_channels = non_varying_frame_metadata[channel_availability]
        available_channels = [available_channels] if not isinstance(available_channels, list) else available_channels
        channel_indices = np.array(available_channels) - 1  # Account for MATLAB indexing
        channel_names = non_varying_frame_metadata["SI.hChannels.channelName"]
        channel_names_available = [channel_names[i] for i in channel_indices]

        return channel_names_available

    def get_dtype(self) -> DtypeType:
        """Get the data type of the video.

        Returns
        -------
        dtype
            Data type of the video.
        """
        return self._dtype

    def get_times(self) -> np.ndarray:
        """Get the timestamps for each frame.

        Returns
        -------
        numpy.ndarray
            Array of timestamps in seconds for each frame.

        Notes
        -----
        This method extracts timestamps from the ScanImage TIFF file(s) for the selected channel.
        It uses the mapping created during initialization to efficiently locate and extract
        timestamps for each sample.
        """
        if self._times is not None:
            return self._times

        # Initialize array to store timestamps
        num_samples = self.get_num_samples()
        num_planes = self.get_num_planes()
        timestamps = np.zeros(num_samples, dtype=np.float64)

        # For each sample, extract its timestamp from the corresponding file and IFD
        for sample_index in range(num_samples):

            # Get the last frame in this sample to get the timestamps
            frame_index = sample_index * num_planes + (num_planes - 1)
            table_row = self._frames_to_ifd_table[frame_index]
            file_index = table_row["file_index"]
            ifd_index = table_row["IFD_index"]

            tiff_reader = self._tiff_readers[file_index]
            image_file_directory = tiff_reader.pages[ifd_index]

            # Extract timestamp from the IFD description
            description = image_file_directory.description
            description_lines = description.split("\n")

            # Use iterator pattern to find frameTimestamps_sec
            timestamp_line = next((line for line in description_lines if "frameTimestamps_sec" in line), None)

            if timestamp_line is not None:
                # Extract the value part after " = "
                _, value_str = timestamp_line.split(" = ", 1)
                try:
                    timestamps[sample_index] = float(value_str.strip())
                except ValueError:
                    # If parsing fails, use sample index / sampling frequency as fallback
                    timestamps[sample_index] = sample_index / self._sampling_frequency
            else:
                # If no timestamp found, throw a warning and use sample index / sampling frequency as fallback
                warnings.warn(
                    f"No frameTimestamps_sec found for sample {sample_index}. Using calculated timestamp instead.",
                    UserWarning,
                )
                timestamps[sample_index] = sample_index / self._sampling_frequency

        # Cache the timestamps
        self._times = timestamps
        return timestamps

    def get_num_planes(self) -> int:
        """Get the number of depth planes.

        For volumetric data, this returns the number of Z-planes in each volume.
        For planar data, this returns 1.

        Returns
        -------
        int
            Number of depth planes.
        """
        return self._num_planes

    def get_plane_extractor(self, plane_index: int) -> ImagingExtractor:
        """Extract a specific depth plane from volumetric data.

        This method allows for extracting a specific depth plane from volumetric imaging data,
        returning a modified version of the extractor that only returns data for the specified plane.

        Parameters
        ----------
        plane_index: int
            Index of the depth plane to extract (0-indexed).

        Returns
        -------
        extractor: ImagingExtractor
            A modified version of the extractor that only returns data for the specified plane.

        Raises
        ------
        ValueError
            If the data is not volumetric (has only one plane).
            If plane_index is out of range.

        Examples
        --------
        >>> extractor = ScanImageImagingExtractor('path/to/volumetric_file.tif')
        >>> # Get only the first plane
        >>> first_plane = extractor.depth_slice(plane_index=0)
        >>> # Get the second plane
        >>> second_plane = extractor.depth_slice(plane_index=1)
        """
        if not self.is_volumetric:
            raise ValueError("Cannot depth slice non-volumetric data. This data has only one plane.")

        # Validate parameters
        if plane_index < 0 or plane_index >= self._num_planes:
            raise ValueError(f"plane_index ({plane_index}) must be between 0 and {self._num_planes - 1}")

        # Create a copy of the current extractor
        import copy

        sliced_extractor = copy.deepcopy(self)

        # Filter the frames_to_ifd_table to only include entries for the specified depth plane
        depth_mask = sliced_extractor._frames_to_ifd_table["depth_index"] == plane_index
        sliced_extractor._frames_to_ifd_table = sliced_extractor._frames_to_ifd_table[depth_mask]

        # Update the number of samples
        sliced_extractor._num_samples_per_channel = len(sliced_extractor._frames_to_ifd_table)

        # Override the is_volumetric flag and num_planes
        sliced_extractor.is_volumetric = False
        sliced_extractor._num_planes = 1

        return sliced_extractor

    @staticmethod
    def get_frames_per_slice(file_path: PathType) -> int:
        """
        Get the number of frames per slice from a ScanImage TIFF file.

        ScanImage can sample mutiple frames per each slice.

        Parameters
        ----------
        file_path : PathType
            Path to the ScanImage TIFF file.

        Returns
        -------
        int
            Number of frames per slice.

        """
        from tifffile import read_scanimage_metadata

        with open(file_path, "rb") as fh:
            all_metadata = read_scanimage_metadata(fh)
            non_varying_frame_metadata = all_metadata[0]

        frames_per_slice = non_varying_frame_metadata.get("SI.hStackManager.framesPerSlice", 1)
        return frames_per_slice

    def __del__(self):
        """Close file handles when the extractor is garbage collected."""
        if hasattr(self, "_tiff_readers"):
            for handle in self._tiff_readers:
                try:
                    handle.close()
                except Exception as e:
                    warnings.warn(f"Error closing TIFF file handle {handle} with error: {e}", UserWarning)
                    pass

class ScanImageImagingInterface(BaseImagingExtractorInterface):

    extractor = ScanImageImagingExtractor

    def __init__(
        self,
        file_path: FilePathType,
        channel_name: Optional[str] = None,
        file_paths: Optional[list[str]] = None,
    ):
        """
        Initialize the ScanImage Imaging Interface.

        Parameters
        ----------
        file_path : PathType
            Path to the TIFF file. If this is part of a multi-file series, this should be the first file.
        channel_name : str, optional
            Name of the channel to extract. If None and multiple channels are available, the first channel will be used.
        file_paths : List[PathType], optional
            List of file paths to use. If provided, this overrides the automatic
            file detection heuristics. Use this if automatic detection does not work correctly and you know
            exactly which files should be included.  The file paths should be provided in an order that
            reflects the temporal order of the frames in the dataset.
        """

        self.channel_name = channel_name
        super().__init__(file_path=file_path, channel_name=channel_name, file_paths=file_paths)

    def get_metadata(self) -> DeepDict:
        """
        Get metadata for the ScanImage imaging data.

        Returns
        -------
        DeepDict
            The metadata dictionary containing imaging metadata from the ScanImage files.
        """
        metadata = super().get_metadata()

        session_start_time = self._get_session_start_time()
        if session_start_time:
            metadata["NWBFile"]["session_start_time"] = session_start_time

        # Extract ScanImage-specific metadata
        if hasattr(self.imaging_extractor, "_general_metadata"):
            # Add general metadata to a custom field
            scanimage_metadata = self.imaging_extractor._general_metadata

            # Update device information
            device_name = "Microscope"
            metadata["Ophys"]["Device"][0].update(name=device_name, description=f"Microscope controlled by ScanImage")

            # Update imaging plane metadata
            imaging_plane_metadata = metadata["Ophys"]["ImagingPlane"][0]
            imaging_plane_metadata.update(
                device=device_name,
                imaging_rate=self.imaging_extractor.get_sampling_frequency(),
                description="Imaging plane from ScanImage acquisition",
            )

            # Update photon series metadata
            photon_series_key = self.photon_series_type  # "TwoPhotonSeries" or "OnePhotonSeries"
            photon_series_metadata = metadata["Ophys"][photon_series_key][0]

            channel_string = self.channel_name.replace(" ", "").capitalize()
            photon_series_name = f"{photon_series_key}{channel_string}"

            photon_series_metadata["name"] = photon_series_name
            photon_series_metadata["description"] = f"Imaging data acquired using ScanImage for {self.channel_name}"

            # Add additional metadata if available
            if "FrameData" in scanimage_metadata:
                frame_data = scanimage_metadata["FrameData"]

                # Calculate scan line rate from line period if available
                if "SI.hRoiManager.linePeriod" in frame_data:
                    scan_line_rate = 1 / float(frame_data["SI.hRoiManager.linePeriod"])
                    photon_series_metadata.update(scan_line_rate=scan_line_rate)
                elif "SI.hScan2D.scannerFrequency" in frame_data:
                    photon_series_metadata.update(
                        scan_line_rate=frame_data["SI.hScan2D.scannerFrequency"]
                    )

                # Add version information to device description if available
                if "SI.VERSION_MAJOR" in frame_data:
                    version = f"{frame_data.get('SI.VERSION_MAJOR', '')}.{frame_data.get('SI.VERSION_MINOR', '')}.{frame_data.get('SI.VERSION_UPDATE', '')}"
                    metadata["Ophys"]["Device"][0][
                        "description"
                    ] = f"Microscope and acquisition data with ScanImage (version {version})"

            # Extract ROI metadata if available
            if "RoiGroups" in scanimage_metadata:
                roi_metadata = scanimage_metadata["RoiGroups"]

                # Extract grid spacing and origin coordinates from scanfields
                grid_spacing = None
                grid_spacing_unit = "n.a"
                origin_coords = None
                origin_coords_unit = "n.a"

                if "imagingRoiGroup" in roi_metadata and "rois" in roi_metadata["imagingRoiGroup"]:
                    rois = roi_metadata["imagingRoiGroup"]["rois"]
                    if isinstance(rois, dict) and "scanfields" in rois:
                        scanfields = rois["scanfields"]
                        if "sizeXY" in scanfields and "pixelResolutionXY" in scanfields:
                            fov_size_in_um = np.array(scanfields["sizeXY"])
                            frame_dimension = np.array(scanfields["pixelResolutionXY"])
                            grid_spacing = fov_size_in_um / frame_dimension
                            grid_spacing_unit = "micrometers"

                        if "centerXY" in scanfields:
                            origin_coords = scanfields["centerXY"]
                            origin_coords_unit = "micrometers"

                # Update imaging plane metadata with grid spacing and origin coordinates
                if grid_spacing is not None:
                    imaging_plane_metadata.update(
                        grid_spacing=grid_spacing.tolist(), grid_spacing_unit=grid_spacing_unit
                    )

                if origin_coords is not None:
                    imaging_plane_metadata.update(origin_coords=origin_coords, origin_coords_unit=origin_coords_unit)

        return metadata

    def _get_session_start_time(self):
            """
            Open a ScanImage TIFF file, read the first frame, extract and parse the 'epoch' metadata
            as the session start time.

            Parameters
            ----------
            tiff_path : str or Path
                Path to the TIFF file.

            Returns
            -------
            datetime
                Parsed datetime from the 'epoch' metadata.

            Raises
            ------
            ValueError
                If 'epoch' metadata is not found or is malformed.
            """
            
            
            from tifffile import TiffReader
            tiff_file_path = self.imaging_extractor.file_path
            with TiffReader(tiff_file_path) as tif:
                image_description = tif.pages[0].tags["ImageDescription"].value

            import re 
            match = re.search(r'epoch\s*=\s*\[([^\]]+)\]', image_description)
            if not match:
                raise ValueError(f"'epoch' field not found in {tiff_file_path}")

            epoch_values = match.group(1).split()
            import warnings
            if len(epoch_values) != 6:
                warnings.warn(
                    f"Expected 6 values in 'epoch' field, found {len(epoch_values)}: \n"
                    f"Epoch field {epoch_values}."
                )
                return None

            year, month, day, hour, minute, seconds = map(float, epoch_values)
            second_int = int(seconds)
            microsecond = int((seconds - second_int) * 1e6)
            
            import datetime
            return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), second_int, microsecond)

class KimLabROIInterface(BaseDataInterface):
    """Data interface for Kim Lab ROI segments data."""

    def __init__(
        self,
        file_path: FilePathType,
        roi_info_file_path: FilePathType,
        timestamps: np.ndarray,
        image_shape: tuple[int, int] = None,
        verbose: bool = False,
    ):
        """Initialize the ROI interface.

        Parameters
        ----------
        file_path : FilePathType
            Path to the df_f.mat file containing ROI segments data.
        roi_info_file_path : FilePathType
            Path to the ROI_info.mat file containing ROI information.
        timestamps : np.ndarray
            Array of timestamps in seconds. Those are used to synch to other data streams.
        image_shape : Optional[tuple[int, int]]
            Shape of the image, can be provided if a PhotonSeries is not found in acquisition
        verbose : bool, default: False
            Whether to print progress information
        """
        super().__init__(file_path=file_path, roi_info_file_path=roi_info_file_path, verbose=verbose)
        self.file_path = Path(file_path)
        self.roi_info_file_path = Path(roi_info_file_path)

        # Validate files exist
        if not self.file_path.is_file():
            raise FileNotFoundError(f"df_f.mat file not found at {self.file_path}")
        if not self.roi_info_file_path.is_file():
            raise FileNotFoundError(f"ROI_info.mat file not found at {self.roi_info_file_path}")
        self.timestamps = timestamps
        self.image_shape = image_shape
        self.verbose = verbose

    def get_metadata(self) -> dict:
        """Get metadata for the ROI segments.

        Returns
        -------
        dict
            The metadata dictionary
        """
        metadata = super().get_metadata()
        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[DeepDict] = None):
        """Add the ROI segments data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add the ROI segments to
        metadata : Optional[DeepDict], optional
            Metadata dictionary
        """

        # Create an ophys processing module if it doesn't exist
        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys", description="Contains optical physiology processed data"
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        # Create device
        if "Microscope" not in nwbfile.devices:
            device = Device(name="Microscope")
            nwbfile.add_device(device)
        else:
            device = nwbfile.devices["Microscope"]

        # Create optical channel with minimum required fields
        optical_channel = OpticalChannel(
            name="OpticalChannel", description="optical channel", emission_lambda=500.0  # TODO: Figure it out
        )

        # Create imaging plane with minimum required fields
        if "SegmentationPlane" not in nwbfile.imaging_planes:
            imaging_plane = nwbfile.create_imaging_plane(
                name="SegmentationPlane",
                optical_channel=optical_channel,
                description="segmentation plane",
                device=device,
                excitation_lambda=600.0,  # TODO: Figure it out
                indicator="unknown",
                location="unknown",
            )
        else:
            imaging_plane = nwbfile.imaging_planes["SegmentationPlane"]

        # Create image segmentation
        img_seg = ImageSegmentation(name="ImageSegmentation")
        ophys_module.add(img_seg)

        # Create plane segmentation
        plane_seg = PlaneSegmentation(
            name="PlaneSegmentation",
            description="Regions of interest from calcium imaging",
            imaging_plane=imaging_plane,
        )
        img_seg.add_plane_segmentation(plane_seg)

        # Load the data
        df_f_data = read_mat(self.file_path)["df_f"]
        roi_info_data = read_mat(self.roi_info_file_path)
        x_coordinates = roi_info_data["x_cor"]
        y_coordinates = roi_info_data["y_cor"]

        # Add ROI entries to the plane segmentation with placeholder image masks
        num_rois = df_f_data.shape[0]  # First dimension is ROIs

        from skimage.draw import polygon

        for roi_index in range(num_rois):
            # Create a small placeholder mask for each ROI
            image_mask = np.zeros(self.image_shape, dtype=bool)
            x = x_coordinates[roi_index]
            y = y_coordinates[roi_index]

            rr, cc = polygon(y, x)
            image_mask[rr, cc] = True

            plane_seg.add_roi(image_mask=image_mask)

        # Create ROI region reference
        roi_table_region = plane_seg.create_roi_table_region(description="all ROIs", region=list(range(num_rois)))
        data = df_f_data.T
        # Add ROI fluorescence data
        roi_response_series = RoiResponseSeries(
            name="RoIResponseSeries",
            data=data,
            rois=roi_table_region,
            unit="a.u.",
            timestamps=self.timestamps,
            description="Change in fluorescence normalized by baseline fluorescence (dF/F)",
        )

        df_over_f_container = DfOverF(roi_response_series=roi_response_series)

        # Add to the ophys processing module
        ophys_module.add(df_over_f_container)
