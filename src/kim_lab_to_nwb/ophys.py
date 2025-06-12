from typing import Optional
from pathlib import Path
from copy import deepcopy

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict, get_json_schema_from_method_signature
from neuroconv.nwbconverter import ConverterPipe
from pydantic import FilePath, validate_call
from pynwb import NWBFile
from pynwb.ophys import PlaneSegmentation, ImageSegmentation, RoiResponseSeries, ImagingPlane, OpticalChannel
from pynwb.device import Device
from pymatreader import read_mat
from pynwb.ophys import DfOverF
from pathlib import Path
from warnings import warn
import numpy as np

from roiextractors import ScanImageImagingExtractor 
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
        self.channel_names = self.get_available_channel_names(file_path=file_path)
        self.photon_series_type = photon_series_type
        if not self.channel_names:
            raise ValueError(f"No channels found in ScanImage file: {file_path}")

        # Create interfaces for each channel
        data_interfaces = {}
        self.interface_channels = {}
        for channel_name in self.channel_names:
            interface_name = f"ScanImage{channel_name.replace(' ', '').capitalize()}ImagingInterface"
            interface = KimScanImageImagingInterface(
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

        return super().add_to_nwbfile(nwbfile=nwbfile, metadata=metadata, conversion_options=conversion_options)


class KIMScanImageImagingExtractor(ScanImageImagingExtractor):
    """
    A wrapper around the ScanImageImagingExtractor for easier access to common functionality.
    and provide a method to get the frame indices that are selected for a specific channel.
    """
    
    def slice_to_valid_samples(self, valid_samples_mask: np.ndarray):
        """
        Slice the extractor to only include valid samples.
        
        This method encapsulates the slicing operation while preserving metadata.
        Typically used when the NIDAQ was turned off before the microscope stopped recording,
        leaving some imaging samples without corresponding sync pulses.
        
        Parameters
        ----------
        valid_samples_mask : np.ndarray
            Boolean mask indicating which samples are valid
            
        Returns
        -------
        KIMScanImageImagingExtractor
            A sliced extractor containing only valid samples
        """
        # Get the last valid sample index
        last_valid_sample_index = np.where(valid_samples_mask)[0][-1]
        
        # Slice the imaging extractor to only include samples with sync pulses
        sliced_extractor = self.slice_samples(
            start_sample=0,
            end_sample=last_valid_sample_index + 1  # +1 because end_sample is exclusive
        )
        
        # Preserve file path and metadata
        sliced_extractor.file_path = self.file_path
        if hasattr(self, '_general_metadata'):
            sliced_extractor._general_metadata = self._general_metadata
        
        return sliced_extractor
    
    def get_original_frame_indices(self, plane_index: Optional[int] = None) -> np.ndarray:
        """
        Get the original frame indices for each sample.

        Returns the index of the original frame for each sample, mapping processed samples
        back to their corresponding frames in the raw microscopy data. This accounts for
        any filtering, subsampling, or exclusions (such as flyback frames) performed by
        the extractor.

        Parameters
        ----------
        plane_index : int, optional
            Which plane to use for frame index calculation in volumetric data.
            If None, plane_index is set to the last plane in the volume. This is because the timestamp of the acquisition of the last plane in a volume is typically set as the timestamp of the volume as a whole. It must be less than the total number of planes.

        Returns
        -------
        np.ndarray
            Array of original frame indices (dtype: int64) with length equal to the
            number of samples. Each element represents the index of the original
            microscopy frame that corresponds to that sample.

        Notes
        -----
        **Frame Index Calculation:**

        - **Planar data**: Frame indices are sequential (0, 1, 2, ...)
        - **Multi-channel data**: Accounts for channel interleaving
        - **Volumetric data**: Uses the specified plane (default: last plane)
        - **Multi-file data**: Includes file offsets for global indexing
        - **Flyback frames**: Automatically excluded from indexing

        **Common Use Cases:**

        - Synchronizing with external timing systems
        - Mapping back to original acquisition timestamps
        - Data provenance and traceability
        - Cross-referencing with raw data files

        **Examples:**

        For a 3-sample volumetric dataset with 5 planes per volume:
        - Default behavior returns indices [4, 9, 14] (last plane of each volumetric sample)
        - With plane_index=0 returns indices [0, 5, 10] (first plane of each volumetric sample)
        """
        num_planes = self.get_num_planes()
        if plane_index is not None:
            assert plane_index < num_planes, f"Plane index {plane_index} exceeds number of planes {num_planes}."
        else:
            plane_index = num_planes - 1

        # Initialize array to store timestamps
        num_samples = self.get_num_samples()
        frame_indices = np.zeros(num_samples, dtype=np.int64)

        # For each sample, extract its timestamp from the corresponding file and IFD
        for sample_index in range(num_samples):

            # Get the last frame in this sample to get the timestamps
            frame_index = sample_index * num_planes + plane_index
            table_row = self._frames_to_ifd_table[frame_index]

            file_index = int(table_row["file_index"])
            ifd_index = int(table_row["IFD_index"])

            # The ifds are local within a file, so we need to add and offset
            # equal to the number of IFDs in the previous files
            file_offset = sum(self._ifds_per_file[:file_index]) if file_index > 0 else 0

            frame_indices[sample_index] = ifd_index + file_offset

        return frame_indices


class KimScanImageImagingInterface(BaseImagingExtractorInterface):

    Extractor = KIMScanImageImagingExtractor

    def __init__(
        self,
        file_path: FilePath,
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
                    origin_coords = [float(x) for x in origin_coords]
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
        file_path: FilePath,
        roi_info_file_path: FilePath,
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
