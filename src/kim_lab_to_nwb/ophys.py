from typing import Optional, Tuple
from pathlib import Path

import numpy as np
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from neuroconv.utils import FolderPathType, DeepDict, FilePathType
from roiextractors import ImagingExtractor
from roiextractors.extraction_tools import PathType
from tifffile import TiffFile
from tqdm import tqdm
from pynwb import NWBFile
from pynwb.ophys import PlaneSegmentation, ImageSegmentation, RoiResponseSeries, ImagingPlane, OpticalChannel
from pynwb.device import Device
from pymatreader import read_mat
from pynwb.ophys import DfOverF

from cohen_u01_nwb_conversion_utils.utils import match_paths

from neuroconv import ConverterPipe


class ScanImageConverter(ConverterPipe):
    """Convert ScanImage data to NWB format."""

    def __init__(self, file_path: FilePathType, sampling_frequency: float = 30.0):
        """
        Initialize the converter.

        Parameters
        ----------
        file_path : str
            Path to the ScanImage file.
        sampling_frequency : float
            Sampling frequency in Hz.
        """

        return None


class MultiTiffMultiPageTiffImagingExtractor(ImagingExtractor):
    """A ImagingExtractor for multiple TIFF files that each have multiple pages."""

    @staticmethod
    def write_imaging(imaging, save_path: PathType, overwrite: bool = False):
        raise NotImplementedError("Writing is not supported for this extractor.")

    extractor_name = "multi-tiff multi-page Imaging Extractor"
    is_writable = False

    def __init__(self, folder_path: PathType, pattern="{}.tif{}", sampling_frequency: float = 30.0):
        """Create a MultiTiffMultiPageImagingExtractor instance.

        Parameters
        ----------
        folder_path : str
            List of path to each TIFF file.
        """

        super().__init__()
        self.folder_path = folder_path

        self.tif_paths = match_paths(folder_path, pattern)

        self.page_tracker = []
        page_counter = 0
        for file_path in tqdm(self.tif_paths, "extracting page lengths"):
            # print(f"{file_path=}")
            with TiffFile(file_path) as tif:
                self.page_tracker.append(page_counter)
                page_counter += len(tif.pages)
                # print(f"num pages: {len(tif.pages)}")
        self.page_tracker = np.array(self.page_tracker)
        page = tif.pages[0]

        self._num_frames = page_counter
        self._num_columns = page.imagewidth
        self._num_rows = page.imagelength

        self._sampling_frequency = sampling_frequency

        self._kwargs = {"folder_path": folder_path}

    def get_video(self, start_frame: int = None, end_frame: int = None, channel: Optional[int] = 0) -> np.ndarray:
        frame_idxs = np.arange(start_frame or 0, end_frame or self._num_frames)
        file_idxs = np.searchsorted(self.page_tracker, frame_idxs + 1) - 1  # index of the file that contains the frame

        file_start_idxs = self.page_tracker[file_idxs]  # index of the first frame of the file
        frame_offset_idxs = frame_idxs - file_start_idxs  # index of the frame in the file
        # dict of file_idx: frame_offset_idxs
        index_dict = {x: frame_offset_idxs[file_idxs == x] for x in np.unique(file_idxs)}

        self.tif_path_list = list(self.tif_paths)
        data = []
        for file_idx, file_frame_offset_idxs in tqdm(index_dict.items()):
            with TiffFile(self.tif_path_list[file_idx]) as tif:
                for frame_offset_idx in file_frame_offset_idxs:
                    page = tif.pages[frame_offset_idx]
                    data.append(page.asarray())

        return np.array(data)

    def get_image_size(self) -> Tuple[int, int]:
        return self._num_rows, self._num_columns

    def get_num_frames(self):
        return self._num_frames

    def get_sampling_frequency(self):
        return self._sampling_frequency

    def get_num_channels(self):
        return 1

    def get_channel_names(self):
        return ["channel_0"]


class MultiTiffMultiPageTiffImagingInterface(BaseImagingExtractorInterface):
    """Interface for multiple multi-page TIFF files."""

    display_name = "MultiTiff MultiPage Imaging Extractor"
    Extractor = MultiTiffMultiPageTiffImagingExtractor

    @classmethod
    def get_source_schema(cls) -> dict:
        source_schema = super().get_source_schema()
        source_schema["properties"]["folder_path"]["description"] = "Directory that contains the TIFF files."
        return source_schema

    def __init__(self, folder_path: FolderPathType, pattern: str, sampling_frequency: float, verbose: bool = False):
        """
        Initialize reading of TIFF files.

        Parameters
        ----------
        folder_path : FilePathType
        verbose : bool, default: False
        """
        super().__init__(
            folder_path=folder_path, pattern=pattern, sampling_frequency=sampling_frequency, verbose=verbose
        )


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
        image_shape : Optional[tuple[int, int]]\
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
