from typing import Optional, Tuple

import numpy as np
from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface
from neuroconv.utils import FolderPathType, DeepDict
from roiextractors import ImagingExtractor
from roiextractors.extraction_tools import PathType
from tifffile import TiffFile
from tqdm import tqdm

from ..utils import match_paths


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
            #print(f"{file_path=}")
            with TiffFile(file_path) as tif:
                self.page_tracker.append(page_counter)
                page_counter += len(tif.pages)
                #print(f"num pages: {len(tif.pages)}")
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
        #print(f"{file_idxs=}")
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
