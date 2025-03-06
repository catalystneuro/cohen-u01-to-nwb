"""Interface for Thor TIFF files with OME metadata."""

from datetime import datetime, timezone
from typing import Optional

import numpy as np
from pydantic import validate_call, FilePath

from neuroconv.datainterfaces.ophys.baseimagingextractorinterface import BaseImagingExtractorInterface 
from neuroconv.utils import DeepDict, get_json_schema_from_method_signature
from neuroconv.datainterfaces import ThorImagingInterface


class ThorImagingInterface(ThorImagingInterface):

    def get_metadata(self) -> DeepDict:
        """
        Retrieve the metadata for the Thor imaging data.

        Returns
        -------
        DeepDict
            Dictionary containing metadata including device information, imaging plane details,
            and photon series configuration.
        """
        metadata = super().get_metadata()

        # Modifications for Dickerson Lab Metadata
        channel_indicators = {"ChanA": "GCaMP", "ChanB": "tdTomato"} 
        excitation_lambda_mapping = {"GCaMP": 495.0, "tdTomato": 554.0}
        emission_lambda_mapping = {"GCaMP": np.nan, "tdTomato": np.nan}  # TODO check out which excitaiton was used
        
        indicator = channel_indicators[self.channel_name]
        excitation_lambda = excitation_lambda_mapping[indicator]
        emission_lambda = emission_lambda_mapping[indicator]
        
        imaging_plane = metadata["Ophys"]["ImagingPlane"][0]
        imaging_plane["excitation_lambda"] = excitation_lambda
        imaging_plane["indicator"] = indicator
        imaging_plane["location"] = "unknown"  # You can add it if you have it

        imaging_plane_name = imaging_plane["name"].replace(f"{self.channel_name}", f"{indicator}")
        imaging_plane["name"] = imaging_plane_name
        
        optical_channel = imaging_plane["optical_channel"][0]
        optical_channel["emission_lambda"] = emission_lambda
        optical_channel_name = optical_channel["name"].replace(f"{self.channel_name}", f"{indicator}")
        optical_channel["name"] = optical_channel_name
        optical_channel["description"] = f"{indicator} channel"    
    
        two_photon_series = metadata["Ophys"]["TwoPhotonSeries"][0]
        two_photon_series["name"] = f"TwoPhotonSeries{indicator}"
        two_photon_series["imaging_plane"] = imaging_plane_name
        return metadata


