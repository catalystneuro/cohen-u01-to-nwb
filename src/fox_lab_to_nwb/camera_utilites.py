from lxml import etree
import os

def extract_phantom_metadata(xml_path):
    # Parse XML using lxml
    tree = etree.parse(xml_path)
    root = tree.getroot()

    # Create metadata dictionary
    metadata = {}

    # Define paths and their corresponding keys/types
    paths = {
        # Camera settings
        "frame_rate": (".//FrameRateDouble", float),
        "total_frames": (".//TotalImageCount", int),
        "first_frame": (".//FirstMovieImage", int),
        "image_count": (".//ImageCount", int),
        # Image properties
        "width": (".//biWidth", int),
        "height": (".//biHeight", int),
        "bit_depth": (".//biBitCount", int),
        "bit_depth_recording": (".//RecBPP", int),
        # Camera info
        "camera_model": (".//CameraModel", str),
        "camera_version": (".//CameraVersion", int),
        "firmware_version": (".//FirmwareVersion", int),
        "software_version": (".//SoftwareVersion", int),
        "serial": (".//Serial", int),
        # Timing
        "shutter_ns": (".//ShutterNs", int),
        "frame_delay_ns": (".//FrameDelayNs", int),
        # Image settings
        "compression": (".//Compression", int),
        "saturation": (".//Saturation", float),
        "brightness": (".//Bright", int),
        "contrast": (".//Contrast", int),
        "gamma": (".//Gamma", float),
        # Trigger settings
        "trigger_frame": (".//TrigFrame", int),
        "post_trigger": (".//PostTrigger", int),
        # Auto exposure
        "auto_exposure": (".//AutoExposure", bool),
        "auto_exp_level": (".//AutoExpLevel", int),
        "auto_exp_speed": (".//AutoExpSpeed", int),
    }

    # Extract all metadata based on paths
    for key, (xpath, type_conv) in paths.items():
        element = root.find(xpath)
        if element is not None and element.text:
            try:
                metadata[key] = type_conv(element.text)
            except (ValueError, TypeError):
                metadata[key] = None
        else:
            metadata[key] = None

    # Special handling for trigger time
    trigger_date = root.find(".//TriggerTime/Date")
    trigger_time = root.find(".//TriggerTime/Time")
    if trigger_date is not None and trigger_time is not None:
        metadata["trigger_time"] = f"{trigger_date.text} {trigger_time.text}"

    # Get image acquisition position and size
    metadata["acquisition"] = {
        "pos_x": int(root.find(".//ImPosXAcq").text) if root.find(".//ImPosXAcq") is not None else None,
        "pos_y": int(root.find(".//ImPosYAcq").text) if root.find(".//ImPosYAcq") is not None else None,
        "width": int(root.find(".//ImWidthAcq").text) if root.find(".//ImWidthAcq") is not None else None,
        "height": int(root.find(".//ImHeightAcq").text) if root.find(".//ImHeightAcq") is not None else None,
    }

    # Get white balance gains
    wb_element = root.find(".//WBGain")
    if wb_element is not None:
        metadata["white_balance"] = {
            "red": float(wb_element.find("Red").text) if wb_element.find("Red") is not None else None,
            "blue": float(wb_element.find("Blue").text) if wb_element.find("Blue") is not None else None,
        }

    return metadata


from typing import Any


def extract_fastec_metadata(file_path: str) -> dict[str, dict[str, Any]]:
    """
    Extract metadata from a Fastec camera metadata file.

    Parameters
    ----------
    file_path : str
        Path to the Fastec metadata file.

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Nested dictionary containing the parsed metadata.
        The top level dictionary has sections as keys ('image', 'camera', 'record', 'normalization').
        Each section contains a dictionary of key-value pairs with automatically converted data types.

    Notes
    -----
    The function automatically converts values to appropriate types:
    - Integers for numeric values
    - Floats for decimal numbers
    - Lists for matrix values [x,y,z]
    - Tuples for bit modes (e.g., "10:3")
    - Strings for text and other values

    Examples
    --------
    >>> metadata = extract_fastec_metadata('metadata.txt')
    >>> frame_rate = metadata['record']['fps']
    >>> resolution = (metadata['image']['width'], metadata['image']['height'])

    Raises
    ------
    FileNotFoundError
        If the metadata file is not found.
    PermissionError
        If there are insufficient permissions to read the file.
    """

    def parse_value(value: str) -> int | float | list[int] | tuple[int, int] | str:
        """
        Parse a string value into its appropriate type.

        Parameters
        ----------
        value : str
            The string value to parse

        Returns
        -------
        int | float | list[int] | tuple[int, int] | str
            Parsed value in its appropriate type
        """
        # Try to convert to int
        try:
            return int(value)
        except ValueError:
            pass

        # Try to convert to float
        try:
            return float(value)
        except ValueError:
            pass

        # Handle matrix values [x,y,z]
        if value.startswith("[") and value.endswith("]"):
            try:
                return [int(x) for x in value[1:-1].split(",")]
            except ValueError:
                return value

        # Handle bit mode (e.g., "10:3")
        if ":" in value and len(value.split(":")) == 2:
            try:
                return tuple(int(x) for x in value.split(":"))
            except ValueError:
                return value

        # Return as string if no other type matches
        return value

    # Check if file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Metadata file not found: {file_path}")

    metadata: Dict[str, Dict[str, Any]] = {}
    current_section: Union[str, None] = None

    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Check if line is a section header
                if line.startswith("[") and line.endswith("]"):
                    current_section = line[1:-1].lower()
                    metadata[current_section] = {}
                    continue

                # Parse key-value pairs
                if "=" in line and current_section is not None:
                    key, value = line.split("=", 1)
                    key = key.strip()
                    value = value.strip()

                    # Parse value to appropriate type
                    parsed_value = parse_value(value)

                    metadata[current_section][key] = parsed_value

        return metadata

    except PermissionError:
        raise PermissionError(f"Insufficient permissions to read file: {file_path}")
    except Exception as e:
        raise Exception(f"Error parsing metadata file: {str(e)}")
