from lxml import etree as ET
from datetime import datetime

fpath = "/Users/bendichter/Downloads/xz_011.xml"


def parse_single_frame_time(date_text, time_text):
    elems = time_text.split(' ')
    time = elems[0]+elems[1][:-3]

    return datetime.strptime(f"{date_text} {time}", '%a %b %d %Y %H:%M:%S.%f')


def parse_start_time(fpath) -> datetime:
    # read xml file
    tree = ET.parse(fpath)
    root = tree.getroot()

    timeblock = root.findall('TIMEBLOCK')[0]
    return parse_single_frame_time(timeblock[0].text, timeblock[1].text).astimezone("UTC")


def parse_times(fpath: str) -> list:
    """
    Parse the time stamps from a Phantom xml file directly from the TIMEBLOCK. This is generally not necessary, because
    the sample times follow the frame rate, but it can be useful for debugging.

    Parameters
    ----------
    fpath : str
    """

    # read xml file
    tree = ET.parse(fpath)
    root = tree.getroot()

    timeblock = root.findall('TIMEBLOCK')[0]

    for i, (date_obj, time_obj) in enumerate(zip(timeblock[::2], timeblock[1::2])):
        dt = parse_single_frame_time(date_obj.text, time_obj.text)
        if not i:
            start_datetime = dt
            tts = [0]
        else:
            tts.append((dt-start_datetime).total_seconds())

    return tts