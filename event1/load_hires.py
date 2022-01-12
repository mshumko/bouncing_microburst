import pathlib
import json

import numpy as np

from event1 import config

def load_hires(sc_id, date):
    """
    Finds and loads the FIREBIRD HiRes data.
    """
    search_str = f'FU{sc_id}_Hires_{date.strftime("%Y-%m-%d")}_L2.txt'
    hr_paths = sorted(pathlib.Path(config.FB_DIR).rglob(search_str))
    assert len(hr_paths) == 1, (
        f'{len(hr_paths)} HiRes paths found in {pathlib.Path(config.FB_DIR)} '
        f'that match {search_str}.'
        )
    hr = spacepy.datamodel.readJSONheadedASCII(str(hr_paths[0]))
    hr['Time'] = pd.to_datetime(hr['Time'])
    return

def readJSONheadedASCII(file_path):
    """
    My basic implementation of spacepy.datamodel.readJSONheadedASCII
    (for some who) have difficulty installing it.
    """
    # Read in the JSON header.
    header = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header.append(line[1:])
            else:
                data = f.readlines()

    parsed_header = json.loads(header)

    return