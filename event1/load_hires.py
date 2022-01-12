import pathlib
import json
from datetime import datetime

import numpy as np
import pandas as pd

from event1 import config

def load_hires(sc_id, date):
    """
    Finds and loads the FIREBIRD HiRes data.
    """
    search_str = f'FU{sc_id}_Hires_{date.strftime("%Y-%m-%d")}_L2.txt'
    hr_paths = sorted(pathlib.Path(config['FB_DIR']).rglob(search_str))
    assert len(hr_paths) == 1, (
        f'{len(hr_paths)} HiRes paths found in {pathlib.Path(config["FB_DIR"])} '
        f'that match {search_str}.'
        )
    hr = readJSONheadedASCII(str(hr_paths[0]))
    hr['Time'] = pd.to_datetime(hr['Time'])
    return

def readJSONheadedASCII(file_path):
    """
    My basic implementation of spacepy.datamodel.readJSONheadedASCII
    (for some who) have difficulty installing it.
    """
    # Read in the JSON header.
    header_list = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_list.append(line[1:])
            else:
                data = f.readlines()
                
    clean_header_str = ''.join(header_list).replace('\n', '')
    parsed_header = json.loads(clean_header_str)

    return

if __name__ == '__main__':
    hr = load_hires(4, datetime(2015, 5, 25))