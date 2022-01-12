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
                raw_data_str = f.readlines()

    # Massage the header
    clean_header_str = ''.join(header_list).replace('\n', '')
    parsed_header = json.loads(clean_header_str)
    # Massage the data
    raw_data_str = [row.replace('\n', '') for row in raw_data_str]
    # Parse times
    times_str = [row.split()[0] for row in raw_data_str]
    # Parse the other columns
    data_converted = np.array([row.split()[1:] for row in raw_data_str]).astype(float)

    data = HiRes()
    data['Time'] = pd.to_datetime(times_str)
    for key in parsed_header:
        if key == 'Time':
            continue
        key_header = parsed_header[key]
        # column attributes
        if isinstance(key_header, dict):
            print(key_header['DIMENSION'])
            start_column = key_header['START_COLUMN']-1
            end_column = key_header['START_COLUMN']-1+key_header['DIMENSION'][0]
            if key_header['DIMENSION'][0] == 1:
                data[key] = data_converted[:, start_column]
            else:
                data[key] = data_converted[:, start_column:end_column]
            
            data.attrs[key] = key_header
        else:
            # Global attributes
            if key in ['CADENCE', 'CAMPAIGN']:
                data.attrs[key] = float(key_header)
            else:
                data.attrs[key] = key_header
        
    return data

class HiRes(dict):
    """
    Expand Python's dict class to include an attr attribute dictionary.
    
    Code credit goes to Matt Anderson:
    https://stackoverflow.com/questions/2390827/how-to-properly-subclass-dict-and-override-getitem-setitem
    (blame him for problems)
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        self.attrs = {}
        return

if __name__ == '__main__':
    hr = load_hires(4, datetime(2015, 5, 25))