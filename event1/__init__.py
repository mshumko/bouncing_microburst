import pathlib
import warnings

import configparser

HERE = pathlib.Path(__file__).parent.resolve()
settings = configparser.ConfigParser()
settings.read(HERE / 'config.ini')

try:
    FB_DIR = settings['Paths'].get('FB_DIR', None)
    FB_DIR = pathlib.Path(FB_DIR)
    config = {'PROJECT_DIR': HERE, 'FB_DIR': FB_DIR}
except KeyError:  # Raised if config.ini does not have Paths.
    warnings.warn(
        "The event1 config.ini file doesn't exist or is missing the FB_DIR. "
        "run 'python3 -m event1 config' to set it up."
        )
