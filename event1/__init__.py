import pathlib

import configparser

HERE = pathlib.Path(__file__).parent.resolve()
settings = configparser.ConfigParser()
settings.read(HERE / 'config.ini')

FB_DIR = settings['Paths'].get('FB_DIR', None)
FB_DIR = pathlib.Path(FB_DIR)

config = {'PROJECT_DIR': HERE, 'FB_DIR': FB_DIR}
