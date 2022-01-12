import sys
import pathlib
import configparser

# Run the configuration script when the user runs
# python3 -m event1 [init, initialize, config, or configure]

here = pathlib.Path(__file__).parent.resolve()


if (len(sys.argv) > 1) and (sys.argv[1] in ['init', 'initialize', 'config', 'configure']):
    print('Running the configuration script.')
    # event1 FIREBIRD dir
    s = (
        f'What is the FIREBIRD data directory?\n'
    )
    FB_DIR = input(s)

    # Create a configparser object and add the user configuration.
    config = configparser.ConfigParser()
    config['Paths'] = {'FB_DIR': FB_DIR}

    with open(here / 'config.ini', 'w') as f:
        config.write(f)

else:
    print(
        'This is a configuration script to set up config.ini file. The config '
        'file contains the FIREBIRD data directory. To configure this package, '
        'run python3 -m event1 config'
    )