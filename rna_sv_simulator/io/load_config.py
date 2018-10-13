import os

from rna_sv_simulator.io.file_ops import file_exists, load_yaml


def load_env_config():
    """Load the env config

    Returns:
        dict
    """

    config_file = get_default_env_config_path()
    file_exists(config_file)
    config = load_env_yaml(config_file)

    return config


def get_default_config_location():
    """Get default location of config dir

    Returns:
        str: abs path to config file
    """

    this_dir = os.path.dirname(
        os.path.abspath(__file__)
    )

    config_dir = os.path.abspath(
        os.path.join(this_dir,
                     '..',
                     '..',
                     'configuration')
    )

    return config_dir


def get_default_env_config_path():
    """Get default location of config file

    Returns:
        str: abs path to config file
    """

    config_dir = get_default_config_location()

    config_file = os.path.join(config_dir,
                               'env.yaml')

    return config_file


def load_env_yaml(config_file):
    """Load the env.yaml file, w/basic checking

    Args:
        config_file:

    Returns:
        dict
    """

    loaded_yaml = load_yaml(config_file)

    # TODO Add error checking on loaded dict

    return loaded_yaml


