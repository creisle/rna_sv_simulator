import os
import yaml


def file_exists(path):
    """Check if file exists and is readable

    Args:
        path (str)

    Returns:
        bool

    Raises:
        IOError
    """

    exists = False

    if os.path.isfile(path):
        if os.access(path, os.R_OK):
            exists = True

        else:
            raise IOError('{file} is not readable'.format(
                file=path
            ))

    else:
        raise IOError('{file} does not exist'.format(
            file=path
        ))

    return exists


def folder_exists(path):
    """Check if a folder exists and is writable

    Args:
        path (str):

    Returns:
        bool
    """

    exists = False

    if os.path.isdir(path):
        if os.access(path, os.W_OK):
            exists = True

    return exists


def make_dir(path):
    """Make a new folder

    Args:
        path (str):

    Returns:
        None
    """

    if not folder_exists():
        os.makedirs(path, 0o777)

    return


def load_yaml(path):
    """Load a yaml file

    Args:
        path (str):

    Returns:
        dict

    Raises:
        IOError
    """

    with open(path, 'r') as stream:
        try:
            loaded_yaml = yaml.load(stream)
        except yaml.YAMLError:
            raise IOError('Invalid YAML file')

    return loaded_yaml


def write_txt_file(path, string):
    """

    Args:
        path (str):     file path
        string (str):   text to write

    Returns:
        None
    """

    with open(path, 'w') as fh:

        try:
            fh.write(string)
        except IOError:
            raise IOError('Failed to write {}'.format(
                path
            ))

    return
