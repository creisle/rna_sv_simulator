import rna_sv_simulator.io.file_ops as io


def make_par_file(env_config, run_config):
    """Create a .par file for Flux

    Args:
        env_config (dict):
        run_config (dict):

    Returns:
        None
    """

    path = run_config.get('paths').get('file_par_name')

    parameters = build_parameter_dict(env_config, run_config)

    test_list_of_items_in_dict(
        env_config.get('flux').get('input').get(
            'required_fields'),
        parameters)

    string = build_par_string(parameters)
    io.write_txt_file(path, string)

    return


def build_parameter_dict(env_config, run_config):
    """Build a parameter dictionary

    Load the default values from the env config and then
    overwrite with anything that is in the run config

    Args:
        env_config (dict):
        run_config (dict):

    Returns:
        dict
    """

    default = env_config.get('flux').get('parameters')
    user = run_config.get('flux').get('parameters')

    parameters = default

    for key, value in user.items():
        parameters[key] = value

    return parameters


def test_list_of_items_in_dict(required, input_dict):
    """Test that the required values are in the parameters dict

    Args:
        required (list):
        input_dict (dict):

    Returns:
        Null

    Raises:
        ValueError
    """

    for item in required:
        if item not in input_dict:
            raise ValueError(
                '{} must be Flux parameters'.format(
                    item)
            )


def build_par_string(parameters, delimiter='\t'):
    """Build the string to write out

    Args:
        parameters (dict)
        delimiter (str): delimiter to separate lines by

    Returns:
        str
    """

    output = ''

    for key, value in parameters.items():

        output += '{key}{delimiter}{value}\n'.format(
            key=key,
            delimiter=delimiter,
            value=value
        )

    return output



