import os
import pytest

import rna_sv_simulator.io.load_config as load_config


class TestLoadEnvConfig:
    """Tests for load_env_config"""

    def test_get_config(self):

        result = load_config.load_env_config()

        assert isinstance(result, dict)
        assert 'flux' in result


class TestGetDefaultConfigLocation:
    """Tests for get_default_config_location"""

    def test_get_config_dir(self):
        """Get default"""

        result = load_config.get_default_config_location()

        assert isinstance(result, str)
        assert os.path.isdir(result) is True


class TestGetDefaultEnvConfigPath:
    """Tests for get_default_env_config_path"""

    def test_get_env_config_path(self):
        """Get default"""

        result = load_config.get_default_env_config_path()

        assert isinstance(result, str)
        assert os.path.isfile(result) is True
        assert result.endswith('env.yaml')


class TestLoadEnvYaml:
    """Tests for load_env_yaml"""

    yaml_file = load_config.get_default_env_config_path()

    def test_valid_yaml(self):
        """Valid YAML"""

        result = load_config.load_env_yaml(self.yaml_file)

        assert isinstance(result, dict)
        assert 'flux' in result
