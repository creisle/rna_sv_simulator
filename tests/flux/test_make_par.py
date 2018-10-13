import os
import pytest
import tempfile

from hypothesis import strategies as st
from hypothesis import given

import rna_sv_simulator.flux.make_par

ENV_CONFIG = {
    'flux': {
        'input': {
            'required_fields': [
                'REF_FILE_NAME',
                'GEN_DIR',
                'NB_MOLECULES',
                'READ_NUMBER',
                'ERR_FILE',
                'FASTA',
                ],
        },
        'paths': {
            'folder_bundle': 'flux-bundle',
            'folder_genome': 'genome',
            'file_par': 'test.par',
            'file_reference': 'annotation.gtf',
        },
        'parameters': {
            'REF_FILE_NAME': 'annotation.gtf',
            'GEN_DIR': 'genome',
            'NB_MOLECULES': 100000,
            'READ_NUMBER': 10000,
            'ERR_FILE': 76,
            'FASTA': 'YES',
        },
    },
}

RUN_CONFIG = {
    'paths': {
      'dir_run': '~/testing/rna_sv_simulator',
      'file_par_name': 'flux_settings.par',
    },
    'flux': {
        'parameters': {
            'number_molecules': 100000,
            'number_reads': 10000,
            'number_read_length': 150,
            'err_file': 76,
            'fasta': 'YES',
        },
    },
}

class TestMakeParFile:
    """Tests for make_par_file"""

    def test_make_par_file(self):

        env_config = ENV_CONFIG
        run_config = RUN_CONFIG

        tmp_dir = tempfile.mktemp()
        run_config['paths']['dir_run'] = tmp_dir

        rna_sv_simulator.flux.make_par.make_par_file(env_config, run_config)


class TestBuildParameterDict:
    """Tests for build_parameter_dict"""

    result = rna_sv_simulator.flux.make_par.build_parameter_dict(ENV_CONFIG, RUN_CONFIG)

    assert isinstance(result, dict)
    assert 'fasta' in result



class TestListOfItemsInDict:
    """Tests for test_list_of_items_in_dict"""

    input_dict = {
        'foo': 'a',
        'bar': 'b',
    }

    def test_item_in_dict(self):

        required = ['foo']
        rna_sv_simulator.flux.make_par.test_list_of_items_in_dict(required, self.input_dict)

    def test_item_not_in_dict(self):

        required = ['baz']

        with pytest.raises(ValueError):
            rna_sv_simulator.flux.make_par.test_list_of_items_in_dict(required, self.input_dict)


class TestBuildParString:

    def test_string_sub(self):

        parameters = {
            'foo': 'bar'
        }
        delimiter = '\t'

        result = rna_sv_simulator.flux.make_par.build_par_string(parameters, delimiter)

        assert delimiter in result
        assert parameters.get('foo') in result



