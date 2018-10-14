# -*- coding: utf-8 -*-

'''Console script for rna_sv_simulator.'''
import argparse
import os


def string_to_boolean(string):
    if str(string).lower() in ['t', 'true', '1', 'y', 'yes']:
        return True
    elif str(string).lower() in ['f', 'false', '0', 'n', 'no']:
        return False
    raise TypeError('{} is not a valid boolean flag value'.format(string))


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required Arguments')
    required.add_argument(
        '-g', '--reference_genome', required=True, nargs='+', metavar='FILE',
        help='The reference genome FASTA reference sequence file(s)')
    required.add_argument(
        '-a', '--annotations', required=True, nargs='+', metavar='FILE',
        help='File(s) with the mavis annotation data'
    )
    flux_args = parser.add_argument_group('Flux Simulator Arguments')
    flux_args.add_argument(
        '-m', '--nb_molecules', dest='NB_MOLECULES', type=int, default=100000,
        help='Number of RNA molecules initially in the experiment.'
    )
    flux_args.add_argument(
        '-r', '--read_number', dest='READ_NUMBER', type=int, default=10000,
        help='Number of reads')
    flux_args.add_argument(
        '-e', '--err_file', dest='ERR_FILE', choices=['35', '76'],
        help='Choose a default error model which are provided by flux for the corresponding read lengths'
    )
    flux_args.add_argument(
        '--unique_ids', dest='UNIQUE_IDS', default=True, type=string_to_boolean,
        help='Create unique read identifiers for paired reads. Information about the relative orientation is left out of the read id and encoded in the pairing information. All /1 reads are sense reads, all /2 reads are anti-sense reads. This option is useful if you want to identify paired reads based on the read ids.'
    )
    args = parser.parse_args()


if __name__ == '__main__':
    main()
