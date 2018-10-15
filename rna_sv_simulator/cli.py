# -*- coding: utf-8 -*-

'''Console script for rna_sv_simulator.'''
import argparse
import os
import subprocess


from mavis.util import read_inputs
from mavis.annotate import file_io


from . import gtf_output
from . import fusion
from .flux import make_par



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
    required.add_argument(
        '-n', '--input', required=True, metavar='FILE',
        help='Input breakpoint pair file for inputting the structural variant'
    )
    required.add_argument(
        '-o', '--output', required=True, metavar='DIR',
        help='Path to the directory where to write the output and intermediary files'
    )
    flux_args = parser.add_argument_group('Flux Simulator Arguments')
    flux_args.add_argument(
        '-m', '--nb_molecules', dest='NB_MOLECULES', type=int, default=100000,
        help='Number of RNA molecules initially in the experiment.'
    )
    flux_args.add_argument(
        '-x', '--expression_profile', metavar='FILE',
        help='Path to the expression profile to be used. Will be randomly generated if not selected'
    )
    flux_args.add_argument(
        '-r', '--read_number', dest='READ_NUMBER', type=int, default=10000,
        help='Number of reads')
    flux_args.add_argument(
        '-e', '--err_file', dest='ERR_FILE', choices=['35', '76'],
        help='Choose a default error model which are provided by flux for the corresponding read lengths'
    )
    flux_args.add_argument(
        '--unique_ids', dest='UNIQUE_IDS', default=True, type=string_to_boolean, choices={'true', 'false'},
        help='Create unique read identifiers for paired reads. Information about the relative orientation is left out of the read id and encoded in the pairing information. All /1 reads are sense reads, all /2 reads are anti-sense reads. This option is useful if you want to identify paired reads based on the read ids.'
    )
    args = parser.parse_args()

    # read the input svs
    sv_list = read_inputs(
        [args.input],
        add_default={
            'stranded': True,
            'protocol': 'transcriptome'
        }
    )
    for sv in sv_list:
        print(sv)
    if len(sv_list) != 1:
        raise NotImplementedError('Currently only supportting input of a single structural variant')
    # load the annotations files
    print('load the reference genome')
    reference_genome = file_io.ReferenceFile('reference_genome', *args.reference_genome, eager_load=True).content
    print('load the annotations')
    annotations = file_io.ReferenceFile('annotations', *args.annotations, eager_load=True).content
    #print('mutate the original genome')
    #mutant_ref_genome, mutant_annotations = fusion.mutate(
    #    reference_genome,
    #    annotations,
    #    sv_list[0]
    #)
    mutant_ref_genome = reference_genome
    mutant_annotations = annotations

    # write the intermediate files for flux simulator
    os.makedirs(os.path.join(args.output, 'genome'), exist_ok=True)
    written = set()
    for chrom, seq in reference_genome.items():
        if chrom in written:
            continue
        written.add(chrom)
        with open(os.path.join(args.output, 'genome', 'chr{}.fa'.format(chrom)), 'w') as fh:
            fh.write('> {}\n{}\n'.format(chrom, str(seq.seq)))

    # write the gtf annotation files
    gtf_file_name = os.path.join(args.output, 'annotations.gtf')
    gtf_output.output_gtf_file(mutant_annotations, gtf_file_name)

    # write the par file
    par_filename = os.path.join(args.output, 'flux.par')
    with open(par_filename, 'w') as fh:
        flux_settings = {
            'NB_MOLECULES': args.NB_MOLECULES,
            'READ_NUMBER': args.READ_NUMBER,
            'ERR_FILE': args.ERR_FILE,
            'UNIQUE_IDS': args.UNIQUE_IDS,
            'GEN_DIR': 'genome',
            'REF_FILE_NAME': os.path.basename(gtf_file_name)
        }
        if args.expression_profile:  # If they selected an expression profile, use it
            flux_settings['PRO_FILE_NAME'] = args.expression_profile
        fh.write(make_par.build_par_string(flux_settings) + '\n')

    # run flux
    subprocess.check_output('flux-simulator -p {}'.format(par_filename))



if __name__ == '__main__':
    main()
