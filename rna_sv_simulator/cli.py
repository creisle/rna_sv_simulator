# -*- coding: utf-8 -*-

"""Console script for rna_sv_simulator."""
import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def process_file(path):
    """
    reads in the genomic and annotation files.
    :return:
    """
    result = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if "i>" in line:
                result.append(line)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genomic", type=str, required=True,
                        help="File with genomic data. /path/to/file")
    parser.add_argument("-a", "--annotation", type=str, required=True,
                        help="File with annotation data. /path/to/file")
    parser.add_argument("-m", "--molecules", nargs=1, type=int, required=True,
                        help="Number of molecules for par files. Defaults to 100,000")
    parser.add_argument("-r", "--read", nargs=1, type=int, required=True,
                        help="Read number. Defaults to 10,000")
    args = parser.parse_args()
    genome_file = process_file(args.genomic)
    annotation_file = process_file(args.annotation)
    # sys.exit(())


if __name__ == "__main__":
    main()
