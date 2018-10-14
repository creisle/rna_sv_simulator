# -*- coding: utf-8 -*-

"""Console script for rna_sv_simulator."""
import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genomic", type=str, required=True,
                        help="File with genomic data")
    parser.add_argument("-a", "--annotation", type=str, required=True,
                        help="File with annotation data")
    parser.add_argument("-p", "--par", type=str, required=True,
                        help="File with par files")
    parser.add_argument("-m", "--molecules", type=int, required=False,
                        help="Number of molecules for par files. Defaults to 100,000")
    parser.add_argument("-r", "--read", type=int, required=False,
                        help="Read number. Defaults to 10,000")
    # sys.exit(())


if __name__ == "__main__":
    main()
