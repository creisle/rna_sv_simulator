# -*- coding: utf-8 -*-

"""Console script for rna_sv_simulator."""
import sys
import click


@click.command()
def main(args=None):
    """Console script for rna_sv_simulator."""

    click.echo("Replace this message by putting your code into "
               "rna_sv_simulator.cli.main")
    click.echo("See click documentation at http://click.pocoo.org/")
    raise SystemExit(1)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
