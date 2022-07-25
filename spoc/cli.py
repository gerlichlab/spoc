"""Console script for spoc."""
import sys
import click
from spoc.contacts import ContactExpander
from spoc.io import FileManager



@click.group()
def main(args=None):
    """Console script for spoc."""
    pass


@click.command()
@click.argument('input')
@click.argument('output')
@click.option(
    "-n",
    "--n_fragments",
    default=3,
    help="Number of fragments per read to expand",
)
def expand(input, output, n_fragments):
    """Script for expanding fragments to contacts"""
    expander = ContactExpander(number_fragments=n_fragments)
    file_manager = FileManager()
    input_fragments = file_manager.load_porec_fragments(input)
    expanded_contacts = expander.expand(input_fragments)
    file_manager.write_multiway_contacts(output, expanded_contacts)


main.add_command(expand)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
