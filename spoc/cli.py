"""Console script for spoc."""
import sys
import click
from spoc.contacts import ContactManipulator
from spoc.fragments import FragmentAnnotator, FragmentExpander
from spoc.io import FileManager
from spoc.pixels import GenomicBinner


@click.group()
def main():
    """Console script for spoc."""


@click.command()
@click.argument("fragments_path")
@click.argument("expanded_contacts_path")
@click.option(
    "-n",
    "--n_fragments",
    default=3,
    help="Number of fragments per read to expand",
)
def expand(fragments_path, expanded_contacts_path, n_fragments):
    """
    Script for expanding labelled fragments to contacts

    Args:
        fragments_path (str): Path to the labelled fragments file.
        expanded_contacts_path (str): Path to the output contacts file.
        n_fragments (int, optional): Number of fragments per read to expand. Defaults to 3.

    """
    expander = FragmentExpander(number_fragments=n_fragments)
    file_manager = FileManager()
    input_fragments = file_manager.load_fragments(fragments_path)
    expanded = expander.expand(input_fragments)
    file_manager.write_multiway_contacts(expanded_contacts_path, expanded)


@click.command()
@click.argument("fragments_path")
@click.argument("label_library_path")
@click.argument("labelled_fragments_path")
def annotate(fragments_path, label_library_path, labelled_fragments_path):
    """Script for annotating porec fragments

    Args:
        fragments_path (str): Path to the input fragments file.
        label_library_path (str): Path to the label library file.
        labelled_fragments_path (str): Path to the output labelled fragments file.

    """
    file_manager = FileManager()
    label_library = file_manager.load_label_library(label_library_path)
    annotator = FragmentAnnotator(label_library)
    input_fragments = file_manager.load_fragments(fragments_path)
    result = annotator.annotate_fragments(input_fragments)
    file_manager.write_fragments(labelled_fragments_path, result)


@click.command()
@click.argument("contact_path")
@click.argument("pixel_path")
@click.option("-b", "--bin_size", default=10_000, type=int)
@click.option("-c", "--same_chromosome", is_flag=True)
def bin_contacts(
    contact_path,
    pixel_path,
    bin_size,
    same_chromosome,
):
    """Script for binning contacts

    Args:
        contact_path (str): Path to the input contact file.
        pixel_path (str): Path to the output pixel file.
        bin_size (int, optional): Size of the bins. Defaults to 10000.
        same_chromosome (bool, optional): Only bin contacts on the same chromosome. Defaults to False.

    """
    # load data from disk
    file_manager = FileManager(use_dask=True)
    contacts = file_manager.load_contacts(contact_path)
    # binning
    binner = GenomicBinner(bin_size=bin_size)
    pixels = binner.bin_contacts(contacts, same_chromosome=same_chromosome)
    # persisting
    file_manager.write_pixels(pixel_path, pixels)


@click.group()
def merge():
    """Functionality to merge files"""


@click.command(name="contacts")
@click.argument("contact_paths", nargs=-1)
@click.option("-o", "--output", help="output path")
def merge_contacts(contact_paths, output):
    """Functionality to merge annotated fragments

    Args:
        contact_paths (tuple): Paths to the input contact files.
        output (str, optional): Path to the output merged contact file.
    """
    file_manager = FileManager(use_dask=True)
    manipulator = ContactManipulator()
    contact_files = [file_manager.load_contacts(path) for path in contact_paths]
    merged = manipulator.merge_contacts(contact_files)
    file_manager.write_multiway_contacts(output, merged)


merge.add_command(merge_contacts)
main.add_command(expand)
main.add_command(annotate)
main.add_command(merge)
main.add_command(bin_contacts)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
