"""Console script for spoc."""
import sys

import click

from spoc.contacts import ContactManipulator
from spoc.contacts import Contacts
from spoc.fragments import FragmentAnnotator
from spoc.fragments import FragmentExpander
from spoc.io import FileManager
from spoc.models.dataframe_models import DataMode
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
    Functionality to expand labelled fragments to contacts. Expands n-way fragments over sequencing reads
    to yield contacts.

    Args:
        fragments_path (str): Path to the labelled fragments file.
        expanded_contacts_path (str): Path to the output contacts file.
        n_fragments (int, optional): Number of fragments per read to expand. Defaults to 3.

    """
    expander = FragmentExpander(number_fragments=n_fragments)
    file_manager = FileManager()
    input_fragments = file_manager.load_fragments(fragments_path)
    expanded = expander.expand(input_fragments)
    file_manager.write_contacts(expanded_contacts_path, expanded)


@click.command()
@click.argument("fragments_path")
@click.argument("label_library_path")
@click.argument("labelled_fragments_path")
def annotate(fragments_path, label_library_path, labelled_fragments_path):
    """
    Functionality to annotate porec fragments. Adds annotating labels and sister identity of mapped read fragments.

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
    """
    Functionality to bin contacts.  Bins higher order contacts into genomic bins of fixed size.
    Contacts path should be an URI.

    Args:
        contact_path (str): Path to the input contact file.
        pixel_path (str): Path to the output pixel file.
        bin_size (int, optional): Size of the bins. Defaults to 10000.
        same_chromosome (bool, optional): Only bin contacts on the same chromosome. Defaults to False.

    """
    # load data from disk
    file_manager = FileManager(DataMode.DASK)
    contacts = Contacts.from_uri(contact_path)
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
    """Functionality to merge annotated fragments. Concatenates contacts with the same
    configuration together and copies contacts with different configurations.

    Args:
        contact_paths (tuple): Paths to the input contact files.
        output (str, optional): Path to the output merged contact file.
    """
    file_manager = FileManager(DataMode.DASK)
    # get list of parameters
    parameters = [file_manager.list_contacts(p) for p in contact_paths]
    # get parameter counts -> if count > 1 then we need to concatenate
    parameter_counter = {}
    for file_index, parameter_list in enumerate(parameters):
        for parameter in parameter_list:
            if parameter not in parameter_counter:
                parameter_counter[parameter] = [file_index]
            else:
                parameter_counter[parameter].append(file_index)
    # iterate over parameters and write
    for parameter, file_indices in parameter_counter.items():
        if len(file_indices) == 1:
            file_index = file_indices[0]
            contacts = file_manager.load_contacts(contact_paths[file_index], parameter)
            file_manager.write_contacts(output, contacts)
        else:
            contact_files = [
                file_manager.load_contacts(contact_paths[i], parameter)
                for i in file_indices
            ]
            manipulator = ContactManipulator()
            merged = manipulator.merge_contacts(contact_files)
            file_manager.write_contacts(output, merged)


merge.add_command(merge_contacts)
main.add_command(expand)
main.add_command(annotate)
main.add_command(merge)
main.add_command(bin_contacts)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
