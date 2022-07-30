"""Console script for spoc."""
import sys
import click
from spoc.contacts import ContactExpander, ContactManipulator
from spoc.labels import FragmentAnnotator
from spoc.io import FileManager



@click.group()
def main(args=None):
    """Console script for spoc."""
    pass


@click.command()
@click.argument('fragments_path')
@click.argument('expanded_contacts_path')
@click.option(
    "-n",
    "--n_fragments",
    default=3,
    help="Number of fragments per read to expand",
)
def expand(fragments_path, expanded_contacts_path, n_fragments):
    """Script for expanding labelled fragments to contacts"""
    expander = ContactExpander(number_fragments=n_fragments)
    file_manager = FileManager(verify_schemas_on_load=True)
    input_fragments = file_manager.load_annotated_fragments(fragments_path)
    expanded = expander.expand(input_fragments)
    file_manager.write_multiway_contacts(expanded_contacts_path, expanded)

@click.command()
@click.argument('fragments_path')
@click.argument('label_library_path')
@click.argument('labelled_fragments_path')
def annotate(fragments_path, label_library_path, labelled_fragments_path):
    """Script for annotating porec fragments"""
    file_manager = FileManager(verify_schemas_on_load=True)
    label_library = file_manager.load_label_library(label_library_path)
    annotator = FragmentAnnotator(label_library)
    input_fragments = file_manager.load_porec_fragments(fragments_path)
    result = annotator.annotate_fragments(input_fragments)
    file_manager.write_annotated_fragments(labelled_fragments_path, result)

@click.group()
def merge():
    """Functionality to merge files"""

@click.command()
@click.argument('contact_paths', nargs=-1)
@click.option('-o', "--output", help="output path")
@click.option(
    "-n",
    "--n_fragments",
    default=3,
    help="Order of contacts",
)
def contacts(contact_paths, n_fragments, output):
    """Functionality to merge annotated fragments"""
    file_manager = FileManager(verify_schemas_on_load=True, use_dask=True)
    manipulator = ContactManipulator(n_fragments, use_dask=True)
    contacts= [file_manager.load_multiway_contacts(path, n_fragments) for path in contact_paths]
    merged = manipulator.merge_contacts(contacts)
    file_manager.write_multiway_contacts(output, merged)


merge.add_command(contacts)
main.add_command(expand)
main.add_command(annotate)
main.add_command(merge)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
