"""This file contains data classes for parameters
of spoc data structures"""
from typing import Optional, List
from pydantic import BaseModel


class ContactsParameters(BaseModel):
    """Parameters for multiway contacts"""

    number_fragments: Optional[int] = None
    metadata_combi: Optional[List[str]] = None
    label_sorted: bool = False
    binary_labels_equal: bool = False
    symmetry_flipped: bool = False

    def __hash__(self) -> int:
        # get metadata hash
        if self.metadata_combi is not None:
            metadata_hash = hash(tuple(self.metadata_combi))
        else:
            metadata_hash = hash(None)
        return hash(
            (
                self.number_fragments,
                metadata_hash,
                self.label_sorted,
                self.binary_labels_equal,
                self.symmetry_flipped,
            )
        )


class PixelParameters(BaseModel):
    """Parameters for genomic pixels"""

    number_fragments: Optional[int] = None
    binsize: Optional[int] = None
    metadata_combi: Optional[List[str]] = None
    label_sorted: bool = False
    binary_labels_equal: bool = False
    symmetry_flipped: bool = False
    same_chromosome: bool = True

    def __hash__(self) -> int:
        # get metadata hash
        if self.metadata_combi is not None:
            metadata_hash = hash(tuple(self.metadata_combi))
        else:
            metadata_hash = hash(None)
        return hash(
            (
                self.number_fragments,
                self.binsize,
                metadata_hash,
                self.label_sorted,
                self.binary_labels_equal,
                self.symmetry_flipped,
                self.same_chromosome,
            )
        )
