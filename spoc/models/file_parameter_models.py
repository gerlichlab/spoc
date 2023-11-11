"""This file contains data classes for parameters
of spoc data structures"""
from typing import Optional, Tuple, List
from pydantic import BaseModel, Field


class GlobalParameters(BaseModel):
    """Base class for global file parameters"""

    number_fragments: Optional[int] = None
    metadata_combi: Optional[Tuple[str, ...]] = None
    label_sorted: bool = False
    binary_labels_equal: bool = False
    symmetry_flipped: bool = False

    @classmethod
    def get_uri_fields(cls) -> List[str]:
        raise NotImplementedError

    def __hash__(self) -> int:
        # get metadata hash
        return hash(
            (
                self.number_fragments,
                self.metadata_combi,
                self.label_sorted,
                self.binary_labels_equal,
                self.symmetry_flipped,
            )
        )


class ContactsParameters(GlobalParameters):
    """Parameters for multiway contacts"""

    @classmethod
    def get_uri_fields(cls) -> List[str]:
        # Specific parameters needed to enforce order
        return [
            "number_fragments",
            "metadata_combi",
            "binary_labels_equal",
            "symmetry_flipped",
            "label_sorted",
        ]


class PixelParameters(GlobalParameters):
    """Parameters for genomic pixels"""

    binsize: Optional[int] = None
    same_chromosome: bool = True

    @classmethod
    def get_uri_fields(cls) -> List[str]:
        # Specific parameters needed to enforce order
        return [
            "number_fragments",
            "binsize",
            "metadata_combi",
            "binary_labels_equal",
            "symmetry_flipped",
            "label_sorted",
            "same_chromosome",
        ]

    def __hash__(self) -> int:
        return hash(
            (
                self.number_fragments,
                self.binsize,
                self.metadata_combi,
                self.label_sorted,
                self.binary_labels_equal,
                self.symmetry_flipped,
                self.same_chromosome,
            )
        )
