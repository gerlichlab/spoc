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


class PixelParameters(BaseModel):
    """Parameters for genomic pixels"""

    number_fragments: Optional[int] = None
    binsize: Optional[int] = None
    metadata_combi: Optional[List[str]] = None
    label_sorted: bool = False
    binary_labels_equal: bool = False
    symmetry_flipped: bool = False
    same_chromosome: bool = True
