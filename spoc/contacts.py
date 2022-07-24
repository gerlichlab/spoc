"""Managing multi-way contacts."""

from itertools import combinations
import pandas as pd
import numpy as np
from .models import annotated_fragment_schema, HigherOrderContactSchema


class ContactExpander:
    """Expands n-way contacts over sequencing reads"""

    def __init__(self, number_fragments: int) -> None:
        self._number_fragments = number_fragments

    def _add_suffix(self, row, suffix):
        """expands contact fields"""
        output = {}
        for key in HigherOrderContactSchema.contact_fields:
            output[key + f"_{suffix}"] = getattr(row, key)
        return output

    def expand(self, fragments: pd.DataFrame) -> pd.DataFrame:
        # check input
        annotated_fragment_schema.validate(fragments)
        # expand fragments
        result = []
        keep_segments = fragments.query("pass_filter == True")
        for (read_name, read_df) in keep_segments.groupby("read_name", as_index=False):
            if len(read_df) < self._number_fragments:
                continue

            rows = list(
                read_df.sort_values(["read_start"], ascending=True)
                .assign(pos_on_read=lambda x: np.arange(len(x)))
                .itertuples()
            )

            read_length = rows[0].read_length
            for alignments in combinations(rows, self._number_fragments):
                contact = {"read_name": read_name, "read_length": read_length}
                # add reads
                for index, align in enumerate(alignments, start=1):
                    contact.update(self._add_suffix(align, index))
                result.append(contact)
        return pd.DataFrame(result)
