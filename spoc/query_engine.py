"""This file contains the classes making up the query engine."""
from typing import Any, Protocol, List, Union, Optional, Dict
import duckdb
import dask.dataframe as dd
import pandas as pd
from pydantic import BaseModel
from spoc.models.dataframe_models import (
    GenomicDataSchema,
    RegionSchema,
    QueryStepDataSchema,
)
from spoc.io import DUCKDB_CONNECTION


class GenomicData(Protocol):
    """Protocol for genomic data
    to be used in the query engine"""

    @property
    def data(self) -> Union[pd.DataFrame, duckdb.DuckDBPyRelation, dd.DataFrame]:
        """Return the data in the object"""
        ...

    def get_schema(self) -> GenomicDataSchema:
        """Return the schema of the underlying data"""
        ...


class QueryStep(Protocol):
    """Protocol for query steps"""

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the query step against the data schema"""

    def __call__(self, *args: Any, **kwds: Any) -> "QueryResult":
        """Apply the query step to the data"""


# TODO: think about allowing anchor composition
class Anchor(BaseModel):
    """Represents an anchor.

    Attributes:
        mode (str): The mode of the anchor.
        anchors (Optional[List[int]]): The list of anchor values (optional).
    """

    mode: str
    anchors: Optional[List[int]] = None

    def __repr__(self) -> str:
        return f"Anchor(mode={self.mode}, anchors={self.anchors})"

    def __str__(self) -> str:
        return self.__repr__()


class Snipper:
    """
    This class represents a snipper used for contact selection.
    It provides methods to validate the filter against a data schema,
    convert data to a duckdb relation, construct a filter string,
    and apply the filter to the data.

    Attributes:
        _regions (List[pd.DataFrame]): The regions to be used for filtering.
        _anchor_mode (Anchor): The anchor mode for filtering.
    """

    def __init__(self, regions: List[pd.DataFrame], anchor_mode: Anchor) -> None:
        # add ids to regions if they don't exist
        if "id" not in regions.columns:
            regions["id"] = range(len(regions))
        self._regions = RegionSchema.validate(regions)
        self._anchor_mode = anchor_mode

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the filter against the data schema"""
        # check whether an anchor is specified that is not in the data
        if self._anchor_mode.anchors is not None:
            if not all(
                anchor in data_schema.get_position_fields().keys()
                for anchor in self._anchor_mode.anchors
            ):
                raise ValueError(
                    "An anchor is specified that is not in the data schema."
                )

    def _convert_to_duckdb(
        self,
        data: Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation],
    ) -> duckdb.DuckDBPyRelation:
        """
        Converts the data to a duckdb relation.

        Parameters:
            data (Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation]): The input data to be converted.

        Returns:
            duckdb.DuckDBPyRelation: The converted duckdb relation.
        """
        if isinstance(data, dd.DataFrame):
            data = data.compute()
        return duckdb.from_df(data, connection=DUCKDB_CONNECTION)

    def _contstruct_filter(self, position_fields: Dict[int, List[str]]) -> str:
        """Constructs the filter string.

        Args:
            position_fields (List[str]): List of position fields.

        Returns:
            str: The constructed filter string.

        Raises:
            NotImplementedError: If the length of fields is not equal to 3.
        """
        query_strings = []
        join_string = " or " if self._anchor_mode.mode == "ANY" else " and "
        # subset on anchor regions
        if self._anchor_mode.anchors is not None:
            subset_positions = [
                position_fields[anchor] for anchor in self._anchor_mode.anchors
            ]
        else:
            subset_positions = list(position_fields.values())
        for fields in subset_positions:
            # check two or three-way definition of fields
            if len(fields) == 3:
                chrom, start, end = fields
                output_string = f"""(data.{chrom} = regions.chrom and
                                     (data.{start} between regions.start and regions.end or 
                                       data.{end} between regions.start and regions.end))"""
            else:
                raise NotImplementedError
            query_strings.append(output_string)
        return join_string.join(query_strings)

    def _get_transformed_schema(
        self, data_schema: GenomicDataSchema
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data.

        Args:
            data_schema (GenomicDataSchema): The input data schema.

        Returns:
            GenomicDataSchema: The schema of the transformed data.
        """
        # get columns of input schema
        input_columns = list(data_schema.get_schema().columns.keys())
        # add region columns
        region_columns = list(self._regions.columns)
        # construct schema
        return QueryStepDataSchema(
            columns=region_columns + input_columns,
            position_fields=data_schema.get_position_fields(),
            contact_order=data_schema.get_contact_order(),
        )

    def __repr__(self) -> str:
        return f"Snipper(anchor_mode={self._anchor_mode})"

    def __call__(self, genomic_data: GenomicData) -> GenomicData:
        """Apply the filter to the data"""
        # get input schema
        input_schema = genomic_data.get_schema()
        # bring input to duckdb dataframe
        if isinstance(genomic_data.data, duckdb.DuckDBPyRelation):
            genomic_data = genomic_data.data
        else:
            genomic_data = self._convert_to_duckdb(genomic_data.data)
        regions = self._convert_to_duckdb(self._regions)
        # get position columns and construct filter
        position_fields = input_schema.get_position_fields()
        # construct query
        df = genomic_data.set_alias("data").join(
            regions.set_alias("regions"), self._contstruct_filter(position_fields)
        )
        return QueryResult(df, self._get_transformed_schema(input_schema))


class MappedRegionFilter:
    """Filter for filtering mapped regions
    of contacts or pixels."""


class Aggregation:
    # TODO: think about how to specify aggregation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the aggregation against the data schema"""


class Transformation:
    # TODO: think about how to specify transformation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the transformation against the data schema"""


class QueryResult:
    """Result of a query"""

    def __init__(
        self,
        data: Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation],
        schema: GenomicDataSchema,
    ) -> None:
        self._data = data
        self._schema = schema

    @property
    def data(self) -> Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation]:
        """Returns the result as a dataframe object, either in memory or as a relation object"""
        return self._data

    def load_result(self) -> pd.DataFrame:
        """Loads the result into memory"""
        if isinstance(self._data, duckdb.DuckDBPyRelation):
            return self._data.to_df()
        elif isinstance(self._data, dd.DataFrame):
            return self._data.compute()
        return self._data

    def get_schema(self) -> GenomicDataSchema:
        """Returns the schema of the result"""
        return self._schema


class BasicQuery:
    def __init__(self, query_plan: List[QueryStep]) -> None:
        self._query_plan = query_plan

    def query(self, input_data: GenomicData) -> QueryResult:
        """Runs the query on the data and returns the result"""
        # instantiate query result
        query_result = QueryResult(input_data.data, input_data.get_schema())
        # run query
        for step in self._query_plan:
            # validate schema
            step.validate(query_result.get_schema())
            # apply step
            query_result = step(query_result)
        return query_result
