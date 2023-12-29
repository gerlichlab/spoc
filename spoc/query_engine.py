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

    def get_schema(self) -> GenomicDataSchema:
        """Return the schema of the underlying data"""


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
        _regions (pd.DataFrame): The regions to be used for filtering.
        _anchor_mode (Anchor): The anchor mode for filtering.
    """

    def __init__(self, regions: pd.DataFrame, anchor_mode: Anchor) -> None:
        # add ids to regions if they don't exist
        if "id" not in regions.columns:
            regions["id"] = range(len(regions))
        self._regions = RegionSchema.validate(regions.add_prefix("region_"))
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
            chrom, start, end = fields
            output_string = f"""(data.{chrom} = regions.region_chrom and
                                    (data.{start} between regions.region_start and regions.region_end or 
                                    data.{end} between regions.region_start and regions.region_end))"""
            query_strings.append(output_string)
        return join_string.join(query_strings)

    def _get_transformed_schema(
        self,
        df: duckdb.DuckDBPyRelation,
        input_schema: GenomicDataSchema,
        position_fields: Dict[int, List[str]],
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data."""
        # construct schema
        return QueryStepDataSchema(
            columns=df.columns,
            position_fields=position_fields,
            contact_order=input_schema.get_contact_order(),
        )

    def _add_end_position(
        self,
        df: duckdb.DuckDBPyRelation,
        bin_size: int,
        position_fields: Dict[int, List[str]],
    ) -> duckdb.DuckDBPyRelation:
        """Adds an end position column to the dataframe"""
        position_columns = set([j for i in position_fields.values() for j in i])
        non_position_columns = [
            column for column in df.columns if column not in position_columns
        ]
        new_position_clause = [
            f"data.chrom as chrom_{position}, data.start_{position} as start_{position}, data.start_{position} + {bin_size} as end_{position}"
            for position in position_fields.keys()
        ]
        # add end position
        return df.set_alias("data").project(
            ",".join(new_position_clause + non_position_columns)
        )

    def __repr__(self) -> str:
        return f"Snipper(anchor_mode={self._anchor_mode})"

    def __call__(self, genomic_data: GenomicData) -> GenomicData:
        """Apply the filter to the data"""
        # get input schema
        input_schema = genomic_data.get_schema()
        # bring input to duckdb dataframe
        if isinstance(genomic_data.data, duckdb.DuckDBPyRelation):
            genomic_df = genomic_data.data
        else:
            genomic_df = self._convert_to_duckdb(genomic_data.data)
        regions = self._convert_to_duckdb(self._regions)
        # get position columns and construct filter
        position_fields = input_schema.get_position_fields()
        # add end position if not present
        if len(position_fields[1]) == 2:
            genomic_df = self._add_end_position(
                genomic_df, input_schema.get_binsize(), position_fields
            )
            position_fields = {
                position: [f"chrom_{position}", f"start_{position}", f"end_{position}"]
                for position in position_fields.keys()
            }
        # construct query
        df = genomic_df.set_alias("data").join(
            regions.set_alias("regions"), self._contstruct_filter(position_fields)
        )
        return QueryResult(
            df, self._get_transformed_schema(df, input_schema, position_fields)
        )


class MappedRegionFilter:
    """Filter for filtering mapped regions
    of contacts or pixels."""


class Aggregation:
    """Aggregation of contacts or pixels"""

    # TODO: think about how to specify aggregation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the aggregation against the data schema"""


class Transformation:
    """Transformation of contacts or pixels"""

    # TODO: think about how to specify transformation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the transformation against the data schema"""


class RegionOffsetTransformation:
    """Adds offset columns for each position field relative
    to required region_columns."""

    def __init__(self, use_mid_point: bool = True) -> None:
        """Initialize the transformation."""
        self._use_mid_point = use_mid_point

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the transformation against the data schema"""
        # check that there are position fields and region columns
        if not data_schema.get_position_fields():
            raise ValueError("No position fields in data schema.")
        schema_columns = data_schema.get_schema().columns
        required_columns = ["region_chrom", "region_start", "region_end"]
        if not all(column in schema_columns for column in required_columns):
            raise ValueError("No region columns in data schema.")

    def _create_transform_columns(
        self, genomic_df: duckdb.DuckDBPyRelation, position_fields: Dict[int, List[str]]
    ) -> duckdb.DuckDBPyRelation:
        """Creates the transform columns for the given position fields"""
        # get existing columns
        transform_strings = [f"data.{column}" for column in genomic_df.columns]
        # create transform columns
        for position_field, fields in position_fields.items():
            _, start, end = fields
            if self._use_mid_point:
                output_string = f"""(FLOOR((data.{start} + data.{end})/2) - FLOOR((data.region_start + data.region_end)/2))
                                         as offset_{position_field}"""
            else:
                output_string = f"""data.{start} - data.region_start as start_offset_{position_field},
                                    data.{end} - data.region_end as end_offset_{position_field}"""
            transform_strings.append(output_string)
        return ",".join(transform_strings)

    def _get_transformed_schema(
        self, df: duckdb.DuckDBPyRelation, input_schema: GenomicDataSchema
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data."""
        # construct schema
        return QueryStepDataSchema(
            columns=df.columns,
            position_fields=input_schema.get_position_fields(),
            contact_order=input_schema.get_contact_order(),
        )

    def __call__(self, genomic_data: GenomicData) -> Any:
        """Apply the transformation to the data"""
        # get input schema
        input_schema = genomic_data.get_schema()
        # bring input to duckdb dataframe
        if isinstance(genomic_data.data, duckdb.DuckDBPyRelation):
            genomic_df = genomic_data.data
        else:
            genomic_df = duckdb.from_df(genomic_data.data, connection=DUCKDB_CONNECTION)
        # get position columns
        position_fields = input_schema.get_position_fields()
        # construct transformation
        df = genomic_df.set_alias("data").project(
            self._create_transform_columns(genomic_df, position_fields)
        )
        return QueryResult(df, self._get_transformed_schema(df, input_schema))


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
        if isinstance(self._data, dd.DataFrame):
            return self._data.compute()
        return self._data

    def get_schema(self) -> GenomicDataSchema:
        """Returns the schema of the result"""
        return self._schema


class BasicQuery:
    """Basic query engine that runs a query plan on the data"""

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
