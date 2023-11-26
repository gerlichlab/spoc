"""This file contains the classes making up the query engine."""
from typing import Any, Protocol, List, Union, Optional
import pandera as pa
import duckdb
import dask.dataframe as dd
import pandas as pd
from pydantic import BaseModel
from spoc.models.dataframe_models import GenomicDataSchema

# Instantiate one duckdb connection to be used for all duckdb relations
DUCKDB_CONNECTION = duckdb.connect(database=":memory:")


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
    """Anchor mode for contact selection"""

    mode: str
    anchors: Optional[List[int]] = None


class Snipper:
    """Snipper for contact selection"""
    # TODO: specify region schema in pandera
    def __init__(self, regions: List[pd.DataFrame], anchor_mode: Anchor) -> None:
        self._regions = regions
        self._anchor_mode = anchor_mode

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the filter against the data schema"""

    def _convert_to_duckdb(
        self, data: Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation],
        table_name: str = "data"
    ) -> duckdb.DuckDBPyRelation:
        """Converts the data to a duckdb relation"""
        if isinstance(data, dd.DataFrame):
            data = data.compute()
        return duckdb.from_df(data, table_name, connection=DUCKDB_CONNECTION)

    def _contstruct_filter(self, position_fields: List[str]) -> str:
        """Constructs the filter string"""
        query_strings = []
        join_string = "or" if self._anchor_mode.mode == "ANY" else "and"
        # subset on anchor regions
        subset_positions = [
             position_fields[anchor] for anchor in self._anchor_mode.anchors
        ]
        for fields in subset_positions:
            # check two or three-way definition of fields
            if len(fields) == 3:
                chrom, start, end = fields
                output_string += f"({chrom} = regions.chrom and {start} >= regions.start and {end} <= regions.end)"
            else:
                raise NotImplementedError
            query_strings.append(output_string)
        return join_string.join(query_strings)


    def __call__(self, data: GenomicData) -> GenomicData:
        """Apply the filter to the data"""
        # bring input to duckdb dataframe
        if isinstance(data.data, duckdb.DuckDBPyRelation):
            data = data.data
        else:
            data = self._convert_to_duckdb(data.data)
        # TODO: handle ids
        regions = self._convert_to_duckdb(self._regions, table_name="regions")
        # get position columns and construct filter
        position_fields = data.get_schema().get_position_fields()
        # construct query
        df =  (
            data.set_alias("data")
            .join(regions.set_alias("regions"), 
                self._contstruct_filter(position_fields)
            )
        )
        


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
