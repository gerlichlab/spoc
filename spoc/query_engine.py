"""This file contains the classes making up the query engine."""
from typing import Any, Protocol, List, Union
import pandera as pa
import duckdb
import dask.dataframe as dd
import pandas as pd
from abc import ABC, abstractmethod


class GenomicData(Protocol):
    """Protocol for genomic data
    to be used in the query engine"""

    @property
    def data(self) -> Union[pd.DataFrame, duckdb.DuckDBPyRelation, dd.DataFrame]:
        """Return the data in the object"""
        ...

    def get_schema(self) -> pa.DataFrameSchema:
        """Return the schema of the underlying data"""
        ...


class QueryStep(Protocol):
    """Protocol for query steps"""

    def validate(self, data_schema: pa.DataFrameSchema) -> None:
        """Validate the query step against the data schema"""

    def __call__(self, *args: Any, **kwds: Any) -> "QueryResult":
        """Apply the query step to the data"""


class RegionFilter:
    # TODO: specify region schema in pandera
    def __init__(self, regions: List[pd.DataFrame]) -> None:
        self._regions = regions

    def validate(self, data_schema: pa.DataFrameSchema) -> None:
        """Validate the filter against the data schema"""


class Aggregation:
    # TODO: think about how to specify aggregation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: pa.DataFrameSchema) -> None:
        """Validate the aggregation against the data schema"""


class Transformation:
    # TODO: think about how to specify transformation
    def __init__(self, data: GenomicData) -> None:
        self._data = data

    def validate(self, data_schema: pa.DataFrameSchema) -> None:
        """Validate the transformation against the data schema"""


class QueryResult:
    """Result of a query"""

    def __init__(
        self,
        data: Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation],
        schema: pa.DataFrameSchema
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

    def get_schema(self) -> pa.DataFrameSchema:
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
