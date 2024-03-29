"""This file contains the classes making up the query engine."""
from enum import Enum
from itertools import product
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional
from typing import Protocol
from typing import Tuple
from typing import TypeVar
from typing import Union

import dask.dataframe as dd
import duckdb
import numpy as np
import pandas as pd
from pydantic import BaseModel

from spoc.io import DUCKDB_CONNECTION
from spoc.models.dataframe_models import GenomicDataSchema
from spoc.models.dataframe_models import QueryStepDataSchema
from spoc.models.dataframe_models import RegionSchema


T = TypeVar("T")


def convert_string_to_enum(enum_class: Callable[[str], T], string: str) -> T:
    """Converts a string to an enum value"""
    try:
        return enum_class(string.upper())
    except ValueError as exc:
        raise ValueError(f"Invalid value for {enum_class.__name__}: {string}") from exc


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

    def __call__(self, *args: Any, **kwds: Any) -> "QueryPlan":
        """Apply the query step to the data"""


# TODO: think about allowing anchor composition
class Anchor(BaseModel):
    """Represents an anchor.

    Attributes:
        mode (str): The mode of the anchor. (Can be "ANY" or "ALL")
        anchors (Optional[List[int]]): The list of anchor values (optional).
    """

    mode: str
    anchors: Optional[List[int]] = None

    def __repr__(self) -> str:
        return f"Anchor(mode={self.mode}, anchors={self.anchors})"

    def __str__(self) -> str:
        return self.__repr__()


class Overlap:
    """
    This class represents an overlap calculation used for contact and pixel selection.
    It provides methods to validate the filter against a data schema,
    convert data to a duckdb relation, construct a filter string,
    and apply the filter to the data.
    """

    def __init__(
        self,
        regions: pd.DataFrame,
        anchor_mode: Union[Anchor, Tuple[str, List[int]]],
        half_window_size: Optional[int] = None,
    ) -> None:
        """
        Initialize the Overlap object.

        Args:
            regions (pd.DataFrame): A DataFrame containing the regions data.
            anchor_mode (Union[Anchor,Tuple[str,List[int]]]): The anchor mode to be used.
            half_window_size (Optional[int]): The window size the regions should be expanded to. Defaults to None and is inferred from the data.

        Returns:
            None
        """
        # add ids to regions if they don't exist
        if "id" not in regions.columns:
            regions["id"] = range(len(regions))
        if half_window_size is not None:
            expanded_regions = regions.copy()
            # create midpoint
            expanded_regions["midpoint"] = (
                expanded_regions["start"] + expanded_regions["end"]
            ) // 2
            # expand regions
            expanded_regions["start"] = expanded_regions["midpoint"] - half_window_size
            expanded_regions["end"] = expanded_regions["midpoint"] + half_window_size
            # drop midpoint
            expanded_regions = expanded_regions.drop(columns=["midpoint"])
            self._regions = RegionSchema.validate(
                expanded_regions.add_prefix("region_")
            )
            self._half_window_size = half_window_size
        else:
            self._regions = RegionSchema.validate(regions.add_prefix("region_"))
            # infer window size -> variable regions will have largest possible window size
            self._half_window_size = int(
                (self._regions["region_end"] - self._regions["region_start"]).max() // 2
            )
        if isinstance(anchor_mode, tuple):
            self._anchor_mode = Anchor(mode=anchor_mode[0], anchors=anchor_mode[1])
        else:
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
        data: Union[pd.DataFrame, dd.DataFrame],
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
                                    (
                                        data.{start} between regions.region_start and regions.region_end or 
                                        data.{end} between regions.region_start and regions.region_end or
                                        regions.region_start between data.{start} and data.{end}
                                    )
                                )"""
            query_strings.append(output_string)
        return join_string.join(query_strings)

    def _get_transformed_schema(
        self,
        data_frame: duckdb.DuckDBPyRelation,
        input_schema: GenomicDataSchema,
        position_fields: Dict[int, List[str]],
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data."""
        # construct schema
        return QueryStepDataSchema(
            columns=data_frame.columns,
            position_fields=position_fields,
            contact_order=input_schema.get_contact_order(),
            binsize=input_schema.get_binsize(),
            region_number=len(self._regions),
            half_window_size=self._half_window_size,
        )

    def _add_end_position(
        self,
        data_frame: duckdb.DuckDBPyRelation,
        bin_size: Optional[int],
        position_fields: Dict[int, List[str]],
    ) -> duckdb.DuckDBPyRelation:
        """Adds an end position column to the dataframe"""
        position_columns = {j for i in position_fields.values() for j in i}
        non_position_columns = [
            column for column in data_frame.columns if column not in position_columns
        ]
        new_position_clause = [
            f"data.chrom as chrom_{position}, data.start_{position} as start_{position}, data.start_{position} + {bin_size} as end_{position}"
            for position in position_fields.keys()
        ]
        # add end position
        return data_frame.set_alias("data").project(
            ",".join(new_position_clause + non_position_columns)
        )

    def __repr__(self) -> str:
        return f"Overlap(anchor_mode={self._anchor_mode})"

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
        snipped_df = genomic_df.set_alias("data").join(
            regions.set_alias("regions"), self._contstruct_filter(position_fields)
        )
        return QueryPlan(
            snipped_df,
            self._get_transformed_schema(snipped_df, input_schema, position_fields),
        )


class AggregationFunction(Enum):
    """Enum for aggregation functions.
    Options are:
        SUM: Sum of values.
        AVG_WITH_EMPTY: Average of values, empty values are counted as 0.
        AVG: Average of values, empty values are not counted.
        COUNT: Number of values.
    """

    SUM: str = "SUM"
    AVG_WITH_EMPTY: str = "AVG_WITH_EMPTY"
    AVG: str = "AVG"
    COUNT: str = "COUNT"


class DistanceAggregation:
    """Aggregation based on distances from a region. Uses all available distances."""

    def __init__(
        self,
        value_column: str,
        function: Union[AggregationFunction, str] = AggregationFunction.AVG,
        densify_output: bool = True,
        position_list: Optional[List[int]] = None,
    ) -> None:
        """Initialize the aggregation.

        Args:
            value_column (str): The name of the column to be aggregated.
            function (Union[AggregationFunction,str]): The aggregation function to be applied. Defaults to AggregationFunction.AVG.
            densify_output (bool, optional): Whether to densify the output. Defaults to True.
                                             This requires a binsize value to be set in the data schema.
            position_list (Optional[List[int]]): The list of positions to use for aggregations, starting with 1. Defaults to using all positions.
        """
        if isinstance(function, str):
            parsed_function = convert_string_to_enum(AggregationFunction, function)
        else:
            parsed_function = function
        self._function = parsed_function
        self._value_column = value_column
        self._densify_output = densify_output
        self._position_list = position_list

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the aggregation against the data schema"""
        # check that at leastl one distance field is present
        if "distance_1" not in data_schema.get_schema().columns:
            raise ValueError("No distance fields in data schema.")
        # check that all position fields are present
        if self._position_list is not None:
            for position in self._position_list:
                if position not in data_schema.get_position_fields():
                    raise ValueError(f"Position {position} not in data schema.")
        # check that value column is present
        if self._value_column not in data_schema.get_schema().columns:
            raise ValueError("Value column not in data schema.")
        # check for binsize -> only pixels have that
        if data_schema.get_binsize() is None:
            raise ValueError("No binsize specified in data schema.")
        # check for window size
        if data_schema.get_half_window_size() is None:
            raise ValueError("No window size specified in data schema.")

    def _get_transformed_schema(
        self,
        data_frame: duckdb.DuckDBPyRelation,
        input_schema: GenomicDataSchema,
        position_fields: Dict[int, List[str]],
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data."""
        # construct schema
        return QueryStepDataSchema(
            columns=data_frame.columns,
            position_fields=position_fields,
            contact_order=len(position_fields),
            binsize=input_schema.get_binsize(),
        )

    def _aggregate_distances(
        self,
        data_frame: duckdb.DuckDBPyRelation,
        input_schema: GenomicDataSchema,
        position_fields: Dict[int, List[str]],
    ) -> duckdb.DuckDBPyRelation:
        """Aggregates the distances."""
        # get distance columns
        distance_columns = [
            f"distance_{position}" for position in position_fields.keys()
        ]
        # construct aggregation
        if self._function == AggregationFunction.COUNT:
            aggregation_string = (
                f"COUNT(*) as {self._value_column}_{self._function.name.lower()}"
            )
        elif self._function == AggregationFunction.AVG_WITH_EMPTY:
            # For average, we need to sum up the values and divide by the number of regions
            aggregation_string = f"SUM({self._value_column})::float/{input_schema.get_region_number()} as {self._value_column}_{self._function.name.lower()}"
        else:
            aggregation_string = f"{self._function.name}({self._value_column}) as {self._value_column}_{self._function.name.lower()}"
        data_frame = (
            data_frame.set_alias("data")
            .aggregate(
                ",".join(distance_columns + [aggregation_string]),
            )
            .order(",".join(distance_columns))
        )
        return data_frame

    def _create_empty_dense_output(
        self,
        input_schema: GenomicDataSchema,
        position_fields: Dict[int, List[str]],
    ) -> duckdb.DuckDBPyRelation:
        """Create dense value columns for all distances."""
        binsize: Optional[int] = input_schema.get_binsize()
        if binsize is None:
            raise ValueError("No binsize specified in data schema.")
        int_binsize: int = binsize
        windowsize: Optional[int] = input_schema.get_half_window_size()
        if windowsize is None:
            raise ValueError("No window size specified in data schema.")
        int_windowsize: int = windowsize
        # create combinations of distances
        distance_combinations = pd.DataFrame(
            product(
                np.arange(
                    -(np.floor(int_windowsize / int_binsize) * int_binsize),
                    (np.floor(int_windowsize / int_binsize) * int_binsize) + 1,
                    int_binsize,
                ),
                repeat=len(position_fields.keys()),
            ),
            columns=[f"distance_{i}" for i in position_fields.keys()],
        )
        # fill value
        if self._function in (
            AggregationFunction.COUNT,
            AggregationFunction.SUM,
            AggregationFunction.AVG_WITH_EMPTY,
        ):
            distance_combinations["fill_value"] = 0
        else:
            distance_combinations["fill_value"] = np.nan
        return duckdb.from_df(distance_combinations, connection=DUCKDB_CONNECTION)

    def _fill_empty_output(
        self,
        data_frame: duckdb.DuckDBPyRelation,
        empty_dense_output: duckdb.DuckDBPyRelation,
        position_fields: Dict[int, List[str]],
    ) -> duckdb.DuckDBPyRelation:
        """Fill empty output with values from dense output."""
        # get distance columns
        distance_columns = [f"distance_{i}" for i in position_fields.keys()]
        # construct join and coalesce output
        data_frame = (
            data_frame.set_alias("data")
            .join(
                empty_dense_output.set_alias("empty_dense_output"),
                ",".join(distance_columns),
                how="right",
            )
            .project(
                ",".join(distance_columns)
                + f", COALESCE(data.{self._value_column}_{self._function.name.lower()}, empty_dense_output.fill_value) as {self._value_column}"
            )
            .set_alias("filled")
            .order(",".join([f"filled.{col}" for col in distance_columns]))
        )
        return data_frame

    def __call__(self, genomic_data: GenomicData) -> GenomicData:
        """Apply the aggregation to the data"""
        # get input schema
        input_schema = genomic_data.get_schema()
        # bring input to duckdb dataframe
        if isinstance(genomic_data.data, duckdb.DuckDBPyRelation):
            genomic_df = genomic_data.data
        else:
            genomic_df = duckdb.from_df(genomic_data.data, connection=DUCKDB_CONNECTION)
        # get position columns
        position_fields = input_schema.get_position_fields()
        if self._position_list is not None:
            position_fields = {
                position: position_fields[position] for position in self._position_list
            }
        # construct transformation
        aggregated_data = self._aggregate_distances(
            genomic_df, input_schema, position_fields
        )
        if self._densify_output:
            empty_dense_output = self._create_empty_dense_output(
                input_schema, position_fields
            )
            aggregated_data = self._fill_empty_output(
                aggregated_data, empty_dense_output, position_fields
            )
        return QueryPlan(
            aggregated_data,
            self._get_transformed_schema(
                aggregated_data, input_schema, position_fields
            ),
        )


class DistanceMode(Enum):
    """Enum for distance modes."""

    LEFT: str = "LEFT"
    RIGHT: str = "RIGHT"
    BOTH: str = "BOTH"
    MIDPOINT: str = "MIDPOINT"


class DistanceTransformation:
    """Adds distance columns for each position field relative
    to required region_columns."""

    def __init__(
        self, distance_mode: Union[DistanceMode, str] = DistanceMode.LEFT
    ) -> None:
        """Initialize the transformation.

        Args:
            distance_mode (Union[DistanceMode,str]): The distance mode to be used. Defaults to DistanceMode.MIDPOINT.
                                      Specifies how the distance is calculated relative to the region midpoint.
                                      Note that the distance is always calculated relative to the midpoint of the region.
                                      If a binsize is specificed in the data schema, this needs to be set to LEFT.
        """
        if isinstance(distance_mode, str):
            distance_mode = convert_string_to_enum(DistanceMode, distance_mode)
        self._distance_mode = distance_mode

    def validate(self, data_schema: GenomicDataSchema) -> None:
        """Validate the transformation against the data schema"""
        # check that there are position fields and region columns
        if not data_schema.get_position_fields():
            raise ValueError("No position fields in data schema.")
        schema_columns = data_schema.get_schema().columns
        required_columns = ["region_chrom", "region_start", "region_end"]
        if not all(column in schema_columns for column in required_columns):
            raise ValueError("No region columns in data schema.")
        if (
            self._distance_mode != DistanceMode.LEFT
            and data_schema.get_binsize() is not None
        ):
            raise ValueError(
                "Binsize specified in data schema, but distance mode is not set to LEFT."
            )

    def _create_transform_columns(
        self, genomic_df: duckdb.DuckDBPyRelation, input_schema: GenomicDataSchema
    ) -> duckdb.DuckDBPyRelation:
        """Creates the transform columns for the given position fields"""
        # position fields
        position_fields = input_schema.get_position_fields()
        # get existing columns
        transform_strings = [f"data.{column}" for column in genomic_df.columns]
        # check whether binsize is specified
        if input_schema.get_binsize() is not None:
            binsize = input_schema.get_binsize()
        else:
            binsize = 1
        # create transform columns
        for position_field, fields in position_fields.items():
            _, start, end = fields
            if self._distance_mode == DistanceMode.MIDPOINT:
                output_string = f"""(FLOOR((data.{start} + data.{end})/2) - FLOOR((data.region_start + data.region_end)/2))
                                         as distance_{position_field}"""
            if self._distance_mode == DistanceMode.LEFT:
                output_string = f"""data.{start} - FLOOR((FLOOR(data.region_start/{binsize}) * {binsize}
                                    + FLOOR(data.region_end/{binsize}) * {binsize})/2) as distance_{position_field}"""
            if self._distance_mode == DistanceMode.RIGHT:
                output_string = f"""data.{end} - FLOOR((data.region_start + data.region_end)/2) as distance_{position_field}"""
            if self._distance_mode == DistanceMode.BOTH:
                output_string = f"""data.{start} - FLOOR((data.region_start + data.region_end)/2) as start_distance_{position_field},
                                    data.{end} - FLOOR((data.region_start + data.region_end)/2) as end_distance_{position_field}"""
            transform_strings.append(output_string)
        return ",".join(transform_strings)

    def _get_transformed_schema(
        self, data_frame: duckdb.DuckDBPyRelation, input_schema: GenomicDataSchema
    ) -> GenomicDataSchema:
        """Returns the schema of the transformed data."""
        # construct schema
        return QueryStepDataSchema(
            columns=data_frame.columns,
            position_fields=input_schema.get_position_fields(),
            contact_order=input_schema.get_contact_order(),
            binsize=input_schema.get_binsize(),
            region_number=input_schema.get_region_number(),
            half_window_size=input_schema.get_half_window_size(),
        )

    def __call__(self, genomic_data: GenomicData) -> GenomicData:
        """Apply the transformation to the data"""
        # get input schema
        input_schema = genomic_data.get_schema()
        # bring input to duckdb dataframe
        if isinstance(genomic_data.data, duckdb.DuckDBPyRelation):
            genomic_df = genomic_data.data
        else:
            genomic_df = duckdb.from_df(genomic_data.data, connection=DUCKDB_CONNECTION)
        # construct transformation
        transformed_df = genomic_df.set_alias("data").project(
            self._create_transform_columns(genomic_df, input_schema)
        )
        return QueryPlan(
            transformed_df, self._get_transformed_schema(transformed_df, input_schema)
        )


class QueryPlan:
    """Result of a query"""

    def __init__(
        self,
        data: Union[pd.DataFrame, duckdb.DuckDBPyRelation],
        schema: GenomicDataSchema,
    ) -> None:
        self._data = data
        self._schema = schema

    @property
    def data(self) -> Union[pd.DataFrame, duckdb.DuckDBPyRelation]:
        """Returns the result as a dataframe object, either in memory or as a relation object"""
        return self._data

    def compute(self) -> pd.DataFrame:
        """Loads the result into memory"""
        if isinstance(self._data, duckdb.DuckDBPyRelation):
            return self._data.to_df()
        return self._data

    def get_schema(self) -> GenomicDataSchema:
        """Returns the schema of the result"""
        return self._schema


# pylint: disable=too-few-public-methods
# this is a wrapper with one task, so it only has one method
class Query:
    """Basic query engine that runs a query plan on the data"""

    def __init__(self, query_steps: List[QueryStep]) -> None:
        self._query_steps = query_steps

    def build(self, input_data: GenomicData) -> QueryPlan:
        """Runs the query on the data and returns the result"""
        # instantiate query result
        query_plan = QueryPlan(input_data.data, input_data.get_schema())
        # run query
        for step in self._query_steps:
            # validate schema
            step.validate(query_plan.get_schema())
            # apply step
            query_plan = step(query_plan)
        return query_plan
