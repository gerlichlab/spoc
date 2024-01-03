"""Tests for the aggregation functions in the query engine."""
import pytest

from spoc.query_engine import AggregationFunction
from spoc.query_engine import BasicQuery
from spoc.query_engine import OffsetAggregation
from spoc.query_engine import RegionOffsetTransformation


@pytest.fixture(name="pixels_without_value_column")
def pixels_without_value_column_fixture(pixels_with_single_region):
    """Pixels with single region"""
    return RegionOffsetTransformation()(pixels_with_single_region)


@pytest.mark.parametrize(
    "genomic_data_fixture",
    [
        "contacts_without_regions",
        "pixels_without_regions",
        "contacts_with_single_region",
        "contacts_with_multiple_regions",
        "pixels_with_single_region",
        "pixels_with_multiple_regions",
    ],
)
def test_input_wo_offset_rejected(genomic_data_fixture, request):
    """Test that the validation fails for incorrect inputs."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    with pytest.raises(ValueError):
        query = BasicQuery(
            query_plan=[
                OffsetAggregation(
                    value_column="value", function=AggregationFunction.COUNT
                ),
            ],
        )
        query.query(genomic_data)


def test_input_wo_data_column_rejected(pixels_without_value_column):
    """Test that the validation fails for incorrect inputs."""
    with pytest.raises(ValueError):
        query = BasicQuery(
            query_plan=[
                OffsetAggregation(
                    value_column="test_column_does_no_exist",
                    function=AggregationFunction.AVG,
                ),
            ],
        )
        query.query(pixels_without_value_column)


def test_aggregation_succeds_on_correct_inputs():
    """Test that the aggregation succeeds on correct inputs."""
    pass
