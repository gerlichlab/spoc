# Interface query

This document defines the high-level interface of the spoc query engine.

## Class relationships

```mermaid
classDiagram
namespace genomicData {
    class Pixels
    class Contacts
    } 

    class gDataProtocol{
        <<Protocol>>
        +DataFrame data
        get_schema(): pandera.Schema
    }

    class BasicQuery
        BasicQuery: +List~Filter|Aggregation|Transformation~ query_plan
        BasicQuery: +compose_with(BasicQuery) BasicQuery
        BasicQuery: +query(gDataProtocol)

    BasicQuery --> QueryResult : query result
    BasicQuery "1" --* "*" Filter
    BasicQuery "1" --* "*" Aggregation
    BasicQuery "1" --* "*" Transformation
    gDataProtocol --> BasicQuery : takes input
    gDataProtocol ..|> Pixels : realization
    gDataProtocol ..|> Contacts : realization
    gDataProtocol ..|> QueryResult : realization
    Filter --|> RegionFilter

    class Filter {
        <<Interface>>
        +gDataProtocol data
        +get_row_subset() pd.Series~gDataIds~
    }

    class RegionFilter {
        +gDataProtocol data
        +pd.DataFrame~Regions~
        +String anchor_mode
        +get_row_subset() pd.Series~gDataIds~
    }
        

    class Aggregation

    class Transformation

    class QueryResult
        QueryResult: +DataFrame data
        QueryResult: get_schema() pandera.Schema
        QueryResult: compute() gDataProtocol
```

## Description

- __gDataProtocol__: Protocol class that defines the interface of genomic data that can be accepted by `BasicQuery`. Implements a method to get it's schema as well as a parameter to get the underlying data
- __BasicQuery__: Central query class that encapsulates querying an object that implements the `gDataProtocol`. Holds references to a query plan, which is a list of filters, aggregations and transformations that are executed in order and specify the filtering, aggregation and transformation operations. Is composable with other basic query instances to capture more complex queries. Performs checks on the proposed operations based on the `get_schema()` method and the requested filters and aggregations.
- __QueryResult__: Result of a BasicQuery that implements the `gDataProtocol` and can either be computed, which manifests the query in memory, or passed to basic query again.
- __Filter__: Interface of a filter that is accepted by `BasicQuery` and encapsulates filtering along rows of genomic data.
- __RegionFilter__: A filter that filters for overlap with specific genomic regions that are passed to the constructor. Anchormode refers to the way that the genomic regions are overlapped (e.g. at least one, exactly one, all, the first contact etc.)