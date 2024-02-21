# Query engine
This technical document describes the spoc query engine, a set of classes that implements spoc's interface for querying multi-dimensional genomic data.

## Principles

### Composable pieces
Spoc's query engine consists of composable pieces that can be combined to produce an expressive query language. These pieces represent basic operations on genomic data that are easily implemented and understood on their own. This allows a great degree of flexibility, while also allowing predefined recipes that less experienced users can get started with.

### Lazy evaluation
The spoc query engine is designed with lazy evaluation as a guiding principle. This means that data queries are only executed when they are needed to minimize loading data into memory and computational overhead. To enable this, spoc queries have a construction phase, which specifies the operations to be executed and an execution phase, that actually executes the query.

## Query plans and query steps

The most important ingredient in this query language is a class that implements the `QueryStep` protocol. This protocol serves two purposes:

- It exposes a way to validate the data schema during query building
- It implements adding itself to a query

This way, query steps can be combined into a query plan that specifies the analysis to be executed. In generaly, spoc supports the following types of query steps:

- **Overlap**: Overlaps genomic data with a set of genomic intervals
- **Aggregation**: Aggregates the data passed into the query using an aggregation function
- **Transform**: Adds columns to query based on a fixed or custom computation

Specific examples of query steps are:

- **DistanceTransformation**: Example of a transfomration. Adds distance of genomic positions to regions added by Overlap
- **DistanceAggregation**: Exaple of an aggregation. Aggregates the distances to genomic regions using an aggregation function.

### Input and output of query steps

A query step takes as input a class that implements the `GenomicData` protocol. This protocol allows retrievel of the data schema (a thin wrapper over a pandera dataframe schema) as well as the data itself. The output of a query step is again a class that ipmlements the `GenomicData` protocol to allow composition. Specific examples of possible inputs are:

- **Pixels**: Represents input pixels
- **Contacts**: Represents input contacts
- **QueryPlan**: The result of a query step

### Composition of query steps

To allow specifying complex queries, query steps need to be combined. This is done using the `Query` class. It takes a query plan (a list of `QueryStep` instances) as input, exposes the `build` method, which takes input data, validates all query steps and adds them to the resulting `QueryPlan` instance that is returned.

### Manifestation of results

So far, we have only talked about specifying the query to be executed, but not how to actually execute it. A `QueryPlan` has a `compute()` method that returns the manifested dataframe as a `pd.DataFrame` instance. This is the step that actually executes the specified query.

## Examples

### Selecting a subset of contacts at a single genomic position
In this example, we want to select a subset of genomic contacts at a single location. For this, we first load the required input data:


```python
from spoc.query_engine import Overlap, Anchor, Query
from spoc.contacts import Contacts
import pandas as pd

contacts = Contacts.from_uri("../tests/test_files/contacts_unlabelled_2d_v2.parquet::2")
```

Then we specify a target region


```python
target_region = pd.DataFrame({
    "chrom": ['chr1'],
    "start": [100],
    "end": [400],
})
```

First, we want to select all contacts where any of the fragments constituting the contact overlaps the target region. To perform this action, we use the Overlap class and pass the target region as well as an instance of the `Anchor` class. The `Anchor` dataclass allows us to specify how we want to filter contacts for region overlap. It has two attributes `mode` and `anchors`. `Anchors` indicates the positions we want to filter on (default is all positions) and `mode` specifies whether we require all positions to overlap or any position to overlap. So for example, if we want all of our two-way contacts for which any of the positions overlap, we would use `Anchor(mode='ANY', anchors=[1,2])`.


```python
query_steps = [
    Overlap(target_region, anchor_mode=Anchor(mode="ANY", anchors=[1,2]))
]
```

A query plan is a list of qury steps that can be used in the basic query class


```python
query = Query(query_steps=query_steps)
```

The `.build` method executes the query plan and retuns a `QueryPlan` object


```python
result = query.build(contacts)
result
```




    <spoc.query_engine.QueryPlan at 0x1804b711a00>



The `.load_result` method of the `QueryResult` object can be executed using `.load_result`, which returns a `pd.DataFrame`. The resulting dataframe has additional columns that represent the regions, with which the input contacts overlapped.


```python
df = result.compute()
print(type(df))
df.filter(regex=r"chrom|start|end|id")
```

    <class 'pandas.core.frame.DataFrame'>
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>region_chrom</th>
      <th>region_start</th>
      <th>region_end</th>
      <th>region_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100</td>
      <td>200</td>
      <td>chr1</td>
      <td>1000</td>
      <td>2000</td>
      <td>chr1</td>
      <td>100</td>
      <td>400</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>2000</td>
      <td>3000</td>
      <td>chr1</td>
      <td>200</td>
      <td>300</td>
      <td>chr1</td>
      <td>100</td>
      <td>400</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>3000</td>
      <td>4000</td>
      <td>chr1</td>
      <td>300</td>
      <td>400</td>
      <td>chr1</td>
      <td>100</td>
      <td>400</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



We can also restrict the positions to filter on, by passing different anchor parameters. For example, we can filter for contacts, where the first position overlaps with our target:


```python
query_steps = [
    Overlap(target_region, anchor_mode=Anchor(mode="ANY", anchors=[1]))
]
Query(query_steps=query_steps)\
    .build(contacts)\
    .compute()\
    .filter(regex=r"chrom|start|end|id")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>region_chrom</th>
      <th>region_start</th>
      <th>region_end</th>
      <th>region_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100</td>
      <td>200</td>
      <td>chr1</td>
      <td>1000</td>
      <td>2000</td>
      <td>chr1</td>
      <td>100</td>
      <td>400</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



This time, only the first contact overlaps.

The same functionality is implemented also for the Pixels class

## Selecting a subset of contacts at multiple genomic regions
The Overlap class is also capable of selecting contacts at multiple genomic regions. Here, the behavior of `Overlap` deviates from a simple filter, because if a given contact overlaps with multiple regions, it will be returned multiple times.

Specify target regions


```python
target_regions = pd.DataFrame({
    "chrom": ['chr1', 'chr1'],
    "start": [100, 150],
    "end": [400, 200],
})
```


```python
query_steps = [
    Overlap(target_regions, anchor_mode=Anchor(mode="ANY", anchors=[1]))
]
Query(query_steps=query_steps)\
    .build(contacts)\
    .compute()\
    .filter(regex=r"chrom|start|end|id")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>region_chrom</th>
      <th>region_start</th>
      <th>region_end</th>
      <th>region_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100</td>
      <td>200</td>
      <td>chr1</td>
      <td>1000</td>
      <td>2000</td>
      <td>chr1</td>
      <td>100</td>
      <td>400</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>100</td>
      <td>200</td>
      <td>chr1</td>
      <td>1000</td>
      <td>2000</td>
      <td>chr1</td>
      <td>150</td>
      <td>200</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



In this example, the contact overlapping both regions is duplicated.

The same functionality is implemented also for the pixels class.

## Calculating the distance to a target region and aggregating the result
In this example, we calculate the distance of pixels to target regions and aggregate based on the distances. This is a very common use case in so-called pileup analyses, where we want to investigate the average behavior around regions of interest.


```python
from spoc.pixels import Pixels
from spoc.query_engine import DistanceTransformation, DistanceMode
import pandas as pd
import numpy as np
from itertools import product
```

First we define a set of target pixels


```python
def complete_synthetic_pixels():
    """Pixels that span two regions densely"""
    np.random.seed(42)
    # genomic region_1
    pixels_1 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr1"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    # genomic region_2
    pixels_2 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr2"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    return pd.concat((pd.DataFrame(pixels_1), pd.DataFrame(pixels_2)))
```


```python
pixels = Pixels(complete_synthetic_pixels(), number_fragments=3, binsize=50_000)
```

Then we define the target regions we are interested in.


```python
target_regions = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [900_000, 900_000],
            "end": [1_100_000, 1_100_000],
        }
    )
```

We are then interested in selecting all contacts that are contained within these pixels and then calculate the distance to them. The selection step can be done with the `Overlap` class that we described above. The distance transformation can be done with the `DistanceTransformation` query step. This query step takes an instance of genomic data that contains regions (as defined by it's schema) and calculates the distance to all position columns. All distances are calculated with regards to the center of each assigned region. Since genomic positions are defined by a start and end,the `DistanceTransformation` query step has a `DistanceMode` parameter that defines whether we would like to calculate the distance with regard to the start of a genomic position, the end or it's center.


```python
query_steps = [
    Overlap(target_regions, anchor_mode=Anchor(mode="ANY")),
    DistanceTransformation(
        distance_mode=DistanceMode.LEFT,
    ),
]
```

We can then execute this query plan using the Query class. This well add an distance column to the genomic dataset returned.


```python
Query(query_steps=query_steps)\
    .build(pixels)\
    .compute()\
    .filter(regex=r"chrom|distance")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom_1</th>
      <th>chrom_2</th>
      <th>chrom_3</th>
      <th>region_chrom</th>
      <th>distance_1</th>
      <th>distance_2</th>
      <th>distance_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>-50000.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>50000.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>chr1</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>100000.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>245</th>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>-100000.0</td>
    </tr>
    <tr>
      <th>246</th>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>-50000.0</td>
    </tr>
    <tr>
      <th>247</th>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>248</th>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>50000.0</td>
    </tr>
    <tr>
      <th>249</th>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>chr2</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>100000.0</td>
    </tr>
  </tbody>
</table>
<p>250 rows × 7 columns</p>
</div>



## Aggregating genomic data based on it's distance to a target region
In this example, we extend the above use-case to aggregate the results based on the distance columns added. This is a common use-case to calculate aggregate statistics for different distance levels. To achieve this, we employ the same query plan as above and extend it using the `DistanceAggregation` query step.


```python
from spoc.query_engine import DistanceAggregation, AggregationFunction
```

The `DistanceAggregation` class requires the following parameters:
- `value_columns`: Thie specifies the value to aggregate
- `function`: The aggregation function to use. This is the enumerated type `AggregationFunction`
- `densify_output`: Whether missing distance values should be filled with empty values (specific empty value depends on the aggregation function)

Note that there are two different average functions available, `AVG` and `AVG_WITH_EMPTY`. `AVG` performs and average over all available columns, where as `AVG_WITH_EMPTY` counts missing distances per regions as 0.


```python
query_steps = [
    Overlap(target_regions, anchor_mode=Anchor(mode="ALL")),
    DistanceTransformation(),
    DistanceAggregation(
        value_column='count',
        function=AggregationFunction.AVG,
    ),
]
```


```python
Query(query_steps=query_steps)\
    .build(pixels)\
    .compute()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>distance_1</th>
      <th>distance_2</th>
      <th>distance_3</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>-50000.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>0.0</td>
      <td>5.5</td>
    </tr>
    <tr>
      <th>3</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>50000.0</td>
      <td>5.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>100000.0</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>120</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>-100000.0</td>
      <td>8.0</td>
    </tr>
    <tr>
      <th>121</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>-50000.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>122</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>0.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>123</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>50000.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>124</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>125 rows × 4 columns</p>
</div>



In addition, we can also aggregate on a subset of distance positions, using the `position_list` parameter:


```python
query_steps = [
    Overlap(target_regions, anchor_mode=Anchor(mode="ALL")),
    DistanceTransformation(),
    DistanceAggregation(
        value_column='count',
        function=AggregationFunction.AVG,
        position_list=[1,2]
    ),
]
```


```python
Query(query_steps=query_steps)\
    .build(pixels)\
    .compute()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>distance_1</th>
      <th>distance_2</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-100000.0</td>
      <td>-100000.0</td>
      <td>4.8</td>
    </tr>
    <tr>
      <th>1</th>
      <td>-100000.0</td>
      <td>-50000.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-100000.0</td>
      <td>0.0</td>
      <td>5.3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>-100000.0</td>
      <td>50000.0</td>
      <td>4.7</td>
    </tr>
    <tr>
      <th>4</th>
      <td>-100000.0</td>
      <td>100000.0</td>
      <td>5.3</td>
    </tr>
    <tr>
      <th>5</th>
      <td>-50000.0</td>
      <td>-100000.0</td>
      <td>4.8</td>
    </tr>
    <tr>
      <th>6</th>
      <td>-50000.0</td>
      <td>-50000.0</td>
      <td>4.4</td>
    </tr>
    <tr>
      <th>7</th>
      <td>-50000.0</td>
      <td>0.0</td>
      <td>4.5</td>
    </tr>
    <tr>
      <th>8</th>
      <td>-50000.0</td>
      <td>50000.0</td>
      <td>5.4</td>
    </tr>
    <tr>
      <th>9</th>
      <td>-50000.0</td>
      <td>100000.0</td>
      <td>3.4</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.0</td>
      <td>-100000.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>0.0</td>
      <td>-50000.0</td>
      <td>3.5</td>
    </tr>
    <tr>
      <th>12</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>4.4</td>
    </tr>
    <tr>
      <th>13</th>
      <td>0.0</td>
      <td>50000.0</td>
      <td>5.4</td>
    </tr>
    <tr>
      <th>14</th>
      <td>0.0</td>
      <td>100000.0</td>
      <td>4.3</td>
    </tr>
    <tr>
      <th>15</th>
      <td>50000.0</td>
      <td>-100000.0</td>
      <td>5.3</td>
    </tr>
    <tr>
      <th>16</th>
      <td>50000.0</td>
      <td>-50000.0</td>
      <td>4.7</td>
    </tr>
    <tr>
      <th>17</th>
      <td>50000.0</td>
      <td>0.0</td>
      <td>4.0</td>
    </tr>
    <tr>
      <th>18</th>
      <td>50000.0</td>
      <td>50000.0</td>
      <td>4.2</td>
    </tr>
    <tr>
      <th>19</th>
      <td>50000.0</td>
      <td>100000.0</td>
      <td>6.1</td>
    </tr>
    <tr>
      <th>20</th>
      <td>100000.0</td>
      <td>-100000.0</td>
      <td>5.4</td>
    </tr>
    <tr>
      <th>21</th>
      <td>100000.0</td>
      <td>-50000.0</td>
      <td>2.8</td>
    </tr>
    <tr>
      <th>22</th>
      <td>100000.0</td>
      <td>0.0</td>
      <td>3.6</td>
    </tr>
    <tr>
      <th>23</th>
      <td>100000.0</td>
      <td>50000.0</td>
      <td>5.2</td>
    </tr>
    <tr>
      <th>24</th>
      <td>100000.0</td>
      <td>100000.0</td>
      <td>4.3</td>
    </tr>
  </tbody>
</table>
</div>


