# Query engine
This technical document describes the spoc query engine, a set of classes that implements spoc's interface for querying multi-dimensional genomic data.

## Principles

### Composable pieces
Spoc's query engine consists of composable pieces that can be combined to produce an expressive query language. These pieces represent basic operations on genomic data that are easily implemented and understood on their own. This allows a great degree of flexibility, while also allowing predefined recipes that less experienced users can get started with.

### Lazy evaluation
The spoc query engine is designed with lazy evaluation as a guiding principle. This means that data queries are only executed when they are needed to minimize loading data into memory and computational overhead. To enable this, spoc queries have a construction phase, which specifies the operations to be executed and an exection phase, that actually executes the query.

## Query plans and query steps

The most important ingredient in this query language is a class that implements the `QueryStep` protocol. This protocol serves two purposes:

- It exposes a way to validate the data schema during query building
- It implements adding itself to a query

This way, query steps can be combined into a query plan that specifies the analysis to be executed. Specific examples of query steps are:

- **Snipper**: Implements selecting overlapping contacts or pixels for a set of genomic regions.
- **Transformation**: Transforms one or more columns to add additional columns
- **Aggregation**: Aggregation of data such as counting contacts per region

### Input and output of query steps

A query step takes as input a class that implements the `GenomicData` protocol. This protocol allows retrievel of the data schema (a thin wrapper over a pandera dataframe schema) as well as the data itself. The output of a query step is again a class that ipmlements the `GenomicData` protocol to allow composition. Specific examples of possible inputs are:

- **Pixels**: Represents input pixels
- **Contacts**: Represents input contacts
- **QueryResult**: The result of a query step

### Composition of query steps

To allow specifying complex queries, query steps need to be combined. This is done using the `BasicQuery` class. It takes a query plan (a list of `QueryStep` instances) as input, exposes the `query` method, which takes input data, validates all query steps and adds them to the resulting `QueryResult` instance that is returned.

### Manifestation of results

So far, we have only talked about specifying the query to be executed, but not how to actually execute it. A `QueryResult` has a `load_result()` method that returns the manifested dataframe as a `pd.DataFrame` instance. This is the step that actually executes the specified query.

## Examples

### Selecting a subset of contacts at a single genomic position
In this example, we want to select a subset of genomic contacts at a single location. For this, we first load the required input data:


```python
from spoc.query_engine import Snipper, Anchor, BasicQuery
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

First, we want to select all contacts where any of the fragments constituting the contact overlaps the target region. To perform this action, we use the Snipper class and pass the target region as well as an instance of the `Anchor` class. The `Anchor` dataclass allows us to specify how we want to filter contacts for region overlap. It has two attributes `mode` and `anchors`. `Anchors` indicates the positions we want to filter on (default is all positions) and `mode` specifies whether we require all positions to overlap or any position to overlap. So for example, if we want all of our two-way contacts for which any of the positions overlap, we would use `Anchor(mode='ANY', anchors=[1,2])`.


```python
query_plan = [
    Snipper(target_region, anchor_mode=Anchor(mode="ANY", anchors=[1,2]))
]
```

A query plan is a list of qury steps that can be used in the basic query class


```python
query = BasicQuery(query_plan=query_plan)
```

The `.query` method executes the query plan and retuns a `QueryResult` object


```python
result = query.query(contacts)
result
```




    <spoc.query_engine.QueryResult at 0x23d0367eaf0>



The `.load_result` method of the `QueryResult` object can be executed using `.load_result`, which returns a `pd.DataFrame`. The resulting dataframe has additional columns that represent the regions, with which the input contacts overlapped.


```python
df = result.load_result()
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>id</th>
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
query_plan = [
    Snipper(target_region, anchor_mode=Anchor(mode="ANY", anchors=[1]))
]
BasicQuery(query_plan=query_plan)\
    .query(contacts)\
    .load_result()\
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>id</th>
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
The Snipper class is also capable of selecting contacts at multiple genomic regions. Here, the behavior of `Snipper` deviates from a simple filter, because if a given contact overlaps with multiple regions, it will be returned multiple times.

Specify target regions


```python
target_regions = pd.DataFrame({
    "chrom": ['chr1', 'chr1'],
    "start": [100, 150],
    "end": [400, 200],
})
```


```python
query_plan = [
    Snipper(target_regions, anchor_mode=Anchor(mode="ANY", anchors=[1]))
]
BasicQuery(query_plan=query_plan)\
    .query(contacts)\
    .load_result()\
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>id</th>
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
