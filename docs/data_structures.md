# Spoc data structures
Let us have a look what data structures are available within spoc and see how they relate to each other. On a high level, spoc provides data structures for all parts of the transformation pipeline, from raw reads to aggregated pixels. 

Often, these data structures (except the pixels class) will not be used within everyday analysis tasks but rather within analysis pipelines.

## Data frame schemas
Spoc data structures are wrappers around tabular data containers such as `panda.DataFrame` or `dask.dataframe.DataFrame`. To ensure that the underlying data complies with the format that spoc expects, spoc implements dataframe validation using `pandera`. The underlying schemas reside in the `spoc.dataframe_models` file.

## I/O
Reading and writing of spoc data structures is managed by the `spoc.io` package, specifically by the `FileManager` class. Examples of using the FileManager can be found with the specific data structure.

## Fragments
Fragments encapsulate a data structure that can hold a dynamic number of aligned fragments per sequencing unit. In a Pore-C experiment, a sequencing unit is the sequencing read that holds multiple fragments per read. In theory, this structure can also be used for other experiment types that generate aligned fragments that are grouped together by an id, for example SPRITE

Reading fragments using `FileManager`


```python
from spoc.io import FileManager
fragments = FileManager().load_fragments("good_porec.parquet")
```

Fragments class has data accessor for fragments.


```python
fragments.data.head()
```




<div style="overflow-x:scroll;">
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
<table border="1" class="dataframe" style="overflow-x:scroll;">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>strand</th>
      <th>read_name</th>
      <th>read_start</th>
      <th>read_end</th>
      <th>read_length</th>
      <th>mapping_quality</th>
      <th>align_score</th>
      <th>align_base_qscore</th>
      <th>pass_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>True</td>
      <td>dummy</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>True</td>
      <td>dummy</td>
      <td>2</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>True</td>
      <td>dummy</td>
      <td>3</td>
      <td>6</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



The fragments class constructor validates the underlying data structure using pandera and the dataframe schemas in `spoc.dataframe_models`


```python
from pandera.errors import SchemaError

try:
    FileManager().load_fragments("bad_porec.parquet")
except SchemaError as e:
    print(str(e).split("\n")[0])

#> column 'chrom' not in dataframe
```

Fragments class also supports reading as dask dataframe


```python
fragments = FileManager(use_dask=True).load_fragments("good_porec.parquet")
fragments.data
```




<div><strong>Dask DataFrame Structure:</strong></div>
<div style="overflow-x:scroll;">
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>strand</th>
      <th>read_name</th>
      <th>read_start</th>
      <th>read_end</th>
      <th>read_length</th>
      <th>mapping_quality</th>
      <th>align_score</th>
      <th>align_base_qscore</th>
      <th>pass_filter</th>
    </tr>
    <tr>
      <th>npartitions=1</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th></th>
      <td>object</td>
      <td>int64</td>
      <td>int64</td>
      <td>bool</td>
      <td>object</td>
      <td>int64</td>
      <td>int64</td>
      <td>int64</td>
      <td>int64</td>
      <td>int64</td>
      <td>int64</td>
      <td>bool</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>
</div>
<div>Dask Name: validate, 2 graph layers</div>



Note that if reading from dask dataframes, schema evaluation is deferred until the dask taskgraph is evaluated.


```python
fragments = FileManager(use_dask=True).load_fragments("bad_porec.parquet")

try:
    fragments.data.compute()
except SchemaError as e:
    print(str(e).split("\n")[0])

#> column 'chrom' not in dataframe
```

### Annotating fragments
Fragments can carry metadata that adds additional information, which can be propagated in the analysis pipeline. `FragmentAnnotator` uses a dictionary called label library that contains compound fragment ids and metainformation to annotate fragments. These ids are concatenations of the read_id, chromosome, start and end of the mapping.


```python
from spoc.fragments import FragmentAnnotator
fragments = FileManager().load_fragments("good_porec.parquet")
label_library = FileManager().load_label_library("ll1.pickle")
label_library
#> {'dummy_chr1_1_4': True, 'dummy_chr1_2_5': False}

annotated_fragments = FragmentAnnotator(label_library)\
                                        .annotate_fragments(fragments)
annotated_fragments.data.head()
```




<div style="overflow-x:scroll;">
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>strand</th>
      <th>read_name</th>
      <th>read_start</th>
      <th>read_end</th>
      <th>read_length</th>
      <th>mapping_quality</th>
      <th>align_score</th>
      <th>align_base_qscore</th>
      <th>pass_filter</th>
      <th>metadata</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>True</td>
      <td>dummy</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
      <td>SisterB</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>True</td>
      <td>dummy</td>
      <td>2</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>True</td>
      <td>SisterA</td>
    </tr>
  </tbody>
</table>
</div>



## Contacts
While the fragment representation retains flexibility, it is often not practical to have contacts of multiple orders and types in different rows of the same file. To this end, we employ the contact representation, where each row contains one contact of a defined order, e.g. a duplet or a triplet. The `Contact` class is a wrapper around the data structure that holds this representation.
The `Contacts` class is a generic interface that can represent different orders.
The class that creates contacts from fragments is called `FragmentExpander`, which can be used to generate contacts of arbitrary order.


```python
import pandas as pd
from spoc.fragments import FragmentExpander

fragments = FileManager().load_fragments("fragments_unlabelled.parquet")
fragments.data.head()
```




<div style='overflow-x:scroll;'>
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
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>strand</th>
      <th>read_name</th>
      <th>read_start</th>
      <th>read_end</th>
      <th>read_length</th>
      <th>mapping_quality</th>
      <th>align_score</th>
      <th>align_base_qscore</th>
      <th>pass_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>True</td>
      <td>dummy</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>True</td>
      <td>dummy</td>
      <td>2</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>True</td>
      <td>dummy</td>
      <td>3</td>
      <td>6</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>4</td>
      <td>7</td>
      <td>True</td>
      <td>dummy</td>
      <td>4</td>
      <td>7</td>
      <td>1</td>
      <td>4</td>
      <td>4</td>
      <td>4</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>5</td>
      <td>8</td>
      <td>True</td>
      <td>dummy2</td>
      <td>5</td>
      <td>8</td>
      <td>1</td>
      <td>5</td>
      <td>5</td>
      <td>5</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>




```python
contacts = FragmentExpander(number_fragments=2).expand(fragments)
contacts.data.head()
```




<div style="overflow-x:scroll;">
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
      <th>read_name</th>
      <th>read_length</th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>mapping_quality_1</th>
      <th>align_score_1</th>
      <th>align_base_qscore_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>mapping_quality_2</th>
      <th>align_score_2</th>
      <th>align_base_qscore_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>chr1</td>
      <td>4</td>
      <td>7</td>
      <td>4</td>
      <td>4</td>
      <td>4</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>chr1</td>
      <td>4</td>
      <td>7</td>
      <td>4</td>
      <td>4</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



Fragment expander also allows us to deal with metadata that is associated with fragments.


```python
fragments_labelled = FileManager().load_fragments("fragments_labelled.parquet")
contacts_labelled = FragmentExpander(number_fragments=2)\
                                      .expand(fragments_labelled)
contacts_labelled.data.head()
```




<div style='overflow-x:scroll;'>
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
      <th>read_name</th>
      <th>read_length</th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>mapping_quality_1</th>
      <th>align_score_1</th>
      <th>align_base_qscore_1</th>
      <th>metadata_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>mapping_quality_2</th>
      <th>align_score_2</th>
      <th>align_base_qscore_2</th>
      <th>metadata_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>SisterB</td>
    </tr>
    <tr>
      <th>1</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>2</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>4</td>
      <td>7</td>
      <td>4</td>
      <td>4</td>
      <td>4</td>
      <td>SisterB</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>3</td>
      <td>6</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>4</th>
      <td>dummy</td>
      <td>1</td>
      <td>chr1</td>
      <td>2</td>
      <td>5</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>4</td>
      <td>7</td>
      <td>4</td>
      <td>4</td>
      <td>4</td>
      <td>SisterB</td>
    </tr>
  </tbody>
</table>
</div>



The contact class retains the information as to whether the expanded contacts contain metadata


```python
contacts_labelled.contains_metadata
#> True
```

### Symmetry


#### Unlabelled contacts
The quantification of genomic interactions in conventional (2-way) Hi-C assumes no difference in the order of interactions. This means that whether a genomic location is measured in the first read or second read of a paired-end sequencing experiment carries the same information. This means that during preprocessing, conventional Hi-C data is flipped based on some convention (often, the first read has a smaller genomic location based on some sort of order), and then only the upper triangular interaction matrix is stored.

When we talk about higher genomic order, a similar reasoning can apply (except for special use cases), and we thus can flip genomic contacts such that genomic coordinates are monotonically increasing from lower to higher order (we mean this order if we refer to flipping below). This produces a symmetric, high-dimensional tensor, meaning that every permutation of dimensions does not change the associated value.

In `spoc`, this logic is implemented in the `ContactManipulator` class, in the `.flip_symmetric_contacts` method.



```python
contacts = FileManager().load_contacts("contacts_unlabelled_2d.parquet")
contacts.data.head().filter(regex="(read_name|start|end)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>start_2</th>
      <th>end_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>1000</td>
      <td>2000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>2000</td>
      <td>3000</td>
      <td>200</td>
      <td>300</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>3000</td>
      <td>4000</td>
      <td>300</td>
      <td>400</td>
    </tr>
  </tbody>
</table>
</div>



This particular contacts dataframe has one contact that conforms with the convention that the first contact should be smaller than the second (`read1`), whereas the other two contacts don't conform to that convention. Using the `.flip_symmetric_contacts` method, we can fix this:


```python
from spoc.contacts import ContactManipulator
flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
flipped_contacts.data.head().filter(regex="(read_name|start|end)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>start_2</th>
      <th>end_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>1000</td>
      <td>2000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>200</td>
      <td>300</td>
      <td>2000</td>
      <td>3000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>300</td>
      <td>400</td>
      <td>3000</td>
      <td>4000</td>
    </tr>
  </tbody>
</table>
</div>



Symmetry flipped contacts have the flag `symmetry_flipped` set to true.


```python
flipped_contacts.symmetry_flipped
#> True
```

These operations are available for arbitrary contact cardinalities:


```python
contacts = FileManager().load_contacts("contacts_unlabelled_3d.parquet")
contacts.data.head().filter(regex="(read_name|start|end)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>start_3</th>
      <th>end_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>1000</td>
      <td>2000</td>
      <td>250</td>
      <td>300</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>2000</td>
      <td>3000</td>
      <td>200</td>
      <td>300</td>
      <td>400</td>
      <td>500</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>3000</td>
      <td>4000</td>
      <td>300</td>
      <td>400</td>
      <td>100</td>
      <td>200</td>
    </tr>
  </tbody>
</table>
</div>




```python
flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
flipped_contacts.data.head().filter(regex="(read_name|start|end)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>start_3</th>
      <th>end_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>250</td>
      <td>300</td>
      <td>1000</td>
      <td>2000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>200</td>
      <td>300</td>
      <td>400</td>
      <td>500</td>
      <td>2000</td>
      <td>3000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>100</td>
      <td>200</td>
      <td>300</td>
      <td>400</td>
      <td>3000</td>
      <td>4000</td>
    </tr>
  </tbody>
</table>
</div>



#### Labelled contacts

For labelled contacts, the situation is more complex as we have to deal with different flavours of symmetry. If we take the example of triplets, that is genomic contacts of order 3, and binary labels (denoted as A or B), there are 8 possible contact orders:
```
 (AAA, BBB, BAA, ABA, AAB, ABB, BAB, BBA)
``` 

If we extend the argument that order of interactions is unimportant, this reduces to 4 possible arrangements of labels:

```
 (AAA, BBB, ABB, BAA)
```
 
It is often the case that for binary contacts, the specific label type is unimportant, the only important information is whether the labels were different (e.g. for sister-specific labels or homologous chromosome labels). In such a situation, the possible arrangements are reduced further to 2 possible label types:

```
 (AAA, AAB)
```

For those contact types, two different “rules” for symmetry apply; for the situation of all similar labels (AAA or BBB), the same rules apply for unlabelled contacts as there is no difference in order. For the other situation (ABB or BBA, which we denote as ABB from now on), only permutations within one label type produce the same value, meaning that if we have the label state ABB, and denote permutations as tuples of length three with (0,1,2) being the identity permutation then only the permutations (0,1,2) and (0,2,1) are identical. Practically, this means that we can flip contacts that are related by these permutations such that their genomic coordinates are monotonically increasing, but we cannot do this for contacts that are not related through a symmetry relation. 

This reasoning can be generalized to higher dimensions, where contacts can be flipped if they can be related via a symmetry relation. Practically speaking, this means that we have a higher-order contact with two possible label states of order n with k labels of type A and (n-k) labels of type B; we can flip the labels within type A and within type B. For example, if we have a contact of order 4 with the configuration AABB, we can flip within A and within B based on genomic coordinates, but not within them. This reasoning also applies to the situation where we have more than one possible label state. Also, here, we can flip within one label state but not between them.

This logic is implemented in `spoc` in the `ContactManipulator` class, with to methods:

- The `.equate_binary_labels` method can be used to specify whether in a binary label situation, the labels shold be different or not (e.g. whether AA is equivaltent to BB) or not. This method is optional and can be used prior to the flipping procedure.
- The `.flip_symmetric_contacts` method flips symmetric contacts base don the rules specified above

Note that all symmetry operations require the contact metadata to be alphabetically sorted; this can be either done explicitly via the `.sort_labels` method, or is performed automatically within the other methods. For example, the metadata order of `ABA` will be converted to `AAB` by the `.sort_labels` operation.


Let's look at an example! Here, we have 3D contacts that have binary labels.


```python
contacts = FileManager().load_contacts("contacts_labelled_3d.parquet")
contacts.data.head().filter(regex="(read_name|start|end|metadata)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>B</td>
      <td>1000</td>
      <td>2000</td>
      <td>B</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>5000</td>
      <td>5500</td>
      <td>A</td>
      <td>2000</td>
      <td>3000</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>B</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>800</td>
      <td>900</td>
      <td>B</td>
      <td>3000</td>
      <td>3200</td>
      <td>B</td>
      <td>3000</td>
      <td>4000</td>
      <td>B</td>
    </tr>
  </tbody>
</table>
</div>




```python
contacts.number_fragments, contacts.get_label_values()
#> (3, {'A', 'B'})
```

This dataframe contains three contacts with the following label state:
```
(ABB)
(AAB)
(BBB)
```
In this analysis use case, we want to equate the binary labels since we don't have a biological reason to believe that there is any difference between the labels. The only information that is important for us is whether the contacts happened between different label states or the same label state. Therefore, we use the `.equate_binary_labels` method to replace all occurrences of the same contact combination with their alphabetically first example. For example, the label state `(ABB)` is a contact, where two parts come from one label state and one part comes from another. In our logic, this is equivalent to `(AAB)`, which is the alphabetically first label, which we therefore use to replace it. Following the same logic, `(BBB)` will be replaced by `(AAA)`.


```python
equated_contacts = ContactManipulator().equate_binary_labels(contacts)
equated_contacts.data.head().filter(regex="(read_name|start|end|metadata)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>A</td>
      <td>1000</td>
      <td>2000</td>
      <td>B</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>5000</td>
      <td>5500</td>
      <td>A</td>
      <td>2000</td>
      <td>3000</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>B</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>800</td>
      <td>900</td>
      <td>A</td>
      <td>3000</td>
      <td>3200</td>
      <td>A</td>
      <td>3000</td>
      <td>4000</td>
      <td>A</td>
    </tr>
  </tbody>
</table>
</div>



As you can see, the occurrence of `(ABB)` has been replaced by `(AAB)`, and the occurrence of `(BBB)` has been replaced by `(AAA)`. We can now reduce the symmetry of these contacts based on the logic explained above using the `.flip_symmetric_contacts` method.


```python
flipped_labelled_contacts = ContactManipulator()\
                .flip_symmetric_contacts(equated_contacts)
flipped_labelled_contacts.data.head()\
          .filter(regex="(read_name|start|end|metadata)")
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
      <th>read_name</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>read1</td>
      <td>100</td>
      <td>200</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>A</td>
      <td>1000</td>
      <td>2000</td>
      <td>B</td>
    </tr>
    <tr>
      <th>1</th>
      <td>read2</td>
      <td>2000</td>
      <td>3000</td>
      <td>A</td>
      <td>5000</td>
      <td>5500</td>
      <td>A</td>
      <td>200</td>
      <td>300</td>
      <td>B</td>
    </tr>
    <tr>
      <th>2</th>
      <td>read3</td>
      <td>800</td>
      <td>900</td>
      <td>A</td>
      <td>3000</td>
      <td>3200</td>
      <td>A</td>
      <td>3000</td>
      <td>4000</td>
      <td>A</td>
    </tr>
  </tbody>
</table>
</div>



### Persisting contacts
Contacts can be persisted using the file manager, which writes the data as well as the global parameters to a parquet file. Loading the file restores the global parameters.


```python
flipped_labelled_contacts.get_global_parameters()
#> ContactsParameters(number_fragments=3,
#>                     metadata_combi=None,
#>                     label_sorted=True,
#>                     binary_labels_equal=True,
#>                     symmetry_flipped=True)
FileManager().write_multiway_contacts("test.parquet", flipped_labelled_contacts)
contacts = FileManager().load_contacts("test.parquet")
contacts.get_global_parameters()
#> ContactsParameters(number_fragments=3,
#>                     metadata_combi=None,
#>                     label_sorted=True,
#>                     binary_labels_equal=True,
#>                     symmetry_flipped=True)
```

## Pixels

### Concept

Pixels represent aggregated information that counts the number of genomic contacts per genomic bin. In this context, a genomic bin is a genomic interval with a fixed size (for example, 10 kb), and a genomic contact is defined above as an interaction between the corresponding bins. Pixels represent a single contact order, binsize and labelling state (if metadata is provided) and are thought to be the central data structure that an analyst will use to interact with multiway genomic data to answer questions about genomic structure.

### Implementation

Within `spoc`, pixel instances can be generated from genomic contacts through the `GenomicBinner` class.


```python
from spoc.pixels import GenomicBinner
contacts = FileManager().load_contacts("contacts_for_pixels_3d.parquet")
contacts.data.head().filter(regex="(read_name|chrom|start|end|metadata)")
```




<div style="overflow-x:scroll;">
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
      <th>read_name</th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>chrom_3</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>a</td>
      <td>chr1</td>
      <td>100010</td>
      <td>100015</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>500010</td>
      <td>500050</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>600100</td>
      <td>600200</td>
      <td>SisterB</td>
    </tr>
    <tr>
      <th>1</th>
      <td>b</td>
      <td>chr1</td>
      <td>5000010</td>
      <td>5000050</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>7000050</td>
      <td>7000070</td>
      <td>SisterA</td>
      <td>chr4</td>
      <td>2000300</td>
      <td>2000400</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>2</th>
      <td>c</td>
      <td>chr1</td>
      <td>10000010</td>
      <td>10000050</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>25000800</td>
      <td>25000900</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>6000050</td>
      <td>6000600</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>3</th>
      <td>d</td>
      <td>chr1</td>
      <td>10000010</td>
      <td>10000050</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>25001000</td>
      <td>25002000</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>6000010</td>
      <td>6000700</td>
      <td>SisterB</td>
    </tr>
  </tbody>
</table>
</div>



`GenomicBinner` takes a binsize and aggregates contacts by counting the number of contacts that fall within a genomic bin. As with all the other preprocessing functionalities, `GenomicBinner` can take arguments as either pandas dataframes or dask dataframes. The output is consistent, meaning that when passing a pandas dataframe, a pandas dataframe is returned, and if passing a dask dataframe, a dask dataframe is returned.


```python
pixels = GenomicBinner(bin_size=100_000).bin_contacts(contacts)
```

The default behaviour of genomic binner is to filter for contacts that are on the same chromosome and produce an intrachromosomal pixels instance. This behaviour allows us to create a compact representation of pixels, which only store one chromosome field per pixel.


```python
pixels.data.head()
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
      <th>chrom</th>
      <th>start_1</th>
      <th>start_2</th>
      <th>start_3</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100000</td>
      <td>500000</td>
      <td>600000</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>10000000</td>
      <td>25000000</td>
      <td>6000000</td>
      <td>2</td>
    </tr>
  </tbody>
</table>
</div>



The pixels class additionally contains all information associated with the data stored. This data is taken from the contacts object passed to the genomic binner.


```python
(
  pixels.number_fragments,
  pixels.binsize,
  pixels.binary_labels_equal,
  pixels.symmetry_flipped,
  pixels.metadata_combi
)
#> (3, 100000, False, False, None)
```

If interchromosomal contacts are needed, this can be passed to the `.bin_contacts` method of genomic binner and will cause the pixel schema to incorporate chromosome columns for each of the corresponding pixel dimensions.


```python
pixels_w_inter = GenomicBinner(bin_size=100_000)\
                              .bin_contacts(contacts, same_chromosome=False)
pixels_w_inter.data.head()
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
      <th>chrom_2</th>
      <th>start_2</th>
      <th>chrom_3</th>
      <th>start_3</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100000</td>
      <td>chr1</td>
      <td>500000</td>
      <td>chr1</td>
      <td>600000</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>10000000</td>
      <td>chr1</td>
      <td>25000000</td>
      <td>chr1</td>
      <td>6000000</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>5000000</td>
      <td>chr1</td>
      <td>7000000</td>
      <td>chr4</td>
      <td>2000000</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



### Lablled contacts
The concept of pixels representing a single labelling state allows simplifying the interfaces to downstream processing functionality but offloads responsibility to the analyst to ensure that the chosen labelling state adequately addresses the biological question at hand. The pixels class does not perform any labelling state checks and forwards the combination(s) of labelling states present in the contacts class used for their construction. This allows for greater flexibility with regard to the questions that can be answered. 

The `ContactManipulator` class contains functionality to filter contacts for a given labelling state to perform downstream aggregation of pixels.


```python
contacts.data.head().filter(regex="(read_name|chrom|start|end|metadata)")
```




<div style="overflow-x:scroll;">
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
      <th>read_name</th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>chrom_3</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>a</td>
      <td>chr1</td>
      <td>100010</td>
      <td>100015</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>500010</td>
      <td>500050</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>600100</td>
      <td>600200</td>
      <td>SisterB</td>
    </tr>
    <tr>
      <th>1</th>
      <td>b</td>
      <td>chr1</td>
      <td>5000010</td>
      <td>5000050</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>7000050</td>
      <td>7000070</td>
      <td>SisterA</td>
      <td>chr4</td>
      <td>2000300</td>
      <td>2000400</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>2</th>
      <td>c</td>
      <td>chr1</td>
      <td>10000010</td>
      <td>10000050</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>25000800</td>
      <td>25000900</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>6000050</td>
      <td>6000600</td>
      <td>SisterA</td>
    </tr>
    <tr>
      <th>3</th>
      <td>d</td>
      <td>chr1</td>
      <td>10000010</td>
      <td>10000050</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>25001000</td>
      <td>25002000</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>6000010</td>
      <td>6000700</td>
      <td>SisterB</td>
    </tr>
  </tbody>
</table>
</div>




```python
filtered_contacts = ContactManipulator()\
                      .subset_on_metadata(contacts,
                                          metadata_combi=['SisterA',
                                                          'SisterB',
                                                          'SisterB'])
filtered_contacts.data.head()\
                .filter(regex="(read_name|chrom|start|end|metadata)")
```




<div style='overflow-x:scroll;'>
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
      <th>read_name</th>
      <th>chrom_1</th>
      <th>start_1</th>
      <th>end_1</th>
      <th>metadata_1</th>
      <th>chrom_2</th>
      <th>start_2</th>
      <th>end_2</th>
      <th>metadata_2</th>
      <th>chrom_3</th>
      <th>start_3</th>
      <th>end_3</th>
      <th>metadata_3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>a</td>
      <td>chr1</td>
      <td>100010</td>
      <td>100015</td>
      <td>SisterA</td>
      <td>chr1</td>
      <td>500010</td>
      <td>500050</td>
      <td>SisterB</td>
      <td>chr1</td>
      <td>600100</td>
      <td>600200</td>
      <td>SisterB</td>
    </tr>
  </tbody>
</table>
</div>



The resulting contacts carry the filtered metadata combination as an attribute:


```python
filtered_contacts.metadata_combi
#> ['SisterA', 'SisterB', 'SisterB']
```
`GenomicBinner` carries this information forward and adds it to the pixels object


```python
pixels_filtered = GenomicBinner(bin_size=100_000)\
                                .bin_contacts(filtered_contacts)
pixels_filtered.data.head()
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
      <th>chrom</th>
      <th>start_1</th>
      <th>start_2</th>
      <th>start_3</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>100000</td>
      <td>500000</td>
      <td>600000</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
pixels_filtered.metadata_combi
#> ['SisterA', 'SisterB', 'SisterB']
```

### Persisting pixels
Pixels come in a multitude of different flavours, and it would thus be cumbersome to have to keep track of pixels of different binsize, labelling states or contact order that belong to the same source fragments. We thus developed a pixel file format that exposes a unified name to pixels that belong together, as well as functionality to load the exact pixels that we want.

#### The pixel file format
The pixel file format is a directory that contains a number of parquet files as well as a config file that holds information about what specific pixels are available. As with the other io-operations, the `FileManager` object is responsible for saving pixels.


```python
import os
FileManager().write_pixels("test_pixels.parquet", pixels)
os.listdir('test_pixels.parquet')
#> ['5f0e4d6a6c0e5afb82c3a6ec4bce1635.parquet', 'metadata.json']
```

The pixel foler at `test.parquet` contains a single specific pixel file as well as a metadata file, which can be read by `FileManager`:

```python
FileManager().list_pixels("test_pixels.parquet")
#> [
#>     PixelParameters(
#>                     number_fragments=3,
#>                     binsize=100000,
#>                     metadata_combi=None,
#>                     label_sorted=False,
#>                     binary_labels_equal=False,
#>                     symmetry_flipped=False,
#>                     same_chromosome=True
#                     )
#> ]
```

- The `list_pixels` method lists all available pixels in the form of `PixelParameters`, a pydantic dataclass that holds the parameters that describe a specific pixel file

Writing additional pixels updates the metadata file


```python
pixels_filtered._metadata_combi = ['A', 'B', 'B'] # make metadata shorter
FileManager().write_pixels("test_pixels.parquet", pixels_filtered)
FileManager().list_pixels("test_pixels.parquet")
#>     [
#>       PixelParameters(...),
#>       PixelParameters(...)
#>     ] 
```

#### Loading pixels

Pixel files can be loaded using the parameter pydantic class


```python
from spoc.file_parameter_models import PixelParameters
loaded_pixels = FileManager()\
                          .load_pixels("test_pixels.parquet",
                                        PixelParameters(number_fragments=3,
                                                        binsize=100_000,
                                                        same_chromosome=True)
                                      )
```

Additionally, pixels can be loaded using a URI that specifies the parameters seperated by `::`. This method of loading requires the parameters path, number_fragments and binsize, but will match the rest to the available pixels.


```python
from spoc.pixels import Pixels
loaded_pixels = Pixels.from_uri('test_pixels.parquet::3::100000::ABB')
loaded_pixels.get_global_parameters()

#>    PixelParameters(
#>                    number_fragments=3,
#>                    binsize=100000,
#>                    metadata_combi=['A', 'B', 'B'],
#>                    label_sorted=False,
#>                    binary_labels_equal=False,
#>                    symmetry_flipped=False,
#>                    same_chromosome=True
#>                   )
```

If the specified URI is not specific enough and multiple pixels match, an error is raised.
