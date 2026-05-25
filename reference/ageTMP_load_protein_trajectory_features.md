# Load manuscript protein trajectory feature universe

Load manuscript protein trajectory feature universe

## Usage

``` r
ageTMP_load_protein_trajectory_features(
  rdata_path,
  list_index = 2,
  matrix_name = "hope.mat0.adj"
)
```

## Arguments

- rdata_path:

  Path to a `protein_rev_tadj62_list.RData`-style object.

- list_index:

  List element containing the manuscript trajectory matrices.

- matrix_name:

  Matrix whose row names define the protein feature universe.

## Value

Character vector of protein features.
