# Cluster samples from temporal molecular profiles

Cluster samples using the stacked TMP matrix produced by
[`ageTMP_combine_tmp_matrices()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_combine_tmp_matrices.md).
When `ConsensusClusterPlus` is available, this uses the consensus
k-means workflow used by the manuscript clustering scripts. Otherwise,
it falls back to deterministic k-means for lightweight development and
testing.

## Usage

``` r
ageTMP_cluster_tmp_samples(
  tmp_matrix,
  k_max = 10,
  reps = 300,
  p_item = 1,
  p_feature = 1,
  algorithm = "km",
  distance = "euclidean",
  seed = 1234,
  title = "ageTMP_consensus_tmp",
  plot = NULL,
  use_consensus = TRUE
)
```

## Arguments

- tmp_matrix:

  Numeric feature-by-sample TMP matrix.

- k_max:

  Maximum number of clusters to evaluate.

- reps:

  Number of consensus resampling repetitions.

- p_item:

  Sample resampling proportion for consensus clustering.

- p_feature:

  Feature resampling proportion for consensus clustering.

- algorithm:

  Clustering algorithm passed to `ConsensusClusterPlus`.

- distance:

  Distance metric passed to `ConsensusClusterPlus`.

- seed:

  Random seed.

- title:

  Output title/directory prefix used by `ConsensusClusterPlus`.

- plot:

  Plot option passed to `ConsensusClusterPlus`.

- use_consensus:

  Whether to use `ConsensusClusterPlus` when available.

## Value

A list with the input `clustme`, cluster assignments by k, and the raw
consensus result when applicable.
