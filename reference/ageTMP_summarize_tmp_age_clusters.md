# Summarize age cutpoints from TMP sample clusters

Order clustered samples by age, relabel clusters by first appearance
along the age axis, and summarize each cluster's age range. This mirrors
the final interpretation step in the manuscript age-class clustering
workflow.

## Usage

``` r
ageTMP_summarize_tmp_age_clusters(
  clusters,
  clinical,
  sample_col = "id",
  age_col = "age"
)
```

## Arguments

- clusters:

  Named cluster-assignment vector.

- clinical:

  Clinical metadata.

- sample_col:

  Clinical sample ID column.

- age_col:

  Clinical age column.

## Value

A list containing relabeled clusters and an age-range summary.
