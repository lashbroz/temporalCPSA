# Select contiguous age classes from an AD-TMP matrix

Partition age-ordered samples into contiguous age segments using
molecular structure in an age-dependent temporal molecular profile
(AD-TMP) matrix. Unlike unconstrained clustering, this procedure only
considers segmentations that can be written as ordered age intervals. It
therefore guarantees age-contiguous classes by construction.

## Usage

``` r
ageTMP_segment_age_classes(
  tmp_matrix,
  sample_age,
  k = 5,
  min_n = 10,
  labels = NULL,
  scale_rows = FALSE
)
```

## Arguments

- tmp_matrix:

  Numeric AD-TMP matrix with features in rows and samples in columns.

- sample_age:

  Numeric sample ages. If named, names are matched to
  `colnames(tmp_matrix)` after sample-ID normalization. If unnamed, the
  vector must be in the same order as `colnames(tmp_matrix)`.

- k:

  Number of contiguous age segments to select.

- min_n:

  Minimum number of samples allowed in each segment.

- labels:

  Optional labels for the resulting age classes. If `NULL`, labels are
  `segment1`, `segment2`, ...

- scale_rows:

  Whether to center and scale AD-TMP rows before scoring segments. Set
  this to `FALSE` when `tmp_matrix` has already been row-scaled.

## Value

An object of class `ageTMP_age_segmentation`, a list containing the
ordered sample metadata, named class assignments, selected cutpoints,
segment summary, dynamic-programming score, and the fitted cost matrix.

## Details

The method is a reusable formalization inspired by the manuscript
age-class workflow, where AD-TMP clustering and an age-ordered green
diagnostic heatmap were used to judge molecularly supported age strata.
The original manuscript step was interpretive: age-contiguous structure
was assessed visually from the ordered clustering diagnostics rather
than solved as a constrained optimization problem. In the manuscript,
the final age classes balanced data-driven AD-TMP structure with
developmental precedent, with particular attention to retaining
adolescence as a distinct biology-informed interval. This function makes
the contiguity principle explicit for new cohorts by optimizing over
possible age cutpoints directly. It is intended as an exploratory
diagnostic for age-class assessment and interpretation, not as a
decisive automated step in the CPSA pipeline.
