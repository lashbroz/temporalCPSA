# Extract sex-specific glycopeptide trajectory score matrices

This is the package analogue of the manuscript `get_terms(..., 2)`
helper used by the glyco ADO survival scripts. It preserves the three
returned matrix types: discovery tumor score, reference trajectory
score, and discovery trajectory score.

## Usage

``` r
ageTMP_extract_glyco_trajectory_terms(
  female_trajectory,
  male_trajectory,
  discovery_clinical,
  reference_clinical,
  sex = c("Male", "Female"),
  slot = 2
)
```

## Arguments

- female_trajectory:

  Female trajectory list, usually `f.gene.mat.adj.list`.

- male_trajectory:

  Male trajectory list, usually `m.gene.mat.adj.list`.

- discovery_clinical:

  Discovery clinical data with `id` and `Gender`.

- reference_clinical:

  Reference clinical data with `id` and `Gender`.

- sex:

  One of `"Male"` or `"Female"`.

- slot:

  Trajectory list element to use. The manuscript glyco SA scripts use
  adaptive slot `2`.

## Value

A list with `score.adj`, `dtt.v`, and `dtt.h` matrices.
