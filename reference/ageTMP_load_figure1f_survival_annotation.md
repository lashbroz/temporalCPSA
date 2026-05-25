# Load Figure 1F survival-days annotation data

Reads an external archived Figure 1F survival-days annotation TSV file.
This helper exists so reproduction scripts can load the Figure 1F
descriptive annotation explicitly, without treating it as a general CPSA
survival modeling input or bundling the manuscript-specific object as
package data.

When `quiet = FALSE`, the loader prints a short note explaining that
this object reflects the manuscript figure-preparation cohort rather
than a final source-derived survival modeling table.

## Usage

``` r
ageTMP_load_figure1f_survival_annotation(path, quiet = FALSE)
```

## Arguments

- path:

  Path to `figure1f_survival_annotation_data.tsv`.

- quiet:

  Logical; if `TRUE`, suppress the contextual message printed when the
  packaged annotation is loaded.

## Value

A data frame with the Figure 1F survival-days annotation data.
