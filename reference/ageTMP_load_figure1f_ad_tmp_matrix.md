# Load Figure 1F archived AD-TMP heatmap matrix

Reads an external archived Figure 1F AD-TMP heatmap matrix from a TSV
file. This helper exists so manuscript reproduction scripts can load the
archived Figure 1F heatmap input with an explicit contextual note, while
keeping the matrix as a manuscript-specific reproduction input rather
than bundled package data.

When `quiet = FALSE`, the loader prints a short note explaining that
this object reflects the manuscript figure-preparation matrix rather
than a source-derived matrix regenerated from final repository data
tables.

## Usage

``` r
ageTMP_load_figure1f_ad_tmp_matrix(path, quiet = FALSE)
```

## Arguments

- path:

  Path to `figure1f_ad_tmp_legacy_clustme.tsv`, with feature names in
  the first column and sample IDs in the remaining columns.

- quiet:

  Logical; if `TRUE`, suppress the contextual message printed when the
  matrix is loaded.

## Value

A numeric matrix with AD-TMP features in rows and samples in columns.
