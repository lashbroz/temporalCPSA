# Plot an ordered age-segmentation diagnostic heatmap

Draw the green ordered-segmentation diagnostic used to evaluate
candidate age-class resolutions. Rows show different values of `k`;
columns are samples sorted by age; colors encode contiguous segment
membership. The selected solution is overlaid with vertical cutpoint
guides so the proposed age classes can be read directly from the
diagnostic display.

## Usage

``` r
ageTMP_plot_segment_diagnostic(
  diagnostic,
  selected_k = NULL,
  class_labels = NULL,
  palette = NULL,
  main = "AD-TMP contiguous age-segmentation diagnostic",
  show_age_axis = TRUE,
  suggested_cutpoints = NULL,
  suggested_labels = NULL,
  ...
)
```

## Arguments

- diagnostic:

  Output of
  [`ageTMP_segment_diagnostic_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_segment_diagnostic_matrix.md).

- selected_k:

  Which fitted `k` to overlay. If `NULL`, defaults to the fitted `k`
  with the largest mean silhouette when available, otherwise the largest
  fitted `k`.

- class_labels:

  Optional labels for the selected classes. If `NULL`, the labels stored
  in the selected fit are used.

- palette:

  Optional color palette. Defaults to a light-to-dark green palette.

- main:

  Plot title.

- show_age_axis:

  Whether to show an age axis below the heatmap.

- suggested_cutpoints:

  Optional numeric age cutpoints to overlay across the diagnostic
  heatmap. Use this to draw proposed or manuscript-selected age classes
  on top of the ordered segmentation landscape.

- suggested_labels:

  Optional labels for the age intervals defined by
  `suggested_cutpoints`.

- ...:

  Additional arguments passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

Invisibly returns `diagnostic`.
