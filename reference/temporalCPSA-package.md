# temporalCPSA: Cross-Population Survival Analysis via Temporal Molecular Profiles

`temporalCPSA` provides tools for both manuscript reproduction and
reusable cross-population survival analysis via temporal molecular
profiles. The package scope includes TMP generation, normal/reference
versus tumor trajectory comparison, trajectory visualization and
divergence analysis, multi-omic age-class clustering, and downstream
association or survival modeling.

## Getting started

Start with
[`ageTMP_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_cpsa_spec.md)
to define a survival model specification and
[`ageTMP_fit_reference_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_cpsa.md)
to fit feature-wise Cox models in a reference cohort. Users may supply
their own age classes, use biologically motivated clinical age classes,
explore age classes from AD-TMP diagnostics, or omit age-class structure
entirely.

## Core functions

- [`ageTMP_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_cpsa_spec.md)
  defines a CPSA model.

- [`ageTMP_fit_reference_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_cpsa.md)
  fits feature-wise reference-cohort survival models.

- [`ageTMP_predict_tumor_trajectory_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_predict_tumor_trajectory_matrix.md)
  derives temporal molecular profile matrices for downstream modeling.

- [`ageTMP_rank_trajectory_sd()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_rank_trajectory_sd.md)
  and
  [`ageTMP_filter_trajectory_sd()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_filter_trajectory_sd.md)
  screen low-dynamic-range trajectory features.

- [`ageTMP_segment_age_classes()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_segment_age_classes.md)
  and
  [`ageTMP_plot_segment_diagnostic()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_plot_segment_diagnostic.md)
  provide optional exploratory age-class diagnostics.

## Installed quick start

A short installed Markdown guide is available at:
`system.file("doc", "temporalCPSA-quick-start.md", package = "temporalCPSA")`.
The GitHub README is available at
<https://github.com/lashbroz/temporalCPSA>.

## See also

Useful links:

- <https://lashbroz.github.io/temporalCPSA/>

- <https://github.com/lashbroz/temporalCPSA>

- Report bugs at <https://github.com/lashbroz/temporalCPSA/issues>

## Author

**Maintainer**: Nicole Tignor <nicole.tignor@gmail.com>

## Examples

``` r
help(package = "temporalCPSA")
?temporalCPSA
?ageTMP_fit_reference_cpsa
```
