# Manuscript glyco ADO Cox model specification

This returns the exact covariate and interaction specification used by
the glyco PED/ADO survival scripts. It is useful when callers want the
package CPSA engine but need the manuscript glyco ADO settings to be
explicit rather than implied by defaults.

## Usage

``` r
ageTMP_glyco_ado_cpsa_spec(
  cohort = c("discovery", "reference"),
  scale_features = TRUE
)
```

## Arguments

- cohort:

  One of `"discovery"` or `"reference"`.

- scale_features:

  Whether features should be row-scaled by the fit engine. The legacy
  scripts scaled rows immediately before model fitting.

## Value

A CPSA specification list.
