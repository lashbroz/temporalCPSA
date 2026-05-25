# Load normal/reference molecular data used for AD-TMP trajectory analysis

The normal/reference data used by the paper's age-dependent tumor
molecular phenotype (AD-TMP) analyses are distributed with `ageTMP` as a
documented RDS object. These data were originally represented in the
paper-analysis code by `for_tadj.RData`.

## Usage

``` r
ageTMP_load_normal_reference(path = NULL)
```

## Arguments

- path:

  Optional path to a normal-reference RDS file. If `NULL`, the package
  copy at `inst/extdata/normal_reference.rds` is used.

## Value

A named list containing normal/reference matrices, sample metadata, and
provenance metadata.

## Details

The normal developmental reference data come from the DEveLopmental
Trajectory Atlas (DELTA) in dorsolateral prefrontal cortex (DLPFC),
PMID: 30518843. DELTA is available at <http://amp.pharm.mssm.edu/DELTA>.

The package copy stores the normal/reference matrices, sample metadata,
and when available, the manuscript protein N-TMP consensus cluster
assignments needed to reproduce the manuscript age-dependent trajectory
analyses without loading the original private `for_tadj.RData` object.

For protein Figure 2 style heatmaps, the stored
`protein$clusters$consensus_k4` vector follows the manuscript remapping
of the original normal consensus classes: raw cluster 3 is displayed as
N-TMP 1, raw cluster 4 as N-TMP 2, raw cluster 1 as N-TMP 3, and raw
cluster 2 as N-TMP 4. These labels are part of the heatmap ordering
logic, not just cosmetic colors.

In the manuscript protein trajectory workflow, these normal/reference
data are adjusted with `score ~ pH + PMI + Ethnicity`. Missing pH is
imputed from the full normal-reference metadata before sex
stratification, ethnicity label `H` is combined with `C`, and `C` is
used as the reference level. These details are implemented in
[`ageTMP_compare_normal_tumor_trajectory()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_compare_normal_tumor_trajectory.md)
for reproducibility.
