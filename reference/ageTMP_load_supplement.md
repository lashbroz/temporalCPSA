# Read a public supplementary workbook sheet

Read a public supplementary workbook sheet

## Usage

``` r
ageTMP_load_supplement(data_dir = "data", table, sheet, ...)
```

## Arguments

- data_dir:

  Path to the public data directory.

- table:

  Supplementary table name, such as `"STable4"`.

- sheet:

  Sheet name.

- ...:

  Additional arguments passed to
  [`readxl::read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html).

## Value

A tibble containing the requested sheet.
