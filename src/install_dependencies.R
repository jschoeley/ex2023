libs <- c(
  "ggplot2",
  "ISOweek",
  "purrr",
  "dplyr",
  "tidyr",
  "yaml",
  "glue",
  "httr",
  "readr",
  "eurostat",
  "gridExtra",
  "stringr",
  "ungroup",
  "lubridate",
  "openxlsx",
  "abind",
  "StMoMo",
  "tidyverse",
  "grImport2",
  #"ggflags",
  "scales",
  "devtools",
  "Cairo",
  "cowplot"
)

install.packages(
  libs,
  repos = c(getOption('repos')),
  dep = TRUE
)

# apt install libcairo2-dev libspectre-dev librsvg2-dev libpoppler-glib-dev r-base-dev
#devtools::install_github("sjp/grConvert")
