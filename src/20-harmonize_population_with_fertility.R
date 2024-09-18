# Harmonize population and fertility data
#
# (1) Jan 1st and midyear population estimates for all nation states
#     from World Population Prospects 22 and 24
# (2) Age specific fertility rates from WPP 24
# (3) Mid-year population estimates for England & Wales,
#     Northern Ireland and Scotland from HMD/HFD
#     - years <= 2021 from HMD midyear pop estimates
#     - years 2022-2024 projected based on HMD and HFD population,
#       mortality and fertility data

# Init ------------------------------------------------------------

library(yaml)
library(readr); library(dplyr); library(tidyr)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  config = './cfg/config.yaml',
  global = './src/00-global.R',
  region_meta = './cfg/region_metadata.csv',
  skeleton = './tmp/10-harmonized_skeleton.rds',
  wpp_population = './dat/wppjoint/11-wppjoint.rds',
  hmdhfdgb = './dat/hmdhfd/11-hmdhfdgb.rds'
)
paths$output <- list(
  harmonized_population = './tmp/20-harmonized_population.rds',
  out = './out'
)

# global configuration
config <- read_yaml(paths$input$config)

# meta data on regions
region_meta <- read_csv(paths$input$region_meta, na = '.')

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  # translation of wpp sex code to harmonized sex code
  code_sex_wpp =
    c(`male` = config$skeleton$sex$Male,
      `female` = config$skeleton$sex$Female)
  # translation of hmdhfd sex code to harmonized sex code
  code_sex_hmdhfd =
    c(`Male` = config$skeleton$sex$Male,
      `Female` = config$skeleton$sex$Female)
  # lookup table for region codes
  # only countries defined in skeleton
  region_lookup_wpp =
    region_meta |>
    filter(region_code_iso3166_2 %in% config$skeleton$region) |>
    select(region_code_iso3166_2, region_code_wpp) |>
    drop_na()
  region_lookup_hmd =
    region_meta |>
    filter(region_code_iso3166_2 %in% config$skeleton$region) |>
    select(region_code_iso3166_2, region_code_hmd)

  # first year in harmonized data set
  skeleton_first_year = config$skeleton$year$start
  # last year in harmonized data set
  skeleton_final_year = config$skeleton$year$end

  # population projection jump-off year
  uk_projection_jumpoff_year = 2021
  # number of years to project GB population data from 2018
  uk_projection_n_years = skeleton_final_year-uk_projection_jumpoff_year

})

# list containers for analysis artifacts
dat <- list()
fig <- list()

# Functions -------------------------------------------------------

# global functions
source(paths$input$global)

#' Project a Population via Stable Assumption
StableProjection1x1 <- function (
    pop_m, pop_f, mx_m, mx_f, fx_f, n = 1, srb = 1.04
) {

  nage = length(pop_m)

  # first project the female population for n steps
  # Leslie projection matrix
  A_f <- matrix(0, nrow = nage, ncol = nage)
  diag(A_f[-1,-nage]) <- head(exp(-mx_f), -1)
  A_f[1,] <- fx_f * 1/(1+srb)
  # population matrix, first column is jump off population
  P_f <- matrix(0, nrow = nage, ncol = n+1)
  P_f[,1] <- pop_f
  # project population
  for (i in 1:n+1) {
    P_f[,i] <- A_f%*%P_f[,i-1]
  }

  # now project the male population with male births derived from
  # projected female births via sex ratio
  A_m <- matrix(0, nage, ncol = nage)
  diag(A_m[-1,-nage]) <- head(exp(-mx_m), -1)
  P_m <- matrix(0, nrow = nage, ncol = n+1)
  P_m[,1] <- pop_m
  for (i in 1:n+1) {
    P_m[,i] <- A_m%*%P_m[,i-1]
    P_m[1,i] <- P_f[1,i]*srb
  }

  projection <-
    expand_grid(n = 1:n, age_start = 1:nage-1) |>
    mutate(Female = c(P_f[,-1]), Male = c(P_m[,-1])) |>
    pivot_longer(
      cols = c(Female, Male),
      names_to = 'sex',
      values_to = 'population_midyear'
    )

  return(projection)

}

#' Make a Grid of Population Pyramids
PopPyramids <- function(
    dat, population, age, sex, year, highlight, facet, title
) {
  require(ggplot2); require(dplyr)
  dat |>
    transmute(
      pop = ifelse(sex == 'Male', -{{population}}, {{population}}),
      age = {{age}}, sex = {{sex}}, year = {{year}}, hl = {{highlight}},
      fct = {{facet}}
    ) |>
    ggplot(
      aes(x = age, y = pop, color = hl,
          group = interaction(sex, year, hl))
    ) +
    geom_line(show.legend = FALSE) +
    annotate('text', x = 90, y = -Inf, label = '\u2642',
             hjust = -1, size = 6) +
    annotate('text', x = 90, y = Inf, label = '\u2640',
             hjust = 1, size = 6) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      labels = function (x) ifelse(x == 100, '100+', x)
    ) +
    scale_y_continuous(
      labels = function(x) {formatC(abs(x*1e-3), format = 'd')}
    ) +
    coord_flip() +
    facet_wrap(~fct, scales = 'free_x') +
    labs(x = 'Age', y = 'Population in 1000s')
}

# Load data -------------------------------------------------------

dat$skeleton <- readRDS(paths$input$skeleton)

dat$wpp <- readRDS(paths$input$wpp_population)

dat$hmdhfdgb <- readRDS(paths$input$hmdhfdgb)

# Harmonize WPP population ----------------------------------------

dat$wpp_clean <-
  dat$wpp |>
  # select columns of interest
  select(
    population_source,
    region = ISO2_code,
    iso_year = Time, age = AgeGrpStart,
    male = Male, female = Female
  ) |>
  # sex to long format
  pivot_longer(
    cols = c(female, male),
    names_to = 'sex',
    values_to = 'value'
  ) |>
  # ensure proper names of factor variables
  mutate(
    sex =
      factor(sex, levels = names(cnst$code_sex_wpp),
             labels = cnst$code_sex_wpp) |>
      as.character()
  ) |>
  # add row id
  mutate(id = GenerateRowID(region, sex, iso_year, age)) |>
  # make popjan1st, popmidyear, fertility different columns
  separate(
    population_source, into = c('population_source', 'measure', 'timeframe'),
    sep = '_'
  ) |>
  mutate(
    measure = paste(measure, population_source, sep = '_')
  ) |>
  select(-timeframe, -population_source) |>
  pivot_wider(
    names_from = measure, values_from = value
  ) |>
  mutate(
    # WPP scales popnumbers in 1000's, so we scale back
    popjan1st_wpp24 = popjan1st_wpp24*1000,
    popmidyear_wpp24 = popmidyear_wpp24*1000,
    popjan1st_wpp22 = popjan1st_wpp22*1000,
    popmidyear_wpp22 = popmidyear_wpp22*1000,
    # set fertility to 0 where we know it to be essentially 0
    # also remove the x1000 scaling
    fertility = ifelse(sex == 'Male' | age < 15 | age > 49,
                       0, fertility_wpp24/1000),
    fertility_source = 'wpp24'
  )

# Check WPP data --------------------------------------------------

# check if data makes sense
fig$wpp_population <-
  PopPyramids(
    dat = dat$wpp_clean,
    population = popmidyear_wpp24,
    age = age, sex = sex, year = iso_year, highlight = iso_year,
    facet = region
  ) +
  scale_color_viridis_c() +
  labs(title = 'WPP midyear population estimates 1990-2024') +
  MyGGplotTheme(scaler = 0.8, panel_border = FALSE) +
  theme(panel.background = element_rect(color = NA, fill = 'grey95'))
fig$wpp_population

fig$wpp24vs22 <-
  dat$wpp_clean |>
  filter(sex == 'Female') |>
  filter(iso_year %in% c(2020:2023)) |>
  mutate(delta = popmidyear_wpp24/popmidyear_wpp22) |>
  ggplot() +
  aes(x = age, y = delta, color = as.factor(iso_year)) +
  geom_hline(yintercept = 1, color = 'grey40') +
  geom_line() +
  scale_color_viridis_d() +
  facet_wrap(~region) +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.9, 1.1)) +
  MyGGplotTheme(panel_border = FALSE) +
  theme(panel.background = element_rect(color = NA, fill = 'grey95')) +
  labs(
    title = 'Ratio of female midyear population WPP 24 vs. 22',
    y = 'Ratio', x = 'Age', color = 'Year'
  )
fig$wpp24vs22

# contrast wpp24 numbers with wpp22
dat$wpp_clean |>
  filter(iso_year == 2023, sex == 'Male', region == 'BG') |>
  select(region, iso_year, age, sex, popmidyear_wpp22, popmidyear_wpp24) |>
  mutate(ratio = popmidyear_wpp24/popmidyear_wpp22)
dat$wpp_clean |>
  filter(iso_year == 2023, sex == 'Male', region == 'BG') |>
  select(region, iso_year, age, sex, popmidyear_wpp22, popmidyear_wpp24) |>
  mutate(ratio = popmidyear_wpp24/popmidyear_wpp22)

# Harmonize pop data for UK regions -------------------------------

# for years not yet in the data we get population estimates via a
# Leslie-Matrix projection of midyear population assuming a stable
# population, i.e. no migration and constant fertility and mortality
# rates. we use the 'female dominant' projection of the male population.

# prepare a data set with required variables for
# Leslie matrix population projection for GB sub-regions
dat$gb_pop_estimates <-
  # skeleton
  expand_grid(
    year = 2010:cnst$uk_projection_jumpoff_year,
    region_code_hmd = c('GBRTENW', 'GBR_NIR', 'GBR_SCO'),
    sex = c('Female', 'Male'),
    age_start = 0:110
  ) |>
  # midyear population
  left_join(
    dat$hmdhfdgb$hmd_pop |>
      pivot_longer(
        c(Female, Male),
        names_to = 'sex', values_to = 'population_midyear'
      ) |>
      select(region_code_hmd, year = Year, age_start = Age,
             sex, population_midyear)
  ) |>
  # death rates
  left_join(
    dat$hmdhfdgb$hmd_mort |>
      pivot_longer(
        c(Female, Male),
        names_to = 'sex', values_to = 'death_rate'
      ) |>
      select(region_code_hmd, year = Year, age_start = Age, sex, death_rate)
  ) |>
  # female fertility rates
  left_join(
    dat$hmdhfdgb$hfd_fert |>
      select(region_code_hmd, year = Year, age_start = Age,
             female_fertility_rate = ASFR)
  ) |>
  mutate(
    # fertility rates are NA outside the age range [12, 55],
    # replace with 0
    female_fertility_rate =
      ifelse(is.na(female_fertility_rate), 0, female_fertility_rate),
    # same with death rates in highest ages
    death_rate =
      ifelse(is.na(death_rate), 0, death_rate)
  )

# replace missing fertility rates in jump-off year 2021
# with rates from 2019
dat$gb_pop_estimates[
  dat$gb_pop_estimates$year == 2021 &
    dat$gb_pop_estimates$sex == 'Female',
  'female_fertility_rate'
] <-
  dat$gb_pop_estimates |>
  filter(year == 2019, sex == 'Female') |>
  pull(female_fertility_rate)

# perform the projections by sub-region
# start 2021 and project 2022, 23, 24
dat$gb_pop_projections <-
  dat$gb_pop_estimates |>
  filter(year == cnst$uk_projection_jumpoff_year) |>
  group_by(region_code_hmd) |>
  group_modify(~{

    male <- filter(.x, sex == 'Male')
    female <- filter(.x, sex == 'Female')

    projection <- StableProjection1x1(
      pop_m = male$population_midyear, pop_f = female$population_midyear,
      mx_m = male$death_rate, mx_f = female$death_rate,
      fx_f = female$female_fertility_rate,
      n = cnst$uk_projection_n_years, srb = 1.04
    )

    return(projection)

  }) |>
  ungroup() |>
  mutate(year = cnst$uk_projection_jumpoff_year+n)

# bind and harmonize estimates and projections
dat$gb_clean <-
  bind_rows(
    hmd_estimates = select(
      dat$gb_pop_estimates,
      region_code_hmd, year, sex, age_start, population_midyear
    ),
    hmd_projections = select(
      dat$gb_pop_projections,
      region_code_hmd, year, age_start, sex, population_midyear
    ),
    .id = 'population_source'
  ) |>
  # make age 100 an open age group
  mutate(
    age_start = ifelse(age_start > 100, 100, age_start)
  ) |>
  group_by(population_source, region_code_hmd, year, sex, age_start) |>
  summarise(population_midyear = sum(population_midyear)) |>
  ungroup() |>
  # derive jan1st population from midyear population
  group_by(region_code_hmd, sex, age_start) |>
  arrange(region_code_hmd, age_start, sex, year) |>
  mutate(
    population_jan1st = (population_midyear+lag(population_midyear))/2
  ) |>
  ungroup() |>
  # ensure proper names of factor variables
  mutate(
    sex =
      factor(sex, levels = names(cnst$code_sex_hmdhfd),
             labels = cnst$code_sex_hmdhfd) |>
      as.character(),
    region_iso = factor(
      region_code_hmd,
      levels = cnst$region_lookup_hmd$region_code_hmd,
      labels = cnst$region_lookup_hmd$region_code_iso3166_2
    ) |> as.character()
  ) |>
  # add row id
  mutate(id = GenerateRowID(region_iso, sex, year, age_start))

# check if projection went well
fig$gb_midyear_population1 <-
  dat$gb_clean |>
  filter(year %in% c(2018:2024)) |>
  PopPyramids(
    population = population_midyear,
    age = age_start, sex = sex, year = year,
    highlight = population_source, facet = region_iso
  ) +
  scale_color_manual(values = c('grey50', 'red')) +
  labs(title = 'HMD population estimates 2018-21 (grey) and stable population projection 2022-24 (red) for U.K. regions') +
  MyGGplotTheme(scaler = 0.8) +
  theme(panel.background = element_rect(color = NA, fill = 'grey95'))
fig$gb_midyear_population1

fig$gb_midyear_population2 <-
  dat$gb_clean |>
  group_by(year, population_source, region_iso) |>
  summarise(
    N = sum(population_midyear)/1000
  ) |>
  ggplot(aes(x = year, y = N, color = population_source)) +
  geom_point() +
  scale_x_continuous(breaks = seq(2010, 2024, 2)) +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_color_manual(values = c('grey50', 'red')) +
  facet_wrap(~region_iso, scales = 'free_y') +
  MyGGplotTheme(scaler = 0.8) +
  theme(
    panel.background = element_rect(color = NA, fill = 'grey95'),
    legend.position = c(0.9, 0.2)
  ) +
  labs(
    x = NULL, y = 'Population in 1000s',
    color = 'Population source'
  )
fig$gb_midyear_population2

# Join population data with skeleton ------------------------------

# join the different sources of population count data
# with the skeleton
dat$pop_joined <-
  dat$skeleton |>
  left_join(
    select(
      dat$wpp_clean,
      id,
      fertility_rate_wpp24 = fertility,
      population_jan1st_wpp22 = popjan1st_wpp22,
      population_jan1st_wpp24 = popjan1st_wpp24,
      population_midyear_wpp22 = popmidyear_wpp22,
      population_midyear_wpp24 = popmidyear_wpp24

    ),
    by = 'id'
  ) |>
  left_join(
    select(
      dat$gb_clean,
      id,
      population_midyear_hmd = population_midyear,
      population_jan1st_hmd = population_jan1st,
      population_source_hmd = population_source
    ),
    by = 'id'
  )

# Choose default sources ------------------------------------------

dat$pop_joined <-
  dat$pop_joined |>
  mutate(
    # midyear population: use wpp estimates unless they are missing (as for GB regions)
    population_midyear = case_when(
      !is.na(population_midyear_wpp24) ~ round(population_midyear_wpp24, 1),
      !is.na(population_midyear_hmd) ~ round(population_midyear_hmd, 1),
      TRUE ~ as.numeric(NA)
    ),
    population_midyear_source = case_when(
      !is.na(population_midyear_wpp24) ~ 'wpp24',
      !is.na(population_midyear_hmd) ~ population_source_hmd,
      TRUE ~ as.character(NA)
    ),
    # jan 1st population: use wpp estimates unless they are missing (as for GB regions)
    population_jan1st = case_when(
      !is.na(population_jan1st_wpp24) ~ round(population_jan1st_wpp24, 1),
      !is.na(population_jan1st_hmd) ~ round(population_jan1st_hmd, 1),
      TRUE ~ as.numeric(NA)
    ),
    population_jan1st_source = case_when(
      !is.na(population_jan1st_wpp24) ~ 'wpp24',
      !is.na(population_jan1st_hmd) ~ population_source_hmd,
      TRUE ~ as.character(NA)
    )
  )

# Sum female and male counts to totals ----------------------------

# # populate sex category "total"
# dat$row_female <- which(dat$pop_joined$sex == 'Female')
# dat$row_male <- which(dat$pop_joined$sex == 'Male')
# dat$row_total <- which(dat$pop_joined$sex == 'Total')
#
# dat$population_varnames <-
#   c(
#     'population_midyear', 'population_jan1st',
#     'population_midyear_wpp22', 'population_midyear_wpp24', 'population_midyear_hmd',
#     'population_jan1st_wpp22', 'population_jan1st_wpp24', 'population_jan1st_hmd'
#   )
#
# for (var in dat$population_varnames) {
#   dat$pop_joined[[var]][dat$row_total] <-
#     dat$pop_joined[[var]][dat$row_female] +
#     dat$pop_joined[[var]][dat$row_male]
# }

# Select variables of interest for export -------------------------

dat$export <-
  dat$pop_joined |>
  # select output variables of interest
  select(
    id,
    population_midyear, population_midyear_source,
    population_jan1st, population_jan1st_source,
    fertility_rate = fertility_rate_wpp24,
    population_midyear_wpp22, population_midyear_wpp24, population_midyear_hmd,
    population_jan1st_wpp22, population_jan1st_wpp24, population_jan1st_hmd
  )

# Export ----------------------------------------------------------

saveRDS(dat$export, file = paths$output$harmonized_population)

ExportFigure(
  fig$wpp_population, path = paths$output$out, device = 'pdf',
  filename = '20-wpp_population',
  width = config$figspec$dimensions$width,
  height = 1.3*config$figspec$dimensions$width
)

ExportFigure(
  fig$wpp24vs22, path = paths$output$out, device = 'pdf',
  filename = '20-wpp24vs22',
  width = config$figspec$dimensions$width,
  height = 1.3*config$figspec$dimensions$width
)

ExportFigure(
  fig$gb_midyear_population1, path = paths$output$out, device = 'pdf',
  filename = '20-gb_midyear_population1',
  width = config$figspec$dimensions$width,
  height = 0.5*config$figspec$dimensions$width
)

ExportFigure(
  fig$gb_midyear_population2, path = paths$output$out, device = 'pdf',
  filename = '20-gb_midyear_population2',
  width = config$figspec$dimensions$width,
  height = 0.5*config$figspec$dimensions$width
)
