# Plot pandemic population deficit plus

# Init ------------------------------------------------------------

library(yaml)
library(readr); library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(ggflags)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  global = './src/00-global.R',
  region = './cfg/region_metadata.csv',
  projections = './out/40-projections.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig = './out/',
  rds_popdeficit = './out/50-deficit.rds',
  xlsx_popdeficit = './out/50-deficit.xlsx'
)

# global configuration
config <- read_yaml(paths$input$config)

# global objects and functions
global <- source(paths$input$global)

# constants specific to this analysis
cnst <- within(list(), {
  region = filter(
    read_csv(paths$input$region),
    region_code_iso3166_2 %in% config$skeleton$regions
  )
})

dat <- list()

# Load data -------------------------------------------------------

dat$proj <- readRDS(paths$input$projections)

# Plot population deficit -----------------------------------------

popdeficit <- list()

popdeficit$data <-
  dat$proj |>
  group_by(year, region, nsim) |>
  summarise(
    expected_total = sum(population_dec31st_noncovidmortality),
    observed_total = sum(population_dec31st_covidmortality),
    delta_absolute = -(observed_total-expected_total),
    delta_relative = delta_absolute/expected_total
  ) |>
  group_by(year, region) |>
  summarise(
    delta_relative_Q05 = quantile(delta_relative, p = 0.05, na.rm = T),
    delta_relative_Q50 = quantile(delta_relative, p = 0.50, na.rm = T),
    delta_relative_Q95 = quantile(delta_relative, p = 0.95, na.rm = T),
    delta_absolute_Q05 = quantile(delta_absolute, p = 0.05, na.rm = T),
    delta_absolute_Q50 = quantile(delta_absolute, p = 0.50, na.rm = T),
    delta_absolute_Q95 = quantile(delta_absolute, p = 0.95, na.rm = T),
    expected_total_Q50 = quantile(expected_total, p = 0.50, na.rm = T),
    observed_total_Q50 = quantile(observed_total, p = 0.50, na.rm = T)
  ) |>
  ungroup()

popdeficit$fig <-
  popdeficit$data |>
  filter(year == 2023) |>
  mutate(
    region_ggflag = tolower(region),
    region_rank = rank(delta_relative_Q50)
  ) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) |>
  ggplot(aes(y = region_rank, yend = region_rank)) +
  geom_vline(aes(xintercept = 0), size = 0.5, color = 'grey80') +
  geom_segment(aes(x = delta_relative_Q05, xend = delta_relative_Q95),
               size = 1, color = 'grey70') +
  geom_text(
    aes(
      x = delta_relative_Q50,
      label = region_name_en
    ),
    position = position_nudge(y = -0.25, x = -0.0015), hjust = 1,
    size = 1.7, color = 'grey60'
  ) +
  geom_text(
    aes(
      x = delta_relative_Q50,
      label = paste0(
        formatC(delta_relative_Q50*100, digits = 2, format = 'f',
                decimal.mark = ','), '%'
      )
    ),
    position = position_nudge(y = +0.25, x = -0.0015), hjust = 1,
    size = 1.7, color = 'grey60'
  ) +
  geom_text(
    aes(
      x = delta_relative_Q50,
      label = paste0(
        formatC(delta_absolute_Q50, format = 'd',
                big.mark = ' ')
      )
    ),
    position = position_nudge(y = +0.25, x = +0.0015), hjust = 0,
    size = 1.7, color = 'grey60'
  ) +
  geom_point(aes(x = delta_relative_Q50), size = 5.5) +
  geom_flag(
    aes(x = delta_relative_Q50, country = region_ggflag), size = 5
  ) +
  scale_x_continuous(labels = ~scales::percent(.x, decimal.mark = '.'),
                     breaks = seq(0, 0.05, 0.0025)) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(1, 1))) +
  MyGGplotTheme(grid = 'x', axis = 'x') +
  labs(
    y = NULL,
    x = 'Pandemic population deficit December 31st 2023'
  )

popdeficit$fig

# Plot population deficit by sex ----------------------------------

popdeficitsex <- list()

popdeficitsex$data <-
  dat$proj |>
  group_by(sex, year, region, nsim) |>
  summarise(
    expected_total = sum(population_dec31st_noncovidmortality),
    observed_total = sum(population_dec31st_covidmortality),
    delta_absolute = -(observed_total-expected_total),
    delta_relative = delta_absolute/expected_total
  ) |>
  group_by(sex, year, region) |>
  summarise(
    delta_relative_Q05 = quantile(delta_relative, p = 0.05, na.rm = T),
    delta_relative_Q50 = quantile(delta_relative, p = 0.50, na.rm = T),
    delta_relative_Q95 = quantile(delta_relative, p = 0.95, na.rm = T),
    delta_absolute_Q05 = quantile(delta_absolute, p = 0.05, na.rm = T),
    delta_absolute_Q50 = quantile(delta_absolute, p = 0.50, na.rm = T),
    delta_absolute_Q95 = quantile(delta_absolute, p = 0.95, na.rm = T),
    expected_total_Q50 = quantile(expected_total, p = 0.50, na.rm = T),
    observed_total_Q50 = quantile(observed_total, p = 0.50, na.rm = T)
  ) |>
  ungroup()

popdeficitsex$fig <- list()
for (stratum in c('Female', 'Male')) {
  popdeficitsex$fig[[stratum]] <-
    popdeficitsex$data |>
    filter(sex == stratum) |>
    filter(year == 2023) |>
    mutate(
      region_ggflag = tolower(region),
      region_rank = rank(delta_relative_Q50)
    ) |>
    left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) |>
    ggplot(aes(y = region_rank, yend = region_rank)) +
    geom_vline(aes(xintercept = 0), size = 0.5, color = 'grey80') +
    geom_segment(aes(x = delta_relative_Q05, xend = delta_relative_Q95),
                 size = 1, color = 'grey70') +
    geom_text(
      aes(
        x = delta_relative_Q50,
        label = region_name_en
      ),
      position = position_nudge(y = -0.25, x = -0.0015), hjust = 1,
      size = 1.7, color = 'grey60'
    ) +
    geom_text(
      aes(
        x = delta_relative_Q50,
        label = paste0(
          formatC(delta_relative_Q50*100, digits = 1, format = 'f',
                  decimal.mark = ','), '%'
        )
      ),
      position = position_nudge(y = +0.25, x = -0.0015), hjust = 1,
      size = 1.7, color = 'grey60'
    ) +
    geom_text(
      aes(
        x = delta_relative_Q50,
        label = paste0(
          formatC(delta_absolute_Q50, format = 'd',
                  big.mark = ' ')
        )
      ),
      position = position_nudge(y = +0.25, x = +0.0015), hjust = 0,
      size = 1.7, color = 'grey60'
    ) +
    geom_point(aes(x = delta_relative_Q50), size = 5.5) +
    geom_flag(
      aes(x = delta_relative_Q50, country = region_ggflag), size = 5
    ) +
    scale_x_continuous(labels = ~scales::percent(.x, decimal.mark = ','),
                       breaks = seq(0, 0.05, 0.005)) +
    scale_y_continuous(breaks = NULL, expand = expansion(add = c(1, 1))) +
    MyGGplotTheme(grid = 'x', axis = 'x') +
    labs(
      y = NULL,
      x = paste0(stratum, ' pandemic population deficit December 31st 2023')
    )
}

popdeficitsex$fig$Female

popdeficitsex$fig$combined <-
  cowplot::plot_grid(popdeficitsex$fig$Male,
                     popdeficitsex$fig$Female, ncol = 2)
popdeficitsex$fig$combined

# Plot population deficit by age ----------------------------------

popdeficitage <- list()

popdeficitage$cnst <- list(
  agegroups = c(seq(0, 85, 5), 101)
)

popdeficitage$data <-
  dat$proj |>
  filter(year == 2023) |>
  mutate(
    agegroup =
           cut(age, popdeficitage$cnst$agegroups, include.lowest = TRUE)
  ) |>
  group_by(year, sex, region, agegroup, nsim) |>
  mutate(
    popcovid = sum(population_dec31st_covidmortality),
    popnoncovid = sum(population_dec31st_noncovidmortality),
    delta_absolute = -pmin(popcovid-popnoncovid, 0),
    delta_relative = delta_absolute/popnoncovid
  ) |>
  group_by(agegroup, year, sex, region) |>
  summarise(
    delta_relative_Q05 = quantile(delta_relative, p = 0.05, na.rm = T),
    delta_relative_Q50 = quantile(delta_relative, p = 0.50, na.rm = T),
    delta_relative_Q95 = quantile(delta_relative, p = 0.95, na.rm = T),
    delta_absolute_Q05 = quantile(delta_absolute, p = 0.05, na.rm = T),
    delta_absolute_Q50 = quantile(delta_absolute, p = 0.50, na.rm = T),
    delta_absolute_Q95 = quantile(delta_absolute, p = 0.95, na.rm = T),
    expected_Q50   = quantile(population_dec31st_noncovidmortality, p = 0.50, na.rm = T),
    observed_Q50   = quantile(population_dec31st_covidmortality, p = 0.50, na.rm = T)
  )

popdeficitage$data |>
  ungroup() |>
  select(agegroup, sex, region, delta_relative_Q50) |>
  pivot_wider(names_from = sex, values_from = delta_relative_Q50) |>
  mutate(Male = -Male) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) |>
  ggplot() +
  aes(y = agegroup) +
  geom_col(aes(x = Male), fill = 'blue') +
  geom_col(aes(x = Female), fill = 'orange') +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_discrete(labels = c(seq(0, 80, 5), '85+')) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  facet_wrap(~region_name_en, scales = 'free_x')

popdeficitage$data |>
  ungroup() |>
  select(agegroup, sex, region, delta_absolute_Q50) |>
  pivot_wider(names_from = sex, values_from = delta_absolute_Q50) |>
  mutate(Male = -Male) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) |>
  ggplot() +
  aes(y = agegroup) +
  geom_col(aes(x = Male), fill = 'blue') +
  geom_col(aes(x = Female), fill = 'orange') +
  scale_y_discrete(labels = c(seq(0, 80, 5), '85+')) +
  facet_wrap(~region_name_en, scales = 'free_x')


# Export ----------------------------------------------------------

ExportFigure(
  popdeficit$fig, paths$output$fig, filename = '41-popdeficit',
  width = 170, height = 170, dpi = 300, device = 'pdf'
)

write.xlsx(
  popdeficit$data, file = paths$output$xlsx_popdeficit,
  keepNA = TRUE, na.string = '.',
  firstRow = TRUE, firstCol = TRUE, overwrite = TRUE
)
