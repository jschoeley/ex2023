# Plot age decomposition of e0 deficit

# Init ------------------------------------------------------------

library(yaml)
library(readr); library(dplyr); library(openxlsx)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  global = './src/00-global.R',
  region = './cfg/region_metadata.csv',
  e0_deficit_decomp = './out/50-arriaga_cntfc.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig = './out/',
  xlsx_e0 = './out/53-e0deficitage.xlsx'
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

dat$lt <- readRDS(paths$input$e0_deficit_decomp)

# Plot e0 deficits by age -----------------------------------------

e0 <- list()

e0$cnst <- list(
  colors = c('2020' = '#CCD800', '2021' = '#D89E00', '2022' = '#C05B00', '2023' = '#C02C00', '2024' = 'black')
)

e0$data$lt <-
  dat$lt |>
  filter(sex == 'Total') |>
  left_join(cnst$region, by = c('region_iso' = 'region_code_iso3166_2')) |>
  filter(year %in% c('2020', '2021', '2022', '2023', '2024')) |>
  mutate(age = as.integer(age), age_group = (age %/% 20)*20) |>
  ungroup() |>
  select(region_iso, region_name_en, year, age_group, e0_cntrb_t_mean) |>
  group_by(region_iso, region_name_en, year, age_group) |>
  summarise(
    e0_cntrb_t_mean = sum(e0_cntrb_t_mean)
  )

e0$fig <-
  e0$data$lt |>
  ggplot() +
  aes(x = e0_cntrb_t_mean, y = age_group) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_path(aes(color = year)) +
  facet_wrap(~region_name_en, ncol = 5) +
  scale_color_manual(values = e0$cnst$colors) +
  scale_y_continuous(
    breaks = c(0, 20, 40, 60, 80, 100), labels = c('0-19', '20-39', '40-59', '60-79', '80-99', '100+')
  ) +
  scale_x_continuous(breaks = seq(-2, 0.5, 0.5), labels = c('-2', '-1.5', '-1', '-0.5', '0', '+0.5')) +
  labs(
    x = 'Years of contribution to total LE deficit', y = 'Age group', color = 'Year',
    title = 'Age-specific contributions to total LE deficit since 2020 by year'
  ) +
  MyGGplotTheme() +
  theme(panel.background = element_rect(fill = 'grey97', color = NA))

e0$fig

# Export ----------------------------------------------------------

ExportFigure(
  e0$fig, paths$output$fig, filename = '53-e0deficitage',
  width = 170, height = 140, device = 'pdf', scale = 1.4
)

write.xlsx(e0$data, file = paths$output$xlsx_e0,
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE, overwrite = TRUE)
