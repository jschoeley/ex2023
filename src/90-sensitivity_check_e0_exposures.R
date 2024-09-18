# Plot life expectancy

# Init ------------------------------------------------------------

library(yaml)
library(readr); library(dplyr); library(openxlsx)
library(purrr)
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
  projectioninput = './out/30-projectioninput.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig = './out/',
  xlsx_e0 = './out/51-e0.xlsx'
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

# Function --------------------------------------------------------

CalculateLifeTable <-
  function (df, x, nx, Dx, Ex) {

    require(dplyr)

    df %>%
      transmute(
        x = {{x}},
        nx = {{nx}},
        mx = {{Dx}}/{{Ex}},
        px = exp(-mx*{{nx}}),
        qx = 1-px,
        lx = head(cumprod(c(1, px)), -1),
        dx = c(-diff(lx), tail(lx, 1)),
        Lx = ifelse(mx==0, lx*nx, dx/mx),
        Tx = rev(cumsum(rev(Lx))),
        ex = Tx/lx
      )

  }

# Load data -------------------------------------------------------

dat$lt <- readRDS(paths$input$projectioninput)

# Calculate and plot e0 -------------------------------------------

e0wppsens <- list()

e0wppsens$data <-
  dat$lt |>
  group_by(year, sex, region) |>
  group_modify(~{
    e0_wpp_22 <-
      CalculateLifeTable(.x, x = age_start, nx = age_width, Dx = death, Ex = population_py_wpp22) |>
      filter(x == 0) |>
      pull(ex)
    e0_wpp_24 <-
      CalculateLifeTable(.x, x = age_start, nx = age_width, Dx = death, Ex = population_py_wpp24) |>
      filter(x == 0) |>
      pull(ex)
    tibble(e0_wpp_22, e0_wpp_24)
  }) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) |>
  filter(year %in% 2000:2023)

e0wppsens$fig <-
  e0wppsens$data %>%
  ggplot(aes(x = year, group = sex, color = sex, fill = sex)) +
  geom_vline(xintercept = 2020, color = 'grey') +
  geom_point(aes(y = e0_wpp_22), shape = 21, fill = 'white') +
  geom_point(aes(y = e0_wpp_24), shape = 21, size = 0.5) +
  geom_point(aes(y = lifeexpectancy_eurostat),
             data = dat$lt %>% filter(age_start == 0),
             shape = 3) +
  scale_x_continuous(
    breaks = seq(2000, 2023, 1),
    labels = c('2000', rep('', 22), '2023'),
    limits = c(2000, 2023),
  ) +
  scale_y_continuous(breaks = seq(70, 90, 2)) +
  scale_color_manual(values = unlist(config$figspec$colors$sex)) +
  scale_fill_manual(values = unlist(config$figspec$colors$sex)) +
  MyGGplotTheme(grid = 'y', axis = 'x', panel_border = TRUE) +
  labs(
    y = 'Period life expectancy', x = NULL, color = NULL, fill = NULL,
    title = 'Life expectancy under WPP 22 (open circle) and WPP 24 (closed) population estimates',
    subtitle = 'Crosses mark the Eurostat reported life expectancy'
  ) +
  facet_wrap(~region, scales = 'free_y', ncol = 4)

e0wppsens$fig

# Export ----------------------------------------------------------

ExportFigure(
  e0wppsens$fig, paths$output$fig, filename = '90-e0wppsens',
  width = 170, height = 140, device = 'pdf', scale = 1.4
)
