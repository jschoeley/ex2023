# Plot life expectancy

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
dat$lc <- readRDS('out/40-leecarter.rds')

# Calculate and plot e0 -------------------------------------------

# First wave countries: those with a pronounced period shock in 2020,
# followed by a recovery
# BE, CH, ES, FR, IT, SE
# Second wave countries: those with a pronounced period shock in 2021,
# followed by a recovery
# BG, CZ, EE, HR, HU, LT, PL, SK
# Delayed shock countries: those with little excess mortality in 2020/21
# but elevated mortality in 22/23
# DE, DK, FI, IS, NO
# New regime countries: those with near constant annual life expectancy
# deficits since 2020
# NL, PT, AT, SI

e0 <- list()

e0$data$lt <-
  dat$lt |>
  group_by(year, sex, region) |>
  transmute(
    x  = age_start,
    nx = age_width,
    mx = death/population_py,
    px = exp(-mx*nx),
    qx = 1-px,
    lx = head(cumprod(c(1, px)), -1),
    dx = c(-diff(lx), tail(lx, 1)),
    Lx = ifelse(mx==0, lx*nx, dx/mx),
    Tx = rev(cumsum(rev(Lx))),
    ex = Tx/lx
  ) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2')) %>%
  filter(year %in% 2000:2023) |>
  filter(x == 0) |>
  ungroup()

e0$data$lc <-
  dat$lc$predicted |>
  filter(age == '0') |>
  mutate(region = as.character(region),
         sex = as.character(sex)) |>
  group_by(region, sex, year) |>
  summarise(
    mean = mean(ex),
    sd = sd(ex)
  ) |>
  ungroup() |>
  filter(year %in% 2000:2023) |>
  left_join(e0$data$lt) |>
  mutate(
    ex_deficit = ex-mean,
    ex_deficit_zscore = ex_deficit/sd
  )

e0$data$lc |>
  group_by(region, sex) |>
  summarise(
    max_ex_deficit_year = which.min(ex_deficit)-1 + 2020,
    z_score_23 = ex_deficit[year == 2023],
    max_z_score = ex_deficit[year == max_ex_deficit_year]
  ) |>
  group_by(max_ex_deficit_year) |>
  summarise(
    n = n(),
    mean_z_max = mean(max_z_score),
    mean_z_2023 = mean(z_score_23)
  )

e0$fig <-
  e0$data$lt %>%
  ggplot(aes(x = year, group = sex, color = sex)) +
  geom_vline(xintercept = 2020, color = 'grey') +
  geom_ribbon(
    aes(x = year, ymin = mean-sd, ymax = mean+sd),
    color = NA, fill = 'grey85',
    data = e0$data$lc
  ) +
  geom_line(
    aes(x = year, y = mean),
    data = e0$data$lc
  ) +
  geom_line(aes(y = ex), data = . %>% filter(year < 2020)) +
  geom_point(aes(y = ex), data = . %>% filter(year < 2020), shape = 21, fill = 'white') +
  geom_point(aes(y = ex), data = . %>% filter(year >= 2020)) +
  scale_x_continuous(
    breaks = seq(2000, 2023, 1),
    labels = c('2000', rep('', 22), '2023'),
    limits = c(2010, 2023),
  ) +
  scale_y_continuous(breaks = seq(70, 90, 2)) +
  scale_color_manual(values = unlist(config$figspec$colors$sex)) +
  MyGGplotTheme(grid = 'y', axis = 'x', panel_border = TRUE) +
  labs(
    y = 'Period life expectancy', x = NULL
  ) +
  facet_wrap(~region, scales = 'free_y')

e0$fig

# Export ----------------------------------------------------------

ExportFigure(
  e0$fig, paths$output$fig, filename = '51-e0',
  width = 170, height = 140, device = 'pdf', scale = 1.4
)

write.xlsx(e0$data, file = paths$output$xlsx_e0,
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE, overwrite = TRUE)
