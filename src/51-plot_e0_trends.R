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

# Load data -------------------------------------------------------

dat$lt <- readRDS('out/50-lifetables.rds')
dat$lc <- readRDS('out/50-arriaga_cntfc.rds')

dat$lc_sim <- readRDS('tmp/50-arriaga_cntfc_sim.rds')

dimnames(dat$lc_sim)


# Calculate p-value of e0 deficit ---------------------------------

pval <- list()

# simulations of actual and expected e0 by region, sex, and year
# sim position 1 reserved for means over simulations
pval$sim <-
  dat$lc_sim['0',as.character(2020:2024),,,,c('ex_actual', 'ex_expected')]
# actual and expected mean
pval$mean <- apply(pval$sim[,,,-1,], c(1:3, 5), mean)

# calculate distribution of test statistic (life expectancy deficit)
# under H0: continuation of pre-pandemic mortality trends.
pval$test <- pval$sim[,,,,'ex_expected']
for (sim in 2:251) {
  pval$test[,,,sim] <-
    abs(pmin(
      0, pval$sim[,,,sim,'ex_expected'] - pval$mean[,,,'ex_expected']
    ))
}
# add observed test statistic
pval$test[,,,1] <- abs(pmin(0,
  pval$mean[,,,'ex_actual'] - pval$mean[,,,'ex_expected']
))

# calculate probability of having at least the observed effect under H0
pval$p <- apply(pval$test, 1:3, function (x) {
  if (any(is.na(x)) | any(is.nan(x))) {
    p <- NA
  } else {
    p <- 1-ecdf(x[-1])(x[1])
    p[x[1]<1e-12] <- 1
  }
  return(p)
})

pval$df <- as.data.frame.table(pval$p, stringsAsFactors = FALSE)
names(pval$df) <- c('year', 'sex', 'region', 'e0_deficit_pval')
pval$df$year <- as.integer(pval$df$year)
ecdf(pval$test['2023','Total','LU',]) |> plot(verticals = TRUE)

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
  filter(scenario == 'actual', sex == 'Total', age == 0, year %in% 2010:2024) |>
  select(year, sex, region = region_iso, e0_actual = ex_mean)

e0$data$lc <-
  dat$lc |> filter(age == '0') |>
  mutate(year = as.integer(year)) |>
  select(
    region = region_iso, sex, year,
    e0_expected_avg = ex_expected_mean,
    e0_expected_q050 = ex_expected_q0.05,
    e0_expected_q950 = ex_expected_q0.95,
    e0_deficit_avg = ex_actual_minus_expected_mean,
    e0_deficit_q050 = ex_actual_minus_expected_q0.05,
    e0_deficit_q950 = ex_actual_minus_expected_q0.95
  )

e0$data$combine <-
  left_join(e0$data$lt, e0$data$lc, by = c('year', 'sex', 'region')) |>
  left_join(pval$df, by = c('year', 'sex', 'region')) |>
  left_join(cnst$region, by = c('region' = 'region_code_iso3166_2'))

e0$fig <-
  e0$data$combine |>
  ggplot(aes(x = year)) +
  geom_vline(xintercept = 2019.5, color = 'grey') +
  geom_ribbon(
    aes(x = year, ymin = e0_expected_q050, ymax = e0_expected_q950),
    color = NA, fill = 'grey80'
  ) +
  geom_line(
    aes(x = year, y = e0_expected_avg),
  ) +
  geom_line(
    aes(y = e0_actual),
    data = . %>% filter(year < 2020)
  ) +
  geom_point(aes(y = e0_actual), data = . %>% filter(year < 2020),
             shape = 21, fill = 'white') +
  geom_point(
    aes(y = e0_actual),
    data = . %>% filter(year >= 2020)
  ) +
  scale_x_continuous(
    breaks = seq(2010, 2024, 1),
    labels = c('2010', rep('', 9), "'20", rep('', 3), "'24"),
    limits = c(2010, 2024),
  ) +
  scale_y_continuous(breaks = seq(70, 90, 1)) +
  MyGGplotTheme(grid = 'y', axis = 'x', panel_border = FALSE) +
  labs(
    y = 'Period life expectancy', x = NULL
  ) +
  facet_wrap(~region_name_en, ncol = 5, scales = 'free_y')

e0$fig

# Export ----------------------------------------------------------

ExportFigure(
  e0$fig, paths$output$fig, filename = '51-e0',
  width = 170, height = 140, device = 'pdf', scale = 1.4
)

write.xlsx(e0$data$combine, file = paths$output$xlsx_e0,
           keepNA = TRUE, na.string = '.',
           firstRow = TRUE, firstCol = TRUE, overwrite = TRUE)
