# Estimate life tables and associated statistics with CIs
#
# Simulation based life table inference based on Poisson samples of
# observed death counts. Implemented in a huge array.

# Init ------------------------------------------------------------

library(glue); library(yaml)
library(dplyr); library(tidyr); library(readr)

# Constants -------------------------------------------------------

# randomness in Poisson sampling
set.seed(1987)

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  projectioninput = './out/30-projectioninput.rds',
  leecarter = './out/40-leecarter.rds',
  figspecs = './src/00-global.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  lifetables = './out/50-lifetables.rds',
  lifetables_csv = './tmp/50-lifetables.csv',
  lifetables_sim = './tmp/50-lifetables_sim.rds',
  arriaga_cntfc = './out/50-arriaga_cntfc.rds',
  arriaga_cntfc_csv = './tmp/50-arriaga_cntfc.csv',
  arriaga_cntfc_sim = './tmp/50-arriaga_cntfc_sim.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- list(); cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
  # number of Poisson life-table replicates
  n_sim = 250
  # quantiles for CI's
  quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
  forecast_period = seq(
    config$forecast$jumpoff,
    config$forecast$jumpoff+config$forecast$h-1
  )
})

tmp <- list()

# Function --------------------------------------------------------

source(paths$input$figspecs)

# this function returns TRUE wherever elements are the same,
# including NA's, and FALSE everywhere else
compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# return a vector of mean and quantiles
QuantileWithMean <- function (x, prob = cnst$quantiles) {
  x <- x[!(is.na(x)|is.nan(x)|is.infinite(x))]
  q <- quantile(x, prob = prob, names = FALSE, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  result <- c(m, q)
  names(result) <- c('mean', paste0('q', prob))
  return(result)
}

# interpolate a range of values with its mean
InterpolateWithMean <- function (x) {
  rep(mean(x), length(x))
}

# Data ------------------------------------------------------------

lt_input <- list()

# input data for life-table calculation
# harmonized death counts and population exposures with open age group 100+
lt_input$openage_100 <- readRDS(paths$input$projectioninput)

# harmonized counterfactual Lee-Carter forecasts
lt_input$leecarter <- readRDS(paths$input$leecarter)

# Create Poisson replicates of counts -----------------------------

# create life table replicates by region, sex, and year
# based on repeatedly sampling death counts from a Poisson
# distribution with mean equal to estimated mean from PCLM

lifetables <- list()
lifetables$cnst <- list(
  nage = length(unique(lt_input$openage_100$age_start)),
  nyears = config$skeleton$year$end-config$skeleton$year$start+1,
  nregions = length(config$skeleton$region)
)

lifetables$input <-
  lt_input$openage_100 |>
  select(
    id, region, sex, year, age_start, age_width,
    population_py, death_observed = death,
    nweeks_year
  )
# vector ordered by sex, region, year, age distributed into
# 7D array [age, year, sex, region_id, sim_id, var_id, scenario]
# DESCRIPTION OF VAR ID'S
# LIFE TABLE
# (1)  <population_py>    person years exposure
# (2)  <death_total>      total death count
# (3)
# (4)  <nmx>              death rate
# (5)  <npx>              conditional probability of surviving age x
# (6)  <nqx>              conditional probability of dying within age x
# (7)  <lx>               probability of surviving to age x
# (8)  <ndx>              probability of dying within age x
# (9)  <nLx>              life table person years of exposure in age x
# (10) <Tx>               life table person years of exposure above age x
# (11) <ex>               life expectancy at age x
# ARRIAGA DECOMPOSITION LE CHANGE
# (12) <ex_lag>           1 year lag in ex
# (13) <ex_diff>          1 year difference in ex
# (14) <ex_diff_lag>      1 year lag in ex difference
# (15) <lx_lag>           1 year lag in lx
# (16) <Lx_lag>           1 year lag in Lx
# (17) <Tx_lag>           1 year lag in Tx
# (18) <e0_cntrb_d>       direct contribution of nmx changes to e0 changes
# (19) <e0_cntrb_i>       indirect contribution of nmx changes to e0 changes
# (20) <e0_cntrb_t>       total contrib. of nmx changes to e0 changes
lifetables$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      sex = 3,
      region_iso = lifetables$cnst$nregions,
      sim_id = cnst$n_sim+1, # first dimension is mean
      var_id = 20,
      scenario = 2
    ),
    dimnames = list(
      0:(lifetables$cnst$nage-1),
      config$skeleton$year$start:config$skeleton$year$end,
      c(unlist(config$skeleton$sex), 'Total'),
      config$skeleton$region,
      1:(cnst$n_sim+1),
      c('population_py', 'death_total', 'empty',
        'nmx', 'npx', 'nqx', 'lx', 'ndx', 'nLx', 'Tx', 'ex',
        'ex_lag', 'ex_diff', 'ex_diff_lag', 'lx_lag',
        'nLx_lag', 'Tx_lag',
        'e0_cntrb_d', 'e0_cntrb_i', 'e0_cntrb_t'),
      c('actual', 'projected')
    )
  )

# actual observables
lifetables$simulation[,,-3,,,'population_py','actual'] <-
  lifetables$input$population_py
lifetables$simulation[,,-3,,1,'death_total','actual'] <-
  lifetables$input$death_observed

# projected observables based on 5 year average nmx change
# here we store all-cause death counts how we would
# expect them under continuing pre-pandemic trends
lifetables$simulation[,,-3,,,'population_py','projected'] <-
  lifetables$simulation[,,-3,,,'population_py','actual']
lifetables$simulation[,as.character(cnst$forecast_period),
                      -3,,-1,'death_total','projected'] <-
  (
    lt_input$leecarter$harmonized_forecast_sim |>
      filter(year %in% cnst$forecast_period) |> pull(mx)
  ) *
  (
    lifetables$input |> filter(year %in% cnst$forecast_period) |>
      pull(population_py)
  )


# simulate observed total death counts
lifetables$simulation[,,-3,,-1,'death_total','actual'] <-
  apply(lifetables$simulation[,,-3,,1,'death_total','actual'],
        MARGIN = 1:4, function (lambda) rpois(n = cnst$n_sim, lambda),
        simplify = TRUE) |>
  aperm(c(2,3,4,5,1))

# check if actuals and projected somewhat line-up
lifetables$simulation[,'2020','Male','AT',2,'death_total',c('actual', 'projected')]

# Add sex-specific counts to total --------------------------------

# [age, year, sex, region_id, sim_id, var_id, scenario]
lifetables$simulation[,,'Total',,,c('death_total', 'population_py'),] <-
  lifetables$simulation[,,'Female',,,c('death_total', 'population_py'),] +
  lifetables$simulation[,,'Male',,,c('death_total', 'population_py'),]

# Calculate lifetables over simulated counts ----------------------

# nmx
lifetables$simulation[,,,,,'nmx',] <-
  lifetables$simulation[,,,,,'death_total',] /
  lifetables$simulation[,,,,,'population_py',]

# npx, using constant hazard assumption,
# i.e. npx = exp(-nmx)) for single year age groups
lifetables$simulation[,,,,,'npx',] <-
  exp(-lifetables$simulation[,,,,,'nmx',])
lifetables$simulation[lifetables$cnst$nage,,,,,'npx',] <- 0

# no need for fancy nax adjustment, I checked, virtually
# same results as with PWE
# nax <- lifetables$simulation[,,,,,'nmx',]
# nax[1:101,,,,,] <- 0.5
# I <- lifetables$simulation[1,,,'F',,'nmx',]<0.107
# I[is.na(I)] <- FALSE
# nax[1,,,'F',,][I] <-
#   0.053+2.8*lifetables$simulation[1,,,'F',,'nmx',][I]
# nax[1,,,'F',,][!I] <- 0.350
# I <- lifetables$simulation[1,,,'M',,'nmx',]<0.107
# I[is.na(I)] <- FALSE
# nax[1,,,'M',,][I] <-
#   0.045+2.684*lifetables$simulation[1,,,'M',,'nmx',][I]
# nax[1,,,'M',,][!I] <- 0.330
# lifetables$simulation[,,,,,'nqx',] <-
#   lifetables$simulation[,,,,,'nmx',] /
#   (1+(1-nax)*lifetables$simulation[,,,,,'nmx',])
# lifetables$simulation[lifetables$cnst$nage,,,,,'nqx',] <- 1
# lifetables$simulation[,,,,,'npx',] <-
#   1-lifetables$simulation[,,,,,'nqx',]

# nqx
lifetables$simulation[,,,,,'nqx',] <-
  1-lifetables$simulation[,,,,,'npx',]
# lx
lifetables$simulation[,,,,,'lx',] <-
  apply(
    lifetables$simulation[,,,,,'npx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (npx) head(cumprod(c(1, npx)), -1)
  )
# ndx
lifetables$simulation[,,,,,'ndx',] <-
  apply(
    lifetables$simulation[,,,,,'lx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (lx) c(-diff(lx), tail(lx, 1))
  )
# nLx = ifelse(mx==0, lx*nx, ndx/nmx)
lifetables$simulation[,,,,,'nLx',] <-
  lifetables$simulation[,,,,,'ndx',]/lifetables$simulation[,,,,,'nmx',]
tmp$I <- compareNA(lifetables$simulation[,,,,,'nmx',],0)
lifetables$simulation[,,,,,'nLx',][tmp$I] <-
  lifetables$simulation[,,,,,'lx',][tmp$I]
# Tx = rev(cumsum(rev(nLx)))
lifetables$simulation[,,,,,'Tx',] <-
  apply(
    lifetables$simulation[,,,,,'nLx',],
    # apply function to vector of data by age
    MARGIN = 2:6, function (nLx) rev(cumsum(rev(nLx)))
  )
# ex = Tx/lx
lifetables$simulation[,,,,,'ex',] <-
  lifetables$simulation[,,,,,'Tx',] /
  lifetables$simulation[,,,,,'lx',]

# Calculate annual ex change --------------------------------------

lifetables$simulation[,,,,,'ex_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex',]
lifetables$simulation[,,,,,'ex_diff',] <-
  lifetables$simulation[,,,,,'ex',]-lifetables$simulation[,,,,,'ex_lag',]
lifetables$simulation[,,,,,'ex_diff_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'ex_diff',]

# Calculate Arriaga decomposition of annual e0 changes ------------

# decompose annual changes in e0 into age specific mortality changes
# see Arriaga (1984)
# Measuring and explaining the change in life expectancies
# DOI 10.2307/2061029

lifetables$simulation[,,,,,'lx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'lx',]
lifetables$simulation[,,,,,'nLx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'nLx',]
lifetables$simulation[,,,,,'Tx_lag',] <-
  lifetables$simulation[,c(NA,1:(lifetables$cnst$nyears-1)),,,,'Tx',]

lifetables$simulation[,,,,,'e0_cntrb_d',] <-
  (
    lifetables$simulation[,,,,,'nLx',]/lifetables$simulation[,,,,,'lx',]-
      lifetables$simulation[,,,,,'nLx_lag',]/lifetables$simulation[,,,,,'lx_lag',]
  ) * lifetables$simulation[,,,,,'lx',]

lifetables$simulation[,,,,,'e0_cntrb_i',] <-
  (
    lifetables$simulation[,,,,,'lx_lag',]/
      lifetables$simulation[,,,,,'lx',]-
      apply(lifetables$simulation[,,,,,'lx_lag',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))/
      apply(lifetables$simulation[,,,,,'lx',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
  ) * apply(lifetables$simulation[,,,,,'Tx',],
            # apply function to vector of data by age
            2:6, function (x) c(x[-1], 0))
lifetables$simulation[lifetables$cnst$nage,,,,,'e0_cntrb_i',] <- 0

lifetables$simulation[,,,,,'e0_cntrb_t',] <-
  lifetables$simulation[,,,,,'e0_cntrb_d',] +
  lifetables$simulation[,,,,,'e0_cntrb_i',]

# test if age decomposition sums to observed e0 difference
arriage_e0diff_test <- abs(
  apply(
    lifetables$simulation[,,,,1,'e0_cntrb_t',], 2:5, sum
  ) -
    lifetables$simulation[1,,,,1,'ex_diff',]
)

# maximum decomposition error less than 0.1 years?
max(arriage_e0diff_test, na.rm = TRUE) < 0.1
# distribution of decomposition errors
hist(log(arriage_e0diff_test), breaks = 50)

# Counterfactual Arriaga decomposition ----------------------------

# decompose life expectancy deficit (observed minus expected e0)
# into age specific mortality changes
# see Arriaga (1984)
# Measuring and explaining the change in life expectancies
# DOI 10.2307/2061029

arriaga_cntfc <- list()
arriaga_cntfc$simulation <-
  array(
    dim = c(
      age = lifetables$cnst$nage,
      year = lifetables$cnst$nyears,
      sex = 3,
      region_iso = lifetables$cnst$nregions,
      sim_id = cnst$n_sim+1,
      var_id = 16
    ),
    dimnames = list(
      0:(lifetables$cnst$nage-1),
      config$skeleton$year$start:config$skeleton$year$end,
      c(unlist(config$skeleton$sex), 'Total'),
      config$skeleton$region,
      1:(cnst$n_sim+1),
      c('death_total_actual', 'death_total_expected',
        'nmx_actual', 'nmx_expected',
        'ex_actual_minus_expected',
        'ex_actual', 'ex_expected',
        'lx_actual', 'lx_expected',
        'nLx_actual', 'nLx_expected',
        'Tx_actual', 'Tx_expected',
        'e0_cntrb_d', 'e0_cntrb_i', 'e0_cntrb_t')
    )
  )

arriaga_cntfc$simulation[,,,,,'death_total_actual'] <-
  lifetables$simulation[,,,,,'death_total','actual']
arriaga_cntfc$simulation[,,,,,'death_total_expected'] <-
  lifetables$simulation[,,,,,'death_total','projected']

arriaga_cntfc$simulation[,,,,,'nmx_actual'] <-
  lifetables$simulation[,,,,,'nmx','actual']
arriaga_cntfc$simulation[,,,,,'nmx_expected'] <-
  lifetables$simulation[,,,,,'nmx','projected']

arriaga_cntfc$simulation[,,,,,'lx_actual'] <-
  lifetables$simulation[,,,,,'lx','actual']
arriaga_cntfc$simulation[,,,,,'lx_expected'] <-
  lifetables$simulation[,,,,,'lx','projected']

arriaga_cntfc$simulation[,,,,,'nLx_actual'] <-
  lifetables$simulation[,,,,,'nLx','actual']
arriaga_cntfc$simulation[,,,,,'nLx_expected'] <-
  lifetables$simulation[,,,,,'nLx','projected']

arriaga_cntfc$simulation[,,,,,'Tx_actual'] <-
  lifetables$simulation[,,,,,'Tx','actual']
arriaga_cntfc$simulation[,,,,,'Tx_expected'] <-
  lifetables$simulation[,,,,,'Tx','projected']

arriaga_cntfc$simulation[,,,,,'ex_actual'] <-
  lifetables$simulation[,,,,,'ex','actual']
arriaga_cntfc$simulation[,,,,,'ex_expected'] <-
  lifetables$simulation[,,,,,'ex','projected']

arriaga_cntfc$simulation[,,,,,'ex_actual_minus_expected'] <-
  arriaga_cntfc$simulation[,,,,,'ex_actual'] -
  arriaga_cntfc$simulation[,,,,,'ex_expected']

arriaga_cntfc$simulation[,,,,,'e0_cntrb_d'] <-
  (
    arriaga_cntfc$simulation[,,,,,'nLx_actual'] /
      arriaga_cntfc$simulation[,,,,,'lx_actual'] -
      arriaga_cntfc$simulation[,,,,,'nLx_expected'] /
      arriaga_cntfc$simulation[,,,,,'lx_expected']
  ) * arriaga_cntfc$simulation[,,,,,'lx_actual']

arriaga_cntfc$simulation[,,,,,'e0_cntrb_i'] <-
  (
    arriaga_cntfc$simulation[,,,,,'lx_expected'] /
      arriaga_cntfc$simulation[,,,,,'lx_actual'] -
      apply(arriaga_cntfc$simulation[,,,,,'lx_expected'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0)) /
      apply(arriaga_cntfc$simulation[,,,,,'lx_actual'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0))
  ) * apply(arriaga_cntfc$simulation[,,,,,'Tx_actual'],
            # apply function to vector of data by age
            2:5, function (x) c(x[-1], 0))
arriaga_cntfc$simulation[lifetables$cnst$nage,,,,,'e0_cntrb_i'] <- 0

arriaga_cntfc$simulation[,,,,,'e0_cntrb_t'] <-
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_d'] +
  arriaga_cntfc$simulation[,,,,,'e0_cntrb_i']

# test if age decomposition sums to observed e0 deficit
arriaga_cntfc$test <- log(
  apply(
    arriaga_cntfc$simulation[,,,,2,'e0_cntrb_t'], 2:4, sum
  ) /
    arriaga_cntfc$simulation[1,,,,2,'ex_actual_minus_expected']
)

# maximum decomposition error less than 10%?
max(abs(arriaga_cntfc$test), na.rm = TRUE) < 0.1
# distribution of decomposition errors
ecdf(arriaga_cntfc$test) |> plot()
hist(arriaga_cntfc$test, breaks = 500, xlim = c(-0.2, 0.2))

# Calculates CI over simulations ----------------------------------

# ci's for life tables
lifetables$ci <-
  apply(
    lifetables$simulation[
      ,,,,-1,
      c('nmx', 'npx', 'nqx', 'lx', 'ex', 'ex_diff',
        'e0_cntrb_t'),
    ],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) |>
  aperm(c(2:7,1))
V <- dimnames(lifetables$ci)
names(attr(lifetables$ci, 'dim'))[7] <- 'quantile'
dimnames(lifetables$ci) <- V

# ci's for age specific contributions to e0 deviations from expectation
arriaga_cntfc$ci <-
  apply(
    arriaga_cntfc$simulation[,,,,-1,],
    -5,
    QuantileWithMean,
    simplify = TRUE
  ) |>
  aperm(c(2:6,1))
V <- dimnames(arriaga_cntfc$ci)
names(attr(arriaga_cntfc$ci, 'dim'))[6] <- 'quantile'
dimnames(arriaga_cntfc$ci) <- V

# Transform to data frame -----------------------------------------

lifetables$ci_df <-
  as.data.frame.table(lifetables$ci, stringsAsFactors = FALSE)
names(lifetables$ci_df) <-
  c(names(attr(lifetables$ci, 'dim')), 'value')
lifetables$ci_df <-
  lifetables$ci_df |>
  as_tibble() |>
  pivot_wider(id_cols = c(age, year, sex, region_iso, scenario),
              names_from = c(var_id, quantile),
              values_from = value) |>
  mutate(across(c(age, year), ~as.integer(.x)))

arriaga_cntfc$ci_df <-
  as.data.frame.table(arriaga_cntfc$ci, stringsAsFactors = FALSE)
names(arriaga_cntfc$ci_df) <-
  c(names(attr(arriaga_cntfc$ci, 'dim')), 'value')
arriaga_cntfc$ci_df <-
  arriaga_cntfc$ci_df |>
  filter(year %in% cnst$forecast_period) |>
  as_tibble() |>
  pivot_wider(id_cols = c(region_iso, sex, year, age),
              names_from = c(var_id, quantile),
              values_from = value)

# Test ------------------------------------------------------------

(
  lifetables$simulation['0','2021','Female',,1,'population_py','actual'] ==
    lt_input$openage_100 |>
    filter(age_start == 0, year == 2021, sex == 'Female') |>
    pull(population_py)
) |> all()

cbind(
  lifetables$simulation['0','2021','Female',,2,'nmx','projected'],
  lt_input$leecarter$harmonized_forecast_sim |>
    filter(age == 0, year == 2021, sex == 'Female', nsim == 1) |>
    pull(mx)
)

# Export ----------------------------------------------------------

saveRDS(lifetables$ci_df, paths$output$lifetables)
lifetables$ci_df |>
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) |>
  write_csv(paths$output$lifetables_csv)
saveRDS(lifetables$simulation, paths$output$lifetables_sim)

saveRDS(arriaga_cntfc$ci_df, paths$output$arriaga_cntfc)
arriaga_cntfc$ci_df |>
  mutate(across(.cols = where(is.numeric), .fns = ~round(.x,6))) |>
  write_csv(paths$output$arriaga_cntfc_csv)
saveRDS(arriaga_cntfc$simulation, paths$output$arriaga_cntfc_sim)
