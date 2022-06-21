# functions to support simulations of population dynamics

# function to set initial conditions based on assumed constant proportional
#   survival over all age classes
set_initial <- function(adult = 500, nsim = 1, age = 1:40, adults = 3:max(age), rate = 0.3, matrix = NULL) {
  
  # work out age frequency from steady state if matrix provided,
  #   from constant survival assumption otherwise
  if (is.null(matrix)) {
    initial_age_frequency <- dexp(age, rate = rate)
  } else {
    initial_age_frequency <- Re(eigen(matrix)$vectors[, 1])
  }
  initial_age_frequency <- initial_age_frequency / sum(initial_age_frequency)
  
  # set total abundance from yoy abundance
  total_abundance <- adult / sum(initial_age_frequency[adults])
  
  # return
  matrix(
    rpois(
      nsim * length(initial_age_frequency),
      lambda = total_abundance * initial_age_frequency
    ),
    nrow = nsim,
    byrow = TRUE
  )
  
}

# function to set simulate settings based on current actions
# effects of actions:
#   gene mixing: based on Lutz et al. 2020, 80 % increase in early life survival
#   stocking: add x fingerlings in a given year
get_settings <- function(
  actions,
  k,
  discharge,
  water_temperature,
  climate,
  p_capture = 0.05
) {
  
  # set up parameters based on actions
  p_capture <- ifelse(any(grepl("fishing_regulations", actions)), 0, p_capture)
  n_stocked <- 0
  if(any(grepl("stocking", actions))) {
    idx <- which(grepl("stocking", actions))
    n_stocked <- as.numeric(paste0(gsub("stocking", "", actions[idx]), "000"))
  }
  
  # pull out discharge for the climate and add environmental flows if required
  discharge_set <- discharge[[paste0("value_", climate)]]
  if (any(grepl("env_flows", actions))) {
    freq_set <- 1
    if (any(grepl("_in", actions))) {
      idx <- which(grepl("_in", actions))
      freq_set <- as.numeric(gsub("env_flows_high3t_in", "", actions[idx]))
    }
    discharge_set <- add_envflows(
      x = discharge_set, 
      date = discharge$date_formatted,
      size = ifelse(any(grepl("_high", actions)), 15000, 10000),
      n = ifelse(any(grepl("3t", actions)), 3, 2),
      frequency = freq_set
    )
  }
  
  # and calculate flow metrics from discharge time series
  covars <- data.frame(
    discharge = calculate_discharge_effects(
      x = discharge_set,
      date = discharge$date_formatted
    ),
    water_temperature = calculate(
      value = water_temperature[[paste0("value_", climate)]],
      date = water_temperature$date_formatted,
      resolution = annual(season = 8:11, subset = 1970:2019)
    )$metric
  )
  
  # return named list of settings
  list(
    p_capture = p_capture,
    covars = covars,
    carrying_capacity = k,
    n_stocked = 0.1226 * 0.5 * n_stocked # account for fingerling survival and sex ratio
  )
  
}

# calculate risk of falling below any threshold prior
#   to some time point
calculate_risk <- function(sims, year = NULL, start = 2, threshold = NULL, adults = 3:30, n= 100) {
  
  # how many years are there?
  nyear <- dim(sims)[3]
  
  # default to final year if year not specified
  if (is.null(year))
    year <- nyear - 1
  
  # how many adults are there?
  adult_abundance <- apply(sims[, adults, ], c(1, 3), sum)
  
  # what is our sequence of thresholds?
  if (is.null(threshold))
    threshold <- c(0, exp(seq(0, log(max(adult_abundance)), length = n)))
  
  # calculate probability of being below threshold prior to or in year
  out <- matrix(NA, nrow = length(year), ncol = length(threshold))
  for (i in seq_along(year)) {
    out[i, ] <- sapply(
      threshold,
      function(x, y, z) mean(
        apply(y[, start:(z + 1)], 1, function(.x, .y) any(.x < .y), .y = x)
      ),
      y = adult_abundance,
      z = year[i]
    )
  }
  
  # return
  out
  
}

# function to calculate expected minimum population from year
#   start to the final year
calculate_emps <- function(sims, start = 1, adults = 3:30) {
  
  # which years do we care about?
  years <- lapply(start, function(x, end) x:end, end = dim(sims)[3])
  
  # pull out adult abundances
  adult_abundance <- apply(sims[, adults, ], c(1, 3), sum)
  
  # subset to those years and get min pop 
  sapply(years, get_emps, abund = adult_abundance)
  
}

# internal function to get EMPS from a single time period
get_emps <- function(years, abund) {
  
  # what is the minimum pop size per trajectory?
  if (length(years) > 1) {
    
    # look over multiple years
    emps <- apply(
      abund[, years],
      1,
      min
    )
    
  } else {
    
    # otherwise only one year, so min = observed
    emps <- abund[, years]
    
  }
  
  # return
  mean(emps)
  
}

# function to do all of the summaries above from a file name and path
summarise_sims <- function(
  file,
  path = "outputs/simulations/", 
  threshold = seq_len(1000),
  start = 5
) {
  
  # load single simulation  
  tmp <- qread(paste0(path, file))
  
  # calculate min pop size and risk curve
  emps <- calculate_emps(tmp, adults = 3:30, start = start)
  risk <- calculate_risk(tmp, adults = 3:30, start = start, threshold = threshold)
  
  # return
  list(
    emps = emps,
    risk = risk
  )
  
}
