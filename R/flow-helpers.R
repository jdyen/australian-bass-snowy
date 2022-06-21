# functions to load and format flow scenarios

# internal function to add a season ID to data.frames
add_season <- function(x, cool_season = 4:10) {
  x %>% mutate(
    season = ifelse(month(date_formatted) %in% cool_season, "cool", "warm")
  )
}

# flow threshold = 15000 (10000, 12500 are options too)
calculate_discharge_effects <- function(
  x, 
  date,
  threshold = 15000,
  nday = 32,
  params = list(), 
  spawning = list()
) {
  
  # set default parameters
  param_set <- list(
    shape1 = 1.25,
    shape2 = 2,
    alpha = 0.95, 
    beta = 0,
    tau = 4000,
    shift = 40000
  )
  param_set[names(params)] <- params
  
  # and set default spawning period
  spawning_set <- list(
    month = 6:12
  )

  # subset to relevant time periods (defaults to Jun-Dec)
  idx <- month(date) %in% spawning_set$month
  x <- x[idx]
  date <- date[idx]
  
  # pull out a single year to help track months (leap years don't matter 
  #   because this drops Feb)
  date_once <- date[year(date) %in% min(year(date))]
  
  # identify consecutive days above threshold and calculate probability of spawning
  #   and survival over first 15 days
  metric <- tapply(
    x,
    year(date), 
    calculate_discharge_effect_annual,
    threshold = threshold, 
    nday = nday, 
    param = param_set,
    date = date_once
  )
  
  # drop final year and return
  metric[-length(metric)]
  
}

# helper function to calculate spawning and survival probability for a single year
calculate_discharge_effect_annual <- function(x, threshold, nday, param, survival, date) {
  
  # set up subsets of data because only Aug-Dec are relevant to migration
  #   and June-Nov to spawning
  xspawn <- x[month(date) %in% c(6:11)]
  xspawn_sub <- x[month(date) %in% c(8:11)]
  xsurv <- x[month(date) %in% c(8:12)]
  
  # calculate probability of spawning from Aug 1 to end of spawning period (Dec 31)
  #   accounting for late starters (minus cumulative dbeta for absent days if arriving late)
  daily_spawning_prob <- dbeta(seq(0, 1, length = length(xsurv)), param$shape1, param$shape2)
  daily_spawning_prob <- daily_spawning_prob / sum(daily_spawning_prob)

  # proportion of eggs estimated to be flushed out to sea under high discharge
  egg_mortality <- param$alpha / (1 + exp(-(xsurv - param$shift) / param$tau)) + param$beta

  # a day is productive for larvae for a set number of days based on discharge
  #   - calculate whether each day is productive so we can multiply it by 
  #     daily spawning probs
  discharge_level <- c(10000, 15000, 18000, 27000, 32000, 36000, 40000, 44000, 49000, 61000, 80000)
  days_productive <- c(20, 30:39)
  productive <- lapply(discharge_level, function(x, y) y >= x, y = xsurv)
  productive <- mapply(expand_productivity, productive, days_productive, SIMPLIFY = FALSE)
  productive <- apply(
    do.call(rbind, productive),
    2,
    any
  )
  
  # now we need to check each day for spawning probability, egg mortality, and productivity
  spawning <- rep(0, length(xspawn_sub))
  for (i in seq_along(xspawn_sub)) {
    spawning[i] <- daily_spawning_prob[i] * prod(1 - egg_mortality[i:(i + 2)]) * all(productive[(i + 2):(i + nday - 1)])
  }
  
  # check arrival dates (one fifth of pop moving per day above threshold) 
  #   and work out proportion moving
  move_days <- xspawn > threshold
  nmove <- min(5, sum(move_days))

  # calculate total spawning probability for each moving cohort (20% of adult pop)
  pr_spawn <- rep(0.2 * sum(spawning), nmove)
  
  # account for late arrivals (after Aug 1) if required
  move_dates <- date[move_days][seq_len(nmove)]
  if (any(month(move_dates) >= 8)) {
    
    late_dates <- move_dates[month(move_dates) >= 8]
    
    for (i in seq_along(late_dates)) {
      days_late <- time_length(
        interval(
          dmy(paste0("01-08-", year(move_dates)[1])),
          late_dates[i]
        ),
        unit = "day"
      )
      pr_spawn[i] <- 0.2 * sum(spawning[(days_late + 1):length(spawning)])
      
    }
    
  }
  
  # calculate and return annual values
  sum(pr_spawn)
  
}

# function to add days following initial productivity hit
#   where productivity remains sufficient for larvae
expand_productivity <- function(thresh, days) {

  # how many days are we considering?  
  n <- length(thresh)
  
  # which days hit the threshold?
  hits <- which(thresh)
  
  # add the days following the initial hit
  prod <- lapply(hits, function(x, day) x:(x + day), day = days)
  
  # collapse this down to a single vector for each hit, removing
  #   duplicate/overlapping hits
  idx <- unique(unlist(prod))
  
  # remove any hits that exceed the length of thresh
  idx <- idx[idx <= n]
  
  # update the threshold check
  thresh[idx] <- TRUE

  # return  
  thresh
  
}

# function to add climate change impacts to discharge sequences
add_climate_change_scenarios <- function(x, catchment, scenario, reference, type = "discharge", variable = "value") {
  
  # check catchment and scenario
  if (!catchment %in% c(
    "upper_murray", 
    "goulburn",
    "yarra", 
    "lower_murray", 
    "campaspe",
    "snowy",
    "werribee",
    "ovens",
    "wimmera"
  )) {
    stop("catchment not implemented; using default Victorian values", call. = FALSE)
  }
  
  # catchments (values for 2040 and 2065 relative to 1995)
  if (type == "discharge") {
    
    effect <- list(
      upper_murray = list(
        rcp85 = list(low = c(17.2, 13.5), medium = c(-8.4, -16.6), high = c(-23.3, -39.4)),
        rcp45 = list(low = c(14.5, 16.3), medium = c(-5.2, -5.6), high = c(-26.4, -37.4))
      ),
      lower_murray = list(
        rcp85 = list(low = c(32.8, 27.1), medium = c(-4.6, -11.4), high = c(-37.5, -47.0)),
        rcp45 = list(low = c(21.6, 29.3), medium = c(1.2, -5.4), high = c(-27.6, -39.3))
      ),
      goulburn = list(
        rcp85 = list(low = c(9.9, 1.3), medium = c(-9.5, -13.7), high = c(-29.1, -41.9)),
        rcp45 = list(low = c(12.2, 8.1), medium = c(-3.8, -11.7), high = c(-28.7, -33.1))
      ),
      campaspe = list(
        rcp85 = list(low = c(10.5, 1.0), medium = c(-12.3, -20.7), high = c(-37.1, -57.0)),
        rcp45 = list(low = c(20.8, 9.1), medium = c(-6.4, -12.3), high = c(-43.9, -43.6))
      ),
      yarra = list(
        rcp85 = list(low = c(10, 0.8), medium = c(-11, -16.4), high = c(-29.2, -44.3)),
        rcp45 = list(low = c(10.6, 8.7), medium = c(-3.1, -11.4), high = c(-30.0, -34.0))
      ),
      snowy = list(
        rcp85 = list(low = c(22.5, 21.0), medium = c(-7.1, -17.9), high = c(-25.3, -36.1)),
        rcp45 = list(low = c(13.8, 17.3), medium = c(-3.3, -9.3), high = c(-29.9, -31.8))
      ),
      werribee = list(
        rcp85 = list(low = c(11.8, 7.5), medium = c(-7.7, -18.1), high = c(-28.9, -45.5)),
        rcp45 = list(low = c(16.5, 10.7), medium = c(-3.8, -7.6), high = c(-35.0, -36.5))
      ),
      ovens = list(
        rcp85 = list(low = c(11.7, 1.2), medium = c(-10.8, -6.0), high = c(-23.3, -43.9)),
        rcp45 = list(low = c(12.8, 9.6), medium = c(-6.0, -15.9), high = c(-31.0, -34.4))
      ),
      wimmera = list(
        rcp85 = list(low = c(12.1, 12.3), medium = c(-6.5, -14.4), high = c(-32.3, -53.1)),
        rcp45 = list(low = c(21.0, 11.5), medium = c(-4.4, -12.0), high = c(-33.8, -38.6))
      ),
      default = list(
        rcp85 = list(low = c(8.7, 1.5), medium = c(-8.5, -15.9), high = c(-24.7, -43.8)),
        rcp45 = list(low = c(14.0, 9.4), medium = c(-1.6, -11.1), high = c(-29.1, -33.3))
      )
    )
    
  } else {
    
    if (type == "water_temperature") {
      
      # use air temperature multiplied by 0.65 based on linear model of observed
      #   air temperature against water temperature in the lower Murray, Goulburn,
      #   and Campaspe systems
      
      effect <- list(
        upper_murray = list(
          rcp85 = list(low = c(1.1, 1.9), medium = c(1.4, 2.6), high = c(1.7, 3.0)),
          rcp45 = list(low = c(0.8, 1.2), medium = c(1.1, 1.6), high = c(1.5, 2.1))
        ),
        lower_murray = list(
          rcp85 = list(low = c(1.1, 2.0), medium = c(1.5, 2.5), high = c(1.7, 3.1)),
          rcp45 = list(low = c(0.7, 1.2), medium = c(1.1, 1.6), high = c(1.4, 2.0))
        ),
        goulburn = list(
          rcp85 = list(low = c(1.0, 2.0), medium = c(1.4, 2.4), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.4, 1.9))
        ),
        campaspe = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.4), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.4), high = c(1.4, 1.8))
        ),
        yarra = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.8))
        ),
        snowy = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.4, 2.5), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.8, 1.2), medium = c(1.1, 1.5), high = c(1.4, 2.0))
        ),
        
        werribee = list(
          rcp85 = list(low = c(1.0, 1.8), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.4), high = c(1.2, 1.8))
        ),
        ovens = list(
          rcp85 = list(low = c(1.0, 2.0), medium = c(1.4, 2.5), high = c(1.6, 3.0)),
          rcp45 = list(low = c(0.7, 1.2), medium = c(1.1, 1.6), high = c(1.5, 2.0))
        ),
        wimmera = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.9))
        ),
        default = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.8))
        )
      )
      
      # apply 0.65 multiplier to all
      effect <- lapply(effect, function(x) lapply(x, function(y) lapply(y, function(z) z * 0.65)))
      
    } else {
      stop("only discharge and water temperature are currently implemented", call. = FALSE)
    }
  }
  
  # pull out relevant effect
  scaling_factor <- effect[[catchment]][[scenario]]
  
  # calculate scaled discharge
  if (reference > 2075)
    stop("scaling will not extend beyond 2075", call. = FALSE)
  lowz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$low[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$low[1] +
        (scaling_factor$low[2] - scaling_factor$low[1]) * (reference - 2040) / (2065 - 2040)
    )
  )
  mediumz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$medium[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$medium[1] +
        (scaling_factor$medium[2] - scaling_factor$medium[1]) * (reference - 2040) / (2065 - 2040)
    )
  )
  highz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$high[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$high[1] +
        (scaling_factor$high[2] - scaling_factor$high[1]) * (reference - 2040) / (2065 - 2040)
    )
  )

  # rescale appropriately for different data types
  if (type == "discharge") {
    x[[paste0(variable, "_", scenario, "low")]] <-
      x[[variable]] * (1 + lowz / 100)
    x[[paste0(variable, "_", scenario, "med")]] <-
      x[[variable]] * (1 + mediumz / 100)
    x[[paste0(variable, "_", scenario, "high")]] <-
      x[[variable]] * (1 + highz / 100)
  } else {
    x[[paste0(variable, "_", scenario, "low")]] <-
      x[[variable]] + lowz
    x[[paste0(variable, "_", scenario, "med")]] <-
      x[[variable]] + mediumz
    x[[paste0(variable, "_", scenario, "high")]] <-
      x[[variable]] + highz
  }

  # and return
  x
  
}

# internal function to rescale discharge
rescale_discharge <- function(data, variable = "value") {
  
  # add a season ID to data.frames
  data <- data %>% mutate(
    season = ifelse(month(date_formatted) %in% 4:10, "cool", "warm")
  )
  
  # decile rescaling for years prior to 1975
  data[[paste0(variable, "_1975")]] <- rescale_segment(
    x = data[[variable]],
    idx = data$date_formatted > ymd("1975-06-30"),
    season = data$season
  )
  
  # repeat decile rescaling for years prior to 1997
  data[[paste0(variable, "_1997")]] <- rescale_segment(
    x = data[[variable]],
    idx = data$date_formatted > ymd("1997-06-30"),
    season = data$season
  )
  
  # remove added columns
  data <- data %>% select(-season)

  # return
  data
  
}

# function to rescale flows based on seasonal quantiles
seasonal_quantile_rescale <- function(
  x,
  reference,
  x_season = NULL,
  reference_season = NULL,
  ...
) {
  
  # only one season if NULL
  if (is.null(x_season)) {
    x_season <- rep(1, length(x))
    reference_season <- rep(1, length(reference))
  }
  
  # need the reference season
  if (is.null(reference_season))
    stop("reference_season must be provided", call. = FALSE)
  
  # use quantile_rescale in each season
  seasons <- unique(x_season)
  for (i in seq_along(seasons)) {
    idx <- x_season == seasons[i]
    idy <- reference_season == seasons[i]
    x[idx] <- quantile_rescale(
      x = x[idx], reference = reference[idy], ...
    )
  }
  
  # and return
  x
  
}

# function to rescale by quantile
quantile_rescale <- function(
  x, 
  reference, 
  probs = seq(0, 1, by = 0.1)
) {
  
  # work out quantiles in reference (e.g. future climate)
  reference_quantiles <- get_quantile(
    x = reference, probs = probs
  )
  
  # work out quantiles in observed
  x_quantiles <- get_quantile(
    x = x, probs = probs
  )
  
  # simplified change ratio: proportional change in each quantile
  change_ratio <-
    reference_quantiles$quantiles / x_quantiles$quantiles

  # apply this change to observed values  
  x <- x * change_ratio[x_quantiles$bins]
  
  # return, ensure positive
  ifelse(x < 0, 0, x)
  
}

# internal function to calculate quantiles and bins with 
#   closed intervals
get_quantile <- function(x, probs) {
  
  # calculate quantiles of x  
  breaks <- quantile(x, probs = probs)
  
  # reduce first quantile to close interval at both ends
  breaks[1] <- breaks[1] - 1e-3
  
  # bin x into its quantiles
  bins <- cut(
    x, breaks = breaks, labels = FALSE
  )
  
  # and return mean of each bin and bins
  list(
    quantiles = tapply(
      x, bins, mean
    ),
    bins = bins
  )
  
}

# function to rescale part of a sequence based on the remaining
#   years in that sequence
rescale_segment <- function(x, idx, season = NULL, ...) {

  # all one season if not specified
  if (is.null(season))
    season <- rep(1, length(x))
  
  # calculate and return
  x_segment <- seasonal_quantile_rescale(
    x = x[!idx],
    reference = x[idx],
    x_season = season[!idx],
    reference_season = season[idx] 
  )

  # update rescaled segment
  x[!idx] <- x_segment
  
  # and return
  x
  
}

# check years of data by site
check_available <- function(x) {
  years_available <- tapply(
    x$value, 
    list(x$site_code, year(x$date_formatted), x$variable_name), 
    function(x) ifelse(sum(!is.na(x)) > 300, 1, 0)
  )
  years_available[is.na(years_available)] <- 0
  years_available
}

# fill gaps in discharge by resampling
resample_discharge <- function(x, available, variable = "value") {
  
  year_set <- as.numeric(colnames(available))
  for (i in seq_len(nrow(available))) {
    
    site_id <- rownames(available)[i]
    
    target <- year_set[available[site_id, , 1] == 0]
    source <- year_set[available[site_id, , 1] == 1]
    
    if (length(target) > 0) {
      
      x[[variable]] <- 
        resample_discharge_internal(
          x[[variable]],
          x$date_formatted,
          target = target,
          source = source
        )
      
    }
    
  }
  
  # return
  x
  
}

# function to resample years of discharge, accounting for leap years
resample_discharge_internal <- function(value, date, target, source) {
  
  # create a data.frame to resample
  data <- data.frame(
    date_formatted = date,
    value = value
  )

  # sample years from source to replace target years,
  #   accounting for leap years
  idx <- leap_year(source)
  idy <- leap_year(target)
  filled_years <- rep(NA, length(target))
  filled_years[idy] <- sample(source[idx], size = sum(idy), replace = sum(idy) > sum(idx))
  filled_years[!idy] <- sample(source[!idx], size = sum(!idy), replace = sum(!idy) > sum(!idx))
  
  # resample observed discharge years and collapse into a single data.frame
  resampled <- do.call(
    rbind,
    lapply(
      filled_years, 
      function(x, data) data %>% filter(year(date_formatted) == x),
      data = data
    )
  )
  
  # fix the dates on these years to match target years
  resampled <- resampled %>% 
    mutate(
      date_formatted = ymd(
        paste(
          rep(target, times = sapply(
            filled_years, 
            function(x) ifelse(leap_year(x), 366, 365))
          ),
          month(date_formatted),
          day(date_formatted),
          sep = "-"
        )
      )
    )
  
  # combine with target years and arrange chronologically
  resampled <- resampled %>% 
    bind_rows(data %>% filter(year(date_formatted) %in% source)) %>%
    arrange(date_formatted)
  
  # and return value only
  resampled[, 2]
  
}

# fill remaining gaps with a fill_rolling_na
fill_missing <- function(x) {
  
  idx <- apply(x, 2, function(.x) sum(is.na(.x))) > 0
  to_fill <- names(idx)[idx]
  for (i in seq_along(to_fill)) {
    x[[to_fill[i]]] <- fill_na_rolling(
      x, variable = to_fill[i], recursive = TRUE, max_iter = 20
    )
  }
  
}

# function to fill missing data iteratively based on rolling means of preceding days
fill_na_rolling <- function(
  flow, variable = "value", recursive = FALSE, max_iter = 20
) {
  
  x <- flow[[variable]]
  
  if (any(is.na(x)) & !all(is.na(x))) {
    
    if (recursive) {
      
      counter <- 1
      while(any(is.na(x)) & counter < (max_iter + 1)) {
        
        n <- length(x)
        idx <- sapply(rev(seq_len(5)) - 1, function(x) rep(x, n))
        idx <- sweep(idx, 1, seq_len(n), "+")
        idx <- ifelse(idx > n, NA, idx)
        df <- matrix(x[idx], nrow = n)
        
        mean_val <- apply(df, 1, mean, na.rm = TRUE)    
        
        idx <- which(is.na(x))
        
        flow[[variable]][idx] <- mean_val[idx]
        
        x <- flow[[variable]]
        
        counter <- counter + 1
        
      }
      
    } else {
      
      n <- length(x)
      idx <- sapply(rev(seq_len(5)) - 1, function(x) rep(x, n))
      idx <- sweep(idx, 1, seq_len(n), "+")
      idx <- ifelse(idx > n, NA, idx)
      df <- matrix(x[idx], nrow = n)
      
      mean_val <- apply(df, 1, mean, na.rm = TRUE)    
      
      idx <- which(is.na(x))
      
      flow[[variable]][idx] <- mean_val[idx]
      
    }
    
  }
  
  flow
  
}

filter_qc <- function(x, variable = "value", quality = "quality_code", threshold = 150) {
  
  idx <- x[[quality]] > threshold
  x[[variable]][idx] <- NA
  
  x
  
}

filter_varcode <- function(x, priority = "141.00") {
  
  varcodes <- unique(x$variable_code)
  
  if (length(varcodes) > 1) {

    if (priority %in% varcodes) {
      x <- x[x$variable_code == priority, ]
    } else {
      x <- x[x$variable_code == varcodes[1], ]
      warning(
        "priority variable code not in data, using first unique",
        " variable code, which is ",
        varcodes[1],
        ". Alternative variable codes for this data set: ", 
        varcodes[-1],
        call. = FALSE
      )
    }

  }
  
  x
  
}

# new impute temperature function (simplified)
fit_temp_model <- function(data, formula) {
  mod <- lm(formula, data = data)
  predict(mod, newdata = data)
}

# function to add environmental flows to discharge sequence
#' @param x discharge sequence
#' @param date date corresponding to discharge observations
#' @param size size of release in ML/day
#' @param n number of releases in each year
#' @param frequency how many years to wait before another release
add_envflows <- function(x, date, size = 15000, n = 3, frequency = 1, stochastic = TRUE) {
  
  # which days of the year need releases?
  release_days <- rep(FALSE, length(date))
  for (i in seq_len(n)) {
    release_days <- release_days | (month(date) == (7 + i) & day(date) %in% c(1:3))
  }
  
  # what if frequency is less than annual?
  release_days <- ifelse(
    (year(date) %% frequency) == 0,
    release_days,
    FALSE
  )
  
  # add desired env flow to observed discharge
  if (stochastic)
    size <- rnorm(sum(release_days), mean = size, sd = 1000)
  x[release_days] <- x[release_days] + size
  
  # return
  x
  
}

discharge_plot <- function(x, col = NULL, ...) {
  
  x <- x %>% 
    select(date_formatted, value_1975, value_1997, value_rcp85low, value_rcp85med, value_rcp85high) %>%
    pivot_longer(
      cols = contains("value"),
      names_prefix = "value_",
      names_to = "scenario"
    ) %>% 
    arrange(scenario)
  
  if (is.null(col))
    col = RColorBrewer::brewer.pal(length(unique(x$scenario)), "Set2")
  
  p <- x %>% ggplot(aes(x = date_formatted, y = value, colour = scenario)) +
    geom_line(size = 1.2, position = position_jitter(w = 0.2)) +
    scale_y_log10() +
    scale_linetype_manual(col) +
    theme_bw()
  
  print(p)
  
}

