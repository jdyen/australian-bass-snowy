# Simulation of population trajectories of estuary perch under
#    three management interventions

# set a seed to ensure reproducibility
set.seed(2021-06-21)

# load some packages to read in data sets
library(qs)

# to work with data sets
library(lubridate)
library(dplyr)
library(tidyr)
library(aae.hydro)
library(ggplot2)

# and to simulate population dynamics
library(aae.pop.templates)

# load some helpers
source("R/actions-helpers.R")
source("R/flow-helpers.R")
source("R/simulate-helpers.R")

# define sites
sites <- c("snowy")
catchment <- c("snowy" = "snowy")

# define management actions
climates <- c("1975", "1997", "rcp85low", "rcp85med", "rcp85high", "rcp45low", "rcp45med", "rcp45high")
actions <- list(
  snowy = c("env_flows_low2t",
            "env_flows_high2t", 
            "env_flows_low3t", 
            "env_flows_high3t", 
            "env_flows_high3t_in2",
            "env_flows_high3t_in3",
            "env_flows_high3t_in5",
            "stocking10",
            "stocking15",
            "stocking30",
            "stocking50")
)

# expand management actions
actions <- lapply(actions, expand_combn)

# don't want to keep multiple e-flows actions in a single run
actions <- lapply(actions, remove_conflicts, conflict = c("env_flows_low2t", "env_flows_high2t", "env_flows_low3t", "env_flows_high3t", "env_flows_high3t_in2", "env_flows_high3t_in3", "env_flows_high3t_in5"))
actions <- lapply(actions, remove_conflicts, conflict = c("stocking10", "stocking15", "stocking30", "stocking50", "stocking100"))

# flatten the actions, sites, and climates to a matrix
actions <- expand_actions(actions, climates = climates)

# load historical flows
# Snowy River @ JARRAHMOND 222200
#    Latitude: 37°39'39.9"S, Longitude: 148°21'40.7"E
discharge_gauges <- c("snowy" = "222200")
temp_gauges <- c("snowy" = "222205")
site_names <- c("Jarrahmond")
site_names_temp <- c("Newmerella")
recompile_discharge <- FALSE
if (recompile_discharge) {
  
  # download discharge from WMIS
  discharge <- lapply(
    discharge_gauges,
    fetch_hydro,
    variables = c("discharge"),
    start = "1970-01-01",
    end = "2020-12-31",
    include_missing = TRUE
  )
  
  # remove excess variable codes if provided
  discharge <- lapply(discharge, filter_varcode, priority = "141.00")
  
  # set poor quality code observations to NA
  discharge <- lapply(discharge, filter_qc, threshold = 150)

  # which years of discharge are available?
  years_available <- lapply(discharge, check_available)

  # save this
  qsave(discharge, file = "data/discharge_compiled.qs")
  
  # download water temp from WMIS
  water_temp <- lapply(
    temp_gauges,
    fetch_hydro,
    variables = c("temperature"),
    start = "1970-01-01",
    end = "2020-12-31",
    include_missing = TRUE
  )
  
  # fill gaps with linear model against air temperature (in data/)
  air_temp <- read.csv("data/IDCJAC0010_084030_1800_Data.csv")
  air_temp_later <- read.csv("data/IDCJAC0010_084145_1800_Data.csv")
  air_temp <- air_temp %>% mutate(
    date_formatted = parse_date_time(paste(Year, Month, Day, sep = "-"), orders = c("ymd"))
  )
  air_temp_later <- air_temp_later %>% mutate(
    date_formatted = parse_date_time(paste(Year, Month, Day, sep = "-"), orders = c("ymd"))
  )

  # fit linear model
  water_temp_estimated <- water_temp$snowy %>% 
    left_join(air_temp %>% select("date_formatted", "Maximum.temperature..Degree.C."), by = "date_formatted") %>%
    rename(airtemp = Maximum.temperature..Degree.C.) %>%
  left_join(air_temp_later %>% select("date_formatted", "Maximum.temperature..Degree.C."), by = "date_formatted") %>%
    mutate(
      airtemp = ifelse(is.na(airtemp), Maximum.temperature..Degree.C., airtemp),
      month = factor(month(date_formatted)),
      day_of_year = day(date_formatted)
    ) %>%
    fit_temp_model(value ~ airtemp + month + day_of_year + day_of_year^2)
  water_temp$snowy <- water_temp$snowy %>%
    mutate(value = water_temp_estimated)

  # use rolling mean to fill small number of remaining gaps
  water_temp <- lapply(
    water_temp,
    fill_na_rolling,
    variable = "value",
    recursive = TRUE,
    max_iter = 20
  )
  
  qsave(water_temp, file = "data/water_temp_compiled.qs")
  
} else {
  
  # load pre-compiled versions
  discharge <- qread("data/discharge_compiled.qs")
  water_temp <- qread("data/water_temp_compiled.qs")
  
}

# calculate scaled versions of discharge and water temp (post-1975 and post-1997)
discharge <- lapply(discharge, rescale_discharge)
water_temp <- lapply(water_temp, rescale_discharge)

# add climate change scenarios
discharge <- mapply(
  add_climate_change_scenarios,
  x = discharge, 
  catchment = catchment, 
  MoreArgs = list(
    scenario = "rcp85", reference = 2065, type = "discharge", variable = "value"
  ),
  SIMPLIFY = FALSE
)
discharge <- mapply(
  add_climate_change_scenarios,
  x = discharge, 
  catchment = catchment, 
  MoreArgs = list(
    scenario = "rcp45", reference = 2065, type = "discharge", variable = "value"
  ),
  SIMPLIFY = FALSE
)
water_temp <- mapply(
  add_climate_change_scenarios,
  x = water_temp, 
  catchment = catchment, 
  MoreArgs = list(
    scenario = "rcp85", reference = 2065, type = "water_temperature", variable = "value"
  ),
  SIMPLIFY = FALSE
)
water_temp <- mapply(
  add_climate_change_scenarios,
  x = water_temp, 
  catchment = catchment, 
  MoreArgs = list(
    scenario = "rcp45", reference = 2065, type = "water_temperature", variable = "value"
  ),
  SIMPLIFY = FALSE
)

# plot discharge
for (i in seq_along(discharge))
  discharge_plot(discharge[[i]])
par(mfrow = c(length(discharge), 1))
for (i in seq_along(water_temp))
  discharge_plot(water_temp[[i]])

# simulate population dynamics
nsim <- 1000
qsave(actions, file = "outputs/simulations/actions_list.qs")

# set some initial conditions by site
init_by_site <- c("snowy" = 5000)
k_by_site <- c("snowy" = 10000)
ab <- australian_bass(
  k = 10000
)
initials <- lapply(
  init_by_site,
  set_initial,
  nsim = nsim,
  matrix = ab$dynamics$matrix
)

# calculate discharge metrics for each climate/site
nsite <- length(unique(actions$site))
nclimate <- length(unique(actions$climate))
scenarios_to_plot <- rep(c(1, 2), times = nclimate) + rep(c(0:7) * 8, each = 2 * nsite)
covar_tmp <- vector("list", length = length(scenarios_to_plot))
for (i in seq_along(scenarios_to_plot)) {
  
  # get settings based on actions
  sim_settings <- get_settings(
    actions = actions$actions[i, ],
    k = k_by_site[actions$site[i]],
    discharge = discharge[[actions$site[i]]],
    water_temperature = water_temp[[actions$site[i]]],
    climate = actions$climate[i],
    p_capture = 0.05
  )

  # store covariates
  covar_tmp[[i]] <- sim_settings$covars
  
}

# loop through all actions
for (i in seq_len(actions$n_required)) {
  
  # get settings based on actions
  sim_settings <- get_settings(
    actions = actions$actions[i, ],
    k = k_by_site[actions$site[i]],
    discharge = discharge[[actions$site[i]]],
    water_temperature = water_temp[[actions$site[i]]],
    climate = actions$climate[i],
    p_capture = 0
  )
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  ab <- australian_bass(
    k = sim_settings$carrying_capacity,
    n = list(
      rep(sim_settings$n_stocked, ntime),
      rep(0, ntime),
      rep(0, ntime),
      rep(0, ntime)
    ),
    ntime = ntime, 
    start = rep(1, 4), 
    end = rep(ntime, 4), 
    add = c(TRUE, TRUE, TRUE, TRUE),
    p_capture = sim_settings$p_capture
  )
  
  # extract correct initials
  init_set <- initials[[actions$site[i]]]

  # replace Ricker model with BH density dependence
  ab$dynamics <- update(
    ab$dynamics,
    density_dependence = density_dependence(
      masks = reproduction(ab$dynamics$matrix),
      funs = beverton_holt(k = sim_settings$carrying_capacity, exclude = 1:3)
    )
  )
  
  # simulate population dynamics
  sims <- simulate(
    ab,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = format_covariates(sim_settings$covars),
      density_dependence = list(theta = 1)
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # save outputs
  qsave(sims, file = paste0("outputs/simulations/sims_bh_theta1_", i, ".qs"))
  
}

# we can load and summarise each
actions <- qread("outputs/simulations/actions_list.qs")
sims_list <- paste0("sims_bh_theta1_", seq_len(actions$n_required), ".qs")

# calculate population-level outcomes
sim_sum <- lapply(
  sims_list,
  summarise_sims,
  start = 11,
  threshold = seq_len(5000)
)

# save summary output
qsave(sim_sum, file = "outputs/simulations/sim_sum_bh_theta1.qs")

# and load back in
sim_sum <- qread("outputs/simulations/sim_sum_bh_theta1.qs")

# extract emps and risk curves
emps <- sapply(sim_sum, function(x) x$emps)
risk_curves <- do.call(
  rbind,
  lapply(sim_sum, function(x) x$risk[1, ])
)

# prepare plots of trajectories
naction <- 40
best_action <- 28
sims_to_plot <- list(
  c(best_action, 1),
  c(best_action, 1) + naction,
  c(best_action, 1) + 5 * naction,  # RCP4.5
  c(best_action, 1) + 6 * naction,
  c(best_action, 1) + 7 * naction,
  c(best_action, 1) + 2 * naction,  # RCP8.5
  c(best_action, 1) + 3 * naction,
  c(best_action, 1) + 4 * naction
)
ylim_traj <- 20000
png(file = paste0("outputs/figs/trajectories_bh_theta1.png"),
     height = 1.25 * 6,
     width = 6,
     res = 600,
     units = "in",
     pointsize = 12)
layout(mat = seq_len(2), heights = c(1, 0.25))
par(mar = c(4.1, 5.2, 1.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(8, "Set1")

# initialise plot
sim_plot <- list(
  best = qread(paste0("outputs/simulations/sims_bh_theta1_", sims_to_plot[[1]][1], ".qs")),
  baseline = qread(paste0("outputs/simulations/sims_bh_theta1_", sims_to_plot[[1]][2], ".qs"))
)
xplot <- seq_len(dim(sim_plot[[1]])[3]) + 1969
plot(subset(sim_plot[[1]], subset = 4:40)[1, 1, ] ~ xplot,
     xlim = c(1969, 2019),
     ylim = c(0, ylim_traj),
     type = "n",
     xlab = "", ylab = "", las = 1, bty = "l")
mtext("Abundance", side = 2, line = 4.1, cex = 0.9)
mtext("Year", side = 1, line = 2.5, cex = 0.9)

for (j in seq_along(sims_to_plot)) {
  
  sim_plot <- list(
    best = qread(paste0("outputs/simulations/sims_bh_theta1_", sims_to_plot[[j]][1], ".qs")),
    baseline = qread(paste0("outputs/simulations/sims_bh_theta1_", sims_to_plot[[j]][2], ".qs"))
  )
  traj_sub <- sample(seq_len(nsim), size = min(c(30, nsim)), replace = FALSE)
  
  for (w in seq_along(sim_plot)) {
    
    # plot trajectories
    for (k in traj_sub) {
      yplot <- apply(
        subset(sim_plot[[w]], subset = 4:40)[k, , ],
        2, 
        sum
      )
      lines(yplot ~ xplot, lty = w, col = scales::alpha(col_pal[j], 0.05)) 
    }
    
    # add mean lines
    ymean <- apply(
      apply(
        subset(sim_plot[[w]], subset = 4:40),
        c(1, 3),
        sum
      ),
      2,
      mean
    )
    lines(ymean ~ xplot, col = col_pal[j], lty = w)
    
  }
  
}

# add a legend
x <- seq(0, 1, by = 0.1)
plot(seq(0, 10, by = 1) ~ x,
     bty = "n", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", type = "n")
climate_label <- c(
  "Historical", "Post-1997",
  "RCP4.5 low", "RCP4.5 medium", "RCP4.5 high",
  "RCP8.5 low", "RCP8.5 medium", "RCP8.5 high"
)
legend(x = 0.5, y = 9, 
       ncol = 2,
       xpd = TRUE,
       cex = 0.7,
       legend = c(
           paste(climate_label, "best interventions", sep = ": "),
           paste(climate_label, "no interventions", sep = ": ")
       ),
       col = rep(col_pal, times = 2), 
       lty = rep(seq_along(sim_plot), each = 8), 
       lwd = 2, xjust = 0.5, yjust = 0.8)

# close plotting device
dev.off()

# prepare plots of risk curves
climate_label <- c(
  "Historical", "Post-1997",
  "RCP4.5 low", "RCP4.5 med.", "RCP4.5 high",
  "RCP8.5 low", "RCP8.5 med.", "RCP8.5 high"
)
action_label <- c(
  "none" = "None",
  "env_flows_low2t" = "Env. water (two small)",
  "env_flows_high2t" = "Env. water (two large)",
  "env_flows_low3t" = "Env. water (three small)",
  "env_flows_high3t" = "Env. water (three large)",
  "env_flows_high3t_in2" = "Env. water (three large every two years)",
  "env_flows_high3t_in3" = "Env. water (three large every three years)",
  "env_flows_high3t_in5" = "Env. water (three large every five years)",
  "stocking10" = "Stocking (10000)",
  "stocking15" = "Stocking (15000)",
  "stocking30" = "Stocking (30000)",
  "stocking50" = "Stocking (50000)"
)
threshold <- seq_len(5000)
col_pal <- RColorBrewer::brewer.pal(11, "RdYlGn")
col_pal <- col_pal[c(1, 5:11, 2:4)]
col_pal <- c(col_pal, col_pal[11])
lty_set <- c(rep(1, 11), 2)

# initialise plot
png(
  file = paste0("outputs/figs/risk_curve_bh_theta1.png"),
  units = "in",
  width = 9,
  height = (9 / 4) * 4.4,
  pointsize = 12,
  res = 300
)

layout(mat = matrix(c(1, 9, 2, 10, 17, 3, 11, 4, 12, 17, 5, 13, 6, 14, 17, 7, 15, 8, 16, 17), ncol = 4), heights = c(rep(1, 4), 0.2))
climates_reordered <- climates[c(1:2, 6:8, 3:5)]
par(mar = c(4.1, 4.5, 2.1, 1.1))
for (j in seq_along(climates_reordered)) {
  
  # subset to relevant risk curves
  idx <- actions$climate == climates_reordered[j]
  risk_sub <- risk_curves[idx, ]

  # subset actions  
  actions_sub <- actions$actions[idx, ]

  # filter to individual actions only
  risk_sub_none <- risk_sub[actions_sub[, 2] == "none" & grepl("none", actions_sub[, 1]), , drop = FALSE]
  risk_sub_ef <- risk_sub[actions_sub[, 2] == "none" & grepl("env_flows", actions_sub[, 1]), ]
  risk_sub_stock <- risk_sub[actions_sub[, 2] == "none" & grepl("stocking", actions_sub[, 1]), ]
  
  # grab set of actions for the legend
  actions_sub <- actions_sub[actions_sub[, 2] == "none", ]
  
  # set axis labels
  xlab_set <- ylab_set <- ""
  if (j %in% c(5:8))
    xlab_set <- "Threshold population size (TPS)"
  if (j %in% c(1, 5))
    ylab_set <- "Pr(N < TPS)"  
  
  # set up plot
  plot(
    risk_sub_ef[1, ] ~ threshold, 
    ylim = c(0, 1), 
    las = 1,
    type = "n", 
    bty = "l",
    xlim = c(0, 5000),
    xlab = xlab_set,
    ylab = ylab_set
  )
  
  # add no actions curve
  lines(risk_sub_none[1, ] ~ threshold, col = col_pal[1], lwd = 2, lty = lty_set[1])
  
  # add each curve
  for (k in seq_len(nrow(risk_sub_ef))) {
    lines(risk_sub_ef[k, ] ~ threshold, col = col_pal[k + 1], lwd = 2, lty = lty_set[k + 1])
  }
  
  # add climate label
  mtext(paste0(climate_label[j], " (env. water)"), side = 3, cex = 0.8, line = 0.1, adj = 0)
  
  # set up plot
  plot(
    risk_sub_stock[1, ] ~ threshold, 
    ylim = c(0, 1), 
    las = 1,
    type = "n", 
    bty = "l",
    xlim = c(0, 5000),
    xlab = xlab_set,
    ylab = ylab_set
  )
  
  # add no actions curve
  lines(risk_sub_none[1, ] ~ threshold, col = col_pal[1], lwd = 2, lty = lty_set[1])

  # add each curve
  for (k in seq_len(nrow(risk_sub_stock))) {
    lines(risk_sub_stock[k, ] ~ threshold, col = col_pal[8 + k], lwd = 2, lty = lty_set[8 + k])
  }
  
  # add climate label
  mtext(paste0(climate_label[j], " (stocking)"), side = 3, cex = 0.8, line = 0.1, adj = 0)
  
  # add legend
  if (j == 8) {
    
    par(mar = rep(0, 4))
    plot(seq(0, 1, length = 10) ~ seq(0, 200, length = 10),
         type = "n",
         xaxt = "n",
         yaxt = "n",
         xlab = "",
         ylab = "",
         bty = "n")
    actions_text <- apply(actions_sub, c(1, 2), function(x) action_label[x])
    actions_text <- apply(actions_text, 1, paste0, collapse = "/")
    actions_text <- gsub("/None", "", actions_text)
    actions_text <- gsub("/", " / ", actions_text)
    legend_text <- actions_text
    legend(
      x = 110,
      y = 0.5, 
      xjust = 0.5,
      yjust = 0.5,
      ncol = 4,
      legend = legend_text,
      lty = lty_set,
      lwd = 2, 
      col = col_pal, 
      bty = "n", 
      cex = 1,
      y.intersp = 0.75,
      x.intersp = 1,
      xpd = TRUE
    )
  }
  
  
}

# shut down plot
dev.off()

# extract pr_persist from risk curves, based on 25 % of K by system
pr_persist <- mapply(
  function(x, y) 1 - x$risk[1, y],
  x = sim_sum, 
  y = 1000
)

# and split pr_persist by system
extract_pop_persist <- function(site, actions, persist, emps) {
  data.frame(
    climate = actions$climate[actions$site == site],
    actions = actions$actions[actions$site == site, ],
    persist = persist[actions$site == site], 
    emps = emps[actions$site == site]
  )
}
pop_persist <- lapply(
  sites,
  extract_pop_persist,
  actions = actions,
  persist = pr_persist, 
  emps = emps
)
names(pop_persist) <- sites

# population level tables
reorder_benefit <- function(x) {
  out <- x[order(x$climate, x$persist, x$emps, decreasing = c(FALSE, TRUE, TRUE)), ]
  out$emps <- round(out$emps)
  out
}
write.csv(
  do.call(
    rbind, lapply(
      pop_persist,
      reorder_benefit
    )
  ),
  file = "outputs/tables/population_outcomes_bh_theta1.csv"
)

# add cost info
cost_per_ml <- 75
cost_per_fingerling <- 1.1
costs <- c(
  "none" = 0,
  "env_flows_low2t" = 2 * 10000 * cost_per_ml,
  "env_flows_high2t" = 2 * 15000 * cost_per_ml,
  "env_flows_low3t" = 3 * 10000 * cost_per_ml,
  "env_flows_high3t" = 3 * 15000 * cost_per_ml,
  "env_flows_high3t_in2" = (3 / 2) * 15000 * cost_per_ml,
  "env_flows_high3t_in3" = (3 / 3) * 15000 * cost_per_ml,
  "env_flows_high3t_in5" = (3 / 5) * 15000 * cost_per_ml,
  "stocking10" = 10000 * cost_per_fingerling,
  "stocking15" = 15000 * cost_per_fingerling,
  "stocking30" = 30000 * cost_per_fingerling,
  "stocking50" = 50000 * cost_per_fingerling
)
write.csv(costs, file = "outputs/tables/intervention-costs_bh_theta1.csv")

# barplot of frequency of action inclusion under each climate for each pop
pop_outcomes <- read.csv("outputs/tables/population_outcomes_bh_theta1.csv")
all_actions <- c(
  "env_flows_low2t",
  "env_flows_low3t",
  "env_flows_high2t",
  "env_flows_high3t",
  "env_flows_high3t_in2",
  "env_flows_high3t_in3",
  "env_flows_high3t_in5",
  "stocking10",
  "stocking15",
  "stocking30",
  "stocking50"
)
action_freq <- matrix(0, nrow = length(climates_reordered), ncol = 12)
colnames(action_freq) <- c("n", all_actions)
for (j in seq_along(climates_reordered)) {
  idx <- pop_outcomes$climate == climates_reordered[j]
  x_sub <- pop_outcomes[idx, ]
  idy <- x_sub$emps >= (0.8 * max(x_sub$emps)) & 
    x_sub$persist >= 0.99
  freq_table <- table(unlist(x_sub[idy, grepl("actions", colnames(x_sub))])) / sum(idy)
  freq_table <- c(sum(idy), freq_table)
  names(freq_table)[1] <- "n"
  freq_table <- freq_table[names(freq_table) != "none"]
  action_freq[j, match(names(freq_table), colnames(action_freq))] <- freq_table
}

# plot results
png(file = "outputs/figs/prop-inclusion_bh_theta1.png", width = 6, height = 6 * 1.25, units = "in", res = 600, pointsize = 12)

# set layout
laymat <- matrix(c(1, 2), nrow = 2)
layout(laymat, heights = c(1, 0.25))

# set plot margins
old_mar <- par()$mar
par(mar = c(6.8, 4.1, 2.1, 0.5))

col_pal <- RColorBrewer::brewer.pal(length(all_actions), "Set3")

# sort action_freq
act_tmp <- action_freq[, order(colnames(action_freq))]

# remove dud columns
idx <- !colnames(act_tmp) %in% c("none", "n", "null")
act_tmp <- act_tmp[, idx]
act_tmp <- act_tmp[, match(all_actions, colnames(act_tmp))]

# match actions and colours
idy <- match(colnames(act_tmp), all_actions)  

# plot it
rownames(act_tmp) <- climate_label
barplot(
  t(as.matrix(act_tmp)),
  beside = TRUE,
  space = c(0, 2),
  col = col_pal[idy],
  cex.names = 0.95,
  las = 2,
  ylim = c(0, 1),
  xlab = "",
  ylab = "Proportional inclusion"
)

par(mar = rep(0, 4))
plot(c(0, 1) ~ 1, bty = "n", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend(x = "center", fill = col_pal, legend = action_label[all_actions], cex = 0.8, bty = "n", horiz = FALSE, ncol = 2)

dev.off()
