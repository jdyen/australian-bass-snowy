# function to expand a vector of actions into a full grid of
#   all possible combinations of actions (including no action)
expand_combn <- function(x) {
  out <- lapply(
    seq_along(x), 
    function(i, .x) combn(x = .x, m = i), 
    .x = x
  )
  max_row <- max(sapply(out, nrow))
  out <- lapply(
    out,
    function(x) rbind(
      x, 
      matrix(
        "none",
        nrow = max_row - nrow(x),
        ncol = ncol(x)
      )
    )
  )
  out <- do.call(cbind, out)
  cbind(rep("none", max_row), out)
}

# function to expand list of actions by site to include multiple climates
expand_actions <- function(x, climates) {
  
  # pull out site names
  site_names <- names(x)
  
  # how many actions are there in each site?
  n_actions <- sapply(x, ncol)
  
  # and what's the max at any one site?
  max_actions <- max(sapply(x, nrow))
  
  # add a "no action" option
  x <- lapply(
    x, 
    function(.x, max_x) rbind(.x, rep("none", ncol(.x)))[seq_len(max_x), ],
    max_x = max_actions
  )
  
  # convert to matrix
  x <- t(do.call(cbind, x))
  
  # and expand over sites and climates to give a useful output
  x <- list(
    site = rep(rep(site_names, n_actions), times = length(climates)),
    climate = rep(climates, each = sum(n_actions)),
    actions = do.call(rbind, lapply(seq_along(climates), function(.x) x))
  )
  
  # how many scenarios do we have?
  x$n_required <- nrow(x$actions)
  
  # return
  x
  
}

# function to remove conflicting actions
remove_conflicts <- function(x, conflict) {
  
  # remove conflicts
  x <- x[, apply(x, 2, function(.x) sum(conflict %in% .x) < 2)]
  
  # and clean up empty rows
  x <- x[apply(x, 1, function(.x, n) sum(.x == "none") < n, n = ncol(x)), ]
  
  # return
  x
  
}
