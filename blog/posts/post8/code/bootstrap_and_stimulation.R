
suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
})

# -----------------------------------------------------------------------------
# 0. SETTINGS: change these values and rerun the script to explore the examples
# -----------------------------------------------------------------------------
SEED <- 6400

# Wage population and bootstrap example
POP_MEAN <- 100
POP_SD <- 20
BOOT_SAMPLE_SIZE <- 25
B_BOOT <- 1000
RESAMPLE_COUNTS <- c(250, 1000, 5000)
PARAMETER_MEANS <- c(80, 100, 120)
PARAMETER_SDS <- c(10, 20, 30)
BOOT_CI_LEVELS <- c(0.90, 0.95, 0.99)

# Repeated-sampling coverage experiment
COVER_SAMPLE_SIZES <- c(10, 25, 50, 100)
B_COVER <- 10000
COVER_CI_LEVELS <- c(0.90, 0.95, 0.99)

# Two-sided coin-toss test
N_TOSSES <- 100
P_NULL <- 0.50
OBSERVED_HEADS <- 62
B_NULL <- 10000
ALPHA_LEVELS <- c(0.10, 0.05, 0.01)

# Check settings that would otherwise make the simulations invalid.
stopifnot(
  POP_SD > 0,
  BOOT_SAMPLE_SIZE >= 2,
  B_BOOT >= 100,
  B_COVER >= 100,
  B_NULL >= 100,
  N_TOSSES >= 1,
  P_NULL > 0, P_NULL < 1,
  OBSERVED_HEADS >= 0, OBSERVED_HEADS <= N_TOSSES,
  all(BOOT_CI_LEVELS > 0 & BOOT_CI_LEVELS < 1),
  all(COVER_CI_LEVELS > 0 & COVER_CI_LEVELS < 1),
  all(ALPHA_LEVELS > 0 & ALPHA_LEVELS < 1)
)

# Resolve outputs relative to this script so running it from another working
# directory still writes into the Bootstrap project.
script_path <- NULL
for (frame in sys.frames()) {
  candidate <- tryCatch(frame$ofile, error = function(e) NULL)
  if (!is.null(candidate) && length(candidate) && nzchar(candidate[1])) {
    script_path <- candidate[1]
    break
  }
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (is.null(script_path) && length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1])
}
if (!is.null(script_path)) {
  PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), ".."),
                                winslash = "/", mustWork = FALSE)
} else {
  PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
project_path <- function(...) file.path(PROJECT_ROOT, ...)
# Create output folders if needed.
dir.create(project_path("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(project_path("results"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. BOOTSTRAP DISTRIBUTION FOR THE MEAN WAGE
# -----------------------------------------------------------------------------
# One sample is drawn from the known synthetic population. Bootstrap samples
# are then drawn from this observed sample, with replacement.

set.seed(SEED)
wages <- rnorm(
  n = BOOT_SAMPLE_SIZE,
  mean = POP_MEAN,
  sd = POP_SD
)

sample_mean <- mean(wages)
sample_sd <- sd(wages)

set.seed(SEED + 1)
bootstrap_means <- replicate(
  B_BOOT,
  mean(sample(wages, size = length(wages), replace = TRUE))
)

bootstrap_se <- sd(bootstrap_means)
bootstrap_ci <- do.call(
  rbind,
  lapply(BOOT_CI_LEVELS, function(level) {
    alpha <- 1 - level
    bounds <- quantile(
      bootstrap_means,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )
    data.frame(
      confidence_level = level,
      lower = bounds[1],
      upper = bounds[2],
      width = bounds[2] - bounds[1]
    )
  })
)

bootstrap_summary <- data.frame(
  statistic = c(
    "Population mean (known in this simulation)",
    "Original sample mean",
    "Original sample standard deviation",
    "Bootstrap standard error"
  ),
  estimate = c(POP_MEAN, sample_mean, sample_sd, bootstrap_se)
)

write.csv(
  data.frame(observation = seq_along(wages), wage = wages),
  project_path("results/original_wage_sample.csv"),
  row.names = FALSE
)
write.csv(bootstrap_summary, project_path("results/bootstrap_summary.csv"), row.names = FALSE)
write.csv(bootstrap_ci, project_path("results/bootstrap_confidence_intervals.csv"), row.names = FALSE)
write.csv(
  data.frame(bootstrap_mean = bootstrap_means),
  project_path("results/bootstrap_distribution.csv"),
  row.names = FALSE
)

# Build common histogram bins so the lower tail, confidence interval, and
# upper tail can be colored separately in the interactive figure.
value_range <- diff(range(bootstrap_means))
if (value_range == 0) value_range <- 1
hist_breaks <- seq(
  min(bootstrap_means) - 0.02 * value_range,
  max(bootstrap_means) + 0.02 * value_range,
  length.out = 41
)
boot_hist <- hist(
  bootstrap_means,
  breaks = hist_breaks,
  plot = FALSE,
  right = FALSE,
  include.lowest = TRUE
)
bin_width <- diff(boot_hist$breaks)[1]
bin_centers <- boot_hist$mids
bin_density <- boot_hist$counts / (sum(boot_hist$counts) * bin_width)
density_max <- max(bin_density) * 1.20

boot_plot <- plot_ly()
boot_buttons <- list()
n_boot_traces_per_level <- 5L
default_boot_level_index <- which.min(abs(BOOT_CI_LEVELS - 0.95))

for (i in seq_along(BOOT_CI_LEVELS)) {
  level <- BOOT_CI_LEVELS[i]
  ci_row <- bootstrap_ci[bootstrap_ci$confidence_level == level, ]
  lower <- ci_row$lower
  upper <- ci_row$upper
  in_interval <- bin_centers >= lower & bin_centers <= upper
  in_lower_tail <- bin_centers < lower
  in_upper_tail <- bin_centers > upper
  show_this_level <- (i == default_boot_level_index)
  
  boot_plot <- add_bars(
    boot_plot,
    x = bin_centers,
    y = ifelse(in_lower_tail, bin_density, 0),
    width = bin_width * 0.95,
    name = "Lower tail (outside CI)",
    marker = list(color = "#E69F8A"),
    visible = show_this_level,
    hovertemplate = "Bootstrap mean: %{x:.2f}<br>Density: %{y:.4f}<extra></extra>"
  )
  boot_plot <- add_bars(
    boot_plot,
    x = bin_centers,
    y = ifelse(in_interval, bin_density, 0),
    width = bin_width * 0.95,
    name = paste0(level * 100, "% confidence interval"),
    marker = list(color = "#4C9F9A"),
    visible = show_this_level,
    hovertemplate = "Bootstrap mean: %{x:.2f}<br>Density: %{y:.4f}<extra></extra>"
  )
  boot_plot <- add_bars(
    boot_plot,
    x = bin_centers,
    y = ifelse(in_upper_tail, bin_density, 0),
    width = bin_width * 0.95,
    name = "Upper tail (outside CI)",
    marker = list(color = "#E69F8A"),
    visible = show_this_level,
    hovertemplate = "Bootstrap mean: %{x:.2f}<br>Density: %{y:.4f}<extra></extra>"
  )
  boot_plot <- add_trace(
    boot_plot,
    x = c(lower, lower), y = c(0, density_max),
    type = "scatter", mode = "lines",
    name = paste0("Lower endpoint: ", round(lower, 2)),
    line = list(color = "#187C78", dash = "dash", width = 2),
    visible = show_this_level, hoverinfo = "skip"
  )
  boot_plot <- add_trace(
    boot_plot,
    x = c(upper, upper), y = c(0, density_max),
    type = "scatter", mode = "lines",
    name = paste0("Upper endpoint: ", round(upper, 2)),
    line = list(color = "#187C78", dash = "dash", width = 2),
    visible = show_this_level, hoverinfo = "skip"
  )
  
  visible <- rep(FALSE, n_boot_traces_per_level * length(BOOT_CI_LEVELS) + 2L)
  first_trace <- n_boot_traces_per_level * (i - 1L) + 1L
  visible[first_trace:(first_trace + n_boot_traces_per_level - 1L)] <- TRUE
  visible[(n_boot_traces_per_level * length(BOOT_CI_LEVELS) + 1L):length(visible)] <- TRUE
  
  boot_buttons[[i]] <- list(
    method = "update",
    args = list(
      list(visible = visible),
      list(title = paste0(
        "Bootstrap Distribution of Mean Wage — ", level * 100,
        "% CI (", round(lower, 2), ", ", round(upper, 2), ")"
      ))
    ),
    label = paste0("confidence level = ", level * 100, "%")
  )
}

# Reference lines are shared across all confidence-level choices.
boot_plot <- add_trace(
  boot_plot,
  x = c(sample_mean, sample_mean), y = c(0, density_max),
  type = "scatter", mode = "lines",
  name = paste0("Original sample mean: ", round(sample_mean, 2)),
  line = list(color = "#202020", dash = "dot", width = 2),
  hoverinfo = "skip"
)
boot_plot <- add_trace(
  boot_plot,
  x = c(POP_MEAN, POP_MEAN), y = c(0, density_max),
  type = "scatter", mode = "lines",
  name = paste0("True population mean: ", POP_MEAN),
  line = list(color = "#386CB0", dash = "solid", width = 2),
  hoverinfo = "skip"
)

boot_plot <- layout(
  boot_plot,
  title = paste0(
    "Bootstrap Distribution of Mean Wage — 95% CI (",
    round(bootstrap_ci$lower[default_boot_level_index], 2),
    ", ",
    round(bootstrap_ci$upper[default_boot_level_index], 2),
    ")"
  ),
  xaxis = list(title = "Bootstrap sample mean"),
  yaxis = list(title = "Density", rangemode = "tozero"),
  barmode = "stack",
  updatemenus = list(list(
    type = "dropdown", direction = "down",
    x = 1, xanchor = "right", y = 1.15, yanchor = "top",
    buttons = boot_buttons,
    active = default_boot_level_index - 1L,
    showactive = TRUE
  )),
  margin = list(t = 115, r = 35, b = 100)
)
boot_plot <- config(boot_plot, displayModeBar = TRUE, responsive = TRUE, scrollZoom = TRUE)

saveWidget(
  boot_plot,
  project_path("figures/bootstrap_confidence_interval.html"),
  selfcontained = FALSE
)

# -----------------------------------------------------------------------------
# 1A. INTERACTIVE BOOTSTRAP DISTRIBUTION BY RESAMPLING SIZE AND PARAMETERS
# -----------------------------------------------------------------------------
make_scenario_histogram <- function(states, file, title, filter_prefix,
                                    initial = "All", reference_mean = NULL,
                                    true_mean = NULL) {
  all_values <- unlist(lapply(states, function(s) s$values), use.names = FALSE)
  x_range <- range(all_values, finite = TRUE)
  padding <- max(diff(x_range) * .04, .25)
  x_range <- x_range + c(-padding, padding)
  breaks <- seq(x_range[1], x_range[2], length.out = 41)
  bin_width <- diff(breaks)[1]
  bin_centers <- head(breaks, -1) + bin_width / 2
  colors <- grDevices::hcl.colors(length(states), "Dark 3")
  trace_count <- length(states)
  extra_reference_count <- sum(!is.null(reference_mean), !is.null(true_mean))
  total_traces <- trace_count + extra_reference_count
  visible_all <- rep(TRUE, total_traces)
  widget <- plot_ly(height = 620)
  maximum_density <- 0

  for (i in seq_along(states)) {
    state <- states[[i]]
    counts <- hist(state$values, breaks = breaks, plot = FALSE,
                   right = FALSE, include.lowest = TRUE)$counts
    density <- counts / (sum(counts) * bin_width)
    maximum_density <- max(maximum_density, density)
    widget <- add_bars(
      widget, x = bin_centers, y = density, width = bin_width * .92,
      name = state$label, marker = list(color = colors[i]),
      opacity = .48, visible = if (identical(initial, "All")) TRUE else identical(initial, state$label),
      hovertemplate = "Bootstrap mean: %{x:.3f}<br>Density: %{y:.4f}<extra>%{fullData.name}</extra>"
    )
  }

  if (!is.null(reference_mean)) {
    widget <- add_trace(
      widget, x = c(reference_mean, reference_mean), y = c(0, maximum_density * 1.2),
      type = "scatter", mode = "lines", name = paste0("Sample mean: ", round(reference_mean, 2)),
      line = list(color = "#202020", dash = "dot", width = 2),
      visible = TRUE, hoverinfo = "skip"
    )
  }
  if (!is.null(true_mean)) {
    widget <- add_trace(
      widget, x = c(true_mean, true_mean), y = c(0, maximum_density * 1.2),
      type = "scatter", mode = "lines", name = paste0("True population mean: ", true_mean),
      line = list(color = "#386CB0", width = 2),
      visible = TRUE, hoverinfo = "skip"
    )
  }

  buttons <- list(list(
    method = "update", label = paste0(filter_prefix, "All scenarios"),
    args = list(list(visible = visible_all),
                list(title = paste0(title, " — All scenarios")))
  ))
  for (i in seq_along(states)) {
    visible <- rep(FALSE, total_traces)
    visible[i] <- TRUE
    if (extra_reference_count > 0) {
      visible[(trace_count + 1L):total_traces] <- TRUE
    }
    buttons[[i + 1L]] <- list(
      method = "update", label = paste0(filter_prefix, states[[i]]$label),
      args = list(list(visible = visible),
                  list(title = paste0(title, " — ", filter_prefix, states[[i]]$label)))
    )
  }
  if (identical(initial, "All")) {
    active_button <- 0L
    title_suffix <- "All scenarios"
  } else {
    active_button <- match(initial, vapply(states, function(s) s$label, character(1)))
    if (is.na(active_button)) stop("Unknown initial histogram scenario.", call. = FALSE)
    title_suffix <- initial
  }

  widget <- layout(
    widget,
    title = paste0(title, " — ", title_suffix),
    xaxis = list(title = "Bootstrap sample mean", range = x_range),
    yaxis = list(title = "Density", range = c(0, maximum_density * 1.25)),
    barmode = "overlay",
    legend = list(orientation = "h", x = .5, y = -.18, xanchor = "center"),
    updatemenus = list(list(
      type = "dropdown", direction = "down", x = 1, xanchor = "right",
      y = 1.15, yanchor = "top", buttons = buttons,
      active = active_button, showactive = TRUE
    )),
    margin = list(t = 115, b = 100, r = 35)
  )
  widget <- config(widget, displayModeBar = TRUE, responsive = TRUE, scrollZoom = TRUE)
  saveWidget(widget, file, selfcontained = FALSE)
  invisible(widget)
}

resample_states <- lapply(RESAMPLE_COUNTS, function(B) {
  set.seed(SEED + 1)
  values <- replicate(B, mean(sample(wages, size = length(wages), replace = TRUE)))
  list(label = as.character(B), values = values)
})
make_scenario_histogram(
  resample_states, project_path("figures/bootstrap_distribution.html"),
  "Bootstrap distribution as resampling size changes", "resampling size = ",
  initial = "All", reference_mean = sample_mean, true_mean = POP_MEAN
)

parameter_specs <- expand.grid(mu = PARAMETER_MEANS, sigma = PARAMETER_SDS,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
parameter_states <- lapply(seq_len(nrow(parameter_specs)), function(i) {
  mu <- parameter_specs$mu[i]
  sigma <- parameter_specs$sigma[i]
  set.seed(SEED + 10L)
  sample_values <- rnorm(BOOT_SAMPLE_SIZE, mean = mu, sd = sigma)
  set.seed(SEED + 11L)
  values <- replicate(B_BOOT, mean(sample(sample_values, size = BOOT_SAMPLE_SIZE, replace = TRUE)))
  list(label = paste0("mu = ", mu, ", sigma = ", sigma), values = values)
})
parameter_default <- which(parameter_specs$mu == POP_MEAN & parameter_specs$sigma == POP_SD)
make_scenario_histogram(
  parameter_states, project_path("figures/bootstrap_parameters.html"),
  "Bootstrap distribution under different Normal parameters", "normal parameters = ",
  initial = parameter_states[[parameter_default]]$label
)
# -----------------------------------------------------------------------------
# 2. MONTE CARLO COVERAGE OF CONFIDENCE INTERVALS
# -----------------------------------------------------------------------------
# Each repetition draws a new sample from the same Normal(POP_MEAN, POP_SD^2)
# population. The population standard deviation is treated as known, as in the
# lecture example. This is repeated sampling, not bootstrap resampling.

set.seed(SEED + 2)
coverage_rows <- list()
coverage_row_id <- 1L

for (n in COVER_SAMPLE_SIZES) {
  simulated_samples <- matrix(
    rnorm(B_COVER * n, mean = POP_MEAN, sd = POP_SD),
    nrow = B_COVER,
    ncol = n
  )
  simulated_means <- rowMeans(simulated_samples)
  standard_error <- POP_SD / sqrt(n)
  
  for (level in COVER_CI_LEVELS) {
    z_critical <- qnorm(1 - (1 - level) / 2)
    lower <- simulated_means - z_critical * standard_error
    upper <- simulated_means + z_critical * standard_error
    covered <- (lower <= POP_MEAN) & (POP_MEAN <= upper)
    coverage_rate <- mean(covered)
    
    coverage_rows[[coverage_row_id]] <- data.frame(
      sample_size = n,
      confidence_level = level,
      simulations = B_COVER,
      covered_intervals = sum(covered),
      coverage_rate = coverage_rate,
      coverage_percent = 100 * coverage_rate,
      monte_carlo_se_percent = 100 * sqrt(
        coverage_rate * (1 - coverage_rate) / B_COVER
      ),
      average_interval_width = 2 * z_critical * standard_error
    )
    coverage_row_id <- coverage_row_id + 1L
  }
}

coverage_results <- do.call(rbind, coverage_rows)
coverage_results <- coverage_results[order(
  coverage_results$sample_size,
  coverage_results$confidence_level
), ]
write.csv(coverage_results, project_path("results/coverage_summary.csv"), row.names = FALSE)

# Build a static table for the requested 95% coverage demonstration. The CSV
# above retains results for all nominal levels; this HTML table focuses on 95%.
coverage_95 <- coverage_results[coverage_results$confidence_level == 0.95, ]
coverage_table <- data.frame(
  `Sample size` = coverage_95$sample_size,
  `Nominal confidence` = "95%",
  `Repeated simulations` = format(coverage_95$simulations, big.mark = ","),
  `Intervals covering true mean` = format(coverage_95$covered_intervals, big.mark = ","),
  `Empirical coverage` = sprintf("%.2f%%", coverage_95$coverage_percent),
  `Monte Carlo SE` = sprintf("%.2f percentage points", coverage_95$monte_carlo_se_percent),
  check.names = FALSE
)
html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}
table_headers <- paste0("<th>", html_escape(names(coverage_table)), "</th>", collapse = "")
table_rows <- vapply(seq_len(nrow(coverage_table)), function(i) {
  cells <- paste0("<td>", html_escape(unlist(coverage_table[i, ], use.names = FALSE)), "</td>", collapse = "")
  paste0("<tr>", cells, "</tr>")
}, character(1))
coverage_html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>95% CI coverage simulation</title>",
  "<style>body{font-family:Arial,sans-serif;color:#334155;max-width:1000px;margin:40px auto;padding:0 24px}h1{color:#172B4D}p{color:#64748B}table{border-collapse:collapse;width:100%;margin-top:24px}th,td{padding:12px 14px;border-bottom:1px solid #E2E8F0;text-align:right}th{background:#F1F5F9}th:first-child,td:first-child{text-align:left}</style></head><body>",
  "<h1>Repeated-simulation coverage of 95% confidence intervals</h1>",
  paste0("<p>Each row uses ", format(B_COVER, big.mark = ","),
         " repeated samples from Normal(", POP_MEAN, ", ", POP_SD,
         "<sup>2</sup>). The population SD is known in this simulation.</p>"),
  paste0("<table><thead><tr>", table_headers, "</tr></thead><tbody>"),
  table_rows, "</tbody></table></body></html>"
)
writeLines(coverage_html, project_path("figures/coverage_simulation.html"), useBytes = TRUE)
# -----------------------------------------------------------------------------
# 3. TWO-SIDED HYPOTHESIS TEST USING A SIMULATED NULL DISTRIBUTION
# -----------------------------------------------------------------------------
# Under H0, X = number of heads follows Binomial(N_TOSSES, P_NULL).
# A two-sided result is at least as far from the null expectation as observed.

set.seed(SEED + 3)
null_heads <- rbinom(
  n = B_NULL,
  size = N_TOSSES,
  prob = P_NULL
)

expected_heads <- N_TOSSES * P_NULL
observed_distance <- abs(OBSERVED_HEADS - expected_heads)
heads_grid <- 0:N_TOSSES
as_extreme_as_observed <- abs(heads_grid - expected_heads) >= observed_distance - 1e-12
extreme_count <- sum(abs(null_heads - expected_heads) >= observed_distance - 1e-12)
monte_carlo_p <- (extreme_count + 1) / (B_NULL + 1)
monte_carlo_p_se <- sqrt(monte_carlo_p * (1 - monte_carlo_p) / B_NULL)
exact_two_sided_p <- binom.test(
  x = OBSERVED_HEADS,
  n = N_TOSSES,
  p = P_NULL,
  alternative = "two.sided"
)$p.value

# Exact equal-tail critical values are used to define the discrete rejection
# region. Its actual probability can be smaller than alpha because X is discrete.
critical_values <- do.call(
  rbind,
  lapply(ALPHA_LEVELS, function(alpha) {
    lower_candidates <- 0:N_TOSSES
    lower_probs <- pbinom(
      lower_candidates,
      size = N_TOSSES,
      prob = P_NULL
    )
    valid_lower <- lower_candidates[lower_probs <= alpha / 2]
    lower_critical <- if (length(valid_lower) > 0) max(valid_lower) else NA_integer_
    
    upper_candidates <- 0:N_TOSSES
    upper_probs <- pbinom(
      upper_candidates - 1,
      size = N_TOSSES,
      prob = P_NULL,
      lower.tail = FALSE
    )
    valid_upper <- upper_candidates[upper_probs <= alpha / 2]
    upper_critical <- if (length(valid_upper) > 0) min(valid_upper) else NA_integer_
    
    actual_size <- if (!is.na(lower_critical) && !is.na(upper_critical)) {
      pbinom(lower_critical, size = N_TOSSES, prob = P_NULL) +
        pbinom(
          upper_critical - 1,
          size = N_TOSSES,
          prob = P_NULL,
          lower.tail = FALSE
        )
    } else {
      NA_real_
    }
    
    reject_observed <- if (!is.na(lower_critical) && !is.na(upper_critical)) {
      OBSERVED_HEADS <= lower_critical || OBSERVED_HEADS >= upper_critical
    } else {
      FALSE
    }
    
    data.frame(
      alpha = alpha,
      lower_critical = lower_critical,
      upper_critical = upper_critical,
      actual_null_rejection_probability = actual_size,
      observed_result_rejects = reject_observed
    )
  })
)

decision_at_5 <- if (exact_two_sided_p < 0.05) "Reject H0" else "Fail to reject H0"
hypothesis_summary <- data.frame(
  item = c(
    "Observed number of heads",
    "Null expectation",
    "Monte Carlo repetitions",
    "Monte Carlo two-sided p-value",
    "Exact two-sided binomial p-value",
    "Monte Carlo standard error of p-value",
    "Decision at alpha = 0.05"
  ),
  value = c(
    as.character(OBSERVED_HEADS),
    as.character(expected_heads),
    as.character(B_NULL),
    sprintf("%.5f", monte_carlo_p),
    sprintf("%.5f", exact_two_sided_p),
    sprintf("%.5f", monte_carlo_p_se),
    decision_at_5
  )
)

write.csv(hypothesis_summary, project_path("results/hypothesis_test_summary.csv"), row.names = FALSE)
write.csv(critical_values, project_path("results/rejection_critical_values.csv"), row.names = FALSE)

null_counts <- tabulate(null_heads + 1L, nbins = N_TOSSES + 1L)
null_probabilities <- null_counts / B_NULL
null_distribution <- data.frame(
  heads = heads_grid,
  simulated_count = null_counts,
  simulated_probability = null_probabilities,
  theoretical_probability = dbinom(heads_grid, size = N_TOSSES, prob = P_NULL),
  as_extreme_as_observed = as_extreme_as_observed
)
write.csv(null_distribution, project_path("results/null_distribution.csv"), row.names = FALSE)

# The p-value tails are shown in orange. Red bars are in the selected rejection
# region but not in the p-value tails; the remaining null outcomes are blue.
null_plot <- plot_ly()
null_buttons <- list()
n_null_traces_per_alpha <- 5L
default_alpha_index <- which.min(abs(ALPHA_LEVELS - 0.05))
null_y_max <- max(null_probabilities) * 1.28

pvalue_heads <- heads_grid[as_extreme_as_observed]
lower_pvalue_cutoff <- if (any(pvalue_heads < expected_heads)) {
  max(pvalue_heads[pvalue_heads < expected_heads])
} else {
  NA_real_
}
upper_pvalue_cutoff <- if (any(pvalue_heads > expected_heads)) {
  min(pvalue_heads[pvalue_heads > expected_heads])
} else {
  NA_real_
}

for (i in seq_along(ALPHA_LEVELS)) {
  alpha <- ALPHA_LEVELS[i]
  critical_row <- critical_values[critical_values$alpha == alpha, ]
  lower_critical <- critical_row$lower_critical
  upper_critical <- critical_row$upper_critical
  
  in_rejection_region <- if (!is.na(lower_critical) && !is.na(upper_critical)) {
    heads_grid <= lower_critical | heads_grid >= upper_critical
  } else {
    rep(FALSE, length(heads_grid))
  }
  pvalue_density <- ifelse(as_extreme_as_observed, null_probabilities, 0)
  rejection_only_density <- ifelse(
    in_rejection_region & !as_extreme_as_observed,
    null_probabilities,
    0
  )
  nonrejection_density <- ifelse(
    !in_rejection_region & !as_extreme_as_observed,
    null_probabilities,
    0
  )
  show_this_alpha <- (i == default_alpha_index)
  
  null_plot <- add_bars(
    null_plot,
    x = heads_grid,
    y = nonrejection_density,
    name = "Outside rejection region",
    marker = list(color = "#7298C7"),
    visible = show_this_alpha,
    hovertemplate = "Heads: %{x}<br>Null probability: %{y:.4f}<extra></extra>"
  )
  null_plot <- add_bars(
    null_plot,
    x = heads_grid,
    y = rejection_only_density,
    name = "Rejection region (outside p-value tails)",
    marker = list(color = "#D55E5E"),
    visible = show_this_alpha,
    hovertemplate = "Heads: %{x}<br>Null probability: %{y:.4f}<extra></extra>"
  )
  null_plot <- add_bars(
    null_plot,
    x = heads_grid,
    y = pvalue_density,
    name = "Two-sided p-value tails",
    marker = list(color = "#E7A33E"),
    visible = show_this_alpha,
    hovertemplate = "Heads: %{x}<br>Null probability: %{y:.4f}<extra></extra>"
  )
  
  if (!is.na(lower_critical) && !is.na(upper_critical)) {
    null_plot <- add_trace(
      null_plot,
      x = c(lower_critical + 0.5, lower_critical + 0.5),
      y = c(0, null_y_max),
      type = "scatter", mode = "lines",
      name = paste0("Lower critical boundary: ", lower_critical + 0.5),
      line = list(color = "#B22222", dash = "dash", width = 2),
      visible = show_this_alpha, hoverinfo = "skip"
    )
    null_plot <- add_trace(
      null_plot,
      x = c(upper_critical - 0.5, upper_critical - 0.5),
      y = c(0, null_y_max),
      type = "scatter", mode = "lines",
      name = paste0("Upper critical boundary: ", upper_critical - 0.5),
      line = list(color = "#B22222", dash = "dash", width = 2),
      visible = show_this_alpha, hoverinfo = "skip"
    )
  } else {
    # Keep the same number of traces for each alpha choice.
    null_plot <- add_trace(
      null_plot,
      x = c(NA, NA), y = c(NA, NA),
      type = "scatter", mode = "lines",
      name = "Lower critical boundary",
      visible = show_this_alpha, hoverinfo = "skip"
    )
    null_plot <- add_trace(
      null_plot,
      x = c(NA, NA), y = c(NA, NA),
      type = "scatter", mode = "lines",
      name = "Upper critical boundary",
      visible = show_this_alpha, hoverinfo = "skip"
    )
  }
  
  visible <- rep(FALSE, n_null_traces_per_alpha * length(ALPHA_LEVELS) + 2L)
  first_trace <- n_null_traces_per_alpha * (i - 1L) + 1L
  visible[first_trace:(first_trace + n_null_traces_per_alpha - 1L)] <- TRUE
  visible[(n_null_traces_per_alpha * length(ALPHA_LEVELS) + 1L):length(visible)] <- TRUE
  
  null_buttons[[i]] <- list(
    method = "update",
    args = list(
      list(visible = visible),
      list(title = paste0(
        "Null Distribution of Heads — significance level = ", round(alpha * 100), "%"
      ))
    ),
    label = paste0("significance level = ", round(alpha * 100), "%")
  )
}

# The observed count is marked in the right tail; the mirror cutoff indicates
# the equally extreme outcome in the left tail under the symmetric fair-coin null.
if (!is.na(lower_pvalue_cutoff)) {
  null_plot <- add_trace(
    null_plot,
    x = c(lower_pvalue_cutoff + 0.5, lower_pvalue_cutoff + 0.5),
    y = c(0, null_y_max),
    type = "scatter", mode = "lines",
    name = paste0("Lower p-value cutoff: ", lower_pvalue_cutoff),
    line = list(color = "#7A5195", dash = "dot", width = 2),
    hoverinfo = "skip"
  )
}
if (!is.na(upper_pvalue_cutoff)) {
  null_plot <- add_trace(
    null_plot,
    x = c(upper_pvalue_cutoff - 0.5, upper_pvalue_cutoff - 0.5),
    y = c(0, null_y_max),
    type = "scatter", mode = "lines",
    name = paste0("Observed extreme cutoff: ", upper_pvalue_cutoff),
    line = list(color = "#7A5195", dash = "dot", width = 2),
    hoverinfo = "skip"
  )
}

null_plot <- layout(
  null_plot,
  title = paste0("Null Distribution of Heads — significance level = ", ALPHA_LEVELS[default_alpha_index] * 100, "%"),
  xaxis = list(title = "Number of heads", dtick = 5),
  yaxis = list(title = "Simulated probability", range = c(0, null_y_max)),
  barmode = "stack",
  updatemenus = list(list(
    type = "dropdown", direction = "down",
    x = 1, xanchor = "right", y = 1.15, yanchor = "top",
    buttons = null_buttons,
    active = default_alpha_index - 1L,
    showactive = TRUE
  )),
  margin = list(t = 115, r = 35, b = 100)
)
null_plot <- config(null_plot, displayModeBar = TRUE, responsive = TRUE, scrollZoom = TRUE)

saveWidget(
  null_plot,
  project_path("figures/null_distribution.html"),
  selfcontained = FALSE
)

# -----------------------------------------------------------------------------
# 4. PRINT A SHORT SUMMARY TO THE CONSOLE
# -----------------------------------------------------------------------------
cat("\n--- Bootstrap summary ---\n")
print(bootstrap_summary)
print(bootstrap_ci)

cat("\n--- Coverage results ---\n")
print(coverage_results)

cat("\n--- Two-sided hypothesis test ---\n")
print(hypothesis_summary)
print(critical_values)

cat("\nCSV results and the static 95% coverage table are saved under results/ and figures. Interactive Plotly figures: bootstrap_distribution, bootstrap_parameters, bootstrap_confidence_interval, and null_distribution.\n")






