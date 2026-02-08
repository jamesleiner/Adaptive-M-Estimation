## ------------------------------------------------------------
## Setup simulation parameters
## ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(grid)



# Simulation parameters
sim_params <- expand.grid(
  tau          = c(1.0),
  prob_sample  = c(1.0),
  scenario     = c(1, 2, 3, 4, 5, 6),
  treat_method = c("active", "epsilon_greedy", "UCB", "Thompson", "Thompson-Clip")
)

sim_params$UPDATE_NUM  <- 500
sim_params$BATCH_SIZE  <- 1
sim_params$NUM_BATCHES <- 2000 / sim_params$prob_sample / sim_params$BATCH_SIZE

# Repeat 100 times (100 trials per parameter combo)
sim_params <- do.call("rbind", replicate(n = 100, sim_params, simplify = FALSE))

# Add file names
sim_params$savename <- sapply(
  rownames(sim_params),
  function(x) paste(
    "/ocean/projects/mth230012p/jleiner/adaptive_inference/sims8_",
    x, ".Rdata", sep = ""
  )
)

# Parameter columns (drop savename — utility only)
param_cols <- setdiff(names(sim_params), "savename")

## ------------------------------------------------------------
## Output directory
## ------------------------------------------------------------

outdir <- "/ocean/projects/mth230012p/jleiner/adaptive_inference/aggregated_results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ------------------------------------------------------------
## Helper: process one component at a time (cover, width, etc.)
##          Keeps only every 10th iteration row
## ------------------------------------------------------------

process_component <- function(component_name, prefix, save_raw = FALSE) {
  cat("\n=== Processing component:", component_name, "===\n")
  
  part_list <- list()
  idx <- 1L
  
  for (i in seq_len(nrow(sim_params))) {
    file <- sim_params$savename[i]
    if (!file.exists(file)) next
    
    cat("Loading:", file, "\n")
    
    # Clear any old 'res'
    if (exists("res")) rm(res)
    
    t <- try(load(file), silent = TRUE)  # expects object 'res'
    if (inherits(t, "try-error")) next
    if (!exists("res") || !is.list(res)) next
    
    comp_obj <- res[[component_name]]
    if (is.null(comp_obj)) next
    
    df <- as.data.frame(comp_obj)
    
    # iteration index: 1..nrow(df)
    df$iter <- seq_len(nrow(df))
    
    # <<<<<< THIN TO EVERY 10TH ROW HERE >>>>>>
    df <- df[df$iter %% 10 == 0, , drop = FALSE]
    if (nrow(df) == 0) next
    
    # parameters for this run (no savename)
    params <- sim_params[i, param_cols, drop = FALSE]
    df <- cbind(df, params[rep(1, nrow(df)), , drop = FALSE])
    
    part_list[[idx]] <- df
    idx <- idx + 1L
  }
  
  if (length(part_list) == 0) {
    cat("No data found for component:", component_name, "\n")
    return(invisible(NULL))
  }
  
  # Trim list in case idx didn't use all entries
  if (idx > 1L) {
    part_list <- part_list[seq_len(idx - 1L)]
  }
  
  cat("Binding rows for", component_name, "...\n")
  big <- dplyr::bind_rows(part_list)
  rm(part_list); gc()
  
  cat("Rows in stacked (thinned)", component_name, ":", nrow(big), "\n")
  
  # measurement columns = everything that's not a parameter or iter
  measure_cols <- setdiff(names(big), c(param_cols, "iter"))
  
  cat("Aggregating (averaging over trials) for", component_name, "...\n")
  big_mean <- big %>%
    group_by(across(all_of(c(param_cols, "iter")))) %>%
    summarise(
      across(all_of(measure_cols), mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save results
  if (save_raw) {
    saveRDS(big, file.path(outdir, paste0(prefix, "_thinned_raw.rds")))
  }
  saveRDS(big_mean, file.path(outdir, paste0(prefix, "_thinned_mean.rds")))
  
  # Drop large objects from memory
  rm(big, big_mean); gc()
  
  cat("Done with", component_name, "\n")
  invisible(NULL)
}

## ------------------------------------------------------------
## Run for each component, one at a time
## ------------------------------------------------------------

# Set save_raw = TRUE if you *also* want the stacked, thinned dataset
process_component("cover",  "cover",  save_raw = FALSE)
process_component("width",  "width",  save_raw = FALSE)
process_component("coverc", "coverc", save_raw = FALSE)
process_component("witdthc", "witdthc", save_raw = FALSE)



outdir <- "/ocean/projects/mth230012p/jleiner/adaptive_inference/aggregated_results"

cover_mean  <- readRDS(file.path(outdir, "cover_thinned_mean.rds"))
width_mean  <- readRDS(file.path(outdir, "width_thinned_mean.rds"))
coverc_mean <- readRDS(file.path(outdir, "coverc_thinned_mean.rds"))
widthc_mean <- readRDS(file.path(outdir, "witdthc_thinned_mean.rds"))
colnames(coverc_mean) <- colnames(cover_mean)
colnames(widthc_mean) <- colnames(width_mean)

outdir <- "figures"

# ------------------------------------------------------------------
# Treatment method readable labels
# ------------------------------------------------------------------
treat_labels <- c(
  "active"          = "Uniform",
  "epsilon_greedy"  = "Epsilon-Greedy",
  "UCB"             = "UCB",
  "Thompson"        = "Thompson",
  "Thompson-Clip"   = "Thompson-Clip"
)

# ------------------------------------------------------------------
# Estimator order (Oracle removed)
# ------------------------------------------------------------------
estimator_order <- c(
  "Naive",
  "IPW",
  "SQ-IPW",
  "MAIPW-External",
  "MAIPW-Reuse",
  "MAIPW-Split"
)

est_colors <- c(
  "Naive"          = "#000000",  # black
  "IPW"            = "#E69F00",  # orange
  "SQ-IPW"         = "#56B4E9",  # sky blue
  "MAIPW-External" = "#009E73",  # teal/green
  "MAIPW-Reuse"    = "#CC79A7",  # magenta
  "MAIPW-Split"    = "#0072B2"   # blue
)

est_linetypes <- c(
  "Naive"          = "solid",
  "IPW"            = "longdash",
  "SQ-IPW"         = "dotdash",
  "MAIPW-External" = "solid",
  "MAIPW-Reuse"    = "dashed",
  "MAIPW-Split"    = "dotted"
)

# ==================================================================
# Helper function to plot one metric (coverage or width)
#   file_suffix controls filenames (e.g. "" vs "_linear")
# ==================================================================
plot_metric_by_scenario <- function(df, metric_name, is_coverage = FALSE,
                                    file_suffix = "") {
  
  for (scenario_to_plot in 1:6) {
    cat("Plotting", metric_name, file_suffix,
        "scenario", scenario_to_plot, "...\n")
    
    # -----------------------------
    # Subset by scenario + iter ≤ 1000
    # -----------------------------
    d <- df %>%
      filter(scenario == scenario_to_plot, iter <= 1000)
    
    # -----------------------------
    # Select estimator columns
    # -----------------------------
    base_estimators <- c("Naive")
    projected_cols  <- grep("\\(Projected\\)", names(d), value = TRUE)
    cols_to_keep    <- intersect(c(base_estimators, projected_cols), names(d))
    
    # -----------------------------
    # Long format
    # -----------------------------
    d_long <- d %>%
      pivot_longer(
        cols = all_of(cols_to_keep),
        names_to  = "estimator_raw",
        values_to = "value"
      ) %>%
      mutate(
        estimator = gsub(" \\(Projected\\)", "", estimator_raw),
        estimator = factor(estimator, levels = estimator_order),
        treat_method = factor(treat_method, levels = names(treat_labels))
      ) %>%
      filter(!is.na(estimator))
    
    # -----------------------------
    # Base plot
    # -----------------------------
    p <- ggplot(
      d_long,
      aes(
        x = iter,
        y = value,
        color = estimator,
        linetype = estimator,
        shape = estimator
      )
    ) +
      geom_line(linewidth = 1.5) +
      geom_point(
        data = d_long %>% filter(iter %% 100 == 0),
        size = 2
      ) +
      facet_wrap(
        ~ treat_method,
        nrow = 1,
        scales = "free_y",
        labeller = labeller(treat_method = treat_labels)
      ) +
      labs(
        x = "Iteration",
        y = metric_name,
        color = NULL,
        linetype = NULL,
        shape = NULL
      ) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "grey90"),
        axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 18),
        legend.key.size = grid::unit(2.0, "lines"),
        legend.spacing.x = grid::unit(1.0, "cm"),
        strip.text = element_text(size = 14)
      )
    
    # -----------------------------
    # Manual scales (MAIPWM separation)
    # -----------------------------
    p <- p +
      scale_color_manual(
        values = c(
          "Naive"          = "#000000",
          "IPW"            = "#E69F00",
          "SQ-IPW"         = "#56B4E9",
          "MAIPW-External" = "#009E73",
          "MAIPW-Reuse"    = "#CC79A7",
          "MAIPW-Split"    = "#0072B2"
        ),
        labels = function(x) gsub("^MAIPW", "MAIPWM", x)
      ) +
      scale_linetype_manual(
        values = c(
          "Naive"          = "solid",
          "IPW"            = "longdash",
          "SQ-IPW"         = "dotdash",
          "MAIPW-External" = "solid",
          "MAIPW-Reuse"    = "dashed",
          "MAIPW-Split"    = "dotted"
        ),
        labels = function(x) gsub("^MAIPW", "MAIPWM", x)
      ) +
      scale_shape_manual(
        values = c(
          "Naive"          = NA,
          "IPW"            = NA,
          "SQ-IPW"         = NA,
          "MAIPW-External" = 16,  # circle
          "MAIPW-Reuse"    = 17,  # triangle
          "MAIPW-Split"    = 15   # square
        ),
        labels = function(x) gsub("^MAIPW", "MAIPWM", x)
      ) +
      guides(
        color    = guide_legend(nrow = 1, override.aes = list(linewidth = 1.6)),
        linetype = guide_legend(nrow = 1, override.aes = list(linewidth = 1.6)),
        shape    = guide_legend(nrow = 1, override.aes = list(size = 3))
      )
    
    # -----------------------------
    # Coverage line / legend handling
    # -----------------------------
    if (is_coverage) {
      p <- p + geom_hline(yintercept = 0.9, linetype = "dotted")
    } else {
      p <- p + theme(legend.position = "none")
    }
    
    # -----------------------------
    # Save
    # -----------------------------
    outfile <- file.path(
      outdir,
      paste0(metric_name, file_suffix, "_scenario_", scenario_to_plot, ".pdf")
    )
    
    ggsave(outfile, p, width = 14, height = 4.5)
  }
}
# ==================================================================
# Original coverage/width
# ==================================================================
plot_metric_by_scenario(
  df = cover_mean,
  metric_name = "Coverage",
  is_coverage = TRUE,
  file_suffix = ""
)

plot_metric_by_scenario(
  df = width_mean,
  metric_name = "Width",
  is_coverage = FALSE,
  file_suffix = ""
)

# ==================================================================
# Linear versions: coverc_mean, widthc_mean
# ==================================================================
plot_metric_by_scenario(
  df = coverc_mean,
  metric_name = "Coverage",
  is_coverage = TRUE,
  file_suffix = "_linear"
)

plot_metric_by_scenario(
  df = widthc_mean,   # make sure this matches your object name
  metric_name = "Width",
  is_coverage = FALSE,
  file_suffix = "_linear"
)