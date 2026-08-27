library(lrstat)
library(ggplot2)

method_labels <- lrstat:::lrsim_bmTrtSel_method_labels()

dimnames <- list(c("High Dose", "Low Dose"),
                 c("Nonresponders", "Responders"))

control_hazards <- c(log(2) / 12, log(2) / 24)
scenario_specs <- list(
  null = list(label = "Null",
              responseProbTreatments = c(0.4, 0.4),
              hazardRateControl = control_hazards,
              hazardRateTreatments = structure(
                rbind(control_hazards, control_hazards),
                dimnames = dimnames)),
  dose_response = list(label = "Effective, dose-response",
                       responseProbTreatments = c(0.6, 0.5),
                       hazardRateControl = control_hazards,
                       hazardRateTreatments = structure(
                         rbind(control_hazards * 0.75,
                               control_hazards * 0.75),
                         dimnames = dimnames)),
  no_dose_response = list(label = "Effective, no dose-response",
                          responseProbTreatments = c(0.6, 0.6),
                          hazardRateControl = control_hazards,
                          hazardRateTreatments = structure(
                            rbind(control_hazards * 0.75,
                                  control_hazards * 0.75),
                            dimnames = dimnames)),
  one_effective = list(label = "One dose effective",
                       responseProbTreatments = c(0.6, 0.4),
                       hazardRateControl = control_hazards,
                       hazardRateTreatments = structure(
                         rbind(control_hazards * 0.75,
                               control_hazards),
                         dimnames = dimnames)),
  unequal_hr = list(label = "Unequal HRs",
                    responseProbTreatments = c(0.6, 0.5),
                    hazardRateControl = c(log(2) / 15, log(2) / 15),
                    hazardRateTreatments = structure(
                      rbind(c(log(2) / 15, log(2) / 15) * 0.70,
                            c(log(2) / 15, log(2) / 15) * 0.66),
                      dimnames = dimnames)),
  nonresponder_effect = list(label = "Nonresponder effect",
                             responseProbTreatments = c(0.6, 0.5),
                             hazardRateControl = control_hazards,
                             hazardRateTreatments = structure(
                               rbind(c(log(2) / 24, log(2) / 24),
                                     c(log(2) / 24, log(2) / 24)),
                               dimnames = dimnames))
)

# Explicit (scenario, phase2SampleSizePerArm) combinations to run, rather
# than expand.grid's full cross product. Edit these lists to add/remove
# scenarios or sample sizes for the base and sensitivity analyses.
base_scenarios <- c("null", "dose_response", "no_dose_response",
                    "one_effective")
base_n2 <- 50
sensitivity_specs <- list(
  list(scenario = c("null", "dose_response"), phase2SampleSizePerArm = 20),
  list(scenario = c("unequal_hr", "nonresponder_effect"),
       phase2SampleSizePerArm = 50)
)

make_design_grid <- function(base_scenarios, base_n2, sensitivity_specs) {
  rows <- list(data.frame(scenario = base_scenarios,
                          phase2SampleSizePerArm = base_n2,
                          stringsAsFactors = FALSE))
  for (spec in sensitivity_specs) {
    rows <- c(rows, list(data.frame(
      scenario = spec$scenario,
      phase2SampleSizePerArm = spec$phase2SampleSizePerArm,
      stringsAsFactors = FALSE)))
  }
  unique(do.call(rbind, rows))
}

design_grid <- make_design_grid(base_scenarios, base_n2, sensitivity_specs)

# --- Splitting the simulation grid across multiple machines ---------------
# Usage:
#   Rscript call_lrsim_bmTrtSel.R                 # run everything, plot (original behavior)
#   Rscript call_lrsim_bmTrtSel.R run 1 2          # machine 1 of 2: runs its share, saves RDS
#   Rscript call_lrsim_bmTrtSel.R run 2 2          # machine 2 of 2: runs its share, saves RDS
#   Rscript call_lrsim_bmTrtSel.R combine 2        # after copying both RDS files to one folder,
#                                                   # combine results and produce the plot
# Rows are assigned round-robin (row %% nparts) so each machine gets a mix of
# scenarios/sample sizes rather than one machine getting all the slow ones.
results_file <- function(part, nparts) {
  sprintf("results_part%d_of_%d.rds", part, nparts)
}

run_scenario <- function(scenario, phase2SampleSizePerArm, seed,
                         scenario_index = NA_integer_,
                         scenario_total = NA_integer_,
                         progress = TRUE) {
  spec <- scenario_specs[[scenario]]
  if (progress) {
    message(sprintf("[%d/%d] Starting %s, Phase II n/arm = %d",
                    scenario_index, scenario_total, spec$label,
                    phase2SampleSizePerArm))
  }
  started_at <- Sys.time()
  sim <- lrsim_bmTrtSel(
    phase2SampleSizePerArm = phase2SampleSizePerArm,
    phase3SampleSizePerArmMin = 113, phase3SampleSizePerArmMax = 163,
    responseProbControl = 0.4,
    responseProbTreatments = spec$responseProbTreatments,
    toxicityProbTreatments = c(0, 0), corrEfficacyToxicity = 0,
    hazardRateControl = spec$hazardRateControl,
    hazardRateTreatments = spec$hazardRateTreatments,
    studyDurationPhase3 = 42.1, toxicityWeight = 0, toxicityUpperLimit = 1,
    efficacyThreshold = 0, safetyThreshold = 0,
    methods = c("ctbonferroni", "ctdunnett", "ctsimes", "ctpooled", "cer",
          "tsssd.k", "tsssd.uk", "tsssd.k.rank", "tsssd.uk.rank",
          "tsssd.k.ce", "tsssd.uk.ce",
          "tsssd.k.rank.ce", "tsssd.uk.rank.ce",
          "naive", "ph3only"),
    accrualRatePhase2 = 3, accrualRatePhase3 = 6, followupTimePhase2 = 6,
    maxNumberOfIterations = 10000, seed = seed,
    nthreads = 10)
  if (progress) {
    message(sprintf("[%d/%d] Finished %s, Phase II n/arm = %d (%.1f minutes)",
                    scenario_index, scenario_total, spec$label,
                    phase2SampleSizePerArm,
                    as.numeric(difftime(Sys.time(), started_at, units = "mins"))))
  }
  ndose <- length(sim$selectionProb)
  ngrid <- length(sim$n2)
  dose_names <- paste0("dose.", seq_len(ndose))
  selection_prob <- as.data.frame(
    matrix(rep(sim$selectionProb, each = ngrid), nrow = ngrid))
  names(selection_prob) <- paste0("selection.prob.", dose_names)

  overview <- data.frame(
    n = sim$n2, n1 = sim$n1,
    numberOfIterations = sim$numberOfIterations, trueOBD = sim$trueOBD,
    pcs = sim$pcs, ave.event.stage1 = sim$ave.event[, 1],
    ave.event.stage2 = sim$ave.event[, 2], ave.event.total = sim$ave.event[, 3],
    Scenario = spec$label, phase2SampleSizePerArm = phase2SampleSizePerArm
  )

  do.call(rbind, unname(lapply(names(sim$byMethod), function(method_id) {
    method <- sim$byMethod[[method_id]]
    prob_rej_each <- as.data.frame(
      matrix(method$prob.rej.each, nrow = ngrid, ncol = ndose))
    names(prob_rej_each) <- paste0("prob.rej.each.", dose_names)
    data.frame(
      overview, selection_prob, method = method_id,
      Method = factor(method_labels[method_id], levels = unname(method_labels)),
      gpower = method$gpower, prob = method$prob.rej.any,
      prob.rej.any = method$prob.rej.any, prob_rej_each,
      row.names = NULL, check.names = FALSE
    )
  })))
}

plot_one <- function(results) {
  results$Method <- factor(method_labels[results$method],
                           levels = unname(method_labels))
  scenario <- unique(results$Scenario)
  n2 <- unique(results$phase2SampleSizePerArm)
  ggplot(results, aes(n, prob, color = Method, shape = Method)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
    scale_color_manual(
      values = c("CT-Bonferroni" = "#E69F00",
             "CT-Dunnett" = "#D55E00",
                 "CT-Simes" = "#0072B2",
                 "CT-Pooled" = "#56B4E9",
                 "CER" = "#009E73",
                 "TSSSD-K" = "#CC79A7",
                 "TSSSD-UK" = "#F0E442",
                 "TSSSD-K-Rank" = "#E69F00",
                 "TSSSD-UK-Rank" = "#44AA99",
                 "TSSSD-K-CE" = "#AA3377",
                 "TSSSD-UK-CE" = "#999933",
                 "TSSSD-K-Rank-CE" = "#EE7733",
                 "TSSSD-UK-Rank-CE" = "#117733",
                 "Naive" = "#000000",
                 "Ph3Only" = "#009E73")) +
    scale_shape_manual(
      values = c("CT-Bonferroni" = 13,
             "CT-Dunnett" = 16,
                 "CT-Simes" = 17,
                 "CT-Pooled" = 15,
                 "CER" = 4,
                 "TSSSD-K" = 3,
                 "TSSSD-UK" = 7,
                 "TSSSD-K-Rank" = 8,
                 "TSSSD-UK-Rank" = 5,
                 "TSSSD-K-CE" = 9,
                 "TSSSD-UK-CE" = 10,
                 "TSSSD-K-Rank-CE" = 11,
                 "TSSSD-UK-Rank-CE" = 14,
                 "Naive" = 0,
                 "Ph3Only" = 18)) +
    labs(x = "Sample size per arm at Stage II (Phase III)",
         y = "Probability of rejecting any null hypothesis",
         title = scenario,
         subtitle = paste0("Phase II n/arm = ", n2)) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

plot_results <- function(results) {
  scenario_levels <- vapply(scenario_specs, `[[`, character(1), "label")
  key <- interaction(factor(results$Scenario, levels = scenario_levels),
                     results$phase2SampleSizePerArm,
                     drop = TRUE, lex.order = TRUE)
  lapply(split(results, key), plot_one)
}

save_plots_pdf <- function(results, file = "lrsim_bmTrtSel_results.pdf") {
  plots <- plot_results(results)
  pdf(file, width = 8, height = 6, onefile = TRUE)
  on.exit(dev.off())
  for (p in plots) print(p)
  message(sprintf("Wrote %d pages to %s", length(plots), file))
  invisible(file)
}

run_grid_subset <- function(part, nparts) {
  row_ids <- seq_len(nrow(design_grid))
  keep <- (row_ids - 1) %% nparts == (part - 1)
  grid_subset <- design_grid[keep, , drop = FALSE]
  results <- do.call(rbind, Map(run_scenario, grid_subset$scenario,
                                grid_subset$phase2SampleSizePerArm,
                                seed = 314159 + row_ids[keep] - 1,
                                scenario_index = row_ids[keep],
                                scenario_total = nrow(design_grid)))
  row.names(results) <- NULL
  saveRDS(results, results_file(part, nparts))
  message(sprintf("Saved %d rows to %s", nrow(results),
                  results_file(part, nparts)))
  results
}

combine_and_plot <- function(nparts) {
  files <- vapply(seq_len(nparts), results_file, character(1), nparts = nparts)
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing result file(s), copy them into the working directory first: ",
        paste(missing, collapse = ", "))
  }
  results <- do.call(rbind, lapply(files, readRDS))
  row.names(results) <- NULL
  save_plots_pdf(results)
  results
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "run" else args[1]

if (mode == "combine") {
  nparts <- if (length(args) >= 2) as.integer(args[2]) else 1
  results <- combine_and_plot(nparts)
} else {
  part <- if (length(args) >= 2) as.integer(args[2]) else 1
  nparts <- if (length(args) >= 3) as.integer(args[3]) else 1
  results <- run_grid_subset(part, nparts)
  if (nparts == 1) {
    save_plots_pdf(results)
  } else {
    message(sprintf(
      "Part %d of %d done. Once all parts finish, copy all %s files into one folder and run: Rscript call_lrsim_bmTrtSel.R combine %d",
      part, nparts, "results_part*_of_*.rds", nparts))
  }
}
