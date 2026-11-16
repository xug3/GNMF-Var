library(FAST)
library(dplyr)

# ====================== USER SETTINGS ==========================
INPUT_DIR <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/inputs"
UTILS_FILE <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/simulation_utils.R"
OUT_DIR <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/simulation_results_graph_only"

N_REP <- 1000
K <- 3
BASE_SEED <- 20261027
DETECTION_THRESHOLD <- 0.15  # Strict >, unchanged.
SAVE_EVERY <- 25

GNMF_LAMBDA_1 <- 10
GNMF_LAMBDA_2 <- 10
GNMF_N_ITER <- 100000
GNMF_CONVERGE <- 1e-5
FACTOR_ZERO_CUTOFF <- 1e-5  # Same W/H cutoff as your uploaded script.

MATCH_METRIC <- "correlation"
MIN_NONZERO_SAMPLES <- 3
NA_TO_ZERO <- TRUE          # Existing missing-value convention, unchanged.
SAVE_REPLICATE_DATA <- FALSE
SAVE_ALL_PROPORTIONS <- TRUE

# No NMF_ENGINE, separate NMF optimizer, or method-specific multiple starts.
METHODS <- c("GNMF-Var", "NMF")

# ====================== INITIALIZE =============================
TRUTH_FILE <- file.path(INPUT_DIR, "truth_proportions.csv")
BARCODE_FILE <- file.path(INPUT_DIR, "barcode_matrix.csv")
SETTINGS_FILE <- file.path(INPUT_DIR, "simulation_settings_example.csv")
source(UTILS_FILE)


dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PROP_DIR <- file.path(OUT_DIR, "estimated_proportions")
if (SAVE_ALL_PROPORTIONS) dir.create(PROP_DIR, showWarnings = FALSE)

write_csv <- function(x, name) {
  write.csv(x, file.path(OUT_DIR, name), row.names = FALSE, na = "NA")
}

# ====================== SHARED MODEL CONFIG ====================
common_config <- dmain_config
common_config$r <- K
common_config$lambda_1 <- GNMF_LAMBDA_1
common_config$lambda_2 <- GNMF_LAMBDA_2
common_config$n_iter <- GNMF_N_ITER
common_config$converge <- GNMF_CONVERGE

gnmf_config <- common_config
nmf_config <- common_config
nmf_config$lambda_1 <- 0  # THE ONLY MODEL-CONFIGURATION CHANGE.

same_fields <- setdiff(names(gnmf_config), "lambda_1")
stopifnot(identical(names(gnmf_config), names(nmf_config)),
          identical(gnmf_config[same_fields], nmf_config[same_fields]),
          identical(nmf_config$lambda_1, 0))

saveRDS(list(GNMF_Var = gnmf_config, NMF = nmf_config),
        file.path(OUT_DIR, "fitting_configs.rds"))
write_csv(data.frame(
  method = METHODS, K = K,
  lambda_1 = c(GNMF_LAMBDA_1, 0), lambda_2 = GNMF_LAMBDA_2,
  max_iterations = GNMF_N_ITER, convergence_setting = GNMF_CONVERGE,
  factor_zero_cutoff = FACTOR_ZERO_CUTOFF,
  fit_seed_rule = "BASE_SEED + 50000000 + replicate",
  starts_per_method = 1, FAST_version = as.character(packageVersion("FAST"))
), "fitting_settings.csv")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
capture.output(print(dmain), file = file.path(OUT_DIR, "dmain_function_used.txt"))

# ====================== HELPERS ================================
row_normalize <- function(M) {
  M <- as.matrix(M)
  if (!is.numeric(M) || any(!is.finite(M)) || any(M < 0)) {
    stop("Cannot normalize nonfinite or negative entries.")
  }
  rs <- rowSums(M)
  if (any(rs <= 0)) {
    stop("All-zero row after factor cutoff: ",
         paste(which(rs <= 0), collapse = ", "),
         ". Stopping rather than treating an invalid fit as no detection.")
  }
  sweep(M, 1, rs, "/")
}

read_simulation_settings <- function(path) {
  z <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(z)) ||
      anyDuplicated(z$parameter)) stop("Invalid simulation settings CSV.")
  setNames(as.list(z$value), z$parameter)
}

get_setting_numeric <- function(s, name, default) {
  if (!name %in% names(s)) return(default)
  x <- suppressWarnings(as.numeric(s[[name]]))
  if (length(x) != 1 || !is.finite(x)) default else x
}

make_chain_graph <- function(n) {
  G <- diag(1, nrow = n, ncol = n)
  if (n > 1) {
    i <- seq_len(n - 1)
    G[cbind(i, i + 1)] <- 1
    G[cbind(i + 1, i)] <- 1
  }
  G
}

extract_dmain_factors <- function(res, X, K) {
  if (is.null(res$W) || is.null(res$H)) stop("dmain returned no W/H.")
  W <- as.matrix(res$W)
  H <- as.matrix(res$H)
  target_W <- c(nrow(X), as.integer(K))
  target_H <- c(ncol(X), as.integer(K))
  if (!identical(dim(W), target_W) && identical(dim(t(W)), target_W)) W <- t(W)
  if (!identical(dim(H), target_H) && identical(dim(t(H)), target_H)) H <- t(H)
  if (!identical(dim(W), target_W) || !identical(dim(H), target_H)) {
    stop("Unexpected dmain W/H dimensions.")
  }
  if (any(!is.finite(W)) || any(!is.finite(H)) || any(W < 0) || any(H < 0)) {
    stop("dmain returned nonfinite or negative factors.")
  }
  list(W = W, H = H)
}

# Shared fitting function: no method-specific optimization logic.
fit_dmain_shared <- function(X, G, config, seed) {
  set.seed(seed)  # Reset immediately before EACH dmain call.
  res <- dmain(X, G, config)
  fit <- extract_dmain_factors(res, X, config$r)
  fit$W[fit$W < FACTOR_ZERO_CUTOFF] <- 0
  fit$H[fit$H < FACTOR_ZERO_CUTOFF] <- 0
  rownames(fit$W) <- rownames(X)
  rownames(fit$H) <- colnames(X)
  fit
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = "pearson")
}

cosine_similarity <- function(x, y) {
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x * y) / den
}

all_permutations <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

match_components_to_variants <- function(W, barcode, variants) {
  shared <- intersect(rownames(W), barcode$mutation_id)
  if (length(shared) < 2) stop("Too few shared barcode mutations.")
  if (ncol(W) != length(variants)) stop("K must equal number of truth lineages.")
  Wsub <- W[shared, , drop = FALSE]
  B <- as.matrix(barcode[match(shared, barcode$mutation_id), variants, drop = FALSE])
  storage.mode(B) <- "double"
  S <- matrix(NA_real_, ncol(W), length(variants))
  for (k in seq_len(ncol(W))) {
    for (j in seq_along(variants)) {
      S[k, j] <- if (MATCH_METRIC == "correlation") {
        safe_cor(Wsub[, k], B[, j])
      } else if (MATCH_METRIC == "cosine") {
        cosine_similarity(Wsub[, k], B[, j])
      } else stop("Unknown MATCH_METRIC.")
    }
  }
  perms <- all_permutations(seq_len(ncol(W)))
  scores <- apply(perms, 1, function(p) {
    vals <- S[cbind(p, seq_along(variants))]
    vals[!is.finite(vals)] <- -1e6
    sum(vals)
  })
  best <- perms[which.max(scores), ]
  similarities <- S[cbind(best, seq_along(variants))]
  if (any(!is.finite(similarities))) stop("No fully defined barcode assignment.")
  list(component_for_variant = best, mapping = data.frame(
    variant = variants, component = best, similarity = similarities
  ))
}

relative_reconstruction_error <- function(X, W, H) {
  den <- sqrt(sum(X^2))
  if (den == 0) return(NA_real_)
  sqrt(sum((X - W %*% t(H))^2)) / den
}

calculate_metrics <- function(truth, estimate) {
  stopifnot(identical(dim(truth), dim(estimate)))
  err <- estimate - truth
  vc <- vapply(seq_len(ncol(truth)), function(j) {
    safe_cor(truth[, j], estimate[, j])
  }, numeric(1))
  list(mae = mean(abs(err)), rmse = sqrt(mean(err^2)),
       correlation = safe_cor(as.vector(truth), as.vector(estimate)),
       mean_variant_correlation = if (any(is.finite(vc))) {
         mean(vc[is.finite(vc)])
       } else NA_real_)
}

first_positive <- function(p, times, threshold = DETECTION_THRESHOLD) {
  stopifnot(length(p) == length(times), all(is.finite(p)))
  ii <- which(p > threshold)
  if (length(ii)) as.numeric(times[ii[1]]) else NA_real_
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) unname(quantile(x, p)) else NA_real_
}

# ====================== FIXED INPUTS ============================
truth_template <- read.csv(TRUTH_FILE, check.names = FALSE)
barcode_template <- read.csv(BARCODE_FILE, check.names = FALSE)
settings <- read_simulation_settings(SETTINGS_FILE)
meta_cols <- intersect(c("sample_id", "time_point", "date"), names(truth_template))
variants <- setdiff(names(truth_template), meta_cols)
truth_template <- truth_template[order(truth_template$time_point), , drop = FALSE]
time_points <- truth_template$time_point

sample_order <- as.character(truth_template$sample_id)
P_reference <- row_normalize(as.matrix(truth_template[, variants, drop = FALSE]))
truth_template[, variants] <- P_reference
amp_ids <- barcode_template$amplicon


truth_first <- data.frame(
  variant = variants,
  truth_time = vapply(seq_along(variants), function(j) {
    first_positive(P_reference[, j], time_points)
  }, numeric(1)),
  threshold = DETECTION_THRESHOLD, criterion = ">"
)

write_csv(truth_first, "ground_truth_first_positive.csv")
write_csv(truth_template, "truth_proportions_used.csv")
print(truth_first, row.names = FALSE)

fixed_amp_prior <- make_amplicon_prior(
  n_amplicons = max(as.integer(amp_ids)), seed = BASE_SEED + 777
)
FIXED_AMP_FILE <- file.path(OUT_DIR, "fixed_amplicon_prior.csv")
write.csv(fixed_amp_prior, FIXED_AMP_FILE, row.names = FALSE)
G <- make_chain_graph(nrow(truth_template))
rownames(G) <- colnames(G) <- sample_order


n_reads <- get_setting_numeric(settings, "n_reads", 100000)
c_amplicon <- get_setting_numeric(settings, "c_amplicon", 200)
dirichlet_model <- get_setting_numeric(settings, "dirichlet_model", 1)
seq_error <- get_setting_numeric(settings, "seq_error", 0.001)
amp_dropout_prob <- get_setting_numeric(settings, "amp_dropout_prob", 0)

# ====================== OUTPUT STORAGE ==========================
replicate_metrics <- list()
variant_metrics <- list()
component_matches <- list()
first_positive_records <- list()

save_progress <- function() {
  tables <- list(
    replicate_metrics = replicate_metrics,
    variant_metrics = variant_metrics,
    component_matching = component_matches,
    first_positive_replicates = first_positive_records
  )
  for (nm in names(tables)) {
    if (length(tables[[nm]])) {
      write_csv(do.call(rbind, tables[[nm]]), paste0(nm, ".csv"))
    }
  }
}

# ====================== PAIRED SIMULATION =======================
for (rep_id in seq_len(N_REP)) {
  cat("Replicate", rep_id, "of", N_REP, "\n")
  rep_seed <- BASE_SEED + (rep_id - 1) * 10000  # Stride keeps per-sample seeds from overlapping.
  fit_seed <- BASE_SEED + 50000000 + rep_id  # Identical for the two method calls.
  rep_dir <- if (SAVE_REPLICATE_DATA) {
    file.path(OUT_DIR, "replicate_data", sprintf("rep_%04d", rep_id))
  } else tempfile(pattern = sprintf("gnmf_rep_%04d_", rep_id))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  
  sim <- simulate_swampy_matrix(
    output_dir = rep_dir, truth_csv = TRUTH_FILE, barcode_csv = BARCODE_FILE,
    amplicon_prior_csv = FIXED_AMP_FILE,
    n_reads = n_reads, c_amplicon = c_amplicon,
    dirichlet_model = dirichlet_model, seq_error = seq_error,
    amp_dropout_prob = amp_dropout_prob,
    include_recurrent_noise = TRUE, recurrent_presence_prob = 0.20,
    unique_keep_fraction = 0.02, seed = rep_seed
  )
  F <- sim$frequency[sample_order, , drop = FALSE]
  if (NA_TO_ZERO) {
    F[is.na(F)] <- 0
  } else if (anyNA(F)) stop("Frequency matrix contains NA.")
  if (any(!is.finite(F)) || any(F < 0)) stop("Invalid frequency matrix.")
  keep <- colSums(F > 0) >= MIN_NONZERO_SAMPLES
  F <- F[, keep, drop = FALSE]
  if (ncol(F) < K) stop("Too few mutations after filtering.")
  X <- t(F)  # mutation x sample; shared by both methods.
  
  P_true <- as.matrix(sim$truth[
    match(sample_order, sim$truth$sample_id), variants, drop = FALSE
  ])
  storage.mode(P_true) <- "double"
  P_true <- row_normalize(P_true)
  rownames(P_true) <- sample_order
  colnames(P_true) <- variants
  if (!isTRUE(all.equal(unname(P_true), unname(P_reference), tolerance = 1e-12))) {
    stop("Replicate truth differs from the fixed input trajectory.")
  }
  
  # SAME X, SAME G, SAME FITTING FUNCTION, SAME SEED.
  # Configs differ ONLY in lambda_1.
  gfit <- fit_dmain_shared(X, G, config = gnmf_config, seed = fit_seed)
  nfit <- fit_dmain_shared(X, G, config = nmf_config, seed = fit_seed)
  
  estimates <- list()
  for (method in METHODS) {
    fit <- if (method == "GNMF-Var") gfit else nfit
    mapping <- match_components_to_variants(fit$W, sim$barcode, variants)
    matched_H <- fit$H[, mapping$component_for_variant, drop = FALSE]
    colnames(matched_H) <- variants
    P_hat <- row_normalize(matched_H)[sample_order, variants, drop = FALSE]
    estimates[[method]] <- P_hat
    met <- calculate_metrics(P_true, P_hat)
    
    replicate_metrics[[length(replicate_metrics) + 1]] <- data.frame(
      replicate = rep_id, seed = rep_seed, fit_seed = fit_seed,
      method = method, K = K, n_mutations = nrow(X),
      lambda_1 = if (method == "GNMF-Var") GNMF_LAMBDA_1 else 0,
      lambda_2 = GNMF_LAMBDA_2,
      MAE = met$mae, RMSE = met$rmse, correlation = met$correlation,
      mean_variant_correlation = met$mean_variant_correlation,
      relative_reconstruction_error = relative_reconstruction_error(X, fit$W, fit$H),
      mean_barcode_similarity = mean(mapping$mapping$similarity)
    )
    for (v in variants) {
      e <- P_hat[, v] - P_true[, v]
      variant_metrics[[length(variant_metrics) + 1]] <- data.frame(
        replicate = rep_id, method = method, variant = v,
        MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)),
        correlation = safe_cor(P_true[, v], P_hat[, v])
      )
      t_true <- first_positive(P_true[, v], time_points)
      t_hat <- first_positive(P_hat[, v], time_points)
      first_positive_records[[length(first_positive_records) + 1]] <- data.frame(
        replicate = rep_id, method = method, variant = v,
        truth_time = t_true, estimated_time = t_hat,
        time_difference = t_hat - t_true,
        threshold = DETECTION_THRESHOLD, criterion = ">"
      )
    }
    map <- mapping$mapping
    map$replicate <- rep_id
    map$method <- method
    component_matches[[length(component_matches) + 1]] <- map[
      , c("replicate", "method", "variant", "component", "similarity")
    ]
  }
  
  gP <- estimates[["GNMF-Var"]]
  nP <- estimates[["NMF"]]
  audit <- data.frame(
    sample_id = sample_order, P_true,
    setNames(as.data.frame(gP), paste0("GNMF_", variants)),
    setNames(as.data.frame(nP), paste0("NMF_", variants)),
    check.names = FALSE
  )
  # Preserve the old audit layout for existing plotting scripts.
  if (rep_id == 1) write_csv(audit, "replicate_001_proportion_audit.csv")
  if (SAVE_ALL_PROPORTIONS) {
    full_audit <- data.frame(
      replicate = rep_id, time_point = time_points, audit, check.names = FALSE
    )
    write.csv(full_audit,
              file.path(PROP_DIR, sprintf("replicate_%04d_proportions.csv", rep_id)),
              row.names = FALSE)
  }
  if (!SAVE_REPLICATE_DATA) unlink(rep_dir, recursive = TRUE, force = TRUE)
  if (rep_id %% SAVE_EVERY == 0 || rep_id == N_REP) save_progress()
}

# ====================== ERROR SUMMARY ==========================
rep_df <- do.call(rbind, replicate_metrics)
first_df <- do.call(rbind, first_positive_records)
summary_rows <- list()
metric_names <- c("MAE", "RMSE", "correlation", "mean_variant_correlation",
                  "relative_reconstruction_error", "mean_barcode_similarity")
for (m in METHODS) {
  z <- rep_df[rep_df$method == m, , drop = FALSE]
  for (metric in metric_names) {
    x <- z[[metric]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      method = m, metric = metric,
      mean = safe_mean(x), sd = sd(x[is.finite(x)]), median = safe_median(x),
      q025 = safe_quantile(x, 0.025), q975 = safe_quantile(x, 0.975)
    )
  }
}
write_csv(do.call(rbind, summary_rows), "summary_metrics.csv")

# ====================== FIRST-POSITIVE SUMMARY ==================
first_summary_rows <- list()
for (v in variants) {
  for (m in METHODS) {
    z <- first_df[first_df$variant == v & first_df$method == m, , drop = FALSE]
    detected <- is.finite(z$estimated_time)
    error <- z$time_difference
    valid_error <- is.finite(error)
    n <- nrow(z)
    exact_n <- sum(valid_error & error == 0)
    within1_n <- sum(valid_error & abs(error) <= 1)
    within2_n <- sum(valid_error & abs(error) <= 2)
    ci <- binom.test(exact_n, n)$conf.int
    first_summary_rows[[length(first_summary_rows) + 1]] <- data.frame(
      variant = v, method = m, truth_time = z$truth_time[1],
      n_replicates = n, n_detected = sum(detected), detection_rate = mean(detected),
      mean_time = safe_mean(z$estimated_time), median_time = safe_median(z$estimated_time),
      q025_time = safe_quantile(z$estimated_time, 0.025),
      q975_time = safe_quantile(z$estimated_time, 0.975),
      mean_time_difference = safe_mean(error),
      median_time_difference = safe_median(error),
      mean_absolute_time_error = safe_mean(abs(error)),
      exact_n = exact_n, exact_rate = exact_n / n,
      within1_n = within1_n, within1_rate = within1_n / n,
      within2_n = within2_n, within2_rate = within2_n / n,
      earlier_rate = sum(valid_error & error < 0) / n,
      later_rate = sum(valid_error & error > 0) / n,
      no_detection_rate = mean(!detected),
      exact_CI_lower = ci[1], exact_CI_upper = ci[2]
    )
  }
}
first_summary <- do.call(rbind, first_summary_rows)
write_csv(first_summary, "first_positive_summary.csv")

# count() keeps the NA group for replicates with no threshold crossing.
first_freq <- first_df |>
  dplyr::count(variant, method, estimated_time, name = "count") |>
  dplyr::arrange(variant, method, estimated_time)
write_csv(first_freq, "first_positive_frequencies.csv")

# All accuracy denominators use all replicates. Time errors use detected replicates.
# These are retrospective threshold-crossing times from full-series fits.
cat("\nCompleted", N_REP, "paired replicates. Output:", OUT_DIR, "\n")
print(first_summary, row.names = FALSE)
