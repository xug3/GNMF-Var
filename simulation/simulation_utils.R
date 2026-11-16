rdirichlet1 <- function(alpha) {
  if (any(alpha <= 0)) stop("All Dirichlet parameters must be > 0.")
  z <- rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

normalize_prob <- function(x) {
  x <- pmax(as.numeric(x), 0)
  if (sum(x) <= 0) stop("Probability vector must have positive sum.")
  x / sum(x)
}

safe_rmultinom <- function(size, prob) {
  if (size <= 0) return(rep(0L, length(prob)))
  as.integer(rmultinom(1, size = as.integer(size), prob = normalize_prob(prob))[, 1])
}

# -------------------- 73-time-point pandemic truth --------------------

# This reproduces the STRUCTURE of the SWAMPy use case (73 time points;
# Alpha -> Delta -> Omicron), not the exact supplementary abundance table.
# If you have the exact supplementary abundances, set truth_csv below.
make_73_timepoint_truth <- function(n_time = 73,
                                    delta_mid = 34,
                                    delta_scale = 5,
                                    omicron_mid = 49,
                                    omicron_scale = 3) {
  t <- seq_len(n_time)
  delta_takeover <- plogis((t - delta_mid) / delta_scale)
  omicron_takeover <- plogis((t - omicron_mid) / omicron_scale)

  p_omicron <- omicron_takeover
  p_delta <- (1 - p_omicron) * delta_takeover
  p_alpha <- (1 - p_omicron) * (1 - delta_takeover)

  P <- cbind(Alpha = p_alpha, Delta = p_delta, Omicron = p_omicron)
  P <- P / rowSums(P)
  data.frame(sample_id = sprintf("T%03d", t), time_point = t, P,
             check.names = FALSE)
}

# ---------------------- Synthetic barcode matrix ----------------------

# Builds a mutation architecture with lineage-specific and shared sites.
# For publication, replacing this with your actual UShER/Freyja barcode
# matrix is preferable. Expected external barcode format:
# mutation_id, position, amplicon, Alpha, Delta, Omicron (0/1 entries).
make_synthetic_barcodes <- function(variants = c("Alpha", "Delta", "Omicron"),
                                    genome_length = 29903,
                                    n_amplicons = 98,
                                    n_specific = c(Alpha = 18, Delta = 24, Omicron = 40),
                                    n_shared_pair = 5,
                                    n_shared_all = 4,
                                    seed = 2026) {
  set.seed(seed)
  n_specific <- n_specific[variants]
  if (any(is.na(n_specific))) stop("n_specific must be named for all variants.")

  # number of mutation sites needed
  pairs <- if (length(variants) >= 2) combn(variants, 2, simplify = FALSE) else list()
  n_sites <- sum(n_specific) + length(pairs) * n_shared_pair + n_shared_all

  positions <- sort(sample(seq_len(genome_length), n_sites, replace = FALSE))
  B <- matrix(0L, nrow = n_sites, ncol = length(variants),
              dimnames = list(NULL, variants))

  idx <- 1L
  source <- character(n_sites)

  for (v in variants) {
    nn <- n_specific[[v]]
    rows <- idx:(idx + nn - 1L)
    B[rows, v] <- 1L
    source[rows] <- paste0("specific_", v)
    idx <- idx + nn
  }

  if (length(pairs) > 0 && n_shared_pair > 0) {
    for (pp in pairs) {
      rows <- idx:(idx + n_shared_pair - 1L)
      B[rows, pp] <- 1L
      source[rows] <- paste0("shared_", paste(pp, collapse = "_"))
      idx <- idx + n_shared_pair
    }
  }

  if (n_shared_all > 0) {
    rows <- idx:(idx + n_shared_all - 1L)
    B[rows, ] <- 1L
    source[rows] <- "shared_all"
  }

  # Simple genome-position -> amplicon map for the matrix-level simulator.
  amp <- pmin(n_amplicons,
              floor((positions - 1) / genome_length * n_amplicons) + 1L)

  out <- data.frame(
    mutation_id = paste0("M", positions),
    position = positions,
    amplicon = amp,
    source = source,
    B,
    check.names = FALSE
  )
  out
}

# ------------------------ Amplicon efficiency -------------------------

# SWAMPy uses an experimentally derived amplicon prior. Here, if no prior
# file is supplied, create one non-uniform fixed prior for the whole study.
# The generated prior is exported so every run is auditable.
make_amplicon_prior <- function(n_amplicons = 98, seed = 2027,
                                log_sd = 0.9) {
  set.seed(seed)
  w <- rlnorm(n_amplicons, meanlog = 0, sdlog = log_sd)
  pi <- normalize_prob(w)
  data.frame(amplicon = seq_len(n_amplicons), prior = pi)
}

# ---------------- High-frequency noise mutation pool -----------------

# SWAMPy samples counts of high-frequency errors from Poisson(L * R).
# To keep a rectangular CSV matrix, we pre-generate a pool of error loci.
# "recurrent" here means errors may recur across samples; "unique" loci
# are sample-specific. This is a matrix-level approximation, not FASTQ ART.
make_recurrent_noise_pool <- function(genome_length = 29903,
                                      n_amplicons = 98,
                                      recurrent_sub_rate = 0.003357,
                                      recurrent_keep_fraction = 0.15,
                                      seed = 2028) {
  set.seed(seed)
  n_rec <- rpois(1, lambda = genome_length * recurrent_sub_rate * recurrent_keep_fraction)
  if (n_rec == 0) return(data.frame())
  pos <- sort(sample(seq_len(genome_length), n_rec, replace = FALSE))
  amp <- pmin(n_amplicons,
              floor((pos - 1) / genome_length * n_amplicons) + 1L)
  data.frame(
    mutation_id = paste0("ERR_REC_", pos),
    position = pos,
    amplicon = amp,
    source = "recurrent_error",
    stringsAsFactors = FALSE
  )
}

# ------------------------ One-sample simulator ------------------------

simulate_one_sample <- function(p,
                                barcode,
                                amp_prior,
                                n_reads = 100000,
                                c_amplicon = 200,
                                dirichlet_model = 1,
                                seq_error = 0.001,
                                amp_dropout_prob = 0.00,
                                recurrent_noise = NULL,
                                recurrent_presence_prob = 0.20,
                                recurrent_vaf_alpha = c(0.29, 1.89),
                                unique_sub_rate = 0.002485,
                                unique_keep_fraction = 0.02,
                                genome_length = 29903,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  variants <- names(p)
  p <- normalize_prob(p)
  names(p) <- variants

  n_amp <- nrow(amp_prior)
  pi_a <- normalize_prob(amp_prior$prior)

  # SWAMPy step: genome read counts N_g ~ Multinomial(N, p_g)
  N_g <- safe_rmultinom(n_reads, p)
  names(N_g) <- variants

  # X[g, a]: number of reads from genome g and amplicon a
  X <- matrix(0L, nrow = length(variants), ncol = n_amp,
              dimnames = list(variants, seq_len(n_amp)))

  if (dirichlet_model == 1) {
    # Model 1: one amplicon proportion vector shared across genomes
    p_a <- rdirichlet1(c_amplicon * pi_a)
    for (g in seq_along(variants)) {
      X[g, ] <- safe_rmultinom(N_g[g], p_a)
    }
  } else if (dirichlet_model == 2) {
    # Model 2: independent amplicon proportions for each genome
    for (g in seq_along(variants)) {
      p_ag <- rdirichlet1(c_amplicon * pi_a)
      X[g, ] <- safe_rmultinom(N_g[g], p_ag)
    }
  } else {
    stop("dirichlet_model must be 1 or 2.")
  }

  # Optional whole-amplicon dropout at the sample level.
  if (amp_dropout_prob > 0) {
    drop <- rbinom(n_amp, size = 1, prob = amp_dropout_prob) == 1
    X[, drop] <- 0L
  }

  # True lineage barcode sites
  B <- as.matrix(barcode[, variants, drop = FALSE])
  amp_index <- as.integer(barcode$amplicon)
  n_mut <- nrow(barcode)
  depth <- integer(n_mut)
  alt <- integer(n_mut)

  for (j in seq_len(n_mut)) {
    a <- amp_index[j]
    d <- sum(X[, a])
    depth[j] <- d
    if (d == 0) {
      alt[j] <- 0L
      next
    }

    true_alt <- sum(X[, a] * B[j, ])
    true_ref <- d - true_alt

    # Symmetric per-read sequencing error approximation at matrix level.
    kept_alt <- if (true_alt > 0) rbinom(1, true_alt, 1 - seq_error) else 0L
    false_alt <- if (true_ref > 0) rbinom(1, true_ref, seq_error / 3) else 0L
    alt[j] <- min(d, kept_alt + false_alt)
  }

  freq <- ifelse(depth > 0, alt / depth, NA_real_)
  names(freq) <- barcode$mutation_id
  names(depth) <- barcode$mutation_id
  names(alt) <- barcode$mutation_id

  # Recurrent high-frequency error sites (optional)
  noise_freq <- numeric(0)
  noise_depth <- integer(0)
  noise_alt <- integer(0)

  if (!is.null(recurrent_noise) && nrow(recurrent_noise) > 0) {
    for (r in seq_len(nrow(recurrent_noise))) {
      if (runif(1) > recurrent_presence_prob) next
      a <- recurrent_noise$amplicon[r]
      d <- sum(X[, a])
      if (d <= 0) next
      # Beta == 2-dim Dirichlet for alternative/reference fraction
      f <- rbeta(1, recurrent_vaf_alpha[1], recurrent_vaf_alpha[2])
      aa <- rbinom(1, d, f)
      id <- recurrent_noise$mutation_id[r]
      noise_freq[id] <- aa / d
      noise_depth[id] <- d
      noise_alt[id] <- aa
    }
  }

  # Small number of sample-specific unique high-frequency substitutions.
  # unique_keep_fraction keeps CSV size practical after bypassing FASTQ.
  n_unique <- rpois(1, genome_length * unique_sub_rate * unique_keep_fraction)
  if (n_unique > 0) {
    pos <- sample(seq_len(genome_length), n_unique, replace = FALSE)
    amps <- pmin(n_amp,
                 floor((pos - 1) / genome_length * n_amp) + 1L)
    for (u in seq_len(n_unique)) {
      d <- sum(X[, amps[u]])
      if (d <= 0) next
      # Unique errors are attached to one source genome, sampled by p_g.
      g <- sample(seq_along(variants), 1, prob = p)
      max_source_reads <- X[g, amps[u]]
      if (max_source_reads <= 0) next
      f <- rbeta(1, 0.29, 1.89)
      aa <- rbinom(1, max_source_reads, f)
      id <- paste0("ERR_UNQ_", pos[u], "_", sample.int(1e8, 1))
      noise_freq[id] <- aa / d
      noise_depth[id] <- d
      noise_alt[id] <- aa
    }
  }

  list(
    freq = c(freq, noise_freq),
    depth = c(depth, noise_depth),
    alt = c(alt, noise_alt),
    X = X,
    N_g = N_g
  )
}

# ----------------------- Full dataset simulator -----------------------

simulate_swampy_matrix <- function(output_dir = "swampy_csv_output",
                                   truth_csv = NULL,
                                   barcode_csv = NULL,
                                   amplicon_prior_csv = NULL,
                                   n_reads = 100000,
                                   c_amplicon = 200,
                                   dirichlet_model = 1,
                                   seq_error = 0.001,
                                   amp_dropout_prob = 0.00,
                                   include_recurrent_noise = TRUE,
                                   recurrent_presence_prob = 0.20,
                                   unique_keep_fraction = 0.02,
                                   seed = 20261027) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  # 1) Truth proportions
  if (is.null(truth_csv)) {
    truth <- make_73_timepoint_truth()
  } else {
    truth <- read.csv(truth_csv, check.names = FALSE)
  }

  meta_cols <- intersect(c("sample_id", "time_point", "date"), names(truth))
  variant_cols <- setdiff(names(truth), meta_cols)
  if (length(variant_cols) < 2) stop("Truth CSV must contain at least 2 variant columns.")

  P <- as.matrix(truth[, variant_cols, drop = FALSE])
  P <- P / rowSums(P)
  truth[, variant_cols] <- P

  # 2) Barcode matrix
  if (is.null(barcode_csv)) {
    default_counts <- setNames(rep(25, length(variant_cols)), variant_cols)
    # Give the default 3-lineage use case distinct mutation architectures.
    if (all(c("Alpha", "Delta", "Omicron") %in% variant_cols)) {
      default_counts[c("Alpha", "Delta", "Omicron")] <- c(18, 24, 40)
    }
    barcode <- make_synthetic_barcodes(
      variants = variant_cols,
      n_specific = default_counts,
      seed = seed + 1
    )
  } else {
    barcode <- read.csv(barcode_csv, check.names = FALSE)
  }

  required <- c("mutation_id", "position", "amplicon", variant_cols)
  miss <- setdiff(required, names(barcode))
  if (length(miss) > 0) stop("Barcode CSV missing columns: ", paste(miss, collapse = ", "))

  n_amplicons <- max(as.integer(barcode$amplicon), na.rm = TRUE)

  # 3) Amplicon prior
  if (is.null(amplicon_prior_csv)) {
    amp_prior <- make_amplicon_prior(n_amplicons = n_amplicons, seed = seed + 2)
  } else {
    amp_prior <- read.csv(amplicon_prior_csv, check.names = FALSE)
    if (!all(c("amplicon", "prior") %in% names(amp_prior))) {
      stop("Amplicon prior CSV must contain columns: amplicon, prior")
    }
    amp_prior <- amp_prior[order(amp_prior$amplicon), ]
    amp_prior$prior <- normalize_prob(amp_prior$prior)
  }

  # 4) Recurrent error pool
  recurrent_noise <- NULL
  if (include_recurrent_noise) {
    recurrent_noise <- make_recurrent_noise_pool(
      n_amplicons = n_amplicons,
      seed = seed + 3
    )
  }

  # 5) Simulate each sample
  sample_ids <- if ("sample_id" %in% names(truth)) truth$sample_id else sprintf("S%03d", seq_len(nrow(truth)))
  results <- vector("list", nrow(truth))
  all_mutations <- barcode$mutation_id

  for (i in seq_len(nrow(truth))) {
    p <- as.numeric(P[i, ])
    names(p) <- variant_cols
    results[[i]] <- simulate_one_sample(
      p = p,
      barcode = barcode,
      amp_prior = amp_prior,
      n_reads = n_reads,
      c_amplicon = c_amplicon,
      dirichlet_model = dirichlet_model,
      seq_error = seq_error,
      amp_dropout_prob = amp_dropout_prob,
      recurrent_noise = recurrent_noise,
      recurrent_presence_prob = recurrent_presence_prob,
      unique_keep_fraction = unique_keep_fraction,
      seed = seed + 1000 + i
    )
    all_mutations <- union(all_mutations, names(results[[i]]$freq))
  }

  # Rectangular matrices
  F <- matrix(0, nrow = nrow(truth), ncol = length(all_mutations),
              dimnames = list(sample_ids, all_mutations))
  D <- matrix(0L, nrow = nrow(truth), ncol = length(all_mutations),
              dimnames = list(sample_ids, all_mutations))
  A <- matrix(0L, nrow = nrow(truth), ncol = length(all_mutations),
              dimnames = list(sample_ids, all_mutations))

  # Keep no-coverage sites as NA in frequency matrix, but zero depth/alt count.
  F[,] <- NA_real_

  for (i in seq_len(nrow(truth))) {
    ids <- names(results[[i]]$freq)
    F[i, ids] <- results[[i]]$freq
    D[i, names(results[[i]]$depth)] <- results[[i]]$depth
    A[i, names(results[[i]]$alt)] <- results[[i]]$alt
  }

  # Sites never explicitly generated as a noise mutation but covered are 0 VAF.
  # Barcode sites already have their observed value/NA. For dynamic noise columns,
  # absence in a covered sample is represented as 0 for benchmarking convenience.
  noise_cols <- setdiff(all_mutations, barcode$mutation_id)
  if (length(noise_cols) > 0) {
    tmp_noise <- F[, noise_cols, drop = FALSE]
    tmp_noise[is.na(tmp_noise)] <- 0
    F[, noise_cols] <- tmp_noise
  }

  # 6) Write CSV outputs
  write.csv(truth, file.path(output_dir, "truth_proportions.csv"), row.names = FALSE)
  write.csv(barcode, file.path(output_dir, "barcode_matrix.csv"), row.names = FALSE)
  write.csv(amp_prior, file.path(output_dir, "amplicon_prior.csv"), row.names = FALSE)

  write.csv(data.frame(sample_id = rownames(F), F, check.names = FALSE),
            file.path(output_dir, "mutation_frequency.csv"), row.names = FALSE, na = "NA")
  write.csv(data.frame(sample_id = rownames(D), D, check.names = FALSE),
            file.path(output_dir, "mutation_depth.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = rownames(A), A, check.names = FALSE),
            file.path(output_dir, "mutation_alt_count.csv"), row.names = FALSE)

  # Compact mutation metadata for all columns
  meta <- data.frame(
    mutation_id = all_mutations,
    mutation_class = ifelse(all_mutations %in% barcode$mutation_id, "lineage_barcode",
                            ifelse(grepl("^ERR_REC_", all_mutations), "recurrent_error", "unique_error")),
    stringsAsFactors = FALSE
  )
  meta <- merge(meta,
                barcode[, intersect(c("mutation_id", "position", "amplicon", "source"), names(barcode)), drop = FALSE],
                by = "mutation_id", all.x = TRUE, sort = FALSE)
  write.csv(meta, file.path(output_dir, "mutation_metadata.csv"), row.names = FALSE)

  settings <- data.frame(
    parameter = c(
      "seed", "n_samples", "n_reads", "c_amplicon", "dirichlet_model",
      "seq_error", "amp_dropout_prob", "include_recurrent_noise",
      "recurrent_presence_prob", "unique_keep_fraction"
    ),
    value = c(
      seed, nrow(truth), n_reads, c_amplicon, dirichlet_model,
      seq_error, amp_dropout_prob, include_recurrent_noise,
      recurrent_presence_prob, unique_keep_fraction
    ),
    stringsAsFactors = FALSE
  )
  write.csv(settings, file.path(output_dir, "simulation_settings.csv"), row.names = FALSE)

  invisible(list(
    truth = truth,
    barcode = barcode,
    frequency = F,
    depth = D,
    alt_count = A,
    amplicon_prior = amp_prior,
    settings = settings
  ))
}


