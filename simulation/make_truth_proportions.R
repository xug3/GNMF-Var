PROJECT_DIR <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation"
OUTPUT_FILE <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/inputs/truth_proportions.csv"

make_73_timepoint_truth <- function(
    n_time = 73,
    ba1_mid = 34,
    ba1_scale = 5,
    ba2_mid = 49,
    ba2_scale = 3) {

  t <- seq_len(n_time)

  # Smooth replacement of Delta by BA.1.
  ba1_takeover <- plogis((t - ba1_mid) / ba1_scale)

  # Smooth replacement of the pre-BA.2 population by BA.2.
  ba2_takeover <- plogis((t - ba2_mid) / ba2_scale)

  p_ba2 <- ba2_takeover
  p_ba1 <- (1 - p_ba2) * ba1_takeover
  p_delta <- (1 - p_ba2) * (1 - ba1_takeover)

  P <- cbind(
    Delta = p_delta,
    `BA.1` = p_ba1,
    `BA.2` = p_ba2
  )

  # Numerical safeguard: force every row to sum exactly to 1.
  P <- P / rowSums(P)

  data.frame(
    sample_id = sprintf("T%03d", t),
    time_point = t,
    P,
    check.names = FALSE
  )
}

truth <- make_73_timepoint_truth()
write.csv(truth, OUTPUT_FILE, row.names = FALSE)

cat("Wrote:", OUTPUT_FILE, "\n")
print(truth[truth$time_point %in% c(1, 20, 34, 40, 49, 55, 73), ])


