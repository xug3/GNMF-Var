
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
)

# ============================= USER SETTINGS =============================

RESULT_DIR <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/simulation_results_graph_only"
REP_FILE   <- file.path(RESULT_DIR, "replicate_metrics.csv")

FIG_DIR <- file.path(RESULT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

METHODS_KEEP <- c("GNMF-Var", "NMF")

# ========================== FIGURE 1 =====================================

rep_metrics <- read.csv(REP_FILE, check.names = FALSE)

# Keep all 100 replicate-level values; do NOT summarize before boxplot
plot1_df <- rep_metrics |>
  filter(method %in% METHODS_KEEP) |>
  select(replicate, method, MAE, RMSE) |>
  pivot_longer(
    cols = c(MAE, RMSE),
    names_to = "metric",
    values_to = "value"
  )

plot1_df$method <- factor(plot1_df$method, levels = METHODS_KEEP)
plot1_df$metric <- factor(plot1_df$metric, levels = c("MAE", "RMSE"))

pal_method <- ggsci::pal_d3("category20")(2)
names(pal_method) <- METHODS_KEEP


y_top <- max(tapply(plot1_df$value,
                    interaction(plot1_df$metric, plot1_df$method),
                    function(v) boxplot.stats(v)$stats[5]))

p1 <- ggplot(
  plot1_df,
  aes(x = metric, y = value, fill = method)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    outlier.shape = NA,          # 改：原来是 outlier.size = 1.5
    linewidth = 0.4
  ) +
  scale_fill_manual(values = pal_method) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(ylim = c(0, y_top * 1.05)) +   # 新增：y 轴收到箱体范围
  labs(
    x = NULL,
    y = "Estimation error",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(6, 10, 6, 6)
  )

print(p1)


tiff(
  file.path(FIG_DIR, "Figure1_MAE_RMSE_boxplot.tiff"),
  width = 7.5,
  height = 5.5,
  units = "in",
  res = 300,
  compression = "lzw"
)
print(p1)
dev.off()


# ========================== FIGURE 2 =====================================
# Ground-truth lineage proportions across 73 time points
# Stacked area plot

INPUT_DIR <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/inputs"

TRUTH_FILE <- file.path(INPUT_DIR, "truth_proportions.csv")

truth <- read.csv(
  TRUTH_FILE,
  check.names = FALSE
)

# Three simulated lineages
variants <- c("Delta", "BA.1", "BA.2")

# Convert wide to long format
plot2_df <- truth |>
  select(time_point, all_of(variants)) |>
  pivot_longer(
    cols = all_of(variants),
    names_to = "lineage",
    values_to = "prevalence"
  ) |>
  mutate(
    lineage = factor(lineage, levels = variants)
  ) |>
  arrange(time_point, lineage)

# Colors consistent with previous figures
pal_variant <- ggsci::pal_d3("category20")(length(variants))
names(pal_variant) <- variants

# Plot
p2 <- ggplot(
  plot2_df,
  aes(
    x = time_point,
    y = prevalence,
    fill = lineage
  )
) +
  geom_area(
    position = "stack",
    colour = NA,
    linewidth = 0
  ) +
  scale_fill_manual(
    values = pal_variant,
    breaks = variants
  ) +
  scale_x_continuous(
    breaks = seq(0, 73, by = 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "Time point",
    y = "Ground-truth variant proportion",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )

print(p2)

# Save TIFF
tiff(
  file.path(FIG_DIR, "Figure2_ground_truth_proportions.tiff"),
  width = 9,
  height = 5,
  units = "in",
  res = 300,
  compression = "lzw"
)

print(p2)
dev.off()

# Figure: First-positive time comparison



first_df <- read.csv(
  file.path(RESULT_DIR, "first_positive_replicates.csv")
)

# Check data
stopifnot(
  all(c("replicate", "variant", "method",
        "truth_time", "estimated_time") %in% names(first_df))
)

if (anyNA(first_df$truth_time)) {
  stop("Some ground-truth first-positive times are missing.")
}

if (anyDuplicated(first_df[
  c("replicate", "variant", "method")
])) {
  stop("Duplicate replicate-variant-method records.")
}

# Calculate first-positive accuracy
first_df <- first_df |>
  mutate(
    time_error = estimated_time - truth_time,
    detected = !is.na(estimated_time),
    exact = detected & abs(time_error) == 0,
    within1 = detected & abs(time_error) <= 1,
    within2 = detected & abs(time_error) <= 2,
    early = detected & time_error < 0,
    late = detected & time_error > 0
  )

# Summary for each variant and method
accuracy_table <- first_df |>
  group_by(variant, method) |>
  summarise(
    n_replicates = n(),
    truth_time = first(truth_time),
    
    exact_n = sum(exact),
    exact_accuracy = mean(exact),
    within1_accuracy = mean(within1),
    within2_accuracy = mean(within2),
    
    detection_rate = mean(detected),
    no_detection_rate = mean(!detected),
    
    early_rate = mean(early),
    late_rate = mean(late),
    
    first_time_MAE = if (any(detected)) {
      mean(abs(time_error[detected]))
    } else NA_real_,
    
    median_estimated_time = if (any(detected)) {
      median(estimated_time[detected])
    } else NA_real_,
    
    .groups = "drop"
  )

# 95% binomial confidence interval for exact accuracy
accuracy_table$exact_CI_lower <- vapply(
  seq_len(nrow(accuracy_table)),
  function(i) {
    binom.test(
      accuracy_table$exact_n[i],
      accuracy_table$n_replicates[i]
    )$conf.int[1]
  },
  numeric(1)
)

accuracy_table$exact_CI_upper <- vapply(
  seq_len(nrow(accuracy_table)),
  function(i) {
    binom.test(
      accuracy_table$exact_n[i],
      accuracy_table$n_replicates[i]
    )$conf.int[2]
  },
  numeric(1)
)

accuracy_table <- accuracy_table |>
  mutate(
    variant = factor(
      variant, levels = c("Delta", "BA.1", "BA.2")
    ),
    method = factor(
      method, levels = c("GNMF-Var", "NMF")
    )
  ) |>
  arrange(variant, method)

write.csv(
  accuracy_table,
  file.path(RESULT_DIR, "first_positive_accuracy_table.csv"),
  row.names = FALSE
)

# Publication-friendly table
pct <- function(x) sprintf("%.1f%%", 100 * x)

publication_table <- accuracy_table |>
  transmute(
    Variant = as.character(variant),
    `Truth time` = truth_time,
    Method = as.character(method),
    `Exact accuracy` = pct(exact_accuracy),
    `Exact 95% CI` = sprintf(
      "[%.1f%%, %.1f%%]",
      100 * exact_CI_lower,
      100 * exact_CI_upper
    ),
    `Within +/-1` = pct(within1_accuracy),
    `Within +/-2` = pct(within2_accuracy),
    `No detection` = pct(no_detection_rate),
    `Time MAE` = round(first_time_MAE, 2),
    check.names = FALSE
  )

print(publication_table, row.names = FALSE)

write.csv(
  publication_table,
  file.path(RESULT_DIR, "Table_first_positive_accuracy.csv"),
  row.names = FALSE
)




first_df <- read.csv(
  file.path(RESULT_DIR, "first_positive_replicates.csv"),
  stringsAsFactors = FALSE
)


# ---------------- settings ----------------
VARIANTS   <- c("BA.1", "BA.2")
METHODS    <- c("GNMF-Var", "NMF")   # first one is drawn on top
MAX_BEFORE <- 5    # single time points shown before the true time; anything earlier
# goes into one "or earlier" column (use Inf to show every time point)
MAX_AFTER  <- 3    # same, after the true time ("or later" column)

# ---------------- prepare ----------------
n_m <- length(METHODS)

dat <- first_df |>
  filter(variant %in% VARIANTS, method %in% METHODS) |>
  mutate(
    variant = factor(variant, levels = VARIANTS),
    method  = factor(method,  levels = METHODS)
  )

tot <- dat |> count(variant, method, name = "n_total")   # includes non-detections

# columns to draw: from the first to the last non-zero time point near the true time
bounds_df <- dat |>
  filter(!is.na(estimated_time)) |>
  group_by(variant) |>
  summarise(
    truth_time = first(truth_time),
    lo = min(c(estimated_time[estimated_time >= truth_time - MAX_BEFORE], truth_time)),
    hi = max(c(estimated_time[estimated_time <= truth_time + MAX_AFTER], truth_time)),
    pool_lo = any(estimated_time < lo),   # an "or earlier" column is needed
    pool_hi = any(estimated_time > hi),   # an "or later" column is needed
    .groups = "drop"
  )

grid_df <- do.call(rbind, lapply(seq_len(nrow(bounds_df)), function(i) {
  b <- bounds_df[i, ]
  pos <- c(if (b$pool_lo) b$lo - 1, b$lo:b$hi, if (b$pool_hi) b$hi + 1)
  lab <- as.character(pos)
  if (b$pool_lo) lab[1] <- paste0(b$lo - 1, " or\nearlier")
  if (b$pool_hi) lab[length(lab)] <- paste0(b$hi + 1, " or\nlater")
  data.frame(variant = b$variant, pos = pos, x_lab = lab,
             is_truth = pos == b$truth_time)
}))
x_levels <- unique(grid_df$x_lab[order(grid_df$pos, grid_df$x_lab)])
grid_df$x_fac <- factor(grid_df$x_lab, levels = x_levels)
grid_df <- merge(grid_df, data.frame(method = factor(METHODS, levels = METHODS)))

count_df <- dat |>
  filter(!is.na(estimated_time)) |>
  left_join(bounds_df |> select(variant, lo, hi), by = "variant") |>
  mutate(pos = pmin(pmax(estimated_time, lo - 1), hi + 1)) |>   # pool the tails
  count(variant, method, pos, name = "n_rep")

tile_df <- grid_df |>
  left_join(count_df, by = c("variant", "method", "pos")) |>
  left_join(tot, by = c("variant", "method")) |>
  mutate(
    n_rep    = ifelse(is.na(n_rep), 0, n_rep),
    prop     = n_rep / n_total,
    label    = ifelse(n_rep == 0, "0%",
                      ifelse(prop < 0.01, "<1%", sprintf("%.0f%%", 100 * prop))),
    text_col = ifelse(prop > 0.6, "grey10", "white"),
    y_pos    = n_m + 1 - as.integer(method)      # first method on top
  )

truth_tile <- tile_df |>
  filter(is_truth) |>
  distinct(variant, x_fac)

# check: share of replicates shown in each row (100% unless some had no detection)
check_df <- tile_df |>
  group_by(variant, method) |>
  summarise(pct_shown = 100 * sum(prop), .groups = "drop") |>
  mutate(pct_no_detection = 100 - pct_shown)
print(as.data.frame(check_df), row.names = FALSE)

# ---------------- plot ----------------
p_heat <- ggplot(tile_df, aes(x = x_fac, y = y_pos, fill = prop)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label, colour = text_col), size = 3.6) +
  # true first-positive time: box (white under black), marker and label
  geom_tile(
    data = truth_tile,
    aes(x = x_fac, y = (n_m + 1) / 2),
    inherit.aes = FALSE, height = n_m, width = 1,
    fill = NA, colour = "white", linewidth = 2.4
  ) +
  geom_tile(
    data = truth_tile,
    aes(x = x_fac, y = (n_m + 1) / 2),
    inherit.aes = FALSE, height = n_m, width = 1,
    fill = NA, colour = "black", linewidth = 1.0
  ) +
  geom_point(
    data = truth_tile,
    aes(x = x_fac, y = n_m + 0.68),
    inherit.aes = FALSE, shape = 25, size = 2.6, fill = "black", colour = "black"
  ) +
  geom_text(
    data = truth_tile,
    aes(x = x_fac, y = n_m + 1.0, label = "True time"),
    inherit.aes = FALSE, size = 3.6, fontface = "bold"
  ) +
  facet_grid(~ variant, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
    values  = c(0, 0.25, 0.5, 0.75, 1),
    limits  = c(0, 1), guide = "none"
  ) +
  scale_colour_identity() +
  scale_y_continuous(breaks = seq_len(n_m), labels = rev(METHODS)) +
  coord_cartesian(ylim = c(0.5, n_m + 1.25), expand = FALSE, clip = "off") +
  labs(
    x = "Estimated first-positive time point",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text       = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    axis.line        = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_text(size = 12, colour = "black"),
    axis.text.x      = element_text(lineheight = 0.9),
    panel.spacing    = unit(1.5, "lines"),
    plot.margin      = margin(8, 12, 8, 8)
  )

print(p_heat)

ggsave(
  file.path(FIG_DIR, "Figure_first_positive_heatmap.tiff"),
  p_heat, width = 9, height = 3, dpi = 300, compression = "lzw"
)



# ---------------- numbers for the text ----------------
summ <- dat |>
  group_by(variant, method) |>
  summarise(
    n_replicates  = n(),
    truth_time    = first(truth_time),
    median_est    = median(estimated_time, na.rm = TRUE),
    exact_pct     = 100 * sum(time_difference ==  0, na.rm = TRUE) / n(),
    early_pct     = 100 * sum(time_difference <   0, na.rm = TRUE) / n(),
    late_pct      = 100 * sum(time_difference >   0, na.rm = TRUE) / n(),
    early_gt1_pct = 100 * sum(time_difference <  -1, na.rm = TRUE) / n(),
    nodetect_pct  = 100 * sum(is.na(estimated_time)) / n(),
    median_signed = median(time_difference, na.rm = TRUE),
    time_MAE      = mean(abs(time_difference), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

print(as.data.frame(summ), row.names = FALSE)

write.csv(
  summ,
  file.path(RESULT_DIR, "Table_first_positive_distribution.csv"),
  row.names = FALSE
)


# Read CSV
summary_table <- read.csv(
  file.path(RESULT_DIR, "summary_metrics.csv")
)

# Round numerical results to 3 decimal places
cols <- c("mean", "sd", "median", "q025", "q975")

summary_table[cols] <- lapply(
  summary_table[cols],
  round,
  digits = 3
)

# Save
write.csv(
  summary_table,
  file.path(RESULT_DIR, "summary_metrics_3decimals.csv"),
  row.names = FALSE
)

