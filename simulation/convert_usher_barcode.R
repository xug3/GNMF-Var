# Convert a wide UShER barcode matrix to the format required by
# the SWAMPy-inspired GNMF-Var/NMF simulator.

# ---------------- USER SETTINGS ----------------
USHER_FILE <- "D:/RAP UNLV/research projects/waste water/clustering mutants/data/Usher_barcode_variant_of_interest_10182023.csv"
OUTPUT_FILE <- "D:/RAP UNLV/research projects/waste water/clustering mutants/simulation/inputs/barcode_matrix.csv"

# Use three real lineages present in the supplied file.
# B.1.617.2 is Delta; BA.1 and BA.2 are Omicron sublineages.
USHER_LINEAGES <- c("B.1.617.2", "BA.1", "BA.2")
OUTPUT_NAMES <- c("Delta", "BA.1", "BA.2")

# The UShER file does not contain amplicon coordinates. The simulator needs
# an amplicon label only to create correlated depth/amplicon variation.
# Here mutations are assigned to 98 equal-width pseudo-amplicon bins.
# Replace this mapping with assay-specific coordinates if available.
GENOME_LENGTH <- 29903L
N_AMPLICONS <- 98L
# ------------------------------------------------

usher <- read.csv(USHER_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# -------------------------------------------------------------------------
# 1. First column is lineage ID regardless of its column name.
#    This directly matches the structure shown in your screenshot.
# -------------------------------------------------------------------------
lineage_id <- trimws(as.character(usher[[1]]))
mutation_df <- usher[, -1, drop = FALSE]

if (anyDuplicated(lineage_id)) {
  stop("Duplicate lineage names are present in the first column.")
}

rownames(mutation_df) <- lineage_id

# Check requested lineages exist.
missing_lineages <- setdiff(USHER_LINEAGES, rownames(mutation_df))
if (length(missing_lineages) > 0) {
  stop(
    "Missing lineage(s): ",
    paste(missing_lineages, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 2. Select the three requested lineage rows.
#    Result: lineage x mutation matrix.
# -------------------------------------------------------------------------
B <- as.matrix(mutation_df[USHER_LINEAGES, , drop = FALSE])
storage.mode(B) <- "numeric"

# Treat any nonzero barcode entry as presence = 1.
B <- 1L * (B != 0)
rownames(B) <- OUTPUT_NAMES

# Remove mutation columns absent from all three selected lineages.
keep <- colSums(B) > 0
B <- B[, keep, drop = FALSE]

mutation_id <- colnames(B)

if (length(mutation_id) == 0) {
  stop("No nonzero barcode mutations found for the selected lineages.")
}

# -------------------------------------------------------------------------
# 3. Extract genomic position from mutation name.
#    Examples:
#      A23403G -> 23403
#      C241T   -> 241
#    The regex extracts the first run of digits, which is more robust than
#    requiring an exact single-nucleotide-substitution naming pattern.
# -------------------------------------------------------------------------
position_text <- sub("^[^0-9]*([0-9]+).*$", "\\1", mutation_id)
position <- suppressWarnings(as.integer(position_text))

if (anyNA(position)) {
  bad <- mutation_id[is.na(position)]
  stop(
    "Could not parse genomic positions for mutation(s): ",
    paste(head(bad, 20), collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 4. Create pseudo-amplicon labels.
# -------------------------------------------------------------------------
amplicon <- pmin(
  N_AMPLICONS,
  floor((position - 1) / GENOME_LENGTH * N_AMPLICONS) + 1L
)
amplicon <- as.integer(amplicon)

# -------------------------------------------------------------------------
# 5. Transpose to mutation x lineage, which is what the simulator expects.
# -------------------------------------------------------------------------
B_out <- as.data.frame(t(B), check.names = FALSE)
colnames(B_out) <- OUTPUT_NAMES

# Describe whether each mutation is lineage-specific or shared.
n_support <- rowSums(B_out != 0)

source <- vapply(
  seq_len(nrow(B_out)),
  function(i) {
    present <- OUTPUT_NAMES[as.numeric(B_out[i, ]) != 0]
    if (length(present) == 1) {
      paste0("specific_", present)
    } else {
      paste0("shared_", paste(present, collapse = "_"))
    }
  },
  FUN.VALUE = character(1)
)

barcode_out <- data.frame(
  mutation_id = mutation_id,
  position = position,
  amplicon = amplicon,
  source = source,
  B_out,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Sort by genomic position for readability.
barcode_out <- barcode_out[order(barcode_out$position, barcode_out$mutation_id), ]
rownames(barcode_out) <- NULL

write.csv(
  barcode_out,
  OUTPUT_FILE,
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------------------------------------------------
# 6. Sanity checks / console summary.
# -------------------------------------------------------------------------
cat("\nWrote:", OUTPUT_FILE, "\n\n")
cat("Selected lineage barcode counts:\n")
print(
  data.frame(
    input_lineage = USHER_LINEAGES,
    output_name = OUTPUT_NAMES,
    n_mutations = rowSums(B),
    row.names = NULL
  )
)

cat("\nUnion mutation count:", nrow(barcode_out), "\n")
cat("Expected for your uploaded file: 88 union mutations.\n")

cat("\nFirst few output rows:\n")
print(head(barcode_out, 10))
