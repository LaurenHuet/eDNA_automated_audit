#!/usr/bin/env Rscript
# check_phyloseq.R
# OceanOmics eDNA Output Validator — phyloseq .rds checker
# Called by validate_edna_output.py via singularity
#
# Usage: Rscript check_phyloseq.R <rds_path> <output_json>
#
# Output: JSON array of {status, message} objects written to <output_json>
# Also outputs sample names as a special record so Python can compare with FAIRe xlsx

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: check_phyloseq.R <rds_path> <output_json>")
}

rds_path   <- args[1]
output_json <- args[2]

results <- list()

add <- function(status, message) {
  results[[length(results) + 1]] <<- list(status = status, message = message)
}

# ── load ──────────────────────────────────────────────────────────────────────

tryCatch({
  ps <- readRDS(rds_path)
  add("PASS", "phyloseq .rds loaded successfully")
}, error = function(e) {
  add("FAIL", paste("Failed to load .rds:", conditionMessage(e)))
  writeLines(toJSON(results, auto_unbox = TRUE, pretty = TRUE), output_json)
  quit(status = 1)
})

# ── required tables present ───────────────────────────────────────────────────

has_otu    <- !is.null(tryCatch(otu_table(ps),    error = function(e) NULL))
has_sam    <- !is.null(tryCatch(sample_data(ps),  error = function(e) NULL))
has_tax    <- !is.null(tryCatch(tax_table(ps),    error = function(e) NULL))

if (has_otu) add("PASS", "otu_table present") else add("FAIL", "otu_table is missing")
if (has_sam) add("PASS", "sample_data present") else add("FAIL", "sample_data (sam_data) is missing")
if (has_tax) add("PASS", "tax_table present") else add("FAIL", "tax_table is missing")

# ── otu_table: sample columns unique ─────────────────────────────────────────

if (has_otu) {
  otu <- otu_table(ps)
  samp_names_otu <- if (taxa_are_rows(ps)) colnames(otu) else rownames(otu)

  if (any(duplicated(samp_names_otu))) {
    dupes <- unique(samp_names_otu[duplicated(samp_names_otu)])
    add("FAIL", paste("Duplicate sample names in otu_table:", paste(head(dupes, 5), collapse = ", ")))
  } else {
    add("PASS", "otu_table sample names are unique")
  }
}

# ── sam_data: emit sample names for Python to compare with FAIRe ─────────────

if (has_sam) {
  ps_samples <- sample_names(ps)
  ps_cols    <- colnames(sample_data(ps))
  add("SAMPLE_NAMES", paste(ps_samples, collapse = ","))
  add("SAM_DATA_COLS", paste(ps_cols,    collapse = ","))
}

# ── write output ──────────────────────────────────────────────────────────────

writeLines(toJSON(results, auto_unbox = TRUE, pretty = TRUE), output_json)
