#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(piano))

usage <- paste(
  "Usage:",
  "gxa_run_piano_gsea.R --tsv <analytics.tsv> --pvalue-col <n> --foldchange-col <n>",
  "                    --go <gene-set.tsv> --descr <term-to-accession.tsv>",
  "                    --out <output-prefix> [--pvalue 0.05] [--gs-fdr 0.1]",
  "                    [--minsize 5] [--maxsize 100]",
  sep = "\n"
)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[name]] <- TRUE
      i <- i + 1
    } else {
      parsed[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  parsed
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name, "\n", usage, call. = FALSE)
  }
  value
}

as_numeric_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value)) {
    return(default)
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (is.na(numeric_value)) {
    stop("Argument --", name, " must be numeric", call. = FALSE)
  }
  numeric_value
}

read_tsv <- function(path, header = TRUE) {
  tryCatch(
    read.table(
      path,
      sep = "\t",
      header = header,
      quote = "",
      comment.char = "",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      stop("Failed to read ", path, ": ", conditionMessage(e), call. = FALSE)
    }
  )
}

write_empty_outputs <- function(output_prefix) {
  write.table(
    data.frame(
      Name = character(),
      "p-value" = numeric(),
      "p adj (non-dir.)" = numeric(),
      "Significant (in gene set)" = integer(),
      "Genes (tot)" = integer(),
      "effect.size" = numeric(),
      check.names = FALSE
    ),
    file = paste0(output_prefix, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  writeLines("gene", paste0(output_prefix, "_list.tsv"))
}

choose_gene_set_columns <- function(gene_sets, gene_ids, term_to_accession, accession_to_term) {
  if (ncol(gene_sets) < 2) {
    stop("Gene set file must have at least two tab-separated columns", call. = FALSE)
  }

  candidates <- list(c(gene = 1, term = 2), c(gene = 2, term = 1))
  scores <- vapply(
    candidates,
    function(candidate) {
      gene_values <- as.character(gene_sets[[candidate[["gene"]]]])
      term_values <- as.character(gene_sets[[candidate[["term"]]]])
      sum(gene_values %in% gene_ids, na.rm = TRUE) +
        sum(term_values %in% names(term_to_accession), na.rm = TRUE) +
        sum(term_values %in% names(accession_to_term), na.rm = TRUE)
    },
    numeric(1)
  )

  candidates[[which.max(scores)]]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

analytics_file <- required_arg(args, "tsv")
gene_set_file <- required_arg(args, "go")
description_file <- required_arg(args, "descr")
output_prefix <- required_arg(args, "out")
pvalue_col <- as.integer(required_arg(args, "pvalue-col"))
foldchange_col <- as.integer(required_arg(args, "foldchange-col"))
pvalue_cutoff <- as_numeric_arg(args, "pvalue", 0.05)
gsea_fdr <- as_numeric_arg(args, "gs-fdr", 0.1)
min_size <- as_numeric_arg(args, "minsize", 5)
max_size <- as_numeric_arg(args, "maxsize", 100)

analytics <- read_tsv(analytics_file, header = TRUE)
if (pvalue_col < 1 || pvalue_col > ncol(analytics)) {
  stop("--pvalue-col is outside the analytics table", call. = FALSE)
}
if (foldchange_col < 1 || foldchange_col > ncol(analytics)) {
  stop("--foldchange-col is outside the analytics table", call. = FALSE)
}

gene_column <- names(analytics)[[1]]
gene_ids <- as.character(analytics[[gene_column]])
pvalues <- suppressWarnings(as.numeric(analytics[[pvalue_col]]))
foldchanges <- suppressWarnings(as.numeric(analytics[[foldchange_col]]))

valid <- nzchar(gene_ids) & !is.na(gene_ids) & is.finite(pvalues)
analytics <- data.frame(
  gene = gene_ids[valid],
  pvalue = pvalues[valid],
  foldchange = foldchanges[valid],
  stringsAsFactors = FALSE
)

if (nrow(analytics) == 0) {
  write_empty_outputs(output_prefix)
  quit(save = "no", status = 0)
}

analytics$abs_foldchange <- abs(analytics$foldchange)
analytics$abs_foldchange[!is.finite(analytics$abs_foldchange)] <- -Inf
analytics <- analytics[order(analytics$gene, analytics$pvalue, -analytics$abs_foldchange), ]
analytics <- analytics[!duplicated(analytics$gene), ]

if (!any(analytics$pvalue <= pvalue_cutoff, na.rm = TRUE)) {
  write_empty_outputs(output_prefix)
  quit(save = "no", status = 0)
}

mapping <- read_tsv(description_file, header = FALSE)
if (ncol(mapping) < 2) {
  stop("Description file must have at least two tab-separated columns", call. = FALSE)
}
mapping <- mapping[, 1:2, drop = FALSE]
mapping[[1]] <- as.character(mapping[[1]])
mapping[[2]] <- as.character(mapping[[2]])
mapping <- mapping[nzchar(mapping[[1]]) & nzchar(mapping[[2]]), , drop = FALSE]
term_to_accession <- setNames(mapping[[2]], mapping[[1]])
accession_to_term <- setNames(mapping[[1]], mapping[[2]])

gene_sets <- read_tsv(gene_set_file, header = FALSE)
columns <- choose_gene_set_columns(gene_sets, analytics$gene, term_to_accession, accession_to_term)
gene_set_data <- data.frame(
  gene = as.character(gene_sets[[columns[["gene"]]]]),
  term = as.character(gene_sets[[columns[["term"]]]]),
  stringsAsFactors = FALSE
)
gene_set_data <- gene_set_data[nzchar(gene_set_data$gene) & nzchar(gene_set_data$term), , drop = FALSE]
gene_set_data$term <- ifelse(
  gene_set_data$term %in% names(accession_to_term),
  accession_to_term[gene_set_data$term],
  gene_set_data$term
)
gene_set_data <- unique(gene_set_data)

if (nrow(gene_set_data) == 0) {
  write_empty_outputs(output_prefix)
  quit(save = "no", status = 0)
}

gsc <- loadGSC(gene_set_data, type = "data.frame")
gsa <- tryCatch(
  runGSAhyper(
    genes = analytics$gene,
    pvalues = analytics$pvalue,
    pcutoff = pvalue_cutoff + .Machine$double.eps,
    gsc = gsc,
    gsSizeLim = c(min_size, max_size),
    adjMethod = "fdr"
  ),
  error = function(e) {
    message("No GSEA result: ", conditionMessage(e))
    NULL
  }
)

if (is.null(gsa) || is.null(gsa$resTab) || nrow(gsa$resTab) == 0) {
  write_empty_outputs(output_prefix)
  quit(save = "no", status = 0)
}

result <- as.data.frame(gsa$resTab, stringsAsFactors = FALSE)
result$Name <- rownames(result)
result[] <- lapply(result, function(column) {
  suppressWarnings(type.convert(column, as.is = TRUE))
})

result$genes_total <- result[["Significant (in gene set)"]] + result[["Non-significant (in gene set)"]]
significant_total <- sum(analytics$pvalue <= pvalue_cutoff, na.rm = TRUE)
expected <- result$genes_total * significant_total / nrow(analytics)
result$effect_size <- result[["Significant (in gene set)"]] / expected
result$effect_size[!is.finite(result$effect_size)] <- NA_real_

result <- result[result[["Adjusted p-value"]] <= gsea_fdr & is.finite(result$effect_size), , drop = FALSE]
result <- result[order(result[["Adjusted p-value"]], -result$effect_size, result$Name), , drop = FALSE]

out <- data.frame(
  Name = result$Name,
  "p-value" = result[["p-value"]],
  "p adj (non-dir.)" = result[["Adjusted p-value"]],
  "Significant (in gene set)" = result[["Significant (in gene set)"]],
  "Genes (tot)" = result$genes_total,
  "effect.size" = round(result$effect_size, 4),
  check.names = FALSE
)

write.table(
  out,
  file = paste0(output_prefix, ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(result) == 0) {
  writeLines("gene", paste0(output_prefix, "_list.tsv"))
} else {
  result_terms <- result$Name
  significant_genes <- analytics$gene[analytics$pvalue <= pvalue_cutoff]
  genes_in_significant_sets <- unique(gene_set_data$gene[
    gene_set_data$term %in% result_terms & gene_set_data$gene %in% significant_genes
  ])
  writeLines(c("gene", genes_in_significant_sets), paste0(output_prefix, "_list.tsv"))
}
