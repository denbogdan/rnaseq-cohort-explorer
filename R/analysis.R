# Pure analysis functions (no Shiny dependencies, so they can be unit-tested) ---

#' Stop with a readable message if genes are not in the matrix.
check_genes <- function(mat, genes) {
  missing <- setdiff(genes, rownames(mat))
  if (length(missing) > 0) {
    stop("Gene symbol not found in this cohort: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Spearman correlation between two genes across patients.
correlate_genes <- function(mat, gene1, gene2, log_transform = TRUE) {
  check_genes(mat, c(gene1, gene2))
  d <- data.frame(x = mat[gene1, ], y = mat[gene2, ])
  if (log_transform) d <- log(d + 1)
  ct <- suppressWarnings(stats::cor.test(d$x, d$y, method = "spearman"))
  list(data = stats::setNames(d, c(gene1, gene2)), rho = unname(ct$estimate), p = ct$p.value)
}

format_p <- function(p) if (p < 0.001) "p < 0.001" else paste("p =", round(p, 3))

#' Rank genes by Spearman correlation with `gene` (decreasing); `gene` itself is dropped.
rank_by_correlation <- function(mat, gene, top_genes) {
  check_genes(mat, gene)
  genes <- setdiff(intersect(top_genes, rownames(mat)), gene)
  rho <- stats::cor(t(mat[genes, , drop = FALSE]), mat[gene, ], method = "spearman")[, 1]
  sort(rho[!is.na(rho)], decreasing = TRUE)
}

#' Pre-ranked GSEA of a correlation-ranked gene list. Returns a tidy data frame.
run_gsea <- function(ranked, term2gene, pval = 0.05) {
  res <- clusterProfiler::GSEA(ranked, TERM2GENE = term2gene, verbose = FALSE,
                               pvalueCutoff = pval, minGSSize = 5)
  out <- as.data.frame(res)
  if (nrow(out) == 0) {
    return(data.frame(Pathway = character(), `Set size` = integer(), NES = numeric(),
                      `p-value` = numeric(), FDR = numeric(), check.names = FALSE))
  }
  out <- out[order(out$NES, decreasing = TRUE), ]
  data.frame(Pathway = gsub("_", " ", out$ID), `Set size` = out$setSize,
             NES = signif(out$NES, 3), `p-value` = signif(out$pvalue, 3),
             FDR = signif(out$p.adjust, 3), check.names = FALSE, row.names = NULL)
}

#' Stratify patients into high/low expression groups for survival analysis.
#' Returns list(data, cutpoint) where `data` has os_months, death, group.
stratify_survival <- function(mat, surv, gene, maxstat = TRUE) {
  check_genes(mat, gene)
  d <- merge(surv, data.frame(patient_id = colnames(mat), expression = mat[gene, ]), by = "patient_id")
  d <- d[stats::complete.cases(d), ]
  if (maxstat) {
    rp <- survminer::surv_cutpoint(d, time = "os_months", event = "death", variables = "expression")
    cat <- survminer::surv_categorize(rp)
    d$group <- cat$expression
    list(data = d, cutpoint = rp)
  } else {
    d$group <- ifelse(d$expression > stats::median(d$expression), "high", "low")
    list(data = d, cutpoint = NULL)
  }
}
