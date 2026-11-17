# Data loading helpers ---------------------------------------------------------

#' Read the cohort configuration table (cohort, label, study_id).
read_cohort_config <- function(path = "config/cohorts.csv") {
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Load one cohort: a numeric gene x patient matrix and its survival table.
load_cohort <- function(cohort, data_dir = "data", top_n = 15000) {
  expr_file <- file.path(data_dir, paste0(cohort, "_expression.tsv.gz"))
  surv_file <- file.path(data_dir, paste0(cohort, "_survival.tsv"))
  if (!file.exists(expr_file) || !file.exists(surv_file)) {
    stop("Missing data for ", cohort, ". Run `python data-raw/fetch_cohorts.py` first.", call. = FALSE)
  }
  df <- data.table::fread(expr_file, data.table = FALSE)
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df[[1]]
  surv <- utils::read.delim(surv_file, stringsAsFactors = FALSE)
  list(matrix = mat, survival = surv, top_genes = top_expressed(mat, top_n))
}

#' Names of the `n` genes with the highest median expression across patients.
top_expressed <- function(mat, n = 15000) {
  med <- apply(mat, 1, stats::median, na.rm = TRUE)
  names(sort(med, decreasing = TRUE))[seq_len(min(n, length(med)))]
}

#' Gene sets (MSigDB Hallmark, KEGG legacy, Reactome) as TERM2GENE data frames.
#' Built once with msigdbr and cached to `cache_file`.
load_genesets <- function(cache_file = "data/genesets.rds") {
  if (file.exists(cache_file)) return(readRDS(cache_file))
  pull <- function(collection, subcollection = NULL) {
    d <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    unique(data.frame(term = d$gs_name, gene = d$gene_symbol, stringsAsFactors = FALSE))
  }
  sets <- list(
    Hallmark = pull("H"),
    Reactome = pull("C2", "CP:REACTOME"),
    KEGG     = pull("C2", "CP:KEGG_LEGACY")
  )
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(sets, cache_file)
  sets
}
