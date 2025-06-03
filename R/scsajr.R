#!/usr/bin/env Rscript

# Need to add package description
# Need to figureout how to keep links in the Roxygen documentation

DELIMITER <- "$"


#' Filter segments and samples based on coverage and group criteria
#'
#' This function takes a `SummarizedExperiment` containing alternative splicing data
#' (with assays “i” for inclusion counts and “e” for exclusion counts), and filters out:
#'   1. pseudobulk samples with fewer than `sample_min_ncells` cells,
#'   2. cell‐type groups with fewer than `celltype_min_samples` samples,
#'   3. segments (rows) that do not meet minimum coverage across pseudobulks,
#'   4. segments with low variability (standard deviation of PSI).
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param sample_filter Logical vector (length = number of columns in `pbas`) indicating which samples to retain initially.
#'   Defaults to `TRUE` for all samples.
#' @param sites Character vector of splice‐site patterns to keep (e.g., `c("ad", "aa", "dd")`). Only segments whose `rowData(pbas)$sites` is in this set will be retained.
#' @param min_cov Minimum total junction coverage (`i + e`) for a pseudobulk to count toward segment inclusion. Default 10.
#' @param sample_min_ncells Minimum number of cells per pseudobulk; pseudobulk samples with fewer cells are filtered out. Default 20.
#' @param celltype_min_samples Minimum number of pseudobulk samples per group (as defined by `groupby`) to keep that group. Default 2.
#' @param seg_min_samples The absolute minimum number of pseudobulks (after filtering by `min_cov`) that a segment must be observed in. Default 4.
#' @param seg_min_samples_fraq Minimum fraction of total pseudobulk samples that must meet `min_cov` for a segment to be retained. Default 0.2.
#'   The effective minimum samples per segment is `max(seg_min_samples, ncol(pbas) * seg_min_samples_fraq)`.
#' @param seg_min_sd Minimum standard deviation of PSI (percent spliced‐in) across pseudobulks for a segment to be retained. Default 0.1.
#' @param groupby Either (a) a character vector of length `ncol(pbas)` giving an existing grouping factor for each column, or
#'   (b) a character scalar or vector of column names in `colData(pbas)` whose pasted‐together values define the grouping factor.
#'   See `\link{get_groupby_factor}` for details.
#'
#' @return A filtered `SummarizedExperiment`, containing only those pseudobulks and segments that satisfy all criteria.
#'   The returned object will always have assays “i”, “e” converted to dense matrices, and will include a “psi” assay (percent spliced‐in) for each retained segment/samples.
#'   The `rowData` will also contain a column “nna” (number of pseudobulks with `i+e >= min_cov`) and “sd” (standard deviation of PSI).
#'
#' @examples
#' \dontrun{
#' data(pbas)
#' filtered_se <- filter_segments_and_samples(
#'   pbas = pbas,
#'   sample_filter = TRUE,
#'   sites = c("ad", "aa", "dd"),
#'   min_cov = 10,
#'   sample_min_ncells = 20,
#'   celltype_min_samples = 2,
#'   seg_min_samples = 4,
#'   seg_min_samples_fraq = 0.2,
#'   seg_min_sd = 0.1,
#'   groupby = "celltype"
#' )
#' head(filtered_se)
#' }
#' @export
filter_segments_and_samples <- function(
    pbas,
    sample_filter = TRUE,
    sites = c("ad", "aa", "dd"),
    min_cov = 10,
    sample_min_ncells = 20,
    celltype_min_samples = 2,
    seg_min_samples = 4,
    seg_min_samples_fraq = 0.2,
    seg_min_sd = 0.1,
    groupby = "celltype") {
  ## 1. Filter out pseudobulk samples with fewer than sample_min_ncells cells
  # Filter pseudobulk samples with fewer cells
  sample_filter <- sample_filter & (SummarizedExperiment::colData(pbas)$ncells >= sample_min_ncells)

  # Determine how many samples remain in each group (celltype)
  group_factor <- get_groupby_factor(pbas, groupby) # returns a factor vector (length = ncol)
  ct_counts <- table(group_factor[sample_filter]) # Count samples per group (only among those passing sample_filter)

  # Only keep samples whose group has at least 'celltype_min_samples' pseudobulks
  sample_filter <- sample_filter & (group_factor %in% names(ct_counts)[ct_counts >= celltype_min_samples])

  pbas <- pbas[, sample_filter] # Subset the SummarizedExperiment columns (samples)


  ## 2. Filter out cell‐type groups with fewer than celltype_min_samples samples
  # Keep only rows where rowData(pbas)$sites is in the specified 'sites' vector
  seg_sel <- SummarizedExperiment::rowData(pbas)$sites %in% sites
  pbas <- pbas[seg_sel, ]

  # Compute, for each segment, how many pseudobulks have (i + e) >= min_cov
  #   - assay(pbas, "i") and assay(pbas, "e") give inclusion/exclusion counts
  #   - Convert these to dense matrices to simplify row/column sums
  SummarizedExperiment::assay(pbas, "i") <- as.matrix(SummarizedExperiment::assay(pbas, "i"))
  SummarizedExperiment::assay(pbas, "e") <- as.matrix(SummarizedExperiment::assay(pbas, "e"))

  # 'nna' = number of pseudobulks with total coverage >= min_cov
  SummarizedExperiment::rowData(pbas)$nna <- rowSums((SummarizedExperiment::assay(pbas, "i") + SummarizedExperiment::assay(pbas, "e")) >= min_cov)


  ## 3. Filter out segments (rows) that do not meet minimum coverage across pseudobulks
  # Determine per‐segment minimum required pseudobulks:
  # If fractional rule (seg_min_samples_fraq) yields more required samples than seg_min_samples,
  # use the larger of the two.
  min_samples_required <- max(seg_min_samples, ceiling(ncol(pbas) * seg_min_samples_fraq))
  seg_sel2 <- SummarizedExperiment::rowData(pbas)$nna >= min_samples_required
  pbas <- pbas[seg_sel2, ]


  ## 4. Compute PSI (percent spliced‐in) and filter by segment variability
  # Calculate PSI = i / (i + e), setting PSI to NA whenever total coverage < min_cov
  psi_mat <- calc_psi(pbas, min_cov = min_cov)
  SummarizedExperiment::assay(pbas, "psi") <- psi_mat

  # Compute standard deviation of PSI across pseudobulks for each segment
  SummarizedExperiment::rowData(pbas)$sd <- apply(SummarizedExperiment::assay(pbas, "psi"), 1, stats::sd, na.rm = TRUE)

  # Keep only segments with standard deviation above threshold
  seg_sel3 <- SummarizedExperiment::rowData(pbas)$sd > seg_min_sd
  pbas <- pbas[seg_sel3, ]

  # Return the filtered SummarizedExperiment
  return(pbas)
}


#' Fit a quasi‐binomial GLM for alternative splicing per segment
#'
#' For each segment (row) in a `SummarizedExperiment` with assays “i” (inclusion counts) and “e” (exclusion counts),
#' this function fits a quasi‐binomial generalized linear model of the form `i / (i + e) ~ predictors`,
#' adding a pseudocount if desired. It can optionally return p‐values via a likelihood‐ratio test.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param formula An `\link[stats]{formula}` object specifying the model, e.g., `i ~ group`.
#'   Internally, the function treats “i” and “e” as a two‐column response.
#' @param data_terms A named list of covariates: data.frame columns referenced in `formula`.
#'   For example, if `formula = i ~ group`, then `data_terms = list(group = group_factor)` where `group_factor` is a vector of length `ncol(pbas)`.
#' @param pseudocount A numeric fraction of each count total to add to both “i” and “e”, expressed as a fraction of the total counts. Default 0 (no pseudocount).
#' @param parallel Logical; if `TRUE`, uses `\link[plyr]{alply}(…, .parallel = TRUE)` to fit models in parallel (requires registered backend). Default `FALSE`.
#' @param progress Character string passed to `.progress` in `\link[plyr]{alply}`. Default `"none"`.
#' @param return_pv Logical; if `TRUE`, perform a quasi‐likelihood ratio test per segment to return p‐values for each term in `formula`.
#'   Otherwise, return fitted `\link[stats]{glm}` objects. Default `FALSE`.
#' @param overdisp Logical; if `TRUE`, account for overdispersion when computing test statistics. Default `TRUE`.
#' @param disp_param Optional numeric vector of length `nrow(pbas)` providing dispersion parameters per segment.
#'   If `NULL`, dispersion is estimated from each model’s residuals.
#'
#' @return If `return_pv = FALSE`, a named list of length `nrow(pbas)` of fitted `glm` objects
#'   (or `NA` for segments where fitting failed). Names correspond to `rownames(pbas)`.
#'   If `return_pv = TRUE`, a data.frame with one row per segment, columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion (or provided `disp_param`).}
#'     \item{<term>}{Likelihood‐ratio p‐value for each term in `formula`.}
#'   }
#'   Row names correspond to `rownames(pbas)`.
#'
#' @details
#' The function loops over each segment, constructs a two‐column matrix `cbind(i, e)`, adds pseudocounts if requested, and fits a quasi‐binomial `\link[stats]{glm}`.
#' If `return_pv = TRUE`, it calls `\link{qbinom_lrt}` to compute a likelihood‐ratio test for each term in `formula`.
#' Any errors or warnings during fitting are caught; segments with fitting errors return `NA`.
#'
#' @seealso `\link{qbinom_lrt}`, `\link[stats]{glm}`, `\link[plyr]{alply}`
#' @export
fit_as_glm <- function(
    pbas,
    formula,
    data_terms,
    pseudocount = 0,
    parallel = FALSE,
    progress = "none",
    return_pv = FALSE,
    overdisp = TRUE,
    disp_param = NULL) {
  # Ensure plyr is available
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Please install the 'plyr' package to use fit_as_glm().")
  }

  # Extract term names from the formula
  term_labels <- attr(stats::terms(formula), "term.labels")

  # For each segment (row), fit a quasi-binomial GLM
  results_list <- plyr::alply(
    .data = seq_len(nrow(pbas)),
    .margins = 1,
    .fun = function(i) {
      # 1. Build the i/e matrix and add pseudocounts if specified
      inc_vec <- SummarizedExperiment::assay(pbas, "i")[i, ]
      exc_vec <- SummarizedExperiment::assay(pbas, "e")[i, ]
      tot <- inc_vec + exc_vec

      # Add pseudocount proportional to total counts for each sample
      if (pseudocount > 0) {
        adj <- tot * pseudocount
        inc_vec <- inc_vec + adj
        exc_vec <- exc_vec + adj
      }

      # Place adjusted counts into the data frame used by glm()
      #    'x' must be a two-column matrix: cbind(i, e)
      data_terms$x <- cbind(inc_vec, exc_vec)

      # 2. Fit the quasi-binomial GLM
      fit <- tryCatch(
        stats::glm(formula, data = data_terms, family = "quasibinomial"),
        error = function(e) {
          warning(paste0("Segment ", rownames(pbas)[i], ": ", e$message))
          return(NA)
        },
        warning = function(w) {
          warning(paste0("Segment ", rownames(pbas)[i], ": ", w$message))
          return(NA)
        }
      )

      # 3. If return_pv and fit succeeded, compute p-values via LRT
      if (return_pv && !is.na(fit)[1]) {
        fit <- qbinom_lrt(
          fit,
          overdisp = overdisp,
          disp_param = disp_param[i],
          term_labels = term_labels,
          seg_id = rownames(pbas)[i]
        )
      }

      fit
    },
    .parallel = parallel,
    .progress = progress
  )

  # If return_pv, combine p-value results into a data.frame
  if (return_pv) {
    pv_df <- do.call(rbind, results_list)
    rownames(pv_df) <- rownames(pbas)
    colnames(pv_df) <- c("overdispersion", term_labels)
    return(pv_df)
  }

  # Otherwise, return a named list of glm objects (or NA)
  names(results_list) <- rownames(pbas)
  attr(results_list, "term.labels") <- term_labels
  return(results_list)
}


#' Quasi‐likelihood ratio test (Chi-square ANOVA) for a fitted quasi‐binomial GLM
#'
#' Given a fitted `glm(..., family = "quasibinomial")`, this function computes a likelihood‐ratio test for each term in the model, accounting for overdispersion.
#'
#' @param fit_obj A fitted `glm` object (quasi‐binomial family).
#' @param overdisp Logical; if `TRUE`, use estimated dispersion from `fit_obj` (or `disp_param`) when computing test statistics.
#'   If `FALSE`, force dispersion = 1 (i.e., treat as simple binomial).
#' @param disp_param Optional numeric; if provided, use this as the dispersion parameter for the test.
#'   Otherwise, estimate dispersion from `fit_obj`’s residuals.
#' @param term_labels Character vector of term names (excluding the intercept) in the original formula, in the same order returned by `attr(terms(formula), "term.labels")`.
#' @param seg_id Character scalar; segment identifier (row name) used in warning messages.
#'
#' @return A numeric vector of length `length(term_labels) + 1`.
#'   The first element is the dispersion estimate actually used (either `disp_param` or `summary(fit_obj)$dispersion`),
#'    and the remaining elements are p‐values (Chisq tests) for each term in `term_labels`, in the same order.
#'   If `fit_obj` is `NA` (fitting failed), returns a vector of `NA`s of appropriate length.
#'
#' @details
#' 1. If `disp_param` is `NULL`, we take `summary(fit_obj)$dispersion` as the overdispersion estimate.
#'    If the model has zero residual degrees of freedom, dispersion cannot be estimated; in that case, we return `NA` for dispersion and issue a warning (unless `overdisp = FALSE`).
#' 2. We set `disp_used = max(1, <disp_estimate>)` so that the effective dispersion is at least 1. If `overdisp = FALSE`, we force `disp_used = 1`.
#' 3. We then call `\link[stats]{anova}(fit_obj, test = "Chisq", dispersion = disp_used)`.
#'    The p‐values for each term appear in column 5 of the ANOVA table.
#' 4. If any error occurs when computing the ANOVA, we return `rep(NA, length(term_labels) + 1)`.
#'
#' @seealso `\link{fit_as_glm}`, `\link[stats]{anova}`, `\link[stats]{glm}`
#' @export
qbinom_lrt <- function(
    fit_obj,
    overdisp = TRUE,
    disp_param = NULL,
    term_labels,
    seg_id) {
  # Number of output elements: 1 for dispersion, plus one per term
  n_out <- length(term_labels) + 1
  result <- rep(NA_real_, n_out)

  # If fitting failed (NA), return all NAs
  if (is.na(fit_obj)[1]) {
    return(result)
  }

  # Determine dispersion estimate
  if (is.null(disp_param)) {
    disp_est <- summary(fit_obj)$dispersion
    # If no residual DF, cannot estimate dispersion
    if (fit_obj$df.residual == 0) {
      disp_est <- NA_real_
      if (overdisp) {
        warning(paste0(
          "Segment ", seg_id,
          ": cannot estimate overdispersion without replicates; returning NA"
        ))
      }
    }
  } else {
    disp_est <- disp_param
  }

  # Use at least 1 for dispersion; if overdisp = FALSE, force 1
  disp_used <- max(1, disp_est, na.rm = TRUE)
  if (!overdisp) {
    disp_used <- 1
  }
  result[1] <- disp_used

  # Run ANOVA (Chisq) with specified dispersion; catch errors
  a_tbl <- tryCatch(
    stats::anova(fit_obj, test = "Chisq", dispersion = disp_used),
    error = function(e) {
      warning(paste0(
        "Segment ", seg_id,
        ": ANOVA Chi-square test failed: ", e$message
      ))
      return(NULL)
    }
  )
  if (is.null(a_tbl)) {
    return(result)
  }

  # Extract p‐values for each term: they occupy rows 2:(1+length(term_labels)), column 5
  # (anova table has rows: intercept, term1, term2, ..., Residuals)
  pvals <- a_tbl[2:(1 + length(term_labels)), 5]
  result[-1] <- pvals
  return(result)
}


#' Determine grouping factor from a data.frame or SummarizedExperiment
#'
#' If `groupby` is a vector whose length matches the number of rows (or columns), this function treats it as the grouping factor directly.
#' Otherwise, if `groupby` is one or more column names in a data.frame (or in `colData(se)`),
#'  it pastes those columns together (with a delimiter) to form a grouping factor.
#'
#' @param x Either a `SummarizedExperiment` (in which case its `colData()` is used) or a `data.frame` whose rows correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `nrow(x)`: treated as the grouping factor directly.
#'     \item A character scalar or vector of column names in `x`: those columns are pasted (row‐wise) with a delimiter to form the grouping factor.
#'   }
#' @param sep Character string used to separate pasted values when `groupby` is multiple column names. Default is `"\$"`.
#'
#' @return A character vector (or factor) of length `nrow(x)`, representing the grouping factor.
#'   If `groupby` was already the same length as `nrow(x)`, it is returned unchanged.
#'   Otherwise, it pastes together columns of `x`.
#'
#' @examples
#' df <- data.frame(
#'   sample = paste0("S", 1:4),
#'   celltype = c("A", "A", "B", "B"),
#'   batch = c("X", "Y", "X", "Y")
#' )
#'
#' # Case 1: groupby is a factor vector
#' grp1 <- get_groupby_factor(df, groupby = c("A", "A", "B", "B"))
#'
#' # Case 2: groupby is one column name
#' grp2 <- get_groupby_factor(df, groupby = "celltype")
#'
#' # Case 3: groupby is multiple column names
#' grp3 <- get_groupby_factor(df, groupby = c("celltype", "batch"))
#'
#' print(grp1)
#' print(grp2)
#' print(grp3)
#'
#' @seealso `\link{filter_segments_and_samples}`, `\link{test_all_groups_as}`
#' @export
get_groupby_factor <- function(x, groupby, sep = DELIMITER) {
  # If x is SummarizedExperiment, extract its colData as a data.frame
  if ("SummarizedExperiment" %in% class(x)) {
    x <- as.data.frame(SummarizedExperiment::colData(x))
  }
  # Now x must be a data.frame; number of rows = number of samples/pseudobulks
  n <- nrow(x)

  # If groupby is a vector whose length equals n, just return it
  if (length(groupby) == n) {
    return(groupby)
  }

  # Otherwise, check that every element of groupby is a column in x
  if (all(groupby %in% colnames(x))) {
    #    a) Select those columns: x[, groupby, drop = FALSE]
    #    b) Paste them together row-wise with the separator “sep”
    pasted <- do.call(paste, c(x[, groupby, drop = FALSE], sep = sep))
    return(pasted)
  }

  stop(
    "`groupby` must be either a factor/vector of length ", n,
    " or a column name (or names) in the provided data.frame."
  )
}


#' Test alternative splicing across all groups simultaneously
#'
#' For a `SummarizedExperiment` of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests whether PSI differs across multiple groups (e.g., cell types) using a quasi‐binomial GLM.
#' It returns a data.frame of p‐values, adjusted FDR, and delta‐PSI for each segment.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segment IDs; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#'
#' @return A data.frame with one row per segment (rownames = `rownames(pbas)` (segment IDs)), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a single grouping vector `group_factor` via `get_groupby_factor()`.
#' 2. Calls `fit_as_glm()` with formula `x ~ group`, where `x` is the per‐segment `cbind(i,e)`.
#' 3. Extracts raw p‐values for the `group` term, adjusts them (Benjamini-Hochberg) into `group_fdr`.
#' 4. Computes delta‐PSI (`dpsi`), `low_state`, and `high_state` for each segment via `get_dpsi()`.
#' 5. Combine results (endure matching row order)
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, `\link{get_dpsi}`, `\link{test_pair_as}`
#' @export
test_all_groups_as <- function(
    pbas,
    groupby,
    parallel = FALSE) {
  # 1. Construct grouping factor
  group_factor <- get_groupby_factor(pbas, groupby)

  # 2. Fit GLMs per segment, returning a data.frame of (overdispersion + p‐values)
  pv_df <- fit_as_glm(
    pbas = pbas,
    formula = x ~ group,
    data_terms = list(group = group_factor),
    return_pv = TRUE,
    parallel = parallel
  )
  # pv_df has columns: overdispersion, group (p‐value)

  # 3. Adjust p‐values (Benjamini-Hochberg)
  pv_df <- as.data.frame(pv_df, stringsAsFactors = FALSE)
  pv_df$group_fdr <- stats::p.adjust(pv_df$group, method = "BH")

  # 4. Compute delta‐PSI, low_state, high_state
  dpsi_df <- get_dpsi(pbas, group_factor)

  # 5. Combine results: ensure matching row order
  res <- cbind(
    overdispersion = pv_df$overdispersion,
    group = pv_df$group,
    group_fdr = pv_df$group_fdr,
    dpsi_df[c("low_state", "high_state", "dpsi")]
  )
  rownames(res) <- rownames(pbas)
  return(res)
}


#' Test alternative splicing between two conditions (pairwise)
#'
#' For a `SummarizedExperiment` of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests for differential splicing between two specified groups (e.g., two cell types) using a quasi‐binomial GLM per segment.
#' It returns a data.frame of p‐values, FDR, and delta‐PSI between the two conditions.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param conditions A length‐2 character vector specifying the two group labels to compare (e.g., `c("A", "B")`).
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#'
#' @return A data.frame with one row per segment (rownames = `rownames(pbas)`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{pv}{Raw p‐value for the group term (two‐level factor).}
#'     \item{fdr}{Benjamini-Hochberg FDR (across all segments).}
#'     \item{dpsi}{Difference in mean PSI between `conditions[2]` and `conditions[1]`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a single grouping vector `group_factor` via `get_groupby_factor()`.
#' 2. Filters `group_factor` and `pbas` to keep only columns corresponding to the two `conditions`.
#' 3. Calls `fit_as_glm()` with `formula = x ~ group` to get raw p‐values for each segment.
#' 4. Adjusts p‐values (Benjamini-Hochberg) into `fdr`.
#' 5. Computes delta‐PSI (`dpsi`) per segment by pseudobulking and calling `calc_psi()`, then taking
#'    mean PSI per condition and subtracting.
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, \code{\link{calc_psi}}, \link{pseudobulk}, `\link{test_all_groups_as}`
#' @export
test_pair_as <- function(
    pbas,
    groupby,
    conditions,
    parallel = FALSE) {
  # 1. Construct grouping factor
  group_factor <- get_groupby_factor(pbas, groupby)

  # 2. Identify columns (samples) that belong to the two conditions
  keep_cols <- which(group_factor %in% conditions)
  group_factor <- group_factor[keep_cols]
  pbas_sub <- pbas[, keep_cols]

  # 3. Fit GLMs per segment on only the two conditions, returning raw p‐values
  pv_df <- fit_as_glm(
    pbas = pbas_sub,
    formula = x ~ group,
    data_terms = list(group = group_factor),
    return_pv = TRUE,
    parallel = parallel
  )
  # pv_df has columns: overdispersion, group

  # 4. Adjust p‐values (Benjamini-Hochberg)
  pv_df <- as.data.frame(pv_df, stringsAsFactors = FALSE)
  pv_df$fdr <- stats::p.adjust(pv_df$group, method = "BH")

  # 5. Compute PSI matrix for the two‐condition pseudobulk
  pb <- pseudobulk(pbas_sub, group_factor)
  psi_mat <- calc_psi(pb)
  # Extract PSI for the two conditions in the order given by 'conditions'
  psi_vals <- psi_mat[, conditions, drop = FALSE]
  dpsi_vec <- psi_vals[, 2] - psi_vals[, 1]

  # 6. Combine results into final data.frame
  res <- data.frame(
    overdispersion = pv_df$overdispersion,
    pv = pv_df$group,
    fdr = pv_df$fdr,
    dpsi = dpsi_vec,
    row.names = rownames(pbas_sub),
    stringsAsFactors = FALSE
  )

  # Ensure rows are in the same order as the original pbas
  res <- res[rownames(pbas), , drop = FALSE]
  return(res)
}


#' Identify marker segments for each group via one‐vs‐rest tests
#'
#' For a `SummarizedExperiment` of splicing segments with assays “psi” (percent spliced‐in),
#'  this function conducts, for each unique group label, a one‐vs‐rest quasi‐binomial GLM test
#'  (i.e., comparing that group to all other samples).
#' It returns p‐values and delta‐PSI matrices (segments × groups).
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays “i”, “e”, and “psi”.
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#' @param verbose Logical; if `TRUE`, prints each group being tested. Default `FALSE`.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{pv}{Matrix: Raw p‐values, with rows = segments, columns = group labels}
#'     \item{fdr}{Matrix: Benjamini-Hochberg FDR (across all segments)}
#'     \item{dpsi}{Matrix: Delta‐PSI (mean PSI(group) - mean PSI(others)) per segment}
#'   }
#'   Row names of each matrix are `rownames(pbas)`; column names are the unique groups.
#'
#' @details
#' 1. Builds a grouping factor via `get_groupby_factor(pbas, groupby)`.
#' 2. For each unique label `group` in that factor:
#'    \enumerate{
#'      \item Creates a binary factor `f` where `f == TRUE` if sample’s group == `group`, FALSE otherwise.
#'      \item Calls `fit_as_glm()` with `formula = x ~ f` and `data_terms = list(f = f)`, requesting p‐values (`return_pv = TRUE`).
#'      \item Extracts the “group” p‐value column (corresponding to `f`) for each segment and stores in `pv[, group]`.
#'      \item Computes `dpsi[, group]` as `mean(psi[segment, f]) - mean(psi[segment, !f])`.
#'    }
#' 3. Adjusts each column of `pv` by Benjamini-Hochberg into `fdr[, group]`.
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, `\link{calc_psi}`
#' @export
find_marker_as <- function(
    pbas,
    groupby,
    parallel = FALSE,
    verbose = FALSE) {
  # 1. Construct grouping vector (length = ncol)
  group_factor <- get_groupby_factor(pbas, groupby)
  unique_groups <- unique(group_factor)
  n_segs <- nrow(pbas)
  n_groups <- length(unique_groups)

  # 2. Initialise result matrices
  pv_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )
  dpsi_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )

  # 3. For each group, perform one‐vs‐rest test
  for (group in unique_groups) {
    if (verbose) {
      message("Testing group: ", group)
    }
    # Binary factor: TRUE if sample belongs to this group
    f <- group_factor == group

    # Fit GLM per segment; return a data.frame with columns: overdispersion and 'group' p‐value
    glm_res <- fit_as_glm(
      pbas       = pbas,
      formula    = x ~ f,
      data_terms = list(f = f),
      return_pv  = TRUE,
      parallel   = parallel
    )
    # glm_res[data.frame]: rownames = segment IDs, col 'group' = p‐value

    # Store raw p‐values for this group
    pv_mat[, group] <- glm_res[, "group"]

    # Compute dpsi: for each segment, mean PSI in 'group' minus mean PSI in others
    psi_vals <- SummarizedExperiment::assay(pbas, "psi")

    # For rows with all-NA in either subset, result is NA automatically
    mean_ct <- rowMeans(psi_vals[, f, drop = FALSE], na.rm = TRUE)
    mean_other <- rowMeans(psi_vals[, !f, drop = FALSE], na.rm = TRUE)
    dpsi_mat[, group] <- mean_ct - mean_other
  }

  # 4. Adjust p‐values column‐wise (BH)
  fdr_mat <- apply(pv_mat, 2, stats::p.adjust, method = "BH")

  # 5. Return list of matrices
  return(list(
    pv = pv_mat,
    fdr = fdr_mat,
    dpsi = dpsi_mat
  ))
}


#' Select markers from all‐celltype test results
#'
#' From a data.frame of per‐segment results comparing all groups simultaneously (e.g., output of `test_all_groups_as()`),
#'   this function selects only those segments that pass both a group‐level FDR threshold and a minimum delta‐PSI threshold.
#' It returns a data.frame with one row per qualifying segment.
#'
#' @param all_celltype_df A data.frame (rownames = segment IDs) (output of `test_all_groups_as()`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#'
#' @return A data.frame with one row per segment satisfying `group_fdr < fdr_thr` and `dpsi > dpsi_thr`.
#'   Columns in the returned data.frame:
#'   \describe{
#'     \item{`pv`}{Raw p‐value from the likelihood‐ratio test for the group effect (copied from `all_celltype_df$group`).}
#'     \item{`fdr`}{Benjamini‐Hochberg adjusted FDR (copied from `all_celltype_df$group_fdr`).}
#'     \item{`dpsi`}{Delta‐PSI (copied from `all_celltype_df$dpsi`).}
#'     \item{`seg_id`}{Segment ID (rownames of `all_celltype_df`).}
#'     \item{`group`}{Group label with highest mean PSI (copied from `all_celltype_df$high_state`).}
#'   }
#'   Row names of the returned data.frame are the segment IDs.
#'
#' @examples
#' \dontrun{
#' # Assume 'all_test' is output of test_all_groups_as(), with rownames = segment IDs.
#' sig_df <- select_markers_from_all_celltype_test(all_test, fdr_thr = 0.05, dpsi_thr = 0.2)
#' }
#'
#' @export
select_markers_from_all_celltype_test <- function(
    all_celltype_df,
    fdr_thr = 0.05,
    dpsi_thr = 0.1) {
  # Identify segments passing both thresholds
  keep <- (all_celltype_df$group_fdr < fdr_thr) & (all_celltype_df$dpsi > dpsi_thr)
  subset_df <- all_celltype_df[keep, , drop = FALSE]

  # Build result data.frame
  result <- data.frame(
    pv = subset_df$group,
    fdr = subset_df$group_fdr,
    dpsi = subset_df$dpsi,
    seg_id = rownames(subset_df),
    group = subset_df$high_state,
    row.names = rownames(subset_df),
    stringsAsFactors = FALSE
  )
  return(result)
}


#' Select top marker segments per group (celltype)
#'
#' From a list of per‐segment test results (as returned by `find_marker_as()`),
#'  this function selects up to `n` top segments per group
#'  whose FDR is below a threshold and whose absolute delta‐PSI exceeds a threshold.
#' Optionally, it removes duplicate segments so that each segment appears only once
#'  (assigned to the group where it has the largest |dpsi|).
#'
#' @param markers A list containing numeric matrices `pv`, `fdr`, and `dpsi`,
#'   each a matrix with rows = segments and columns = group labels (as returned by `find_marker_as()`).
#'   \describe{
#'     \item{`markers$pv`}{Matrix of raw p‐values, dimensions: ∣segments∣ × ∣groups∣.}
#'     \item{`markers$fdr`}{Matrix of Benjamini‐Hochberg adjusted FDR p‐values, same dimensions.}
#'     \item{`markers$dpsi`}{Matrix of delta‐PSI (mean PSI(group) - mean PSI(others)), same dimensions.}
#'   }
#' @param n Integer: maximum number of markers to return per cell type (default: 5).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#' @param clean_duplicates Logical: if `TRUE`, keep only one entry per segment
#'                         (the one with highest |dpsi|) across all cell types (default: `TRUE`).
#'
#' @return A data.frame with selected marker segments. Columns:
#'   \describe{
#'     \item{`pv`}{Raw p‐value from the likelihood‐ratio test for the group effect in the selected group.}
#'     \item{`fdr`}{Benjamini‐Hochberg adjusted FDR value for the segment in the selected group.}
#'     \item{`dpsi`}{Delta‐PSI (signed) for the segment in the selected group.}
#'     \item{`seg_id`}{Segment identifier (rownames of `markers$pv`).}
#'     \item{`group`}{Group label for which this segment is selected.}
#'   }
#'   Rows are ordered by descending |dpsi| across all selected segments.
#'
#' @details
#' For each group (column) in `markers$fdr`, the function:
#' \enumerate{
#'   \item Identifies segments `seg` satisfying
#'         `markers$fdr[seg, group] <= fdr_thr` and
#'         `|markers$dpsi[seg, group]| >= dpsi_thr`.
#'   \item If at least one segment passes, constructs a temporary data.frame `t` containing
#'         `pv = markers$pv[seg, group]`,
#'         `fdr = markers$fdr[seg, group]`,
#'         `dpsi = markers$dpsi[seg, group]`,
#'         `seg_id = seg`,
#'         `group = group`.
#'   \item Selects up to `n` rows from `t`, ordering by `abs(dpsi)` in descending order.
#'   \item Appends these rows to an aggregate `res` data.frame.
#' }
#' If `clean_duplicates = TRUE`, it then splits `res` by `seg_id` and retains only the row
#'  where `abs(dpsi)` is maximal across groups, removing duplicates.
#' Finally, it re‐orders `res` by `abs(dpsi)` in descending order and returns it.
#'
#' @examples
#' # Load example data
#' data(pbasf)
#'
#' # Select top 10 markers with stricter thresholds:
#' df <- select_markers(pbasf@metadata$markers, n = 10, fdr_thr = 0.01, dpsi_thr = 0.2)
#' head(df)
#'
#' @seealso `\link{find_marker_as}`, `\link{test_all_groups_as}`, `\link{test_pair_as}`
#' @export
select_markers <- function(
    markers,
    n = 5,
    fdr_thr = 0.05,
    dpsi_thr = 0.1,
    clean_duplicates = TRUE) {
  # Initialise an empty result
  result <- NULL

  # Iterate over each group (celltype) (column name) in the fdr matrix
  for (group in colnames(markers$fdr)) {
    # 1. Identify segments passing both FDR and |dpsi| thresholds
    pass_fdr <- markers$fdr[, group] <= fdr_thr
    pass_dpsi <- abs(markers$dpsi[, group]) >= dpsi_thr
    keep_idx <- which(pass_fdr & pass_dpsi)

    # Replace NA flags with FALSE
    keep_idx <- keep_idx[!is.na(keep_idx)]

    # If no segments pass, skip to next group
    if (length(keep_idx) == 0) {
      next
    }

    # 2. Build a temporary data.frame of candidate segments for this group
    temp_df <- data.frame(
      pv = markers$pv[keep_idx, group],
      fdr = markers$fdr[keep_idx, group],
      dpsi = markers$dpsi[keep_idx, group],
      seg_id = rownames(markers$pv)[keep_idx],
      group = group,
      stringsAsFactors = FALSE
    )

    # 3. Order by descending |dpsi| and take up to 'n' rows
    order_idx <- order(abs(temp_df$dpsi), decreasing = TRUE)
    top_n_idx <- utils::head(order_idx, n)
    temp_sel <- temp_df[top_n_idx, , drop = FALSE]

    # 4. Append to the cumulative result
    result <- if (is.null(result)) {
      temp_sel
    } else {
      rbind(result, temp_sel)
    }
  }

  # Reset row names of the aggregated result
  rownames(result) <- NULL

  # 5. If clean_duplicates, ensure each seg_id appears only once
  if (clean_duplicates && nrow(result) > 0) {
    # Split by seg_id; within each group, keep row with largest |dpsi|
    by_seg <- split(result, result$seg_id)
    cleaned_list <- lapply(by_seg, function(df_seg) {
      # Order rows for this segment by descending |dpsi|
      idx <- order(abs(df_seg$dpsi), decreasing = TRUE)
      df_seg[idx[1], , drop = FALSE]
    })
    result <- do.call(rbind, cleaned_list)
    rownames(result) <- NULL
  }

  # 6. Finally, order all rows by descending |dpsi|
  if (nrow(result) > 0) {
    result <- result[order(abs(result$dpsi), decreasing = TRUE), , drop = FALSE]
  }

  return(result)
}


#' Select combined marker segments from per‐group and all‐groups tests
#'
#' This function merges two sets of marker results:
#'   1. Per‐group markers (output of `select_markers`): segments deemed significant in one‐vs‐rest tests.
#'   2. All‐groups markers (output of `select_markers_from_all_celltype_test`): segments significant across all groups.
#'
#' For each group, it prioritizes segments identified by one‐vs‐rest tests; any segment not already selected but significant in the all‐groups test is labeled as a “background” marker.
#' The result contains at most `n` segments per group, with ties broken by largest |dpsi|.
#'
#' @param markers A list containing numeric matrices `pv`, `fdr`, and `dpsi`,
#'   each a matrix with rows = segments and columns = group labels (as returned by `find_marker_as()`).
#'   \describe{
#'     \item{`markers$pv`}{Matrix of raw p‐values, dimensions: ∣segments∣ × ∣groups∣.}
#'     \item{`markers$fdr`}{Matrix of Benjamini‐Hochberg adjusted FDR p‐values, same dimensions.}
#'     \item{`markers$dpsi`}{Matrix of delta‐PSI (mean PSI(group) - mean PSI(others)), same dimensions.}
#'   }
#' @param all_celltype_df A data.frame (rownames = segment IDs) (output of `test_all_groups_as()`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#' @param n Integer: maximum number of markers to return per cell type (after combining) (default: `Inf` (no limit)).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#'
#' @return A data.frame with one row per selected segment, containing:
#'   \describe{
#'     \item{`pv`}{Raw p‐value (from either per‐group or all‐groups test).}
#'     \item{`fdr`}{Adjusted FDR (from the same test).}
#'     \item{`dpsi`}{Delta‐PSI (for the group that nominated this segment).}
#'     \item{`seg_id`}{Segment IDs (rownames).}
#'     \item{`group`}{Group label for which the segment is selected.}
#'     \item{`is_marker`}{Logical: `TRUE` if selected by the per‐group test, `FALSE` if only from the all‐groups test.}
#'   }
#'   Row names are the segment IDs. Rows are ordered within each group by descending |dpsi| and limited to `n` per group.
#'
#' @details
#' 1. Per‐group markers (priority):
#'    Call `select_markers(markers, n = Inf, fdr_thr = fdr_thr, dpsi_thr = dpsi_thr, clean_duplicates = TRUE)`
#'     to gather all segments that pass thresholds in any group. Label these as `is_marker = TRUE`.
#' 2. All‐groups markers (background):
#'    Call `select_markers_from_all_celltype_test(all_celltype_df, fdr_thr = fdr_thr, dpsi_thr = dpsi_thr)`
#'     to gather segments passing the all‐groups thresholds. Exclude any segment already in step 1. Label remaining as `is_marker = FALSE`.
#' 3. Combine and trim:
#'    For each group, among combined rows (priority first), keep at most `n` segments ranked by descending |dpsi|.
#'
#' @seealso `\link{select_markers}`, `\link{select_markers_from_all_celltype_test}`
#' @export
select_all_markers <- function(
    markers,
    all_celltype_df,
    n = Inf,
    fdr_thr = 0.05,
    dpsi_thr = 0.1) {
  # 1. Fetch all per‐group markers without limiting n
  per_group_df <- select_markers(
    markers,
    n = Inf,
    fdr_thr = fdr_thr,
    dpsi_thr = dpsi_thr,
    clean_duplicates = TRUE
  )
  # Add a column to flag these as true markers
  per_group_df$is_marker <- TRUE

  # 2. Fetch all‐groups background markers
  background_df <- select_markers_from_all_celltype_test(
    all_celltype_df,
    fdr_thr = fdr_thr,
    dpsi_thr = dpsi_thr
  )
  # Rename 'seg_id' column consistently
  # (select_markers_from_all_celltype_test already has seg_id and group)
  background_df$is_marker <- FALSE

  # 3. Exclude any segments already in per_group_df
  bg_keep <- setdiff(background_df$seg_id, per_group_df$seg_id)
  background_df <- background_df[background_df$seg_id %in% bg_keep, , drop = FALSE]

  # 4. Combine the two sets
  combined <- rbind(per_group_df, background_df)

  # 5. For each group, limit to top 'n' by |dpsi|
  # Split by group
  split_by_group <- split(combined, combined$group)
  trimmed_list <- lapply(split_by_group, function(df_group) {
    # Order by descending |dpsi|
    df_group <- df_group[order(abs(df_group$dpsi), decreasing = TRUE), , drop = FALSE]
    # Trim to at most n rows
    if (nrow(df_group) > n) {
      df_group <- df_group[seq_len(n), , drop = FALSE]
    }
    return(df_group)
  })
  final_df <- do.call(rbind, trimmed_list)
  rownames(final_df) <- final_df$seg_id

  # 6. Re‐order rows by group (alphabetical) then |dpsi| descending
  final_df <- final_df[order(final_df$group, -abs(final_df$dpsi)), , drop = FALSE]

  return(final_df)
}


###### prepare_reference.R ######

#' Annotate coding status of segments using a GTF file
#'
#' Given a path to an Ensembl GTF (gene annotation) and an AS result list containing segment coordinates,
#'  this function determines whether each segment overlaps coding sequence (CDS) regions.
#' It assigns coding status:
#'   - 'c' if the segment is fully contained within a CDS (i.e., likely constitutive exon),
#'   - 'p' if the segment partially overlaps a CDS,
#'   - 'n' if the segment does not overlap any CDS.
#' It also flags whether the segment’s gene has at least one coding overlap.
#'
#' @param gtf_path Character: path to an Ensembl GTF file
#'   (tab‐delimited, with columns: `chr_id`, `feature` (e.g., "CDS"), `start`, `stop`, `strand`, and gene metadata).
#'   Must be readable by `read.table(..., comment.char = "#")`.
#' @param as_list SAJR::loadSAData() output list containing segment coordinates.
#' List containing at least these components:
#'   \describe{
#'     \item{`seg`}{A data.frame with columns `chr_id`, `start`, `stop`, `strand`, and `gene_id` for each segment. Row names are segment IDs.}
#'   }
#'
#' @return The input `as_list` with additional columns:
#'   \describe{
#'     \item{`seg$cod`}{Character vector of length ∣segments∣.
#'           Each element is `'c'` if the segment is completely within a CDS, `'p'` if partially overlaps, or `'n'` if no overlap.}
#'     \item{`seg$cod.gene`}{Logical vector of length ∣segments∣. `TRUE` if the segment’s gene contains at least one segment with `cod != 'n'`.}
#'   }
#'   The modified `as_list` is returned invisibly.
#'
#' @details
#' 1. Read the GTF file via `read.table(gtf_path, sep = "\t", header = FALSE, comment.char = "#")` and
#'    filter to rows where the second column (`feature`) == "CDS".
#'    Keep columns: `chr_id`, `start`, `stop`, `strand`.
#' 2. Force CDS entries to a `GRanges` using `GenomicRanges::GRanges()`.
#'    - If some `chr_id` values lack the "chr" prefix but segments use it (or vice versa), prefix or strip "chr" to match.
#' 3. Reduce the CDS ranges per chromosome/strand to merge overlapping exons (via `reduce()`).
#' 4. Build a `GRanges` for all segments from `as_list$seg`, using `chr_id`, `start`, `stop`, and `strand`.
#' 5. Use `findOverlaps(seg_gr, cds_gr, type = "any")` to find any overlap, and `findOverlaps(seg_gr, cds_gr, type = "within")` to find full containment.
#' 6. Initialise `cod = 'n'` for all segments. For indices in the "any" overlap set, set `cod = 'p'`. For indices in the "within" set, set `cod = 'c'`.
#' 7. Determine `cod.gene` by marking as `TRUE` any segment whose `gene_id` appears among those with `cod != 'n'`.
#'
#' @note
#' - Requires `GenomicRanges` for `GRanges`, `reduce()`, and `findOverlaps()`.
#' - Uses `IRanges::IRanges()` internally for range construction.
#' - If chromosome naming mismatches occur ("chr1" vs. "1"), this function attempts to harmonize by adding "chr" prefix where needed.
#'
#' @seealso `\link[GenomicRanges]{GRanges}`, `\link[GenomicRanges]{findOverlaps}`, `\link[GenomicRanges]{reduce}`
#' @export
add_is_coding_by_ens_gtf <- function(gtf_path, as_list) {
  # 1. Read GTF and keep only CDS lines
  gtf_df <- utils::read.table(gtf_path,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  # Columns: V1=chr_id, V3=feature, V4=start, V5=stop, V7=strand
  cds_df <- gtf_df[gtf_df[, 3] == "CDS", , drop = FALSE]
  colnames(cds_df)[c(1, 4, 5, 7)] <- c("chr_id", "start", "stop", "strand") # Rename columns for clarity
  cds_df <- cds_df[, c("chr_id", "start", "stop", "strand")]


  # 2. CDS Entries: GRanges
  # Harmonize chromosome naming: if some segment chr_ids have "chr" prefix but CDS lack it, or vice versa
  seg_chr <- as_list$seg$chr_id
  if (any(startsWith(seg_chr, "chr")) && !any(startsWith(cds_df$chr_id, "chr"))) {
    cds_df$chr_id <- paste0("chr", cds_df$chr_id)
  } else if (!any(startsWith(seg_chr, "chr")) && any(startsWith(cds_df$chr_id, "chr"))) {
    cds_df$chr_id <- sub("^chr", "", cds_df$chr_id)
  }

  # Build GRanges for CDS
  cds_gr <- GenomicRanges::GRanges(
    seqnames = cds_df$chr_id,
    ranges   = IRanges::IRanges(start = cds_df$start, end = cds_df$stop),
    strand   = cds_df$strand
  )


  # 3. Merge overlapping CDS ranges
  cds_gr <- GenomicRanges::reduce(cds_gr)


  # 4. Build GRanges for segments
  seg_df <- as_list$seg
  seg_gr <- GenomicRanges::GRanges(
    seqnames = seg_df$chr_id,
    ranges = IRanges::IRanges(start = seg_df$start, end = seg_df$stop),
    strand = ifelse(is.na(seg_df$strand), "*",
      ifelse(seg_df$strand == 1, "+", ifelse(seg_df$strand == -1, "-", "*"))
    )
  )


  # 5. Find overlaps: any vs. within
  ol_any <- GenomicRanges::findOverlaps(seg_gr, cds_gr, type = "any", ignore.strand = FALSE)
  ol_within <- GenomicRanges::findOverlaps(seg_gr, cds_gr, type = "within", ignore.strand = FALSE)


  # 6. Initialise coding status
  n_segs <- length(seg_gr)
  cod_stat <- rep("n", n_segs)

  # Mark partial overlaps as "p"
  ol_any_idx <- unique(S4Vectors::queryHits(ol_any))
  cod_stat[ol_any_idx] <- "p"

  # Mark full containment as "c"
  ol_within_idx <- unique(S4Vectors::queryHits(ol_within))
  cod_stat[ol_within_idx] <- "c"

  # Assign to as_list$seg$cod
  as_list$seg$cod <- cod_stat


  # 7. Determine which genes have at least one coding segment
  coding_genes <- unique(seg_df$gene_id[cod_stat != "n"])
  as_list$seg$cod.gene <- seg_df$gene_id %in% coding_genes

  return(as_list)
}

#################################


#' Aggregate single-cell counts into pseudobulk per group
#'
#' Given a `SummarizedExperiment` of single-cell assays (e.g., inclusion/exclusion counts),
#'  this function sums counts within each group (as defined by `groupby`), removes specified derivative assays,
#'  and returns a new `SummarizedExperiment` where each column corresponds to a group-level pseudobulk.
#'
#' @param se A `SummarizedExperiment` where rows are features (e.g., segments) and columns are individual cells or barcodes.
#'           Assays include `i`, `e`, `counts`, etc.
#' @param groupby Either:
#'   \itemize{
#'     \item A character vector of length `ncol(se)` giving the group label for each column, or
#'     \item One or more column names in `colData(se)`, whose values (pasted together if multiple)
#'           define the grouping factor via `get_groupby_factor()`.
#'   }
#'   After grouping, columns with the same label will be summed together.
#' @param clean_derivatives Character vector of assay names to remove before summation
#'   (e.g., `c("psi", "cpm")`). Default `c("psi", "cpm")`.
#'
#' @return A new `SummarizedExperiment` with:
#'   \describe{
#'     \item{Assays}{Each assay from `se`, except those in `clean_derivatives`, replaced by a matrix of summed values per group.}
#'     \item{rowRanges}{Copied from `rowRanges(se)`.}
#'     \item{colData}{A data.frame of group-level metadata: one row per group, including any columns from `colData(se)`
#'                      that are constant within each group, plus aggregated columns via `pseudobulk_metadata()`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a grouping vector via `get_groupby_factor(se, groupby)`.
#' 2. Drops any assays listed in `clean_derivatives` from `se`.
#' 3. For each remaining assay, uses `visutils::calcColSums()` to sum counts for columns that share the same
#'    group label.
#' 4. Builds a metadata data.frame by calling `pseudobulk_metadata(colData(se), groupby)`, which:
#'    \enumerate{
#'      \item Splits `colData(se)` into sub-data.frames by group.
#'      \item For each group, retains columns that are constant across cells and applies aggregation functions
#'            (e.g., sum of `ncells`) for specified columns.
#'      \item Recombines group-level rows into a single data.frame, with row names matching group labels.
#'    }
#' 5. Constructs a new `SummarizedExperiment` with summed assay matrices, original `rowRanges(se)`, and
#'    the aggregated `colData` (one row per group).
#'
#' @seealso \code{\link{get_groupby_factor}}, \code{\link{pseudobulk_metadata}}, \code{\link[visutils]{calcColSums}}
#' @export
pseudobulk <- function(
    se,
    groupby,
    clean_derivatives = c("psi", "cpm")) {
  # 1. Determine grouping labels for each column
  group_factor <- get_groupby_factor(se, groupby)

  # 2. Remove specified derivative assays before summing
  to_remove <- intersect(SummarizedExperiment::assayNames(se), clean_derivatives)
  if (length(to_remove) > 0) {
    for (assay_name in to_remove) {
      SummarizedExperiment::assay(se, assay_name) <- NULL
    }
  }

  # 3. Sum each assay by group using visutils::calcColSums
  summed_assays <- list()
  for (assay_name in SummarizedExperiment::assayNames(se)) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    summed_assays[[assay_name]] <- visutils::calcColSums(mat, group_factor)
  }

  # 4. Build group-level metadata
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  group_meta <- pseudobulk_metadata(meta, group_factor)
  # Reorder to match the order of columns in summed_assays (rownames of each summed matrix)
  group_meta <- group_meta[colnames(summed_assays[[1]]), , drop = FALSE]

  # 5. Construct new SummarizedExperiment
  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays = summed_assays,
    rowRanges = SummarizedExperiment::rowRanges(se),
    colData = group_meta
  )

  return(new_se)
}


#' Calculate percent spliced‐in (PSI) per segment
#'
#' Given a `SummarizedExperiment` containing inclusion (`i`) and exclusion (`e`) counts,
#'  this function computes, for each segment and sample, the percent spliced‐in: `PSI = i / (i + e)`
#' If the total counts (`i + e`) for a segment in a sample are below `min_cov`, that PSI is set to `NA`.
#'
#' @param se A `SummarizedExperiment` with assays named `i` (inclusion counts) and `e` (exclusion counts).
#'   Rows are segments; columns are pseudobulk samples.
#' @param min_cov Integer scalar; minimum total junction coverage (`i + e`) required to compute a valid PSI.
#'   For any `i + e < min_cov`, the PSI is set to `NA`. Default is `10`.
#'
#' @return A numeric matrix of dimensions `nrow(se)` × `ncol(se)` where each entry is the PSI value
#'   for that segment and sample, or `NA` if coverage is insufficient.
#'
#' @examples
#' \dontrun{
#' # Assume 'pseudobulk_se' is a SummarizedExperiment with assays 'i' and 'e':
#' psi_mat <- calc_psi(pseudobulk_se, min_cov = 10)
#' head(psi_mat)
#' }
#'
#' @seealso \code{\link{pseudobulk}}, \code{\link{calc_cpm}}
#' @export
calc_psi <- function(se, min_cov = 10) {
  # Verify that both 'i' and 'e' assays exist
  if (!all(c("i", "e") %in% SummarizedExperiment::assayNames(se))) {
    warning("Assays 'i' and/or 'e' not found in the SummarizedExperiment; returning NULL.")
    return(NULL)
  }

  # Extract inclusion and exclusion count matrices
  inc_mat <- SummarizedExperiment::assay(se, "i")
  exc_mat <- SummarizedExperiment::assay(se, "e")

  # Compute total counts per segment × sample
  total_mat <- inc_mat + exc_mat

  # Initialize PSI matrix
  psi_mat <- inc_mat / total_mat

  # Where total < min_cov, set PSI to NA
  psi_mat[total_mat < min_cov] <- NA_real_

  return(as.matrix(psi_mat))
}
