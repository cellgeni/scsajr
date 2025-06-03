options(error = function(e) quit("no", 1))

DELIMETER <- "$"


#' Filter segments and samples based on coverage and group criteria
#'
#' This function takes a \code{SummarizedExperiment} containing alternative splicing data
#' (with assays “i” for inclusion counts and “e” for exclusion counts), and filters out:
#'   1. pseudobulk samples with fewer than \code{sample_min_ncells} cells,
#'   2. cell‐type groups with fewer than \code{celltype_min_samples} samples,
#'   3. segments (rows) that do not meet minimum coverage across pseudobulks,
#'   4. segments with low variability (standard deviation of PSI).
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param sample_filter Logical vector (length = number of columns in \code{pbas}) indicating which samples to retain initially.
#'   Defaults to \code{TRUE} for all samples.
#' @param sites Character vector of splice‐site patterns to keep (e.g., \code{c("ad", "aa", "dd")}). Only segments whose \code{rowData(pbas)$sites} is in this set will be retained.
#' @param min_cov Minimum total junction coverage (\code{i + e}) for a pseudobulk to count toward segment inclusion. Default 10.
#' @param sample_min_ncells Minimum number of cells per pseudobulk; pseudobulk samples with fewer cells are filtered out. Default 20.
#' @param celltype_min_samples Minimum number of pseudobulk samples per group (as defined by \code{groupby}) to keep that group. Default 2.
#' @param seg_min_samples The absolute minimum number of pseudobulks (after filtering by \code{min_cov}) that a segment must be observed in. Default 4.
#' @param seg_min_samples_fraq Minimum fraction of total pseudobulk samples that must meet \code{min_cov} for a segment to be retained. Default 0.2.
#'   The effective minimum samples per segment is \code{max(seg_min_samples, ncol(pbas) * seg_min_samples_fraq)}.
#' @param seg_min_sd Minimum standard deviation of PSI (percent spliced‐in) across pseudobulks for a segment to be retained. Default 0.1.
#' @param groupby Either (a) a character vector of length \code{ncol(pbas)} giving an existing grouping factor for each column, or
#'   (b) a character scalar or vector of column names in \code{colData(pbas)} whose pasted‐together values define the grouping factor.
#'   See \code{\link{get_groupby_factor}} for details.
#'
#' @return A filtered \code{SummarizedExperiment}, containing only those pseudobulks and segments that satisfy all criteria.
#'   The returned object will always have assays “i”, “e” converted to dense matrices, and will include a “psi” assay (percent spliced‐in) for each retained segment/samples.
#'   The \code{rowData} will also contain a column “nna” (number of pseudobulks with \code{i+e >= min_cov}) and “sd” (standard deviation of PSI).
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
  # Fleter pseudobulk samples with fewer cells
  sample_filter <- sample_filter & (SummarizedExperiment::colData(pbas)$ncells >= sample_min_ncells)

  # Determine how many samples remain in each group (celltype)
  group_factor <- get_groupby_factor(pbas, groupby) # returns a factor vector (length = ncol)
  ct_counts <- table(group_factor[sample_filter]) # Count samples per group (only among those passing sample_filter)

  # Only keep samples whose group has at least 'celltype_min_samples' pseudobulks
  sample_filter <- sample_filter & (group_factor %in% names(ct_counts)[ct_counts >= celltype_min_samples])

  pbas <- pbas[, sample_filter] # Subset the SummarizedExperiment columns (samples)


  ## 2. Filter out cell‐type groups with fewer than celltype_min_samples} samples
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
#' For each segment (row) in a \code{SummarizedExperiment} with assays “i” (inclusion counts) and “e” (exclusion counts),
#' this function fits a quasi‐binomial generalized linear model of the form \code{i / (i + e) ~ predictors},
#' adding a pseudocount if desired. It can optionally return p‐values via a likelihood‐ratio test.
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param formula An \code{\link[stats]{formula}} object specifying the model, e.g., \code{i ~ group}.
#'   Internally, the function treats “i” and “e” as a two‐column response.
#' @param data_terms A named list of covariates: data.frame columns referenced in \code{formula}.
#'   For example, if \code{formula = i ~ group}, then \code{data_terms = list(group = group_factor)} where \code{group_factor} is a vector of length \code{ncol(pbas)}.
#' @param pseudocount A numeric fraction of each count total to add to both “i” and “e”, expressed as a fraction of the total counts. Default 0 (no pseudocount).
#' @param parallel Logical; if \code{TRUE}, uses \code{\link[plyr]{alply}(…, .parallel = TRUE)} to fit models in parallel (requires registered backend). Default \code{FALSE}.
#' @param progress Character string passed to \code{.progress} in \code{\link[plyr]{alply}}. Default \code{"none"}.
#' @param return_pv Logical; if \code{TRUE}, perform a quasi‐likelihood ratio test per segment to return p‐values for each term in \code{formula}.
#'   Otherwise, return fitted \code{\link[stats]{glm}} objects. Default \code{FALSE}.
#' @param overdisp Logical; if \code{TRUE}, account for overdispersion when computing test statistics. Default \code{TRUE}.
#' @param disp_param Optional numeric vector of length \code{nrow(pbas)} providing dispersion parameters per segment.
#'   If \code{NULL}, dispersion is estimated from each model’s residuals.
#'
#' @return If \code{return_pv = FALSE}, a named list of length \code{nrow(pbas)} of fitted \code{glm} objects
#'   (or \code{NA} for segments where fitting failed). Names correspond to \code{rownames(pbas)}.
#'   If \code{return_pv = TRUE}, a data.frame with one row per segment, columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion (or provided \code{disp_param}).}
#'     \item{<term>}{Likelihood‐ratio p‐value for each term in \code{formula}.}
#'   }
#'   Row names correspond to \code{rownames(pbas)}.
#'
#' @details
#' The function loops over each segment, constructs a two‐column matrix \code{cbind(i, e)}, adds pseudocounts if requested, and fits a quasi‐binomial \code{\link[stats]{glm}}.
#' If \code{return_pv = TRUE}, it calls \code{\link{qbinom_lrt}} to compute a likelihood‐ratio test for each term in \code{formula}.
#' Any errors or warnings during fitting are caught; segments with fitting errors return \code{NA}.
#'
#' @seealso \code{\link{qbinom_lrt}}, \code{\link[stats]{glm}}, \code{\link[plyr]{alply}}
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
#' Given a fitted \code{glm(..., family = "quasibinomial")}, this function computes a likelihood‐ratio test for each term in the model, accounting for overdispersion.
#'
#' @param fit_obj A fitted \code{glm} object (quasi‐binomial family).
#' @param overdisp Logical; if \code{TRUE}, use estimated dispersion from \code{fit_obj} (or \code{disp_param}) when computing test statistics.
#'   If \code{FALSE}, force dispersion = 1 (i.e., treat as simple binomial).
#' @param disp_param Optional numeric; if provided, use this as the dispersion parameter for the test.
#'   Otherwise, estimate dispersion from \code{fit_obj}’s residuals.
#' @param term_labels Character vector of term names (excluding the intercept) in the original formula, in the same order returned by \code{attr(terms(formula), "term.labels")}.
#' @param seg_id Character scalar; segment identifier (row name) used in warning messages.
#'
#' @return A numeric vector of length \code{length(term_labels) + 1}.
#'   The first element is the dispersion estimate actually used (either \code{disp_param} or \code{summary(fit_obj)$dispersion}),
#'    and the remaining elements are p‐values (Chisq tests) for each term in \code{term_labels}, in the same order.
#'   If \code{fit_obj} is \code{NA} (fitting failed), returns a vector of \code{NA}s of appropriate length.
#'
#' @details
#' 1. If \code{disp_param} is \code{NULL}, we take \code{summary(fit_obj)$dispersion} as the overdispersion estimate.
#'    If the model has zero residual degrees of freedom, dispersion cannot be estimated; in that case, we return \code{NA} for dispersion and issue a warning (unless \code{overdisp = FALSE}).
#' 2. We set \code{disp_used = max(1, <disp_estimate>)} so that the effective dispersion is at least 1. If \code{overdisp = FALSE}, we force \code{disp_used = 1}.
#' 3. We then call \code{\link[stats]{anova}(fit_obj, test = "Chisq", dispersion = disp_used)}.
#'    The p‐values for each term appear in column 5 of the ANOVA table.
#' 4. If any error occurs when computing the ANOVA, we return \code{rep(NA, length(term_labels) + 1)}.
#'
#' @seealso \code{\link{fit_as_glm}}, \code{\link[stats]{anova}}, \code{\link[stats]{glm}}
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
#' If \code{groupby} is a vector whose length matches the number of rows (or columns), this function treats it as the grouping factor directly.
#' Otherwise, if \code{groupby} is one or more column names in a data.frame (or in \code{colData(se)}),
#'  it pastes those columns together (with a delimiter) to form a grouping factor.
#'
#' @param x Either a \code{SummarizedExperiment} (in which case its \code{colData()} is used) or a \code{data.frame} whose rows correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length \code{nrow(x)}: treated as the grouping factor directly.
#'     \item A character scalar or vector of column names in \code{x}: those columns are pasted (row‐wise) with a delimiter to form the grouping factor.
#'   }
#' @param sep Character string used to separate pasted values when \code{groupby} is multiple column names. Default is \code{"\$"}.
#'
#' @return A character vector (or factor) of length \code{nrow(x)}, representing the grouping factor.
#'   If \code{groupby} was already the same length as \code{nrow(x)}, it is returned unchanged.
#'   Otherwise, it pastes together columns of \code{x}.
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
#' @seealso \code{\link{filter_segments_and_samples}}, \code{\link{test_all_groups_as}}
#' @export
get_groupby_factor <- function(x, groupby, sep = DELIMETER) {
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
#' For a \code{SummarizedExperiment} of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests whether PSI differs across multiple groups (e.g., cell types) using a quasi‐binomial GLM.
#' It returns a data.frame of p‐values, adjusted FDR, and delta‐PSI for each segment.
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segment IDs; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length \code{ncol(pbas)} that directly gives a group label for each column, or
#'     \item One or more column names in \code{colData(pbas)}, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See \code{\link{get_groupby_factor}} for details.
#' @param parallel Logical; if \code{TRUE}, fit per‐segment GLMs in parallel (requires a registered \code{plyr} backend). Default \code{FALSE}.
#'
#' @return A data.frame with one row per segment (rownames = \code{rownames(pbas)} (segment IDs)), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the \code{group} term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from \code{get_dpsi()}).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between \code{high_state} and \code{low_state}.}
#'   }
#'
#' @details
#' 1. Converts \code{groupby} into a single grouping vector \code{group_factor} via \code{get_groupby_factor()}.
#' 2. Calls \code{fit_as_glm()} with formula \code{x ~ group}, where \code{x} is the per‐segment \code{cbind(i,e)}.
#' 3. Extracts raw p‐values for the \code{group} term, adjusts them (Benjamini–Hochberg) into \code{group_fdr}.
#' 4. Computes delta‐PSI (\code{dpsi}), \code{low_state}, and \code{high_state} for each segment via \code{get_dpsi()}.
#' 5. Combine results (endure matching row order)
#'
#' @seealso \code{\link{fit_as_glm}}, \code{\link{get_groupby_factor}}, \code{\link{get_dpsi}}, \code{\link{test_pair_as}}
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

  # 3. Adjust p‐values (Benjamini–Hochberg)
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
#' For a \code{SummarizedExperiment} of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests for differential splicing between two specified groups (e.g., two cell types) using a quasi‐binomial GLM per segment.
#' It returns a data.frame of p‐values, FDR, and delta‐PSI between the two conditions.
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length \code{ncol(pbas)} that directly gives a group label for each column, or
#'     \item One or more column names in \code{colData(pbas)}, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See \code{\link{get_groupby_factor}} for details.
#' @param conditions A length‐2 character vector specifying the two group labels to compare (e.g., \code{c("A", "B")}).
#' @param parallel Logical; if \code{TRUE}, fit per‐segment GLMs in parallel (requires a registered \code{plyr} backend). Default \code{FALSE}.
#'
#' @return A data.frame with one row per segment (rownames = \code{rownames(pbas)}), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{pv}{Raw p‐value for the group term (two‐level factor).}
#'     \item{fdr}{Benjamini–Hochberg FDR (across all segments).}
#'     \item{dpsi}{Difference in mean PSI between \code{conditions[2]} and \code{conditions[1]}.}
#'   }
#'
#' @details
#' 1. Converts \code{groupby} into a single grouping vector \code{group_factor} via \code{get_groupby_factor()}.
#' 2. Filters \code{group_factor} and \code{pbas} to keep only columns corresponding to the two \code{conditions}.
#' 3. Calls \code{fit_as_glm()} with \code{formula = x ~ group} to get raw p‐values for each segment.
#' 4. Adjusts p‐values (Benjamini–Hochberg) into \code{fdr}.
#' 5. Computes delta‐PSI (\code{dpsi}) per segment by pseudobulking and calling \code{calc_psi()}, then taking
#'    mean PSI per condition and subtracting.
#'
#' @seealso \code{\link{fit_as_glm}}, \code{\link{get_groupby_factor}}, \code{\link{calc_psi}}, \code{\link{pseudobulk}}, \code{\link{test_all_groups_as}}
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

  # 4. Adjust p‐values (Benjamini–Hochberg)
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
#' For a \code{SummarizedExperiment} of splicing segments with assays “psi” (percent spliced‐in),
#'  this function conducts, for each unique group label, a one‐vs‐rest quasi‐binomial GLM test
#'  (i.e., comparing that group to all other samples).
#' It returns p‐values and delta‐PSI matrices (segments × groups).
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays “i”, “e”, and “psi”.
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length \code{ncol(pbas)} that directly gives a group label for each column, or
#'     \item One or more column names in \code{colData(pbas)}, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See \code{\link{get_groupby_factor}} for details.
#' @param parallel Logical; if \code{TRUE}, fit per‐segment GLMs in parallel (requires a registered \code{plyr} backend). Default \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, prints each group being tested. Default \code{FALSE}.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{pv}{Matrix: Raw p‐values, with rows = segments, columns = group labels}
#'     \item{fdr}{Matrix: Benjamini–Hochberg FDR (across all segments)}
#'     \item{dpsi}{Matrix: Delta‐PSI (mean PSI(group) – mean PSI(others)) per segment}
#'   }
#'   Row names of each matrix are \code{rownames(pbas)}; column names are the unique groups.
#'
#' @details
#' 1. Builds a grouping factor via \code{get_groupby_factor(pbas, groupby)}.
#' 2. For each unique label \code{celltype} in that factor:
#'    \enumerate{
#'      \item Creates a binary factor \code{f} where \code{f == TRUE} if sample’s group == \code{celltype}, FALSE otherwise.
#'      \item Calls \code{fit_as_glm()} with \code{formula = x ~ f} and \code{data_terms = list(f = f)}, requesting p‐values (\code{return_pv = TRUE}).
#'      \item Extracts the “group” p‐value column (corresponding to \code{f}) for each segment and stores in \code{pv[, celltype]}.
#'      \item Computes \code{dpsi[, celltype]} as \code{mean(psi[segment, f]) – mean(psi[segment, !f])}.
#'    }
#' 3. Adjusts each column of \code{pv} by Benjamini–Hochberg into \code{fdr[, celltype]}.
#'
#' @seealso \code{\link{fit_as_glm}}, \code{\link{get_groupby_factor}}, \code{\link{calc_psi}}
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

  # 2. Initialize result matrices
  pv_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )
  dpsi_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )

  # 3. For each group, perform one‐vs‐rest test
  for (celltype in unique_groups) {
    if (verbose) {
      message("Testing group: ", celltype)
    }
    # Binary factor: TRUE if sample belongs to this group
    f <- group_factor == celltype

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
    pv_mat[, celltype] <- glm_res[, "group"]

    # Compute dpsi: for each segment, mean PSI in 'celltype' minus mean PSI in others
    psi_vals <- SummarizedExperiment::assay(pbas, "psi")

    # For rows with all-NA in either subset, result is NA automatically
    mean_ct <- rowMeans(psi_vals[, f, drop = FALSE], na.rm = TRUE)
    mean_other <- rowMeans(psi_vals[, !f, drop = FALSE], na.rm = TRUE)
    dpsi_mat[, celltype] <- mean_ct - mean_other
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
#' From a data.frame of per‐segment results comparing all groups simultaneously (e.g., output of \code{test_all_groups_as()}),
#'   this function selects only those segments that pass both a group‐level FDR threshold and a minimum delta‐PSI threshold.
#' It returns a data.frame with one row per qualifying segment.
#'
#' @param all_celltype_df A data.frame (rownames = segment IDs), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the \code{group} term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from \code{get_dpsi()}).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between \code{high_state} and \code{low_state}.}
#'   }
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#'
#' @return A data.frame with one row per segment satisfying \code{group_fdr < fdr_thr} and \code{dpsi > dpsi_thr}.
#'   Columns in the returned data.frame:
#'   \describe{
#'     \item{pv}{Raw p‐value from the likelihood‐ratio test for the group effect (copied from \code{all_celltype_df$group}).}
#'     \item{fdr}{Benjamini‐Hochberg adjusted FDR (copied from \code{all_celltype_df$group_fdr}).}
#'     \item{dpsi}{Delta‐PSI (copied from \code{all_celltype_df$dpsi}).}
#'     \item{seg_id}{Segment ID (rownames of \code{all_celltype_df}).}
#'     \item{group}{Group label with highest mean PSI (copied from \code{all_celltype_df$high_state}).}
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


#' Select marker segments per cell type
#'
#' Filters markers by false discovery rate (FDR) and absolute delta PSI (dPSI),
#' then selects the top n markers per cell type, and removing duplicates
#' (keeping only the strongest marker per segment across all cell_type).
#'
#' @param segment_stats A list containing numeric matrices `pv`, `fdr`, and `dpsi`
#'                      with rows = segments and cols = cell types.
#' @param n Integer: maximum number of markers to return per cell type (default: 5).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#' @param clean_duplicates Logical: if TRUE, keep only one entry per segment
#'                         (the one with highest |dpsi|) across all cell types (default: TRUE).
#' @return A data.frame with columns `pv`, `fdr`, `dpsi`, `seg_id`, and `group`(cell type),
#'         sorted by decreasing |dpsi|.
#' @examples
#' # Load example data
#' data(pbasf)
#'
#' # Select top 10 markers with stricter thresholds:
#' df <- select_markers(pbasf@metadata$markers, n = 10, fdr_thr = 0.01, dpsi_thr = 0.2)
#' head(df)
#' @export
select_markers <- function(segment_stats, n = 5, fdr_thr = 0.05, dpsi_thr = 0.1, clean_duplicates = TRUE) {
  # dataframe to store selected markers
  selected <- NULL

  # Loop over each cell_type (column in the FDR matrix)
  for (cell_type in colnames(segment_stats$fdr)) {
    # Filter segments based on FDR and dPSI thresholds
    is_significant <-
      (segment_stats$fdr[, cell_type] <= fdr_thr) &
        (abs(segment_stats$dpsi[, cell_type]) >= dpsi_thr)
    is_significant[is.na(is_significant)] <- FALSE # Treat NA values as non-significant

    # Skip if no segments pass for this cell_type
    if (!any(is_significant)) next

    # Build dataframe of all significant segments for this cell_type
    top <- data.frame(
      pv = segment_stats$pv[is_significant, cell_type],
      fdr = segment_stats$fdr[is_significant, cell_type],
      dpsi = segment_stats$dpsi[is_significant, cell_type],
      seg_id = rownames(segment_stats$pv)[is_significant],
      group = cell_type
    )

    # Order by decending |dPSI| and select top n segments
    keep_idx <- order(abs(top$dpsi), decreasing = TRUE)[seq_len(min(n, nrow(top)))]
    top <- top[keep_idx, ]

    # Append to output dataframe
    selected <- rbind(selected, top)
  }

  # Reset row name
  rownames(selected) <- NULL

  # (Optional) remove duplicate segment IDs (keep the highest |dpsi|)
  if (clean_duplicates) {
    selected <- lapply(split(selected, selected$seg_id), function(x) x[order(abs(x$dpsi), decreasing = TRUE)[1], ])
    selected <- do.call(rbind, selected)
  }

  # Final sort by descending |dpsi|
  selected <- selected[order(abs(selected$dpsi), decreasing = TRUE), ]

  # Return dataframe
  selected
}
