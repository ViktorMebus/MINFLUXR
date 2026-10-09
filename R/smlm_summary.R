## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_summary.R  -  Per-file summaries for group comparisons
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   One row per file per metric, tagged with the file's dataset group, so
##   condition comparisons use the cell (file) as the unit of replication.
##
## CONTENTS (in file order)
##   per_file_summary
##
## Part of the SMLM helpers; loaded by helpers_smlm.R (see its header for the data
## hierarchy and units). Plain functions, no Shiny: testable in the console.
## ─────────────────────────────────────────────────────────────────────────────

# =============================================================================
# per_file_summary() - ONE ROW PER FILE PER METRIC, FOR BETWEEN-GROUP STATISTICS
# =============================================================================
#
# GOAL
#   Collapse each distribution analysis to per-file summary statistics, tagged
#   with the file's dataset group, so group comparisons can use the CELL (file)
#   as the unit of replication.
#
# WHY IT EXISTS
#   In group mode the distribution panels (Batch NN, trace IDs per Median
#   Cluster, RoG, crowding, within-cluster NN vs size, cross-channel NN) pool
#   every unit from every file in a group into one histogram. Units within a cell
#   are not independent (same labelling, same drift, same biology), so a test run
#   on pooled units is pseudo-replicated: its n is the number of distances or
#   clusters when the number of independent samples is the number of cells. The
#   spatial panels already average per-file curves; this gives the distribution
#   panels the same footing. The app runs no test itself; the export feeds one.
#
# ARGS
#   results   : named list, one entry per panel, each list(df = <stored result
#               data.frame with a `dataset` column>, metrics = <character vector of
#               numeric column names to summarise>, unit = <what one row is>).
#               NULL / empty entries are skipped (panel not run).
#   group_map : named character vector file -> group (file_group_map()); may be empty.
#   files     : loaded file names, used to resolve faceted tags ("file · channel")
#               back to their file so the group lookup works.
#
# RETURNS
#   data.frame(panel, metric, unit, dataset, file, group, n, mean, sd, median, q25, q75)
#   or an empty data.frame when nothing has been run. `dataset` keeps the facet
#   tag (a channel split gives one row per channel per file).
per_file_summary <- function(results, group_map = character(0), files = character(0)) {
  empty <- data.frame(panel = character(), metric = character(), unit = character(),
                      dataset = character(), file = character(), group = character(),
                      n = integer(), mean = numeric(), sd = numeric(), median = numeric(),
                      q25 = numeric(), q75 = numeric(), stringsAsFactors = FALSE)
  if (!length(results)) return(empty)
  # Resolve a dataset tag to its source file: exact match, else the longest loaded
  # file name the tag starts with (facets are named "<file> · <channel>").
  file_of <- function(tag) {
    if (tag %in% files) return(tag)
    hit <- files[startsWith(tag, paste0(files, " "))]
    if (length(hit)) hit[which.max(nchar(hit))] else tag
  }
  rows <- list()
  for (panel in names(results)) {
    ent <- results[[panel]]
    df  <- ent$df
    if (is.null(df) || !is.data.frame(df) || !nrow(df) || !"dataset" %in% names(df)) next
    for (metric in intersect(ent$metrics, names(df))) {
      for (tag in unique(as.character(df$dataset))) {
        v <- suppressWarnings(as.numeric(df[[metric]][as.character(df$dataset) == tag]))
        v <- v[is.finite(v)]                              # summaries on finite values only
        if (!length(v)) next
        f <- file_of(tag)
        g <- if (length(group_map) && f %in% names(group_map)) unname(group_map[[f]]) else NA_character_
        q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(
          panel = panel, metric = metric, unit = ent$unit %||% NA_character_,
          dataset = tag, file = f, group = g, n = length(v),
          mean = mean(v), sd = if (length(v) > 1L) stats::sd(v) else NA_real_,
          median = q[2], q25 = q[1], q75 = q[3], stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out[order(out$panel, out$metric, out$group, out$dataset), , drop = FALSE]
}

