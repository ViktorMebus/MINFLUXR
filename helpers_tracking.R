## ─────────────────────────────────────────────────────────────────────────────
## helpers_tracking.R  —  tracking-mode-only helper functions
## ─────────────────────────────────────────────────────────────────────────────
# NOTE: jsonlite is NOT library()'d — use jsonlite:: explicitly everywhere.
# NOTE: All input/output IDs in the tracking tab are prefixed with trk_
#       to avoid collisions with SMLM tab IDs.

tid_color_map <- function(tids, palette = "Set3", seed = 123) {
  lev <- sort(unique(as.character(tids)))
  n   <- length(lev)
  info    <- RColorBrewer::brewer.pal.info
  if (!palette %in% rownames(info)) palette <- "Set3"
  maxcols <- info[palette, "maxcolors"]
  base    <- RColorBrewer::brewer.pal(max(3, min(maxcols, n)), palette)
  cols    <- if (n > length(base)) grDevices::colorRampPalette(base)(n) else base[seq_len(n)]
  set.seed(seed); cols <- sample(cols, n)
  stats::setNames(as.character(cols), lev)
}

convert_to_qs <- function(in_path, out_path = sub("\\.csv(\\.gz)?$", ".qs", in_path)) {
  dt <- data.table::fread(in_path, showProgress = TRUE)
  dt[, tid := as.character(tid)]
  dt[, c("x","y","tim") := .(as.numeric(x), as.numeric(y), as.numeric(tim))]
  qs::qsave(as.data.frame(dt), out_path, preset = "high")
  out_path
}

split_by_tid <- function(in_path, out_dir, tids_per_chunk = 500, fmt = c("qs","fst")) {
  fmt <- match.arg(fmt)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dt <- data.table::fread(in_path, showProgress = TRUE)
  dt[, tid := as.character(tid)]
  tids   <- sort(unique(dt$tid))
  chunks <- split(tids, ceiling(seq_along(tids)/tids_per_chunk))
  manifest <- data.frame(chunk = seq_along(chunks), n_tids = lengths(chunks), file = NA_character_)
  for (i in seq_along(chunks)) {
    sub <- dt[tid %in% chunks[[i]]]
    f <- file.path(out_dir, sprintf("chunk_%03d.%s", i, fmt))
    if (fmt == "qs") {
      qs::qsave(as.data.frame(sub), f, preset = "high")
    } else {
      fst::write_fst(sub, f, compress = 50)
    }
    manifest$file[i] <- f
  }
  write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
}

prep_tracking <- function(df, normalize_per_tid = TRUE) {
  stopifnot(all(c("x","y","tid","tim") %in% names(df)))
  df %>%
    dplyr::select(tid, x, y, tim, dplyr::everything()) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(tim)) %>%
    dplyr::arrange(tid, tim) %>%
    dplyr::group_by(tid) %>%
    dplyr::mutate(tim0 = if (normalize_per_tid) tim - dplyr::first(tim) else tim) %>%
    dplyr::ungroup()
}

na_lines <- function(x, y) list(x = c(x, NA_real_), y = c(y, NA_real_))

nearest_idx <- function(tt, t) {
  if (!length(tt)) return(NA_integer_)
  i <- base::findInterval(t, tt)
  if (i <= 0) return(1L)
  if (i >= length(tt)) return(length(tt))
  if (abs(tt[i+1] - t) < abs(tt[i] - t)) i + 1L else i
}

decimate_per_tid <- function(df, max_pts = 400, tcol = "tim0", extra_cols = character()) {
  has_z <- "z" %in% names(df)
  extra_cols <- intersect(extra_cols, names(df))     # carry these through unchanged (e.g. per-point colour)
  outx <- outy <- outt <- outid <- NULL
  outz <- NULL
  out_extra <- stats::setNames(vector("list", length(extra_cols)), extra_cols)
  sel_cols <- c("tid","x","y", if (has_z) "z", tcol, extra_cols)
  df_ord <- df[order(df$tid, df[[tcol]]), sel_cols, drop = FALSE]
  sp <- split(df_ord, df_ord$tid)
  for (nm in names(sp)) {
    tb <- sp[[nm]]; n <- nrow(tb)
    keep  <- if (n <= max_pts) seq_len(n) else unique(round(seq(1, n, length.out = max_pts)))
    outid <- c(outid, rep(nm, length(keep)))
    outx  <- c(outx,  tb$x[keep])
    outy  <- c(outy,  tb$y[keep])
    if (has_z) outz <- c(outz, tb$z[keep])
    outt  <- c(outt,  tb[[tcol]][keep])
    for (ec in extra_cols) out_extra[[ec]] <- c(out_extra[[ec]], tb[[ec]][keep])
  }
  out <- data.frame(tid = outid, x = outx, y = outy, t = outt,
                    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  if (has_z) out$z <- outz
  for (ec in extra_cols) out[[ec]] <- out_extra[[ec]]
  out
}

split_by_tid_to_dir <- function(df, out_dir,
                                tids_per_chunk = TIDS_PER_CHUNK_DEFAULT,
                                fmt = CHUNK_FORMAT) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  df$tid <- as.character(df$tid)
  df <- df[order(df$tid, df$tim), , drop = FALSE]
  tids <- sort(unique(df$tid))
  if (!length(tids)) stop("No TIDs to split.")
  chunks <- split(tids, ceiling(seq_along(tids) / max(1L, as.integer(tids_per_chunk))))

  manifest <- data.frame(chunk  = seq_along(chunks),
                         n_tids = lengths(chunks),
                         file   = NA_character_)

  withProgress(message = "Chunking large dataset...", value = 0, {
    for (i in seq_along(chunks)) {
      incProgress(1/length(chunks), detail = paste("Chunk", i, "of", length(chunks)))
      sub <- df[df$tid %in% chunks[[i]], , drop = FALSE]
      f <- file.path(out_dir, sprintf("chunk_%03d.%s", i, fmt))
      if (identical(fmt, "qs")) {
        qs::qsave(sub, f, preset = "high")
      } else if (identical(fmt, "fst")) {
        fst::write_fst(sub, f, compress = 50)
      } else {
        stop("Unsupported chunk format: ", fmt)
      }
      manifest$file[i] <- f
    }
  })
  write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
  invisible(file.path(out_dir, "manifest.csv"))
}

even_idx <- function(n, k) {
  k <- as.integer(k)
  if (!is.finite(k) || k <= 0 || n <= 0) return(integer(0))
  unique(round(seq(1, n, length.out = min(k, n))))
}
