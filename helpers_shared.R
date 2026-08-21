## -----------------------------------------------------------------------------
## helpers_shared.R  -  helpers used by BOTH the SMLM and Tracking tabs
## -----------------------------------------------------------------------------
##
## PURPOSE
##   This is the app's FRONT DOOR. Almost all of it is the import layer: it
##   turns the several file formats an Abberior MINFLUX microscope can emit into
##   ONE canonical data.frame that the rest of the app can rely on. If data looks
##   wrong everywhere in the app, the bug is probably in here.
##
## THE CANONICAL DATA.FRAME (the contract everything downstream assumes)
##   x, y, z ... position in NANOMETRES (files store metres - converted on read)
##   tid ....... trace ID: one tracked molecule / one blinking event. Character.
##   tim ....... timestamp. SECONDS on disk; the tracking tab multiplies by 1000
##               and works in ms from then on. The SMLM tab keeps seconds.
##   itr ....... MINFLUX iteration index (the localization refinement stage)
##   efo ....... effective photon rate at offset (Hz) - equivalent to eco
##   cfr ....... centre-frequency ratio -> localization quality metric
##   eco ....... effective photon count at offset - periphieral photon count
##   dcr ....... detector-channel ratio -> the colour axis for 2-colour work
##   Optional: dwell, fluo, lnc, fbg, ...
##
## WHAT LIVES HERE (in reading order)
##   1. Tiny utilities .......... `%||%`, get_ext()
##   2. Labelling presets ....... minflux_label_presets()  <- SINGLE SOURCE OF
##                                TRUTH shared by the Simulator's forward model
##                                and the SMLM tab's inverse GT estimator
##   3. Readers, one per format:
##        CSV  read_minflux_csv()          (pyMINFLUX export)
##        NPY  read_npy() + read_minflux_npy_as_df()   (pure R, no Python)
##        JSON read_minflux_json_as_df()   (Imspector export)
##      All of them converge on .minflux_finalize_df(), which applies the shared
##      post-processing so every format ends up identical downstream.
##   4. Format normalisation .... canonicalize_minflux_df(),
##                                collapse_minflux_iterations()
##   5. Image overlay ........... make_scale_bar(), apply_background_lut(),
##                                auto_align_overlay()
##   6. Geometry ................ .pnpoly_vec(), filter_by_cell_boundary()
##
## WHY A HAND-WRITTEN .npy READER
##   Abberior also exports .npy, which is far smaller than JSON. Rather than take
##   a Python/reticulate dependency for one file format, read_npy() decodes the
##   NumPy container directly. It handles the structured (record) dtypes those
##   exports use, including the float16 fields R's readBin cannot read natively.
##
## CONVENTION: jsonlite is NOT attached - always call jsonlite:: explicitly, or
## its `validate` shadows shiny::validate.
## -----------------------------------------------------------------------------

# Null-coalescing operator: `a %||% b` returns a unless it is NULL, else b.
# Used constantly to give Shiny inputs a default before they have rendered.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# σ with a real subscript, e.g. sigloc() -> σ_loc, sigloc("loc-wECO") -> σ_loc-wECO.
# Returns an HTML object, so it renders as a subscript anywhere that accepts tags
# (helpText, tags$*, labels, actionButton). Plain character strings are ESCAPED by
# Shiny, which is why the symbol must be passed as its own argument rather than
# embedded in a string: helpText("use ", sigloc(), " for the floor").
sigloc <- function(sub = "loc") htmltools::HTML(sprintf("&sigma;<sub>%s</sub>", sub))

# Collapsible per-method infobox: what the method is, where it is from, how it is
# adapted for MINFLUX (irregular Δt / short tracks / nm precision / blinking), and
# the honest caution. `rating` is a 🟢/🟡/🟠 chip:
#   green  = standard method, well-suited to MINFLUX as used here
#   amber  = standard method, but a MINFLUX property strains an assumption
#   orange = bespoke to this tool, no published MINFLUX equivalent (exploratory)
# Lives here (not in a ui_*.R) because helpers_shared.R is sourced first, so both
# ui_smlm.R and ui_tracking.R can call it at source time.
.method_box <- function(rating, source, adaptation, caution) {
  chip <- switch(rating,
    green  = tags$span("🟢 well-adapted",  style = "color:#1a8a1a;"),
    amber  = tags$span("🟡 adapted, with caution", style = "color:#b26a00;"),
    orange = tags$span("🟠 exploratory for MINFLUX", style = "color:#e67e22;"),
    tags$span(rating))
  tags$details(style = "margin:4px 0 10px; font-size:0.82em;",
    tags$summary(style = "cursor:pointer; color:#2c3e50; font-weight:600;",
                 "ℹ Method & MINFLUX adaptation  ", chip),
    tags$div(style = "padding:6px 4px 2px 10px; color:#555; line-height:1.5; border-left:3px solid #cfe0ee;",
      tags$p(style = "margin:2px 0;", tags$b("Origin: "), source),
      tags$p(style = "margin:2px 0;", tags$b("MINFLUX adaptation: "), adaptation),
      tags$p(style = "margin:2px 0; color:#7a5c00;", tags$b("Caution: "), caution)))
}

# Collapsed detail block. The visible text under a control should be one or two
# sentences saying what it IS; everything else - the reasoning, the caveats, the
# numbers behind a default - goes in here, closed by default.
#
# WHY THIS EXISTS: the Analysis tab had grown ~22,000 characters of always-visible
# helpText. Every sentence was true and most were worth keeping, but a panel that
# opens with eight lines of caveats reads as a warning label rather than a tool,
# and the important sentence gets lost among the merely accurate ones. Collapsing
# is not deleting: nothing here is dropped, it is one click away.
#
# Use .method_box() for the Origin / MINFLUX-adaptation / Caution triple, and
# .detail() for everything else that does not fit in two sentences.
.detail <- function(..., label = "More detail") {
  tags$details(style = "margin:2px 0 8px; font-size:0.8em;",
    tags$summary(style = "cursor:pointer; color:#5a6b7c;", label),
    tags$div(style = "padding:4px 4px 2px 10px; color:#666; line-height:1.5; border-left:2px solid #dfe6ec;",
             ...))
}

# Lower-case file extension of a path/filename ("Sample.NPY" -> "npy").
# Used to dispatch an upload to the right reader.
get_ext <- function(filename) {
  tolower(sub("^.*\\.", "", filename))
}

# -----------------------------------------------------------------------------
# Labeling-method presets - SINGLE SOURCE OF TRUTH.
# Used by the Simulator (forward model: GT -> traces) and by the SMLM tab's
# inverse GT-estimator (traces -> candidate GT). Keeping one table means the
# estimator inverts exactly the numbers the simulator generated with.
#
#   n_fluor       mean fluorophores per labeled molecule (Poisson-drawn, min 1)
#   link_nm       mean distance label-complex centre -> molecule centre (nm).
#                 This is the COMMON-MODE, per-molecule offset shared by every
#                 fluorophore/blink/localization of that molecule - the
#                 irreducible position floor for single-molecule GT estimation.
#   dye_spread_nm Gaussian SD of individual dye positions within the complex
#   footprint_nm  approximate excluded-zone diameter around each labeled protein
#   blink_lambda  mean EXTRA TID appearances per fluorophore (dye re-activation)
# -----------------------------------------------------------------------------
minflux_label_presets <- function() {
  list(
    # -- Covalent / minimal-linkage tags --------------------------------------
    "Direct tag (SNAP / Halo / Click, DOL=1)" = list(
      n_fluor = 1L, link_nm = 3.5, dye_spread_nm = 0,
      footprint_nm = 2,  blink_lambda = 0.2,
      note    = "1 fluorophore / molecule - ~3.5 nm linkage (SNAP/Halo tag domain)"
    ),
    "Nanobody (DOL ≈1)" = list(
      n_fluor = 1L, link_nm = 2.5, dye_spread_nm = 0,
      footprint_nm = 3,  blink_lambda = 0.2,
      note    = "1 fluorophore / molecule - ~2.5 nm (nanobody VHH domain)"
    ),
    # -- Primary monoclonal (1 antibody / epitope) ----------------------------
    "Primary monoclonal IgG (DOL ≈1)" = list(
      n_fluor = 1L, link_nm = 8, dye_spread_nm = 0,
      footprint_nm = 10, blink_lambda = 0.3,
      note    = "1 fluorophore / molecule - ~8 nm (IgG Fab arm to epitope)"
    ),
    "Primary monoclonal IgG (NHS, DOL ≈2–4)" = list(
      n_fluor = 3L, link_nm = 8, dye_spread_nm = 4,
      footprint_nm = 10, blink_lambda = 0.5,
      note    = "~3 fluorophores on one IgG body (~4 nm spread) - ~8 nm from target"
    ),
    # -- Primary polyclonal (3-5 antibodies / antigen) ------------------------
    "Primary polyclonal IgG (NHS, DOL ≈2–4)" = list(
      n_fluor = 9L, link_nm = 8, dye_spread_nm = 9,
      footprint_nm = 12, blink_lambda = 0.5,
      note    = "~9 fluorophores across 3–5 Ab bodies (~9 nm spread) - ~8 nm from target"
    ),
    # -- Monoclonal + secondary -----------------------------------------------
    "Primary monoclonal IgG + secondary nanobody" = list(
      n_fluor = 2L, link_nm = 10.5, dye_spread_nm = 3,
      footprint_nm = 10, blink_lambda = 0.4,
      note    = "~2 fluorophores, nanobody on primary IgG - ~10.5 nm from target"
    ),
    "Primary monoclonal IgG + secondary IgG" = list(
      n_fluor = 8L, link_nm = 20, dye_spread_nm = 12,
      footprint_nm = 16, blink_lambda = 0.6,
      note    = "~8 fluorophores across ~3 secondary IgG bodies (~12 nm spread) - ~20 nm from target"
    ),
    # -- Polyclonal + secondary -----------------------------------------------
    "Primary polyclonal IgG + secondary IgG" = list(
      n_fluor = 18L, link_nm = 20, dye_spread_nm = 18,
      footprint_nm = 20, blink_lambda = 0.8,
      note    = "~18 fluorophores, large complex (~18 nm spread) - ~20 nm from target"
    ),
    # -- Biotin-streptavidin --------------------------------------------------
    "Streptavidin (biotinylated Ab + SA-dye)" = list(
      n_fluor = 3L, link_nm = 15, dye_spread_nm = 3,
      footprint_nm = 8,  blink_lambda = 0.4,
      note    = "~3 fluorophores on SA tetramer (~3 nm spread) - ~15 nm from target"
    ),
    "Custom" = list(
      n_fluor = NA_integer_, link_nm = NA_real_,
      dye_spread_nm = NA_real_, footprint_nm = NA_real_,
      blink_lambda = NA_real_, note = ""
    )
  )
}

# Read MINFLUX CSV (pyMINFLUX format); uses data.table::fread when available
read_minflux_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- data.table::fread(path, data.table = FALSE)
  } else {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  df
}

# -----------------------------------------------------------------------------
# Pure-R NumPy (.npy) reader - no Python/reticulate dependency.
# Supports v1/v2/v3 headers and C-order (fortran_order = FALSE) arrays, including
# *structured* (record) dtypes with named fields and fixed-shape sub-arrays
# (e.g. loc = (3,) f8, dcr = (2,) f8) - exactly the layout Abberior MINFLUX
# "tabular" exports use. Nested struct fields (legacy m2205 itr blocks) are not
# supported and raise an informative error.
# -----------------------------------------------------------------------------

# Map a NumPy dtype token ("<f8", "|b1", "?", ">i4", "<u4", ...) to readBin args.
.npy_dtype_info <- function(tok) {
  endian <- "little"
  ch1 <- substr(tok, 1, 1)
  rest <- tok
  if (ch1 %in% c("<", ">", "=", "|")) {
    endian <- if (ch1 == ">") "big" else "little"
    rest   <- substr(tok, 2, nchar(tok))
  }
  if (rest == "?") return(list(what = "logical", size = 1L, signed = TRUE,  endian = endian))
  kind <- substr(rest, 1, 1)
  size <- suppressWarnings(as.integer(substr(rest, 2, nchar(rest))))
  if (is.na(size)) size <- 1L
  switch(kind,
    b = list(what = "logical", size = 1L,   signed = TRUE,  endian = endian),
    f = list(what = "double",  size = size,  signed = TRUE,  endian = endian),
    i = list(what = "integer", size = size,  signed = TRUE,  endian = endian),
    u = list(what = "integer", size = size,  signed = FALSE, endian = endian),
    stop(sprintf("Unsupported NumPy dtype kind '%s' in token '%s'", kind, tok))
  )
}

# Convert IEEE-754 half-precision (float16) bit patterns (uint16, 0..65535) to
# double. R's readBin has no 2-byte float type, and Abberior MINFLUX stores cfr
# and dcr as <f2, so we decode them by hand.
.half_to_double <- function(h) {
  sign <- bitwShiftR(h, 15L)
  exp  <- bitwAnd(bitwShiftR(h, 10L), 0x1FL)
  mant <- bitwAnd(h, 0x3FFL)
  val  <- numeric(length(h))
  norm <- exp > 0L & exp < 31L
  val[norm] <- (1 + mant[norm] / 1024) * 2^(exp[norm] - 15L)
  sub  <- exp == 0L
  val[sub]  <- (mant[sub] / 1024) * 2^(-14)
  neg <- sign == 1L
  val[neg] <- -val[neg]
  # Abberior writes half-precision Inf/NaN for skipped values (e.g. cfr when the
  # CFR check is disabled). Mirror the JSON reader, which nulls Inf/NaN -> NA.
  val[exp == 31L] <- NA_real_
  val
}

# Decode a column of `n` * `cnt` scalars from a raw vector `bytes` (already the
# tightly-packed bytes for this field across all records, record-major).
.npy_decode_bytes <- function(bytes, info, n_total) {
  if (info$what == "double" && info$size == 2L) {
    h <- readBin(bytes, "integer", n = n_total, size = 2L, signed = FALSE,
                 endian = info$endian)
    return(.half_to_double(h))
  }
  if (info$what == "integer" && info$size == 8L) {
    # readBin has no 8-byte integer; combine two 4-byte words -> double.
    w <- readBin(bytes, "integer", n = n_total * 2L, size = 4L,
                 signed = FALSE, endian = info$endian)
    lo <- w[c(TRUE, FALSE)]; hi <- w[c(FALSE, TRUE)]
    if (info$endian == "big") { tmp <- lo; lo <- hi; hi <- tmp }
    return(hi * 2^32 + lo)
  }
  if (info$what == "integer" && info$size %in% c(1L, 2L)) {
    return(readBin(bytes, "integer", n = n_total, size = info$size,
                   signed = info$signed, endian = info$endian))
  }
  if (info$what == "integer") {  # size 4
    # signed=FALSE only honoured for sizes 1/2; 4-byte unsigned read as signed
    # (fine for trace IDs / counts well under 2^31).
    return(readBin(bytes, "integer", n = n_total, size = info$size, endian = info$endian))
  }
  if (info$what == "double") {
    return(readBin(bytes, "double", n = n_total, size = info$size, endian = info$endian))
  }
  # logical (1 byte)
  as.logical(readBin(bytes, "integer", n = n_total, size = 1L, signed = FALSE,
                     endian = info$endian))
}

# Read an .npy file. Returns, for a structured dtype, a named list of columns
# (vectors for scalar fields, N x k matrices for sub-array fields). For a plain
# dtype, returns a vector or matrix. Attribute "shape" holds the array shape.
read_npy <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", n = 6L)
  if (!identical(as.integer(magic), c(0x93L, 0x4eL, 0x55L, 0x4dL, 0x50L, 0x59L)))
    stop("Not a .npy file (bad magic).")
  ver <- as.integer(readBin(con, "raw", n = 2L))
  if (ver[1] == 1L) {
    hlen <- readBin(con, "integer", n = 1L, size = 2L, signed = FALSE, endian = "little")
  } else {
    hlen <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
  }
  header <- rawToChar(readBin(con, "raw", n = hlen))

  # ---- Parse the Python-dict header (regex; the format is tightly specified) --
  fortran <- grepl("'fortran_order'\\s*:\\s*True", header)
  if (fortran) stop("Fortran-ordered .npy not supported (re-save C-ordered).")
  shp_raw <- sub(".*'shape'\\s*:\\s*\\(([^)]*)\\).*", "\\1", header)
  shape   <- suppressWarnings(as.integer(strsplit(trimws(shp_raw), "\\s*,\\s*")[[1]]))
  shape   <- shape[!is.na(shape)]
  n       <- if (length(shape)) shape[1] else 0L

  descr_blob <- sub(".*'descr'\\s*:\\s*", "", header)
  is_struct  <- grepl("^\\s*\\[", descr_blob)

  if (!is_struct) {
    # Plain dtype: descr is a single quoted token e.g. '<f8'
    tok  <- sub(".*?'([^']+)'.*", "\\1", descr_blob)
    info <- .npy_dtype_info(tok)
    total <- if (length(shape)) prod(shape) else 0L
    raw   <- readBin(con, "raw", n = total * info$size)
    vec   <- .npy_decode_bytes(raw, info, total)
    if (length(shape) >= 2) { dim(vec) <- rev(shape); vec <- t(vec) }  # C-order
    attr(vec, "shape") <- shape
    return(vec)
  }

  # ---- Structured dtype: parse each ('name','dtype'[, (shape)]) tuple ---------
  if (grepl("\\[\\s*\\(\\s*'[^']*'\\s*,\\s*\\[", descr_blob))
    stop("Nested structured .npy dtype (legacy MINFLUX m2205 format) not supported; export the tabular format.")
  rx <- "\\(\\s*'([^']*)'\\s*,\\s*'([^']*)'\\s*(?:,\\s*\\(([^)]*)\\))?\\s*\\)"
  m  <- gregexpr(rx, descr_blob, perl = TRUE)[[1]]
  if (m[1] == -1L) stop("Could not parse structured .npy dtype descriptor.")
  starts <- as.integer(m); lens <- attr(m, "match.length")
  fields <- lapply(seq_along(starts), function(i) {
    piece <- substr(descr_blob, starts[i], starts[i] + lens[i] - 1L)
    g <- regmatches(piece, regexec(rx, piece, perl = TRUE))[[1]]
    sub_shape <- if (nchar(g[4])) {
      v <- suppressWarnings(as.integer(strsplit(trimws(g[4]), "\\s*,\\s*")[[1]]))
      v[!is.na(v)]
    } else integer(0)
    info <- .npy_dtype_info(g[3])
    cnt  <- as.integer(max(1L, prod(sub_shape)))
    list(name = g[2], info = info, cnt = cnt,
         sub_shape = sub_shape, esz = as.integer(info$size * cnt))
  })

  itemsize <- sum(vapply(fields, `[[`, integer(1), "esz"))
  raw <- readBin(con, "raw", n = n * itemsize)
  if (length(raw) < n * itemsize) stop("Truncated .npy data block.")
  # Reshape to (itemsize x n): each column is one record's bytes.
  recs <- matrix(raw, nrow = itemsize, ncol = n)

  out <- vector("list", length(fields)); names(out) <- vapply(fields, `[[`, "", "name")
  off <- 0L
  for (i in seq_along(fields)) {
    f <- fields[[i]]
    rows <- (off + 1L):(off + f$esz)
    fb   <- as.raw(recs[rows, , drop = FALSE])  # column-major -> record-major bytes
    vals <- .npy_decode_bytes(fb, f$info, n * f$cnt)
    if (f$cnt > 1L) { dim(vals) <- c(f$cnt, n); vals <- t(vals) }  # N x cnt
    out[[f$name]] <- vals
    off <- off + f$esz
  }
  attr(out, "shape") <- shape
  out
}

# Canonical JSON reader for Abberior MINFLUX exports.
# Returns a data.frame with x/y/z (nm), tid, itr, tim, efo, eco, cfr, ecc, efc,
# fbg, gri, sta, bot, eot, fnl, sqi, thi, dcr1, dcr2, dcr, dwell.
# The tracking tab only uses x/y/tid/tim/efo/cfr/dcr - extra columns are ignored.
read_minflux_json_as_df <- function(path,
                                    keep_final_only = TRUE,
                                    keep_valid_only = TRUE,
                                    final_itr = NA_integer_,
                                    cfr_strategy = c("last_nonzero", "itr6", "none"),
                                    z_scale = 0.7,
                                    dwell_time_ms = 1.0,
                                    pool_dcr = FALSE) {

  cfr_strategy <- match.arg(cfr_strategy)

  safe_inc <- function(amount, detail = NULL) {
    try(shiny::incProgress(amount, detail = detail), silent = TRUE)
    invisible(NULL)
  }

  # ---- 1. Read file (single bulk read - much faster than readLines+paste) ----
  safe_inc(0.04, "Reading JSON file")
  txt <- tryCatch({
    raw <- readBin(path, "raw", n = file.info(path)$size)
    rawToChar(raw)
  }, error = function(e) {
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })

  # ---- 2. Patch non-standard JSON tokens (NaN / Infinity) -------------------
  # Signed forms first so the bare pass doesn't double-replace already-nulled
  # values.  stringi's fixed replacer is 3-5x faster than gsub on large strings
  # because it skips regex compilation and uses a direct byte-search algorithm.
  safe_inc(0.04, "Patching NaN/Inf tokens")
  if (requireNamespace("stringi", quietly = TRUE)) {
    txt <- stringi::stri_replace_all_fixed(
      txt,
      c("-Infinity", "Infinity", "-NaN", "NaN"),
      c("null",      "null",     "null", "null"),
      vectorize_all = FALSE
    )
  } else {
    txt <- gsub("-Infinity",       "null", txt, fixed = TRUE)
    txt <- gsub("\\bInfinity\\b", "null", txt)
    txt <- gsub("-NaN",            "null", txt, fixed = TRUE)
    txt <- gsub("\\bNaN\\b",      "null", txt)
  }

  # ---- 3. Parse - try simplifyVector = TRUE, fall back without re-reading ----
  safe_inc(0.22, "Parsing JSON")
  parsed <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = TRUE),
    error = function(e) {
      message("[MINFLUXR] simplifyVector=TRUE failed (", conditionMessage(e), ") - using list parse")
      NULL
    }
  )
  # txt is kept alive here so the slow path can reuse it without re-reading.

  # Helper: extract a flat numeric/integer/logical column from the parsed df,
  # guarding against unexpected matrix or list columns.
  .pull_num <- function(src, col, n) {
    v <- src[[col]]
    if (!is.null(v) && !is.list(v) && !is.matrix(v)) suppressWarnings(as.numeric(v))
    else rep(NA_real_, n)
  }
  .pull_int <- function(src, col, n) {
    v <- src[[col]]
    if (!is.null(v) && !is.list(v) && !is.matrix(v)) suppressWarnings(as.integer(v))
    else rep(NA_integer_, n)
  }
  .pull_lgl <- function(src, col, n) {
    v <- src[[col]]
    if (!is.null(v) && is.logical(v)) v else rep(NA, n)
  }

  # Determine whether the simplified parse produced a usable data frame.
  # Note: is.list(data.frame) is TRUE in R, so we cannot use !is.list() to
  # detect matrix/data.frame columns - we must check both explicitly.
  coord_col     <- if ("loc" %in% names(parsed)) "loc" else
                   if ("lnc" %in% names(parsed)) "lnc" else NULL
  # Gate: only require a data.frame with core columns and a coord column.
  # dcr being a list column is fine - we handle both forms in extraction below.
  use_simplified <- is.data.frame(parsed) &&
    all(c("tid", "tim") %in% names(parsed)) &&
    !is.null(coord_col)

  if (use_simplified) {
    # ---- FAST PATH (simplified data frame, all iteration counts) -------------
    safe_inc(0.02, "Creating simplified data frame...")
    src <- parsed; rm(parsed, txt)
    n   <- nrow(src)

    # vld filter
    if (isTRUE(keep_valid_only) && "vld" %in% names(src) && is.logical(src$vld)) {
      src <- src[!is.na(src$vld) & src$vld, , drop = FALSE]
      n   <- nrow(src)
    }
    if (!n) return(data.frame())

    # Coordinates - handle matrix/data.frame or list column
    col_coord <- src[[coord_col]]
    if (is.matrix(col_coord) || is.data.frame(col_coord)) {
      cm <- as.matrix(col_coord)
      x  <- if (ncol(cm) >= 1) cm[, 1] else rep(NA_real_, n)
      y  <- if (ncol(cm) >= 2) cm[, 2] else rep(NA_real_, n)
      z  <- if (ncol(cm) >= 3) cm[, 3] else rep(NA_real_, n)
    } else if (is.list(col_coord)) {
      x <- vapply(col_coord, function(v) if (length(v) >= 1) as.numeric(v[[1]]) else NA_real_, numeric(1))
      y <- vapply(col_coord, function(v) if (length(v) >= 2) as.numeric(v[[2]]) else NA_real_, numeric(1))
      z <- vapply(col_coord, function(v) if (length(v) >= 3) as.numeric(v[[3]]) else NA_real_, numeric(1))
    } else if (all(c("x", "y") %in% names(src))) {
      x <- as.numeric(src$x); y <- as.numeric(src$y); z <- rep(NA_real_, n)
    } else {
      return(data.frame())
    }
    rm(col_coord)

    # DCR - handle matrix/data.frame or list column (jsonlite produces list for
    # fixed-length inner arrays in some versions).
    # dcr[:, 1] == 1 - dcr[:, 0] in Abberior files - redundant, so only dcr1 is kept.
    col_dcr <- src[["dcr"]]
    if (is.matrix(col_dcr) || is.data.frame(col_dcr)) {
      dm   <- as.matrix(col_dcr)
      dcr1 <- if (ncol(dm) >= 1) dm[, 1] else rep(NA_real_, n)
    } else if (is.list(col_dcr)) {
      dcr1 <- vapply(col_dcr, function(v) if (length(v) >= 1) as.numeric(v[[1]]) else NA_real_, numeric(1))
    } else {
      dcr1 <- rep(NA_real_, n)
    }
    rm(col_dcr)

    safe_inc(0.06, "Building table...")
    df_all <- data.frame(
      x    = x,   y    = y,   z    = z,
      tid  = .pull_int(src, "tid", n),
      itr  = .pull_int(src, "itr", n),
      tim  = .pull_num(src, "tim", n),
      efo  = .pull_num(src, "efo", n),
      eco  = .pull_num(src, "eco", n),
      cfr  = .pull_num(src, "cfr", n),
      ecc  = .pull_num(src, "ecc", n),
      efc  = .pull_num(src, "efc", n),
      fbg  = .pull_num(src, "fbg", n),
      gri  = .pull_int(src, "gri", n),
      sta  = .pull_num(src, "sta", n),
      bot  = .pull_lgl(src, "bot", n),
      eot  = .pull_lgl(src, "eot", n),
      fnl  = .pull_lgl(src, "fnl", n),
      sqi  = .pull_num(src, "sqi", n),
      thi  = .pull_num(src, "thi", n),
      dcr1 = dcr1,
      stringsAsFactors = FALSE
    )
    rm(src, x, y, z, dcr1)

  } else {
    # ---- SLOW PATH FALLBACK (list-of-lists) ------------------------------------
    safe_inc(0.00, "List parse fallback: building data frame from list of records...")
    if (!is.null(parsed)) rm(parsed)

    # txt is still in scope - reuse the already-read-and-patched string.
    dat <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    rm(txt)
    if (!length(dat)) return(data.frame())

    if (isTRUE(keep_valid_only) && !is.null(dat[[1]]$vld)) {
      keep <- vapply(dat, function(r) isTRUE(r$vld), logical(1))
      dat  <- dat[keep]
      if (!length(dat)) return(data.frame())
    }

    `%||%` <- function(a, b) if (!is.null(a)) a else b
    pull_num <- function(key) vapply(dat, function(r) suppressWarnings(as.numeric(r[[key]] %||% NA)), numeric(1))
    pull_int <- function(key) vapply(dat, function(r) suppressWarnings(as.integer(r[[key]] %||% NA)), integer(1))
    pull_lgl <- function(key) vapply(dat, function(r) isTRUE(r[[key]]), logical(1))

    coord_key <- if (!is.null(dat[[1]]$loc)) "loc" else "lnc"
    coord_mat <- do.call(rbind, lapply(dat, function(r) {
      l <- r[[coord_key]]
      if (is.null(l)) return(c(NA_real_, NA_real_, NA_real_))
      l <- as.numeric(l)
      c(l, rep(NA_real_, max(0L, 3L - length(l))))[1:3]
    }))
    dcr_vec <- vapply(dat, function(r) {
      d <- r$dcr
      if (is.null(d)) NA_real_ else as.numeric(d[[1]])
    }, numeric(1))

    df_all <- data.frame(
      x    = coord_mat[, 1], y = coord_mat[, 2], z = coord_mat[, 3],
      tid  = pull_int("tid"), itr = pull_int("itr"), tim = pull_num("tim"),
      efo  = pull_num("efo"), eco = pull_num("eco"), cfr = pull_num("cfr"),
      ecc  = pull_num("ecc"), efc = pull_num("efc"), fbg = pull_num("fbg"),
      gri  = pull_int("gri"), sta = pull_int("sta"),
      bot  = pull_lgl("bot"), eot = pull_lgl("eot"), fnl = pull_lgl("fnl"),
      sqi  = pull_int("sqi"), thi = pull_int("thi"),
      dcr1 = dcr_vec,
      stringsAsFactors = FALSE
    )
    rm(dat, coord_mat, dcr_vec)
  }

  # ---- Shared post-processing (JSON + .npy paths converge here) -------------
  .minflux_finalize_df(
    df_all,
    keep_final_only = keep_final_only,
    final_itr       = final_itr,
    cfr_strategy    = cfr_strategy,
    z_scale         = z_scale,
    dwell_time_ms   = dwell_time_ms,
    pool_dcr        = pool_dcr,
    use_fnl         = FALSE,
    safe_inc        = safe_inc
  )
}

# Shared post-processing for a raw per-iteration MINFLUX table (`df_all`, one row
# per iteration with columns x/y/z [m], tid, itr, tim, efo, eco, cfr, ecc, efc,
# dcr1, ...). Used by both the JSON and .npy readers. Performs: cfr fallback,
# x/y NA drop, optional ECO-weighted DCR pooling, per-TID cfr/ecc/efc carry-
# forward, final-iteration selection, dwell computation and m -> nm scaling.
# `use_fnl = TRUE` selects the final localization via the boolean `fnl` flag
# (Abberior .npy) instead of itr == max(itr) (JSON).
.minflux_finalize_df <- function(df_all,
                                 keep_final_only = TRUE,
                                 final_itr       = NA_integer_,
                                 cfr_strategy    = "last_nonzero",
                                 z_scale         = 0.7,
                                 dwell_time_ms   = 1.0,
                                 pool_dcr        = FALSE,
                                 use_fnl         = FALSE,
                                 safe_inc        = function(...) invisible(NULL)) {

  if (!"dcr1" %in% names(df_all)) df_all$dcr1 <- NA_real_
  df_all$dcr <- df_all$dcr1  # dcr[0] in pyMINFLUX (0-indexed)

  # Some acquisitions deliberately skip the CFR (centre-frequency-ratio) check to
  # keep cycle times short, so the file carries no cfr field and it parses as
  # all-NA. Treat a wholly-absent cfr as 0 throughout: the column then exists for
  # every downstream consumer and its filter becomes a harmless no-op (range
  # [0, 0]) instead of an undefined / NA-driven filter.
  if (!("cfr" %in% names(df_all)) || all(!is.finite(df_all$cfr))) {
    df_all$cfr <- 0
  }

  safe_inc(0.02, "Dropping rows missing x/y")
  df_all <- df_all[is.finite(df_all$x) & is.finite(df_all$y), , drop = FALSE]
  if (!nrow(df_all)) return(df_all)

  # ---- Pool DCR across all iterations per TID (optional) ----
  # Matches pyMINFLUX: ECO-weighted mean of dcr1, restricted to non-zero iterations.
  # Computed before keep_final_only so all iterations contribute.
  dcr_pool_tbl <- NULL
  if (isTRUE(pool_dcr) && "dcr1" %in% names(df_all) && "eco" %in% names(df_all)) {
    safe_inc(0.04, "Pooling DCR across iterations (ECO-weighted)")
    dcr_pool_tbl <- df_all |>
      dplyr::filter(is.finite(dcr1), dcr1 > 0, is.finite(eco), eco > 0) |>
      dplyr::group_by(tid) |>
      dplyr::summarise(dcr_pooled = sum(dcr1 * eco) / sum(eco), .groups = "drop")
  }

  # ---- Carry-forward per TID: CFR + ECC + EFC (single sorted pass) ----------
  safe_inc(0.08, "Computing per-TID carry-over (cfr/ecc/efc)")

  carry_cols <- switch(cfr_strategy,
    none   = character(0),
    itr6   = intersect(c("cfr", "ecc", "efc"), names(df_all)),
    intersect(c("cfr", "ecc", "efc"), names(df_all))  # last_nonzero
  )

  carry_tbls <- list()
  if (length(carry_cols)) {
    # Vectorised replacement for  arrange(tid, itr)  +  per-column
    # group_by(tid) |> slice_tail(n = 1): sort the rows ONCE by (tid, itr); for each
    # column the carried value is the last finite, positive entry per tid, i.e. the
    # last row of each tid-run in the sorted valid subset. Grouped slice_tail is one
    # of dplyr's slowest operations and dominated .npy import (~0.7 s on an 11k-loc
    # file); this base-R form gives the identical result an order of magnitude faster.
    base_idx <- if (cfr_strategy == "itr6") which(df_all$itr == 6L) else seq_len(nrow(df_all))
    if (length(base_idx)) {
      o     <- order(df_all$tid[base_idx], df_all$itr[base_idx])   # radix, stable
      idx_s <- base_idx[o]
      tid_s <- df_all$tid[idx_s]
      for (col in carry_cols) {
        v  <- df_all[[col]][idx_s]
        ok <- is.finite(v) & v > 0
        if (!any(ok)) next
        ts <- tid_s[ok]; vs <- v[ok]
        # keep the last row of each tid-run (value at the highest itr for that tid)
        keep <- c(ts[-length(ts)] != ts[-1], TRUE)
        tbl  <- data.frame(tid = ts[keep], vs[keep])
        names(tbl)[2] <- col
        carry_tbls[[col]] <- tbl
      }
    }
  }

  # ---- Keep only final-iteration rows ----------------------------------------
  use_fnl_now <- isTRUE(keep_final_only) && isTRUE(use_fnl) &&
                 "fnl" %in% names(df_all) && any(df_all$fnl %in% TRUE)
  if (isTRUE(keep_final_only) && !use_fnl_now &&
      (is.na(final_itr) || !is.finite(as.numeric(final_itr)))) {
    final_itr <- suppressWarnings(max(df_all$itr, na.rm = TRUE))
  }

  df <- if (!isTRUE(keep_final_only)) {
    df_all
  } else if (use_fnl_now) {
    safe_inc(0.06, "Filtering to final localizations (fnl flag)")
    df_all[df_all$fnl %in% TRUE, , drop = FALSE]
  } else {
    safe_inc(0.06, sprintf("Filtering to final iteration (itr=%s)", as.integer(final_itr)))
    df_all[!is.na(df_all$itr) & df_all$itr == as.integer(final_itr), , drop = FALSE]
  }
  rm(df_all)
  if (!nrow(df)) return(df)

  # ---- Join carried fields and pooled DCR ------------------------------------
  safe_inc(0.04, "Joining per-TID carry-over")
  for (col in names(carry_tbls)) {
    df <- df |> dplyr::select(-dplyr::all_of(col)) |>
      dplyr::left_join(carry_tbls[[col]], by = "tid")
  }
  if (!is.null(dcr_pool_tbl) && nrow(dcr_pool_tbl)) {
    df <- df |> dplyr::left_join(dcr_pool_tbl, by = "tid")
  }

  # ---- Dwell times -----------------------------------------------------------
  safe_inc(0.02, "Computing dwell times")
  df$dwell <- NA_real_
  ok <- is.finite(df$eco) & is.finite(df$efo) & df$efo > 0 &
        is.finite(dwell_time_ms) & dwell_time_ms > 0
  df$dwell[ok] <- round((df$eco[ok] / (df$efo[ok] / 1000)) / dwell_time_ms, 0)

  # ---- Units: meters -> nm; Z scaling ----------------------------------------
  safe_inc(0.04, "Converting units (m -> nm)")
  df$x <- df$x * 1e9
  df$y <- df$y * 1e9
  df$z <- df$z * 1e9
  if (is.finite(z_scale)) df$z <- df$z * z_scale

  safe_inc(0.04, "Done")
  df
}

# Canonical .npy reader for Abberior MINFLUX "tabular" exports.
# Mirrors read_minflux_json_as_df(): returns a data.frame with the same columns
# (x/y/z in nm, tid, itr, tim, efo, eco, cfr, dcr1, dcr, dwell, ...). The binary
# .npy is typically 5-10x smaller than the equivalent JSON and far faster to read.
read_minflux_npy_as_df <- function(path,
                                   keep_final_only = TRUE,
                                   keep_valid_only = TRUE,
                                   final_itr = NA_integer_,
                                   cfr_strategy = c("last_nonzero", "itr6", "none"),
                                   z_scale = 0.7,
                                   dwell_time_ms = 1.0,
                                   pool_dcr = FALSE) {

  cfr_strategy <- match.arg(cfr_strategy)
  safe_inc <- function(amount, detail = NULL) {
    try(shiny::incProgress(amount, detail = detail), silent = TRUE)
    invisible(NULL)
  }

  safe_inc(0.10, "Reading .npy file")
  arr <- read_npy(path)
  if (!is.list(arr))
    stop("Expected a structured MINFLUX .npy (named fields); got a plain array.")

  nm <- names(arr)
  n  <- attr(arr, "shape")[1]
  if (is.null(n) || !length(n)) n <- length(arr[[1]])
  if (!n) return(data.frame())

  # Coordinates: prefer 'loc', fall back to 'lnc' (same convention as JSON).
  coord_key <- if ("loc" %in% nm) "loc" else if ("lnc" %in% nm) "lnc" else NULL
  if (is.null(coord_key)) stop(".npy is missing the loc/lnc coordinate field.")
  cm <- arr[[coord_key]]
  if (is.null(dim(cm))) cm <- matrix(cm, ncol = 1)

  # Column getters tolerant of absent fields.
  gnum <- function(k) if (k %in% nm) {
    v <- arr[[k]]; if (!is.null(dim(v))) v[, 1] else as.numeric(v)
  } else rep(NA_real_, n)
  gint <- function(k) if (k %in% nm) as.integer(arr[[k]]) else rep(NA_integer_, n)
  glgl <- function(k) if (k %in% nm) as.logical(arr[[k]]) else rep(NA, n)

  safe_inc(0.10, "Building table")
  dcr_mat <- if ("dcr" %in% nm) arr[["dcr"]] else NULL
  dcr1 <- if (!is.null(dcr_mat)) {
    if (!is.null(dim(dcr_mat))) dcr_mat[, 1] else as.numeric(dcr_mat)
  } else rep(NA_real_, n)

  # itr may be absent in the tabular format; synthesize a per-TID order from tim
  # so carry-forward / iteration logic still works.
  tid_v <- if ("tid" %in% nm) as.integer(arr[["tid"]]) else rep(NA_integer_, n)
  tim_v <- gnum("tim")
  itr_v <- if ("itr" %in% nm) as.integer(arr[["itr"]]) else {
    ave_rank <- function(x) rank(x, ties.method = "first")
    as.integer(stats::ave(tim_v, tid_v, FUN = ave_rank)) - 1L
  }

  df_all <- data.frame(
    x    = cm[, 1],
    y    = if (ncol(cm) >= 2) cm[, 2] else rep(NA_real_, n),
    z    = if (ncol(cm) >= 3) cm[, 3] else rep(NA_real_, n),
    tid  = tid_v,
    itr  = itr_v,
    tim  = tim_v,
    efo  = gnum("efo"),
    eco  = gnum("eco"),
    cfr  = gnum("cfr"),
    ecc  = gnum("ecc"),
    efc  = gnum("efc"),
    fbg  = gnum("fbg"),
    sta  = gnum("sta"),
    bot  = glgl("bot"),
    eot  = glgl("eot"),
    fnl  = glgl("fnl"),
    dcr1 = dcr1,
    stringsAsFactors = FALSE
  )

  # vld filter (matches JSON path).
  if (isTRUE(keep_valid_only) && "vld" %in% nm) {
    vld <- as.logical(arr[["vld"]])
    df_all <- df_all[!is.na(vld) & vld, , drop = FALSE]
  }
  rm(arr)
  if (!nrow(df_all)) return(data.frame())

  # Final-iteration selection uses the fnl flag when present (Abberior tabular).
  .minflux_finalize_df(
    df_all,
    keep_final_only = keep_final_only,
    final_itr       = final_itr,
    cfr_strategy    = cfr_strategy,
    z_scale         = z_scale,
    dwell_time_ms   = dwell_time_ms,
    pool_dcr        = pool_dcr,
    use_fnl         = TRUE,
    safe_inc        = safe_inc
  )
}

# -- collapse_minflux_iterations() --------------------------------------------
# GOAL: reduce MINFLUX's multi-stage localization to ONE row per localization.
#
# MINFLUX does not measure a position in one shot. It ITERATES: each iteration
# (`itr`) tightens the doughnut pattern around the molecule, so a single
# localization event appears in the raw file as several rows - one per
# iteration, progressively more precise. Downstream we want the final answer,
# not the refinement history.
#
# The grouping key is (tid, gri): one trace, one localization event.
#
# Two different rules are needed, which is the subtlety here:
#   pick_itr()       take the value from the TARGET iteration (the last, most
#                    refined one). Used for position and most metrics - an
#                    earlier, coarser iteration would simply be wrong.
#   latest_finite()  take the value from the highest iteration that actually
#                    has one. Used for `cfr`, which is not reported at every
#                    iteration, so demanding the target iteration would give NA.
#   vld              TRUE if ANY iteration flagged the event valid.
#
# df       : raw long-format table with tid, gri, itr (+ the value columns)
# keep_itr : which iteration is the final one (default 9)
# returns  : one row per (tid, gri)
collapse_minflux_iterations <- function(df, keep_itr = 9L) {
  stopifnot(is.data.frame(df))
  if (!nrow(df)) return(df)

  if (!all(c("tid","gri","itr") %in% names(df))) {
    stop("JSON table missing required columns tid/gri/itr for iteration collapsing.")
  }

  latest_finite <- function(v, itr) {
    ok <- is.finite(v)
    if (!any(ok)) return(NA_real_)
    v[ which.max(itr * ok) ]
  }

  pick_itr <- function(v, itr, target) {
    idx <- which(itr == target)
    if (length(idx)) v[idx[1]] else NA_real_
  }

  df %>%
    dplyr::group_by(tid, gri) %>%
    dplyr::summarise(
      x   = pick_itr(x,   itr, keep_itr),
      y   = pick_itr(y,   itr, keep_itr),
      z   = pick_itr(z,   itr, keep_itr),
      cfr = latest_finite(cfr, itr),
      tim = pick_itr(tim, itr, keep_itr),
      efo = pick_itr(efo, itr, keep_itr),
      eco = pick_itr(eco, itr, keep_itr),
      dcr = pick_itr(dcr, itr, keep_itr),
      sta = pick_itr(sta, itr, keep_itr),
      vld = any(isTRUE(vld)),
      .groups = "drop"
    )
}

# -- canonicalize_minflux_df() ------------------------------------------------
# GOAL: accept the many column-naming conventions found in the wild and rename
# them to the ONE canonical schema documented in the file header, so nothing
# downstream ever has to guess.
#
# The same quantity is spelled differently by different exporters and software
# versions - "x", "X", "x_nm", "x [nm]", "pos_x", "xPosition" are all the same
# column. `map` lists the accepted aliases for each canonical name, in priority
# order, and `pick()` takes the first one actually present.
#
# It also repairs the things that are legitimately ABSENT rather than misnamed:
#   * missing `cfr` -> set to 0. Some acquisitions skip the centre-frequency
#     check to keep cycle times short. Setting 0 makes the CFR filter a harmless
#     no-op instead of leaving NA to propagate through every later filter.
#   * numeric coercion, because CSV readers sometimes infer character.
#
# df      : a data.frame from any reader / any MINFLUX export vintage
# returns : the same data with canonical column names and types
# =============================================================================
# detect_csv_z_calibration() - WHAT AXIAL FACTOR IS ALREADY BAKED INTO A CSV?
# =============================================================================
#
# GOAL
#   A CSV never contains raw instrument z. Something upstream already applied an
#   axial (refractive-index) correction, and applying ours on top would square
#   the factor. This identifies what was already applied so the caller can
#   divide it out and re-apply the app-wide value - making every source respond
#   identically to the axial calibration control.
#
# THE THREE CSV PROVENANCES
#   1. OUR export, current   - carries a `z_scale_applied` column. Exact: the
#                              factor is read from the file, never guessed.
#   2. OUR export, older     - has `dataset_tag` but no `z_scale_applied`. Every
#                              such export was written at 0.7.
#   3. pyMINFLUX export      - neither column; header is
#                              tid,tim,x,y,z,efo,cfr,eco,dcr,dwell,fluo.
#                              pyMINFLUX applies its 0.7 default before writing.
#
# THE ASSUMPTION, STATED PLAINLY
#   Cases 2 and 3 are inferred, not measured - nothing in those files records
#   the factor. Both are assumed to be 0.7 because every known producer of a
#   MINFLUX CSV ships that default. A CSV written with pyMINFLUX's z-scaling
#   turned off, or by some other tool, would be misread, and there is no way to
#   detect that from the data. Hence `certain`: the caller should say out loud
#   which branch it took.
#
# ARGUMENTS  df : the freshly read CSV (before canonicalization)
# RETURNS
#   list(factor, source, certain) - `factor` is what is already applied to z.
# =============================================================================
PYMINFLUX_Z_SCALE_DEFAULT <- 0.7

detect_csv_z_calibration <- function(df) {
  nms <- names(df)

  # 1. Explicit provenance - the only non-guessing branch.
  if ("z_scale_applied" %in% nms) {
    v <- suppressWarnings(as.numeric(df$z_scale_applied))
    v <- v[is.finite(v) & v > 0]
    if (length(v)) {
      return(list(factor = stats::median(v),
                  source = "this app (factor recorded in the file)",
                  certain = TRUE))
    }
  }

  # 2. Our own older export - tagged, but predates the provenance column.
  if (any(c("dataset_tag", ".dataset_tag") %in% nms)) {
    return(list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
                source = "this app (older export, no recorded factor)",
                certain = FALSE))
  }

  # 3. Anything else - pyMINFLUX or an unknown producer.
  py_cols <- c("tid", "tim", "x", "y", "z", "efo", "cfr", "eco", "dcr", "dwell", "fluo")
  looks_py <- all(py_cols %in% nms)
  list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
       source  = if (looks_py) "pyMINFLUX export" else "unrecognized CSV",
       certain = FALSE)
}

# =============================================================================
# has_axial_data() - DOES THIS DATASET CARRY REAL z AT ALL?
# =============================================================================
#
# GOAL
#   One predicate both tabs consult BEFORE doing, or saying, anything about the
#   axial calibration. A 2D acquisition should never see a message about z being
#   detected, removed, rescaled or adopted - those are answers to a question the
#   data does not pose.
#
# WHY A COLUMN CHECK IS NOT ENOUGH
#   A planar dataset usually still has a z column: MINFLUX 2D exports write
#   z = 0, and the loaders synthesise `z <- 0` when it is missing entirely. So
#   presence proves nothing - only VARIATION does. A constant z (whatever its
#   value) carries no axial structure and cannot be calibrated.
#
# ARGUMENTS
#   df     : data.frame that may have a z column
#   min_sd : SD in nm below which z counts as flat. The default is a float-noise
#            floor; the tracking tab passes 1 nm to match .has_real_z(), which
#            already governs its 2D/3D view toggle.
#   min_n  : minimum finite z values before an SD is trustworthy
#
# RETURNS  TRUE only when z varies enough to be worth calibrating.
# =============================================================================
has_axial_data <- function(df, min_sd = 1e-6, min_n = 1L) {
  if (is.null(df) || !is.data.frame(df) || !("z" %in% names(df))) return(FALSE)
  zv <- suppressWarnings(as.numeric(df$z))
  zv <- zv[is.finite(zv)]
  if (length(zv) < max(1L, min_n)) return(FALSE)
  s <- suppressWarnings(stats::sd(zv))
  is.finite(s) && s > min_sd
}

# =============================================================================
# strip_csv_z_calibration() - RETURN A CSV'S z TO RAW, UNCORRECTED nm
# =============================================================================
# GOAL
#   Undo whatever axial factor a CSV arrived with, so the caller can apply the
#   one factor the user controls. Every tab must normalise to raw before
#   applying its calibration, or .npy/.json and .csv drift apart and a CSV gets
#   the correction applied twice.
#
# ARGUMENTS  df : freshly read CSV (before canonicalization)
# RETURNS    df with raw z, the provenance column dropped, and the detection
#            result attached as attr(df, "z_cal") for the caller to report.
# =============================================================================
strip_csv_z_calibration <- function(df) {
  cal <- detect_csv_z_calibration(df)
  # Whether the caller should say anything about this at all. Computed before
  # the division so it describes the file as it arrived.
  cal$has_z <- has_axial_data(df)
  if ("z" %in% names(df) && is.finite(cal$factor) && cal$factor > 0 && cal$factor != 1) {
    df$z <- suppressWarnings(as.numeric(df$z)) / cal$factor
  }
  # Consumed here - a stale value must not survive into the analysis frame and
  # be re-read on a later round-trip.
  df$z_scale_applied <- NULL
  attr(df, "z_cal") <- cal
  df
}

canonicalize_minflux_df <- function(df) {
  stopifnot(is.data.frame(df))

  nms <- names(df)

  pick <- function(options) {
    hit <- options[options %in% nms]
    if (length(hit)) hit[1] else NA_character_
  }

  map <- list(
    x     = c("x", "X", "x_nm", "x [nm]", "pos_x", "xpos", "xPos", "x_position", "xPosition"),
    y     = c("y", "Y", "y_nm", "y [nm]", "pos_y", "ypos", "yPos", "y_position", "yPosition"),
    z     = c("z", "Z", "z_nm", "z [nm]", "pos_z", "zpos", "zPos", "z_position", "zPosition"),
    tid   = c("tid", "trace_id", "traceId", "traceID", "track_id", "trackId", "trackID", "id"),
    tim   = c("tim", "time", "t", "frame", "frame_idx", "frameIndex", "timestamp"),
    efo   = c("efo", "EFO", "photon_count", "photons", "nPhotons", "n_photons"),
    cfr   = c("cfr", "CFR"),
    dcr   = c("dcr", "DCR"),
    eco   = c("eco", "ECO"),
    dwell = c("dwell", "dwell_time", "dwellTime"),
    fluo  = c("fluo", "channel", "ch", "color_channel", "colorChannel")
  )

  for (target in names(map)) {
    src <- pick(map[[target]])
    if (!is.na(src) && src != target) {
      df[[target]] <- df[[src]]
    }
  }

  if (!("x" %in% names(df) && "y" %in% names(df))) {
    stop("Loaded data is missing x/y columns after canonicalization.")
  }

  as_num <- function(v) suppressWarnings(as.numeric(v))
  df$x <- as_num(df$x)
  df$y <- as_num(df$y)
  if ("z" %in% names(df)) df$z <- as_num(df$z)

  if ("tid" %in% names(df)) df$tid <- as.character(df$tid)

  for (nm in c("tim","efo","cfr","dcr","eco","dwell","fluo")) {
    if (nm %in% names(df)) df[[nm]] <- as_num(df[[nm]])
  }

  # CFR is optional: some acquisitions skip the CFR check (it lengthens cycle
  # times), so a CSV may have no cfr column. Set it to 0 throughout when absent or
  # entirely non-finite, so the column always exists and its filter is a no-op.
  if (!("cfr" %in% names(df)) || all(!is.finite(df$cfr))) {
    df$cfr <- 0
  }

  med_abs <- function(v) stats::median(abs(v[is.finite(v)]), na.rm = TRUE)
  mx <- med_abs(df$x); my <- med_abs(df$y)
  if (is.finite(mx) && is.finite(my) && mx < 1e-3 && my < 1e-3) {
    df$x <- df$x * 1e9
    df$y <- df$y * 1e9
    if ("z" %in% names(df)) df$z <- df$z * 1e9
  }

  df
}

# Build scale-bar shapes + annotation for a 2-D Plotly layout.
# Uses paper coordinates (0-1 fractions of the plot area) so the bar always
# sits in the bottom-right corner regardless of zoom or pan.
# thickness_nm controls the visual bar height (scaled to paper units).
# Returns list(shapes = list(...), annotations = list(...)) or NULL.
make_scale_bar <- function(length_nm, x_range_nm, thickness_nm = 10, dark = FALSE, scale_font_size = 18) {
  if (!is.finite(length_nm) || length_nm <= 0) return(NULL)
  if (!is.finite(x_range_nm) || x_range_nm <= 0) return(NULL)
  if (!is.finite(thickness_nm) || thickness_nm <= 0) thickness_nm <- 10
  bar_col <- if (dark) "white" else "black"
  # bar height in paper units: ~0.01 for thickness 10, scales linearly, capped
  bar_h   <- max(0.005, min(0.06, thickness_nm / 1000))
  # bar width in paper units proportional to the actual axis range
  bar_w   <- min(0.85, length_nm / x_range_nm)
  x1 <- 0.97; x0 <- x1 - bar_w
  y0 <- 0.02; y1 <- y0 + bar_h
  # auto-detect µm
  label <- if (length_nm >= 1000) {
    val <- length_nm / 1000
    paste0(if (val == round(val)) as.integer(val) else val, "\u00a0\u00b5m")
  } else {
    paste0(length_nm, "\u00a0nm")
  }
  list(
    shapes = list(list(
      type = "rect",
      xref = "paper", yref = "paper",
      x0 = x0, x1 = x1, y0 = y0, y1 = y1,
      fillcolor = bar_col,
      line = list(color = bar_col, width = 0)
    )),
    annotations = list(list(
      x = (x0 + x1) / 2, y = y1 + 0.015,
      xref = "paper", yref = "paper",
      text = label,
      showarrow = FALSE,
      font = list(color = bar_col, size = scale_font_size),
      xanchor = "center", yanchor = "bottom"
    ))
  )
}

# Single-colour LUT palette for background overlays. Values are the full-intensity
# colour each LUT ramps up to (black -> colour). Names are the UI/setting keys.
lut_overlay_colors <- c(
  green   = "#00FF00",
  magenta = "#FF00FF",
  yellow  = "#FFFF00",
  cyan    = "#00FFFF",
  red     = "#FF0000",
  blue    = "#0000FF",
  orange  = "#FF7F00"
)

# Recolour + adjust a background overlay image through a single-colour LUT.
#
# MINFLUX confocal exports often arrive in a Fire (multi-colour) LUT. To recolour
# one to match the fluorophore, we first collapse it to intensity (grayscale
# luminance) and then linearly map black -> chosen colour with a Multiply
# composite, so intensity v maps to v * colour. This reproduces ImageJ-style
# single-colour LUTs (Green, Magenta, Yellow, ...). A "custom" LUT tints toward
# any hex colour.
#
# Brightness (gain) and contrast (symmetric linear stretch) are applied to the
# intensity before tinting, and also apply in "original" mode (adjusting the Fire
# image in place), so the two sliders behave the same whatever the LUT.
#
# NB: grayscale luminance is only an approximate inverse of the Fire LUT - it
# recovers relative intensity for *display*, not the exact raw photon counts. For
# an exact recolour the user would need to re-export the raw grayscale image.
#
# img          : magick image object.
# lut          : names(lut_overlay_colors), "original" (keep colours), "gray", or
#                "custom" (tint toward custom_color).
# custom_color : hex "#RRGGBB" used when lut == "custom".
# brightness   : gain in % (100 = unchanged).
# contrast     : symmetric linear stretch, 0 (none) .. ~90 (strong).
# Returns a magick image (unchanged when no LUT and both sliders neutral).
apply_background_lut <- function(img, lut = "original", custom_color = NULL,
                                 brightness = 100, contrast = 0) {
  lut <- tolower(as.character(lut %||% "original"))
  if (length(lut) != 1L) lut <- "original"
  brightness <- suppressWarnings(as.numeric(brightness %||% 100))
  if (!is.finite(brightness)) brightness <- 100
  contrast <- suppressWarnings(as.numeric(contrast %||% 0))
  if (!is.finite(contrast)) contrast <- 0

  # Resolve the recolour target: NULL = keep original colours; is_gray = desaturate;
  # otherwise a hex colour to tint toward (preset or custom).
  is_gray <- lut %in% c("gray", "grays", "grey", "greyscale", "grayscale")
  col <- if (lut == "custom") {
           cc <- as.character(custom_color %||% "")[1L]
           if (grepl("^#[0-9A-Fa-f]{6}$", cc)) cc else NULL
         } else if (lut %in% names(lut_overlay_colors)) {
           lut_overlay_colors[[lut]]
         } else NULL
  recolor <- is_gray || !is.null(col)

  neutral_bc <- abs(brightness - 100) <= 1e-6 && contrast <= 0
  if (!recolor && neutral_bc) return(img)            # original, no adjustment -> untouched

  # Collapse to intensity when recolouring (strips the Fire LUT); keep RGB otherwise.
  if (recolor) img <- magick::image_convert(img, colorspace = "Gray")

  if (abs(brightness - 100) > 1e-6)
    img <- magick::image_modulate(img, brightness = brightness)
  if (contrast > 0) {
    c2  <- min(contrast, 95) / 2
    img <- magick::image_level(img, black_point = c2, white_point = 100 - c2)
  }

  if (!recolor) return(img)                                             # original + B/C
  if (is_gray || is.null(col))
    return(magick::image_convert(img, colorspace = "sRGB"))            # grayscale + B/C

  gray_rgb <- magick::image_convert(img, colorspace = "sRGB")
  info <- magick::image_info(gray_rgb)
  tint <- magick::image_blank(info$width[1L], info$height[1L], color = col)
  magick::image_composite(gray_rgb, tint, operator = "Multiply")
}

# Auto-align a confocal overlay image to MINFLUX/SMLM localizations.
#
# img      : magick image object (already rotated/flipped as displayed)
# df_xy    : data.frame with finite x, y columns (nm)
# px_size  : nm per pixel
# mode     : "foci"     - multi-peak voting: top N bright spots vs top N dense
#                         MINFLUX regions; offset with most consensus wins
#            "boundary" - centroid of bright image region to centroid of data
# invert_y : TRUE (default) = image origin top-left, row 0 → high data y
# n_peaks  : number of candidate peaks to extract for foci mode (default 4)
#
# Returns list(x0, y0) for the image origin input, or NULL on failure.
auto_align_overlay <- function(img, df_xy, px_size,
                                mode        = c("foci", "boundary"),
                                invert_y    = TRUE,
                                n_peaks        = 6L,
                                bright_pct     = 0.95,  # top-X% pixels = "bright"
                                border_pct     = 0.10,  # mask this fraction of each edge
                                aspect_max     = 4.0,   # blobs more elongated than this are rejected
                                eps_pct        = 0.03,  # dbscan radius as fraction of image size
                                density_thr    = 0.20,  # MINFLUX density threshold (fraction of peak)
                                density_cell   = 100,   # density grid cell size (nm); 0 = auto
                                density_smooth = 4) {   # Gaussian smoothing sigma (nm; 0 = no smooth)
  mode <- match.arg(mode)

  info <- magick::image_info(img)
  W    <- info$width
  H    <- info$height

  # Convert a gray magick image to a numeric matrix[row, col], values 0-255.
  # magick::image_data returns dim [channels, width, height] filled col-major,
  # so channel=1, col varies fastest → matrix [W, H] then transpose to [H, W].
  img_to_mat <- function(g) {
    gi  <- magick::image_info(g)
    W_  <- as.integer(gi$width[1L])
    H_  <- as.integer(gi$height[1L])
    tmp <- tempfile(fileext = ".gray")
    on.exit(unlink(tmp), add = TRUE)
    magick::image_write(g[1L], path = tmp, format = "gray", depth = 8L)
    fsz <- as.integer(file.size(tmp))
    n_  <- W_ * H_
    if (fsz == n_) {
      v <- as.integer(readBin(tmp, raw(), n = n_))
    } else if (fsz == 2L * n_) {
      r <- readBin(tmp, raw(), n = fsz)
      v16 <- as.integer(r[seq(1L,fsz,2L)]) * 256L + as.integer(r[seq(2L,fsz,2L)])
      mx  <- max(v16); if (mx == 0L) mx <- 1L
      v   <- as.integer(round(v16 / mx * 255L))
    } else {
      v <- as.integer(magick::image_data(g, channels = "gray"))[seq_len(n_)]
    }
    matrix(v, nrow = H_, ncol = W_, byrow = TRUE)
  }

  # Downsample to ≤ 256 px on longest side for fast processing
  max_dim <- 256L
  scale   <- min(1.0, max_dim / max(W, H))
  W2      <- max(1L, round(W * scale))
  H2      <- max(1L, round(H * scale))
  gray    <- magick::image_convert(img, colorspace = "Gray", depth = 8L)
  gray2   <- if (scale < 1) magick::image_resize(gray, paste0(W2, "x", H2, "!")) else gray

  x <- df_xy$x[is.finite(df_xy$x)]
  y <- df_xy$y[is.finite(df_xy$y)]
  if (!length(x)) return(NULL)

  # Helper: pixel centre (0-indexed full-res col/row) → data-space nm coords
  # assuming image origin x0=0, y0=0.
  px_to_nm <- function(col_f, row_f) {
    xnm <- (col_f + 0.5) * px_size
    ynm <- if (invert_y) (H - row_f - 0.5) * px_size else (row_f + 0.5) * px_size
    c(x = xnm, y = ynm)
  }

  img_pks_nm   <- NULL   # filled in foci mode only
  data_pks_nm  <- NULL

  if (mode == "foci") {
    # ---- Foci alignment: find bright blobs in image, dense clusters in MINFLUX,
    # then match using a bounds constraint (MINFLUX extent must land inside the
    # physical image after alignment) with rank preference.

    # -- Robust 1-D Gaussian smooth (kernel clamped to vector length) ----------
    gauss1d <- function(v, sigma) {
      n <- length(v)
      if (n < 3L || !is.finite(sigma) || sigma <= 0) return(v)
      w  <- min(max(1L, as.integer(ceiling(sigma * 3))),
                max(1L, as.integer(floor((n - 1L) / 2L))))
      kv <- exp(-((-w):w)^2 / (2.0 * sigma^2)); kv <- kv / sum(kv)
      fv <- stats::filter(v, kv, sides = 2L)
      ifelse(is.na(fv), v, as.double(fv))
    }
    gauss2d <- function(m, sigma) {
      nr <- nrow(m); nc <- ncol(m)
      m  <- matrix(as.double(m), nr, nc)
      for (i in seq_len(nr)) m[i, ] <- gauss1d(m[i, ], sigma)
      for (j in seq_len(nc)) m[, j] <- gauss1d(m[, j], sigma)
      m
    }

    # -- Peak finding via which.max (linear index → row/col): crash-proof ------
    # Avoids the which(m==max, arr.ind=TRUE)[1,] pattern that fails on edge cases.
    find_peaks <- function(m, n_max, min_sep) {
      nr   <- nrow(m); nc <- ncol(m)
      work <- matrix(as.double(m), nr, nc)
      pks  <- vector("list", n_max)
      nfound <- 0L
      for (k in seq_len(n_max)) {
        idx <- which.max(work)       # single integer: linear index into matrix
        if (!length(idx)) break
        v <- work[[idx]]
        if (!is.finite(v) || v <= 0) break
        ri <- ((idx - 1L) %% nr) + 1L   # row (1-indexed)
        ci <- ((idx - 1L) %/% nr) + 1L  # col (1-indexed)
        nfound <- nfound + 1L
        pks[[nfound]] <- c(row = ri, col = ci)
        r1 <- max(1L, ri - min_sep); r2 <- min(nr, ri + min_sep)
        c1 <- max(1L, ci - min_sep); c2 <- min(nc, ci + min_sep)
        work[r1:r2, c1:c2] <- -Inf
      }
      pks[seq_len(nfound)]
    }

    # -- 1. Confocal: gray-convert → resize → write raw → read back → matrix ----
    # We bypass image_data() entirely and write the grayscale image to a temp
    # "GRAY" raw file (ImageMagick raw format: exactly W×H bytes, row-major,
    # no header).  This is unambiguous regardless of bit depth or colorspace quirks.
    max_dim_c <- 512L
    W_orig    <- as.integer(W[1L]); H_orig <- as.integer(H[1L])  # scalar originals
    sc_c_t    <- min(1.0, max_dim_c / max(W_orig, H_orig))
    gray_c    <- magick::image_convert(img[1L], colorspace = "Gray")
    img_c     <- if (sc_c_t < 1) {
      Wc_t <- max(1L, round(W_orig * sc_c_t))
      Hc_t <- max(1L, round(H_orig * sc_c_t))
      magick::image_resize(gray_c, paste0(Wc_t, "x", Hc_t, "!"))
    } else gray_c
    if (length(img_c) > 1L) img_c <- img_c[1L]

    # Use image_info to get actual post-resize dimensions
    info_c <- magick::image_info(img_c)
    Wc     <- as.integer(info_c$width[1L])
    Hc     <- as.integer(info_c$height[1L])
    sc_c   <- Wc / W_orig

    # Write as raw GRAY (8-bit): exactly Wc*Hc bytes, row-major (no header).
    # Fall back to 16-bit if the file is 2× the expected size.
    tmp_gray <- tempfile(fileext = ".gray")
    on.exit(unlink(tmp_gray), add = TRUE)
    magick::image_write(img_c, path = tmp_gray, format = "gray", depth = 8L)
    n_px    <- Wc * Hc
    fsz     <- as.integer(file.size(tmp_gray))
    raw_len <- fsz

    if (fsz == n_px) {
      bit_mode <- "8-bit"
      data_v   <- as.integer(readBin(tmp_gray, raw(), n = n_px))
    } else if (fsz == 2L * n_px) {
      bit_mode <- "16-bit"
      raw16    <- readBin(tmp_gray, raw(), n = fsz)
      hi       <- as.integer(raw16[seq(1L, fsz, 2L)])
      lo       <- as.integer(raw16[seq(2L, fsz, 2L)])
      v16      <- hi * 256L + lo
      mx16     <- max(v16); if (mx16 == 0L) mx16 <- 1L
      data_v   <- as.integer(round(v16 / mx16 * 255L))
    } else {
      # Unexpected size - try image_data as last resort
      bit_mode <- paste0("fallback(", fsz, "!=", n_px, ")")
      raw_arr  <- magick::image_data(img_c, channels = "gray")
      data_v   <- as.integer(raw_arr)[seq_len(n_px)]
    }

    # Normalize to full 0-255 range so detection is scale-independent
    mn_v <- min(data_v); mx_v <- max(data_v)
    if (mx_v > mn_v)
      data_v <- as.integer(round((data_v - mn_v) / (mx_v - mn_v) * 255L))
    mat_c <- matrix(data_v, nrow = Hc, ncol = Wc, byrow = TRUE)

    # Border mask - fraction set by caller
    brd <- max(1L, round(min(Hc, Wc) * border_pct))
    mat_c[seq_len(brd), ]            <- 0L
    mat_c[seq(Hc - brd + 1L, Hc), ] <- 0L
    mat_c[, seq_len(brd)]            <- 0L
    mat_c[, seq(Wc - brd + 1L, Wc)] <- 0L

    # Small Gaussian smooth: sigma ~0.5 % of image (≈ 2-3 px at 512 px) is enough
    # to merge pixels within the same focus without destroying the signal.
    # A larger sigma averages too much dark background into each bright focus,
    # reducing the peak to near-zero and making detection unreliable.
    sig_c  <- max(2.0, min(Wc, Hc) * 0.005)
    mat_cs <- gauss2d(mat_c, sig_c)

    # eps_pct controls the MINIMUM SEPARATION between two detected peaks
    # (i.e. how close two green circles can be), not the smoothing radius.
    sep_c       <- max(3L, round(min(Wc, Hc) * eps_pct))
    img_pks_all <- find_peaks(mat_cs, n_peaks * 2L, sep_c)
    if (!length(img_pks_all))
      stop("foci: no bright regions found in confocal image")

    # Sensitivity post-filter: drop peaks below rel_frac × brightest smoothed peak.
    smooth_vals <- sapply(img_pks_all, function(pk) mat_cs[pk[["row"]], pk[["col"]]])
    rel_frac    <- ((bright_pct - 0.70) / 0.29) ^ 0.5
    keep        <- smooth_vals >= rel_frac * max(smooth_vals)
    img_pks     <- img_pks_all[keep][seq_len(min(sum(keep), n_peaks))]
    if (!length(img_pks))
      stop("foci: no peaks above sensitivity threshold; lower the sensitivity slider")

    # Record raw grayscale values at surviving peaks for the debug notification
    peak_vals     <- sapply(img_pks, function(pk) mat_c[pk[["row"]], pk[["col"]]])
    img_mat_range <- range(mat_c)

    # Peak (row, col) in downsampled image → 0-indexed full-res pixel → nm
    img_pks_nm <- lapply(img_pks, function(pk) {
      col_f <- (pk[["col"]] - 0.5) / sc_c - 0.5
      row_f <- (pk[["row"]] - 0.5) / sc_c - 0.5
      nm    <- px_to_nm(col_f, row_f)
      list(x = nm[["x"]], y = nm[["y"]])
    })

    # -- 2. MINFLUX: TID-density histogram → smooth → top-N clusters ----------
    fin   <- is.finite(df_xy$x) & is.finite(df_xy$y)
    x_f   <- df_xy$x[fin]; y_f <- df_xy$y[fin]
    if (!length(x_f)) stop("foci: no finite MINFLUX coordinates")
    # Use cluster column when available (groups multiple blinks from one molecule);
    # fall back to tid, then to sequential index if neither column is present.
    id_v <- if (!is.null(df_xy$cluster) &&
                  any(df_xy$cluster[fin] > 0L, na.rm = TRUE))
              df_xy$cluster[fin]
            else if (!is.null(df_xy$tid))
              df_xy$tid[fin]
            else
              seq_along(x_f)

    extent_nm <- max(diff(range(x_f)), diff(range(y_f)))
    # Cell size: user value (nm) or auto (aim for ≤80 bins on longer axis, min 50 nm)
    cell_nm   <- if (density_cell > 0) max(10.0, as.double(density_cell))
                 else max(50.0, extent_nm / 80.0)
    xr   <- range(x_f); yr <- range(y_f)
    nx_d <- max(5L, as.integer(ceiling((xr[2] - xr[1]) / cell_nm)))
    ny_d <- max(5L, as.integer(ceiling((yr[2] - yr[1]) / cell_nm)))
    bx_d <- seq(xr[1], xr[2], length.out = nx_d + 1L)
    by_d <- seq(yr[1], yr[2], length.out = ny_d + 1L)
    ix_d <- pmax(1L, pmin(nx_d, findInterval(x_f, bx_d)))
    iy_d <- pmax(1L, pmin(ny_d, findInterval(y_f, by_d)))

    # Unique molecule count per cell (prevents dwell-heavy tracks from dominating)
    dedup <- !duplicated(cbind(ix_d, iy_d, as.integer(factor(id_v))))
    h_raw <- matrix(
      tabulate(ix_d[dedup] + (iy_d[dedup] - 1L) * nx_d, nbins = nx_d * ny_d),
      nrow = nx_d, ncol = ny_d)

    # Smoothing sigma in bins (density_smooth is in nm; 0 = skip smoothing)
    sig_d <- if (density_smooth > 0) max(0.5, density_smooth / cell_nm) else 0
    h_s   <- if (sig_d >= 0.5) gauss2d(h_raw, sig_d) else matrix(as.double(h_raw), nx_d, ny_d)

    thr_d  <- density_thr * max(h_s, na.rm = TRUE)
    h_st   <- h_s; h_st[h_st < thr_d] <- 0

    sep_d    <- max(5L, round(min(nx_d, ny_d) * 0.15))
    data_pks <- find_peaks(h_st, n_peaks, sep_d)
    if (!length(data_pks))
      stop("foci: no density clusters found in MINFLUX data")

    # h_raw row=ix (x-bin), col=iy (y-bin) → bin centre nm
    data_pks_nm <- lapply(data_pks, function(pk) {
      ix <- pk[["row"]]; iy <- pk[["col"]]
      list(x = (bx_d[ix] + bx_d[ix + 1L]) / 2.0,
           y = (by_d[iy] + by_d[iy + 1L]) / 2.0)
    })

    # -- 3. Geometric permutation matching -------------------------------------
    # Generate every ordered k-length selection of data peaks assigned to the
    # top-k image peaks (image peak rank 1..k always stays in position).
    # For each assignment compute the mean translation (optimal rigid offset for
    # that pairing) and score it by sum-of-squared residuals - how consistently
    # all k pairs agree on the same offset.  A perfect match scores 0; a poor
    # one scores large.  Among valid candidates (MINFLUX extent inside the image
    # with 10 % margin) we keep the lowest-score assignment.
    # Falls back to rank-weighted Hough voting if no candidate passes bounds.
    img_w_nm  <- W * px_size
    img_h_nm  <- H * px_size
    margin_nm <- 0.10 * min(img_w_nm, img_h_nm)

    # All ordered k-selections from 1:n (partial permutations), n! / (n-k)! rows.
    .ordered_sel <- function(n, k) {
      if (k == 0L || n == 0L) return(matrix(integer(0), 1L, 0L))
      if (n < k)              return(matrix(integer(0), 0L, k))
      if (k == 1L)            return(matrix(seq_len(n), ncol = 1L))
      result <- list()
      for (first in seq_len(n)) {
        sub <- .ordered_sel(n - 1L, k - 1L)
        if (!nrow(sub)) next
        adj <- sub; adj[adj >= first] <- adj[adj >= first] + 1L
        result[[length(result) + 1L]] <- cbind(first, adj)
      }
      do.call(rbind, result)
    }

    k_i    <- length(img_pks_nm)
    k_d    <- length(data_pks_nm)
    k_pair <- min(k_i, k_d)   # number of (image peak, data peak) pairs to form

    Px <- sapply(img_pks_nm[seq_len(k_i)],  function(p) p$x)
    Py <- sapply(img_pks_nm[seq_len(k_i)],  function(p) p$y)
    Qx <- sapply(data_pks_nm[seq_len(k_d)], function(p) p$x)
    Qy <- sapply(data_pks_nm[seq_len(k_d)], function(p) p$y)

    # sel_mat: each row is one ordered selection of k_pair data-peak indices
    sel_mat    <- .ordered_sel(k_d, k_pair)
    best_x0    <- NULL; best_y0 <- NULL; best_score <- Inf

    for (ri in seq_len(nrow(sel_mat))) {
      s     <- sel_mat[ri, ]                        # data indices matched to image ranks 1..k_pair
      dx    <- Qx[s] - Px[seq_len(k_pair)]
      dy    <- Qy[s] - Py[seq_len(k_pair)]
      ox    <- mean(dx);  oy <- mean(dy)
      # Score = geometric consistency: how tightly all k_pair pairs agree on (ox, oy)
      score <- sum((dx - ox)^2 + (dy - oy)^2)
      x_ok  <- (xr[1] >= ox - margin_nm) && (xr[2] <= ox + img_w_nm + margin_nm)
      y_ok  <- (yr[1] >= oy - margin_nm) && (yr[2] <= oy + img_h_nm + margin_nm)
      if (x_ok && y_ok && score < best_score) {
        best_score <- score; best_x0 <- ox; best_y0 <- oy
      }
    }

    if (is.null(best_x0)) {
      # No selection satisfied bounds - fall back to rank-weighted Hough voting
      all_ox <- all_oy <- all_wt <- numeric(k_i * k_d)
      idx <- 0L
      for (i in seq_len(k_i)) {
        for (j in seq_len(k_d)) {
          idx <- idx + 1L
          all_ox[idx] <- Qx[j] - Px[i]
          all_oy[idx] <- Qy[j] - Py[i]
          all_wt[idx] <- 1.0 / (i * j)
        }
      }
      bin      <- max(30.0, 0.06 * extent_nm)
      key      <- paste(round(all_ox / bin), round(all_oy / bin), sep = "_")
      wtab     <- tapply(all_wt, key, sum)
      best_key <- names(wtab)[which.max(wtab)]
      sel_h    <- key == best_key
      best_x0  <- stats::weighted.mean(all_ox[sel_h], all_wt[sel_h])
      best_y0  <- stats::weighted.mean(all_oy[sel_h], all_wt[sel_h])
    }

    x0 <- best_x0
    y0 <- best_y0

  } else {
    # ---- Boundary: bright-region centroid → data centroid + polygon extraction
    mat    <- img_to_mat(gray2)
    thresh <- as.integer(stats::quantile(as.integer(mat), 0.50))
    mask   <- mat >= thresh
    if (!any(mask)) return(NULL)

    row_idx <- matrix(rep(seq_len(H2) - 1L,        W2), nrow = H2)
    col_idx <- matrix(rep(seq_len(W2) - 1L, each = H2), nrow = H2)
    wt      <- as.double(mat) * mask
    total   <- sum(wt)
    if (!is.finite(total) || total == 0) return(NULL)

    # Scale centroid (0-indexed downsampled) back to full-res
    img_row <- (sum(row_idx * wt) / total + 0.5) / scale - 0.5
    img_col <- (sum(col_idx * wt) / total + 0.5) / scale - 0.5
    nm      <- px_to_nm(img_col, img_row)

    x0 <- mean(x) - nm[["x"]]
    y0 <- mean(y) - nm[["y"]]

    # -- Extract cell boundary polygon via convex hull of bright (membrane) pixels
    # The hull of membrane ring pixels approximates the outer cell boundary.
    # Vertices are converted to data-nm coordinates (using the alignment offset
    # computed above) so the polygon can be used directly for trace filtering.
    bright_rc        <- which(mask, arr.ind = TRUE)  # [N x 2]: col 1 = row, col 2 = col (1-indexed)
    boundary_polygon <- NULL
    if (nrow(bright_rc) >= 3L) {
      hull_idx <- grDevices::chull(bright_rc[, 2L], bright_rc[, 1L])  # chull(x=col, y=row)
      h_col_0  <- bright_rc[hull_idx, 2L] - 1L   # 0-indexed downsampled col
      h_row_0  <- bright_rc[hull_idx, 1L] - 1L   # 0-indexed downsampled row
      # Convert downsampled pixel centres to full-res 0-indexed centres, then to nm
      h_col_f  <- (h_col_0 + 0.5) / scale - 0.5
      h_row_f  <- (h_row_0 + 0.5) / scale - 0.5
      hull_x   <- sapply(seq_along(h_col_f), function(k) px_to_nm(h_col_f[k], h_row_f[k])[["x"]])
      hull_y   <- sapply(seq_along(h_col_f), function(k) px_to_nm(h_col_f[k], h_row_f[k])[["y"]])
      # Apply alignment offset → data nm coordinates (same space as MINFLUX data)
      boundary_polygon <- data.frame(x = hull_x + x0, y = hull_y + y0)
    }
  }

  pk_to_df <- function(pks) {
    if (!length(pks)) return(NULL)
    data.frame(x = sapply(pks, function(p) p[["x"]]),
               y = sapply(pks, function(p) p[["y"]]))
  }
  list(
    x0               = round(x0, 3),
    y0               = round(y0, 3),
    data_peaks       = pk_to_df(data_pks_nm),
    img_peaks        = pk_to_df(img_pks_nm),
    # boundary_polygon: data.frame(x, y) in data nm coords; non-NULL only in
    # boundary mode.  Pass to filter_by_cell_boundary() to keep inside traces.
    boundary_polygon = if (exists("boundary_polygon")) boundary_polygon else NULL,
    # alignment_score: sum of squared residuals from geometric permutation match
    # (foci mode only; lower = all peaks agree more tightly on the same offset)
    alignment_score  = if (exists("best_score") && is.finite(best_score)) round(best_score, 1) else NULL,
    img_peaks_debug  = if (exists("peak_vals"))     peak_vals     else NULL,
    mat_range_debug  = if (exists("img_mat_range")) img_mat_range else NULL,
    img_mat_debug    = if (exists("mat_c"))         mat_c         else NULL,
    img_mat_smooth   = if (exists("mat_cs"))        mat_cs        else NULL,
    img_dim_debug    = if (exists("Wc"))            c(W_orig, H_orig, Wc, Hc, raw_len, bit_mode) else NULL
  )
}

# -- Vectorised point-in-polygon via ray casting ------------------------------
# Returns logical vector; TRUE where (test_x[i], test_y[i]) is inside the
# polygon defined by vertices (poly_x, poly_y).  Polygon is auto-closed.
.pnpoly_vec <- function(test_x, test_y, poly_x, poly_y) {
  n      <- length(poly_x)
  inside <- rep(FALSE, length(test_x))
  j      <- n
  for (i in seq_len(n)) {
    cross  <- ((poly_y[i] > test_y) != (poly_y[j] > test_y)) &
              (test_x < (poly_x[j] - poly_x[i]) * (test_y - poly_y[i]) /
                        (poly_y[j] - poly_y[i]) + poly_x[i])
    inside <- xor(inside, cross)
    j <- i
  }
  inside
}

# Filter a localisation data frame to keep only traces whose median xy position
# falls inside the cell boundary polygon returned by
#   auto_align_overlay(..., mode = "boundary")$boundary_polygon
# The polygon must already be in data nm coordinates.
# Returns df unchanged if boundary_polygon is NULL or has fewer than 3 vertices.
filter_by_cell_boundary <- function(df, boundary_polygon) {
  if (is.null(boundary_polygon) || !is.data.frame(boundary_polygon) ||
      nrow(boundary_polygon) < 3L)                    return(df)
  if (!all(c("x", "y", "tid") %in% names(df)))       return(df)

  tid_med <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      mx = stats::median(x, na.rm = TRUE),
      my = stats::median(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(mx), is.finite(my))

  if (!nrow(tid_med)) return(df)

  inside    <- .pnpoly_vec(tid_med$mx, tid_med$my,
                           boundary_polygon$x, boundary_polygon$y)
  keep_tids <- tid_med$tid[inside]
  dplyr::filter(df, tid %in% keep_tids)
}
