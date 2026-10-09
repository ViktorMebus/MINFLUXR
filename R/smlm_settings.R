## ─────────────────────────────────────────────────────────────────────────────
## R/smlm_settings.R  -  the settings that decide what the SMLM tab computes
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   One definition of every setting that changes SMLM results: the filters, the
##   channel split, de-blinking and median clustering, and the two opt-in gates
##   that rewrite the data before clustering (fiducial exclusion, density
##   thinning). For each one it records the default, the type and the allowed
##   values. The UI widgets, the helper functions' fallbacks and the reset buttons
##   all read their defaults from here, so the three can no longer disagree.
##
## THE SETTINGS OBJECT
##   A named list whose names are the Shiny input IDs, so it can be passed to any
##   helper that takes `input` (apply_meta_filters_batch, apply_clustering_pipeline,
##   resolve_channels, ...) without changing that helper:
##       s  <- smlm_settings(enable_clustering = TRUE, eps2 = 150)
##       df <- apply_clustering_pipeline(apply_meta_filters_batch(raw, s), s)
##   `$` on it matches names exactly, so a misspelt field returns NULL instead of
##   silently matching a longer name.
##
## FUNCTIONS
##   smlm_settings(...) ............ defaults, overridden by `...`, validated
##   smlm_settings_from_input(input) the same, read from a Shiny `input`
##   as_smlm_settings(x) ........... what a helper calls on its `settings` argument:
##                                   passes settings through, wraps `input` in a
##                                   lazy view (see the function for why lazy)
##   smlm_set(s, ...) .............. a copy with some settings replaced
##   smlm_settings_key(s) .......... hash of the settings, for cache keys
##   smlm_settings_from_file(raw) .. settings from a saved settings file (either layout)
##   smlm_settings_to_legacy() ..... write settings into the older file layout
##   smlm_settings_restore_plan() .. which widget update restores each setting, in order
##   SMLM_DISPLAY_DEFAULTS ......... display options (colour, palette, medians); not settings
##   smlm_display_from_input() ..... display options from the app
##   smlm_channel_colors() ......... channel A/B colours for attach_channel_colors()
##   SMLM_SETTINGS_SPEC ............ default / type / limits / stage per setting
##   SMLM_DEFAULTS ................. the defaults alone, for UI widgets and fallbacks
##
## Display-only inputs (colours, palettes, point sizes) are not here: they change
## the picture, never a number.
## ─────────────────────────────────────────────────────────────────────────────

# GOAL: one row per setting. `type` decides how a value is coerced:
#   flag   TRUE/FALSE                 num    one number      int   one whole number
#   range  NULL (no filter) or c(lo, hi)                     choice one of `choices`
# `min`/`max` bound a number; `na_ok` lets NA through (an unset spread target).
# `stage` is documentation: where in the pipeline the setting acts.
.smlm_min_locs <- if (exists("MIN_LOCS_FOR_SPREAD_DEFAULT")) MIN_LOCS_FOR_SPREAD_DEFAULT else 5L
.ss <- function(default, type, stage, min = -Inf, max = Inf, choices = NULL, na_ok = FALSE)
  list(default = default, type = type, stage = stage, min = min, max = max,
       choices = choices, na_ok = na_ok)

SMLM_SETTINGS_SPEC <- list(
  # ── Filters. Ranges default to NULL = no filter; the app sets them from each
  #    file's data range on upload.
  efo_range        = .ss(NULL, "range", "filter"),
  cfr_range        = .ss(NULL, "range", "filter"),
  dcr_range        = .ss(NULL, "range", "filter"),
  eco_range        = .ss(NULL, "range", "filter"),
  tim_range        = .ss(NULL, "range", "filter"),
  dwell_range      = .ss(NULL, "range", "filter"),
  fluo_range       = .ss(NULL, "range", "filter"),
  tid_length_range = .ss(NULL, "range", "filter"),
  trim_outliers    = .ss(FALSE, "flag", "filter"),
  trim_k           = .ss(4, "num", "filter", min = 2, max = 10),
  sd_min_locs      = .ss(.smlm_min_locs, "int", "filter", min = 2),
  view_SD          = .ss(FALSE, "flag", "filter"),
  target_sd_x      = .ss(NA_real_, "num", "filter", min = 0, na_ok = TRUE),
  target_sd_y      = .ss(NA_real_, "num", "filter", min = 0, na_ok = TRUE),

  # ── Channels: what decides which points may never be grouped together.
  enable_two_color = .ss(FALSE, "flag", "channels"),
  color_by_fluo    = .ss(FALSE, "flag", "channels"),
  use_pooled_dcr   = .ss(FALSE, "flag", "channels"),
  dcr_group1       = .ss(c(0, 0.369), "range", "channels"),
  dcr_group2       = .ss(c(0.3691, 1), "range", "channels"),

  # ── De-blinking.
  enable_clustering     = .ss(FALSE, "flag", "deblink"),
  deblink_method        = .ss("centroid", "choice", "deblink", choices = c("centroid", "vote")),
  eps_centroid          = .ss(10, "num", "deblink", min = 0.1),
  centroid_k            = .ss(2.8, "num", "deblink", min = 0.5, max = 6),
  deblink_anisotropic   = .ss(FALSE, "flag", "deblink"),
  deblink_shared_spread = .ss(FALSE, "flag", "deblink"),
  vote_frac             = .ss(0.3, "num", "deblink", min = 0.02, max = 0.9),
  vote_rho              = .ss(3.0, "num", "deblink", min = 0.5, max = 8),
  vote_link             = .ss("mutual", "choice", "deblink", choices = c("mutual", "either")),

  # ── Median clustering (cluster2).
  enable_cluster2        = .ss(FALSE, "flag", "cluster2"),
  cluster2_method        = .ss("dbscan", "choice", "cluster2", choices = c("dbscan", "hdbscan", "optics")),
  eps2                   = .ss(200, "num", "cluster2", min = 3),
  # min 4: the server raises anything lower, because hulls and shape fits need 4 points
  minpts2                = .ss(15, "int", "cluster2", min = 4),
  cluster2_xi            = .ss(0.05, "num", "cluster2", min = 0.001, max = 0.5),
  cluster2_prefilter     = .ss(FALSE, "flag", "cluster2"),
  cluster2_prefilter_r   = .ss(200, "num", "cluster2", min = 1),
  cluster2_prefilter_min = .ss(3, "int", "cluster2", min = 1),
  cluster2_prune         = .ss(FALSE, "flag", "cluster2"),
  cluster2_prune_tol     = .ss(2, "num", "cluster2", min = 1.1),
  cluster2_core_only     = .ss(FALSE, "flag", "cluster2"),
  cluster2_min_members   = .ss(0, "int", "cluster2", min = 0),
  # "within" = every channel separately; "only:<k>" = one channel, for inspection
  cluster_group_scope    = .ss("within", "choice", "cluster2"),

  # ── Opt-in gates that rewrite the data before clustering.
  fid_exclude      = .ss(FALSE, "flag", "gates"),
  dens_norm_enable = .ss(FALSE, "flag", "gates"),
  dens_thin_enable = .ss(FALSE, "flag", "gates"),
  dens_norm_unit   = .ss("auto", "choice", "gates", choices = c("auto", "cluster", "tid"))
)
rm(.ss)

# GOAL: coerce one value to its setting's type and limits. Returns the default when
# the value is missing or unusable, and says so through `warn` when the caller gave
# a value that had to be replaced or clamped (a script user should hear about it;
# a widget that has not rendered yet should not). `clamp = FALSE` keeps an
# out-of-range number as entered: values read from the app are used as the user
# typed them, exactly as before this file existed; only code-built settings clamp.
.smlm_coerce <- function(id, v, warn = TRUE, clamp = TRUE) {
  sp  <- SMLM_SETTINGS_SPEC[[id]]
  def <- sp$default
  say <- function(msg) if (warn) warning(sprintf("smlm_settings: %s %s", id, msg), call. = FALSE)
  if (is.null(v)) return(def)
  switch(sp$type,
    flag = isTRUE(as.logical(v)[1]),
    num = , int = {
      x <- suppressWarnings(as.numeric(v)[1])
      if (!length(x) || is.na(x)) {
        if (sp$na_ok) return(NA_real_)
        say("is not a number; using the default"); return(def)
      }
      if (clamp && (x < sp$min || x > sp$max)) {
        say(sprintf("= %s is outside [%s, %s]; clamped", format(x), sp$min, sp$max))
        x <- min(max(x, sp$min), sp$max)
      }
      if (sp$type == "int") as.integer(round(x)) else x
    },
    range = {
      x <- suppressWarnings(as.numeric(v))
      if (length(x) != 2L) { say("needs two values; using the default"); return(def) }
      x
    },
    choice = {
      x <- as.character(v)[1]
      if (identical(id, "cluster_group_scope")) {
        if (!identical(x, "within") && !grepl("^only:[0-9]+$", x)) {
          say("must be \"within\" or \"only:<channel>\"; using the default"); return(def)
        }
        return(x)
      }
      if (!x %in% sp$choices) {
        say(sprintf("= \"%s\" is not one of %s; using the default", x,
                    paste(sp$choices, collapse = ", ")))
        return(def)
      }
      x
    })
}

# Defaults go through the same coercion as every value read later, so a default
# and an equal value from the app are identical objects (15L, not 15). A cache key
# built from the settings would otherwise miss when a widget reports a value equal
# to the default it never changed.
for (.id in names(SMLM_SETTINGS_SPEC))
  SMLM_SETTINGS_SPEC[[.id]]["default"] <- list(.smlm_coerce(.id, SMLM_SETTINGS_SPEC[[.id]]$default,
                                                            warn = FALSE, clamp = FALSE))
rm(.id)
SMLM_DEFAULTS <- lapply(SMLM_SETTINGS_SPEC, `[[`, "default")

# =============================================================================
# smlm_settings() - THE SETTINGS FOR ONE SMLM RUN, FROM CODE
# =============================================================================
# GOAL: every data-affecting setting at its app default, with the ones passed in
# `...` replacing them. Values are coerced and checked; an unknown name is an error,
# because a typo would otherwise leave the default in place without a word.
#
# ARGS     ... named settings, e.g. smlm_settings(enable_clustering = TRUE, eps2 = 150)
# RETURNS  a list of class "smlm_settings" with one entry per SMLM_SETTINGS_SPEC row
# =============================================================================
smlm_settings <- function(...) {
  over <- list(...)
  if (length(over) && (is.null(names(over)) || any(!nzchar(names(over)))))
    stop("smlm_settings: every setting must be named", call. = FALSE)
  bad <- setdiff(names(over), names(SMLM_SETTINGS_SPEC))
  if (length(bad))
    stop("smlm_settings: unknown setting(s): ", paste(bad, collapse = ", "), call. = FALSE)
  out <- SMLM_DEFAULTS
  for (id in names(over)) out[id] <- list(.smlm_coerce(id, over[[id]], warn = TRUE))
  structure(out, class = c("smlm_settings", "list"))
}

# =============================================================================
# smlm_settings_from_input() - THE SAME SETTINGS, READ FROM THE APP
# =============================================================================
# GOAL: build the settings object from a Shiny `input` (or any list with the same
# names). A field that is NULL - its widget has not rendered yet, or a saved file
# predates it - takes the default from SMLM_DEFAULTS, the same value the widget
# would show.
#
# REACTIVITY: reading every field creates a dependency on every field. Call it
# inside an observeEvent/eventReactive handler or isolate(), so that a change to,
# say, vote_frac while the centroid method is selected does not invalidate
# anything.
#
# ARGS     input : Shiny input, reactivevalues inside isolate(), or a plain list
# RETURNS  a list of class "smlm_settings"
# =============================================================================
smlm_settings_from_input <- function(input) {
  out <- SMLM_DEFAULTS
  for (id in names(SMLM_SETTINGS_SPEC))
    out[id] <- list(.smlm_coerce(id, input[[id]], warn = FALSE, clamp = FALSE))
  structure(out, class = c("smlm_settings", "list"))
}

# =============================================================================
# as_smlm_settings() - THE ADAPTER THAT LETS A HELPER TAKE EITHER FORM
# =============================================================================
# GOAL: let every data-chain helper take `settings` while its callers still pass
# the Shiny `input`, so the server's call sites can move over one at a time.
#
# HOW: a settings object passes through unchanged. Anything else (Shiny `input`,
# reactivevalues, a plain list) is wrapped in a VIEW that reads a field only when
# the helper asks for it, through the same coercion and defaults as
# smlm_settings_from_input(). Reading lazily matters: a reactive that calls the
# pipeline then depends only on the inputs the pipeline actually used, as it did
# when the helper read `input$...` directly. Reading all 45 up front would make
# every such reactive re-run when any setting changes.
#
# Asking a view or a settings object for a name that is not an SMLM setting is an
# error, so a helper cannot read a field the settings (and any cache key built on
# them) do not know about.
#
# ARGS     x : smlm_settings, Shiny input / reactivevalues, list, or NULL (defaults)
# RETURNS  something `$` and `[[` read settings from
# =============================================================================
as_smlm_settings <- function(x) {
  if (inherits(x, c("smlm_settings", "smlm_settings_view"))) return(x)
  if (is.null(x)) return(smlm_settings())
  structure(list(), class = "smlm_settings_view", src = x, over = list())
}

# One field of a settings object or view; the single read path for both.
.smlm_get <- function(x, name) {
  if (!name %in% names(SMLM_SETTINGS_SPEC))
    stop(sprintf("\"%s\" is not an SMLM setting (see SMLM_SETTINGS_SPEC)", name), call. = FALSE)
  if (inherits(x, "smlm_settings")) return(.subset2(x, name, exact = TRUE))
  ov <- attr(x, "over")
  if (name %in% names(ov)) return(ov[[name]])
  .smlm_coerce(name, attr(x, "src")[[name]], warn = FALSE, clamp = FALSE)
}

# Exact-name access that refuses unknown names. A list's own `$` partially
# matches, so `s$cluster2_prefilter` would quietly return `cluster2_prefilter_r`
# if the flag were ever missing.
`$.smlm_settings`       <- function(x, name) .smlm_get(x, name)
`[[.smlm_settings_view` <- function(x, i, ...) .smlm_get(x, i)
`$.smlm_settings_view`  <- function(x, name) .smlm_get(x, name)

# GOAL: a copy of `s` with some settings replaced, for a helper that needs to force
# one (assign_two_color_group switches the colour split on). Works on settings
# objects and views alike; a view keeps reading everything else lazily.
smlm_set <- function(s, ...) {
  v <- list(...)
  bad <- setdiff(names(v), names(SMLM_SETTINGS_SPEC))
  if (length(bad)) stop("smlm_set: unknown setting(s): ", paste(bad, collapse = ", "), call. = FALSE)
  s <- as_smlm_settings(s)
  if (inherits(s, "smlm_settings_view")) {
    ov <- attr(s, "over")
    for (id in names(v)) ov[id] <- list(.smlm_coerce(id, v[[id]]))
    attr(s, "over") <- ov
  } else {
    for (id in names(v)) s[id] <- list(.smlm_coerce(id, v[[id]]))
  }
  s
}

# =============================================================================
# smlm_settings_key() - ONE HASH FOR "THESE SETTINGS"
# =============================================================================
# GOAL: the part of a cache key that the settings decide. Two settings objects
# with the same values give the same hash, however the values arrived (15, 15L
# and "15" all coerce to 15L), and changing any setting changes it.
#
# ARGS     s : smlm_settings, or anything as_smlm_settings accepts (read in full)
# RETURNS  a short hex string
# =============================================================================
smlm_settings_key <- function(s) {
  if (!inherits(s, "smlm_settings")) s <- smlm_settings_from_input(s)
  digest::digest(unclass(s)[names(SMLM_SETTINGS_SPEC)], algo = "xxhash64")
}

# =============================================================================
# SETTINGS FILES - the smlm_settings block, and the layout older files used
# =============================================================================
# A settings file (Download settings in the SMLM tab) carries the complete
# settings object as `smlm_settings`, so every setting that changes results is
# saved by construction. Files written before that block existed kept most
# settings at the top level under their input ID, and a few in sub-lists; the
# paths below say where. Settings that are absent from such a file are reported
# as not present, so a loader leaves them alone rather than resetting them.
SMLM_SETTINGS_FILE_VERSION <- 2L
.SMLM_LEGACY_PATH <- list(
  enable_two_color = c("two_color", "enable"),
  color_by_fluo    = c("two_color", "use_fluo"),
  dcr_group1       = c("two_color", "dcr_group1"),
  dcr_group2       = c("two_color", "dcr_group2"),
  use_pooled_dcr   = c("analysis", "use_pooled_dcr"),
  dens_norm_enable = c("analysis", "dens_norm_enable"),
  dens_thin_enable = c("analysis", "dens_thin_enable"),
  fid_exclude      = c("fiducial", "fid_exclude")
)
.smlm_legacy_path <- function(id) .SMLM_LEGACY_PATH[[id]] %||% id

# Value at a path of nested names, or NULL when any level is missing.
.smlm_dig <- function(x, path) {
  for (p in path) { if (!is.list(x) || is.null(x[[p]])) return(NULL); x <- x[[p]] }
  x
}

# =============================================================================
# smlm_settings_from_file() - THE SETTINGS IN A LOADED SETTINGS FILE
# =============================================================================
# GOAL: read the settings out of a parsed settings file (jsonlite::fromJSON), in
# either layout, and say which ones the file actually holds.
#
# ARGS     raw : the parsed file (a list)
# RETURNS  list(settings = smlm_settings object, present = character ids,
#               layout = "block" | "legacy")
#          Settings not present take the defaults in `settings`; a loader should
#          apply only `present`.
# =============================================================================
smlm_settings_from_file <- function(raw) {
  if (is.list(raw$smlm_settings)) {
    blk <- raw$smlm_settings
    present <- intersect(names(blk), names(SMLM_SETTINGS_SPEC))
    return(list(settings = smlm_settings_from_input(blk), present = present, layout = "block"))
  }
  src <- list()
  for (id in names(SMLM_SETTINGS_SPEC)) {
    v <- .smlm_dig(raw, .smlm_legacy_path(id))
    if (!is.null(v)) src[[id]] <- v
  }
  list(settings = smlm_settings_from_input(src), present = names(src), layout = "legacy")
}

# GOAL: write settings into the older file layout, so code that restores from
# those places (the filter, crop and QC restore in server_smlm.R) sees the values
# of a new-layout file too.
#   raw : parsed file; s : settings; ids : which settings to write
smlm_settings_to_legacy <- function(raw, s, ids) {
  for (id in ids) {
    path <- .smlm_legacy_path(id)
    v <- .subset2(s, id)
    if (length(path) == 1L) raw[path] <- list(v)
    else {
      sub <- if (is.list(raw[[path[1]]])) raw[[path[1]]] else list()
      sub[path[2]] <- list(v)
      raw[[path[1]]] <- sub
    }
  }
  raw
}

# =============================================================================
# smlm_settings_restore_plan() - HOW TO PUT SETTINGS BACK INTO THE WIDGETS
# =============================================================================
# GOAL: one widget update per setting, in an order that works: every parameter
# before the switch that puts it into effect (enable_clustering,
# enable_cluster2, enable_two_color last), so a restored analysis never runs for
# a moment with the old parameters.
#
# The widget kind follows the setting's type; the two choices shown as drop-downs
# are named here, every other choice is a radio group.
#
# ARGS     s : settings; ids : which settings to restore
# RETURNS  data.frame(id, widget) in application order, with the values in
#          attr "values" (a list, since ranges are length 2)
# =============================================================================
smlm_settings_restore_plan <- function(s, ids) {
  ids <- intersect(names(SMLM_SETTINGS_SPEC), ids)       # spec order, unknowns dropped
  last <- c("enable_two_color", "enable_clustering", "enable_cluster2")
  ids <- c(setdiff(ids, last), intersect(last, ids))
  selects <- c("cluster_group_scope", "dens_norm_unit")
  kind <- vapply(ids, function(id) switch(SMLM_SETTINGS_SPEC[[id]]$type,
    flag = "checkbox", num = "numeric", int = "numeric", range = "range",
    choice = if (id %in% selects) "select" else "radio"), character(1))
  plan <- data.frame(id = ids, widget = unname(kind), stringsAsFactors = FALSE)
  attr(plan, "values") <- lapply(ids, function(id) .subset2(s, id))
  plan
}

# =============================================================================
# DISPLAY OPTIONS - how points are drawn, never what is computed
# =============================================================================
# Kept apart from the settings above on purpose: they change the picture, not a
# number, so they stay out of the cache key and a palette change never
# re-clusters. The plot helpers (build_pts_2d/3d, build_colors_for_pts,
# color_mapping_for_pts, attach_channel_colors) take these as plain arguments,
# and the UI widgets read their initial values from here.
SMLM_DISPLAY_DEFAULTS <- list(
  color_mode              = "default",   # default / multiplex / unicolor / z / tim / locdens
  numeric_palette         = "Spectral",
  numeric_palette_reverse = FALSE,
  unicolor_hex            = "#1f77b4",
  dens_radius_nm          = 500,         # radius for the trace-median density colouring
  show_medians            = FALSE,       # draw trace (or molecule) medians instead of localizations
  color_group1            = "#1f77b4",   # channel A
  color_group2            = "#d62728"    # channel B
)

# GOAL: the display options from a Shiny `input` (or a partial list), each
# missing one at its default. Reads only these eight inputs.
smlm_display_from_input <- function(input) {
  stats::setNames(lapply(names(SMLM_DISPLAY_DEFAULTS), function(id)
    input[[id]] %||% SMLM_DISPLAY_DEFAULTS[[id]]), names(SMLM_DISPLAY_DEFAULTS))
}

# GOAL: what the plot helpers call on their `display` argument: a complete list,
# whether they were given display options, a partial list, or Shiny's input.
.smlm_display <- function(display) {
  if (is.null(display)) return(SMLM_DISPLAY_DEFAULTS)
  if (inherits(display, "reactivevalues")) return(smlm_display_from_input(display))
  utils::modifyList(SMLM_DISPLAY_DEFAULTS, as.list(display)[intersect(names(display),
                                                                       names(SMLM_DISPLAY_DEFAULTS))])
}

# The two channel colours, A then B, as attach_channel_colors() takes them.
smlm_channel_colors <- function(display = NULL) {
  d <- .smlm_display(display)
  c(d$color_group1, d$color_group2)
}

# Print only what differs from the defaults: a settings object is mostly defaults.
print.smlm_settings <- function(x, ...) {
  diff <- names(x)[!vapply(names(x), function(id)
    identical(x[[id]], SMLM_DEFAULTS[[id]]), logical(1))]
  cat(sprintf("<smlm_settings: %d settings, %d changed from the defaults>\n",
              length(x), length(diff)))
  for (id in diff)
    cat(sprintf("  %-24s %s\n", id, paste(format(x[[id]]), collapse = ", ")))
  invisible(x)
}
