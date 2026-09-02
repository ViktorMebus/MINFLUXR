# ── check_dependencies.R ──────────────────────────────────────────────────────
# Run with:  Rscript tools/check_dependencies.R      (from the repo root)
#
# GOAL: catch the failure mode where the app calls pkg::something but packages.R
# never declares `pkg`, so the app installs "successfully" on a clean machine and
# then dies at runtime the first time that feature is opened. That is how `DT`
# escaped: nothing else in the dependency tree pulls it, so the SMLM tables threw
# "there is no package called 'DT'" only when a table first rendered.
#
# It compares two things and reports the difference:
#   referenced  every `pkg::` and `library(pkg)` in the SOURCED app files
#               (deadcode_* / retired_* are excluded - global.R never sources them)
#   declared    the pkgs_* vectors in packages.R, obtained by PARSING the file
#               rather than pattern-matching it, so multi-line vectors with
#               comments in the middle are read correctly
#
# Exits non-zero if anything is referenced but not declared, so it can be wired
# into CI or a pre-commit hook.

files <- c("helpers_shared.R", "helpers_smlm.R", "helpers_tracking.R",
           "ui_smlm.R", "ui_tracking.R", "sim_ui.R", "ui.R",
           "server_smlm.R", "server_tracking.R", "sim_server.R", "server.R", "global.R")
files <- files[file.exists(files)]
if (!length(files)) stop("Run this from the repository root.")

# ── referenced ──────────────────────────────────────────────────────────
# Uses R's own parser rather than regex. getParseData() tags the `pkg` in `pkg::fn`
# as SYMBOL_PACKAGE, so strings can never be mistaken for code - which matters here
# because sim_ui.R embeds CSS containing `summary::before` inside a multi-line
# string, and any line-based regex reads that as a package called `summary`.
refs <- list()
for (f in files) {
  pd <- tryCatch(utils::getParseData(parse(f, keep.source = TRUE)), error = function(e) NULL)
  if (is.null(pd)) { message("could not parse ", f); next }
  hit <- unique(pd$text[pd$token == "SYMBOL_PACKAGE"])
  # library(x) / require(x): the argument is a plain symbol or a string
  lib <- character(0)
  fn  <- which(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% c("library", "require"))
  for (i in fn) {
    sib <- pd[pd$parent == pd$parent[i] & pd$token %in% c("SYMBOL", "STR_CONST"), "text"]
    lib <- c(lib, gsub('^"|"$', "", sib))
  }
  for (h in unique(c(hit, lib))) refs[[h]] <- c(refs[[h]], f)
}

# ── declared ──────────────────────────────────────────────────────────────────
e <- new.env()
for (ex in parse("packages.R")) {
  if (is.call(ex) && identical(as.character(ex[[1]]), "<-") &&
      grepl("^pkgs_", as.character(ex[[2]])))
    tryCatch(eval(ex, e), error = function(err) NULL)
}
declared <- sort(unique(unlist(mget(ls(e), envir = e))))

base <- rownames(installed.packages(priority = c("base", "recommended")))
base <- unique(c(base, "base", "utils", "stats", "grDevices", "graphics", "methods",
                 "tools", "parallel", "compiler", "datasets"))
IGNORE <- c(base, "pkg")   # `pkg` is the installer's own loop variable

used    <- setdiff(names(refs), IGNORE)
missing <- setdiff(used, declared)

cat(sprintf("declared: %d   referenced: %d\n\n", length(declared), length(used)))
if (length(missing)) {
  cat("REFERENCED BUT NOT DECLARED - add these to packages.R:\n")
  for (m in sort(missing))
    cat(sprintf("  %-14s used in: %s\n", m, paste(unique(refs[[m]]), collapse = ", ")))
} else {
  cat("No gaps: everything the app calls is declared.\n")
}

notinst <- declared[!vapply(declared, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (length(notinst))
  cat("\nDeclared but not installed in this library:\n  ", paste(notinst, collapse = ", "), "\n")

quit(status = if (length(missing)) 1L else 0L)
