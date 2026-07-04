# MINFLUXR Plugins — how to build your own

MINFLUXR is a Shiny app made of a few core tabs (SMLM Analysis, Tracking
Analysis, MINFLUX Simulator). The **plugin system** lets you add your own tab
with your own analysis **without editing any core file**. You write a small
self-contained folder; the app discovers and loads it at startup.

This folder (`plugin/`) is the **template**. It is not loaded by the app — it is
the thing you copy. `plugin/plugin.R` is a working example you can run as-is and
then modify.

---

## 1. Quick start (5 minutes)

1. **Copy** this `plugin/` folder to `plugins/myplugin/`
   (the app loads from `plugins/`, plural — see [§4](#4-why-plugin-and-plugins)).

   ```
   MINFLUXR/
   ├── plugins/
   │   └── myplugin/
   │       └── plugin.R      ← your copy
   └── plugin/               ← the untouched template (this folder)
   ```

2. **Run the app** (`app.R`, or the "Run App" button). A new **"My Plugin"** tab
   appears at the end of the navbar.

3. Load a file in the **SMLM Analysis** tab, switch to **My Plugin**, click the
   button — it reports how many localizations are loaded. That round-trip proves
   your plugin can read the core app's live data.

4. **Rename and edit.** Change the title, the `value` id, the `myplugin_*` input
   ids, and the body of the server function. Rebuild it into your own analysis.

To **remove** a plugin, delete its folder under `plugins/`. Nothing else changes.

---

## 2. The three things every plugin does

A plugin is just an R script that, when sourced, defines a **UI tab** and a
**server function**, then **registers** them. Look at `plugin.R` alongside this
list — the sections are numbered to match.

### 2.1 Define a UI tab

A tab is a Shiny `tabPanel()`. It gets appended to the main navbar.

```r
my_ui_tab <- tabPanel(
  title = "My Plugin",       # navbar label
  value = "tab_myplugin",    # UNIQUE internal id (see naming rules below)
  fluidPage( ... )           # your layout
)
```

### 2.2 Define a server function

Your server function receives the live Shiny `input`, `output`, `session`, plus
`smlm` — a bundle of the core SMLM tab's data (see [§3](#3-the-smlm-data-handle)).

```r
my_server <- function(input, output, session, smlm) {
  output$myplugin_summary <- renderPrint({ ... })
}
```

### 2.3 Register both

The single required call. Everything above is inert until you do this:

```r
register_plugin(
  tabs      = list(my_ui_tab),   # one or more tabPanels
  server_fn = my_server          # your server function
)
```

`register_plugin(tabs, server_fn)` is provided by the core app (defined in
`global.R`). It appends your tabs to the navbar and arranges for `server_fn` to
run inside the main Shiny session.

---

## 3. The `smlm` data handle

The core app calls your server function like this (in `server.R`):

```r
.srv(input, output, session, smlm = smlm_exports)
```

So the last argument, `smlm`, is a named list of **reactive handles** exported by
the SMLM tab. Call them like functions (`smlm$plottable_data()`) inside a
reactive context to read the user's current data.

| Handle | Call | What you get |
|---|---|---|
| `smlm$plottable_data` | `smlm$plottable_data()` | Data frame of the **currently displayed** localizations (after filters/crop). Has `x`, `y`, and usually `z` columns in nm. **Start here.** |
| `smlm$uploaded_df` | `smlm$uploaded_df()` | The **raw** uploaded data frame, before filtering. |
| `smlm$sd_filtered_tbl` | `smlm$sd_filtered_tbl()` | Per-trace spread table (one row per trace id) after current filters. |
| `smlm$sd_all_tbl` | `smlm$sd_all_tbl()` | Per-trace spread table for all traces, unfiltered. |
| `smlm$file_store` | `smlm$file_store()` | Named list: `filename → data frame`, every file the user loaded. |
| `smlm$sd_filt` | `smlm$sd_filt$x`, `$y`, `$z` | Committed spread-filter ranges (each `c(min, max)`). |
| `smlm$sd_display_ranges` | `smlm$sd_display_ranges()` | Effective spread ranges with `Inf` replaced by the data max. |
| `smlm$meta_full_rng` | `smlm$meta_full_rng()` | Full metadata value ranges. |

**Always guard against "no data yet."** Reactives error before a file is loaded,
so wrap reads in `shiny::validate(shiny::need(...))` (shows a friendly message)
or `req(...)` (silently waits):

```r
df <- smlm$plottable_data()
shiny::validate(shiny::need(is.data.frame(df) && nrow(df) > 0,
                            "Load a file in the SMLM Analysis tab first."))
```

> These handles are **read-only** for your purposes — treat them as inputs. Don't
> try to write back into the core app's state; build your own outputs instead.

---

## 4. Why `plugin/` and `plugins/`?

- **`plugins/`** (plural) is where **active** plugins live. `global.R` loops over
  `plugins/*/plugin.R` at startup and sources each one.
- **`plugin/`** (singular, this folder) is the **template**. It is deliberately
  *not* under `plugins/`, so the app never loads it. It exists only to be copied.

This means you can keep the template around for reference and copy it as many
times as you like into `plugins/`.

---

## 5. Naming rules (avoid clashes)

All plugins and the core app share **one Shiny namespace**. Two things must be
globally unique or the app will misbehave:

1. **Tab `value` ids** — don't reuse the core ids (`tab_smlm`, `tab_tracking`,
   `tab_sim`). Use something like `tab_myplugin`.
2. **Every `input$...` and `output$...` id you create** — prefix them all with
   your plugin name, e.g. `myplugin_go`, `myplugin_summary`. This is the single
   most common source of plugin bugs.

R object names in your file (`my_ui_tab`, `my_server`) also live in the global
environment, so give them plugin-specific names too if you add several plugins.

---

## 6. Growing beyond one file

For anything non-trivial, keep `plugin.R` as a thin entry point and split the
rest into sibling files. Create them, then uncomment the matching `source()`
lines near the top of `plugin.R` (`pd` is your plugin folder path):

```
plugins/myplugin/
├── plugin.R      ← entry point: source()s the others, then register_plugin()
├── packages.R    ← install/load any extra CRAN packages your plugin needs
├── helpers.R     ← pure functions (no Shiny), easy to test in isolation
├── ui.R          ← defines your tabPanel(s)
└── server.R      ← defines your server function
```

```r
# in plugin.R
pd <- .plugin_dir
source(file.path(pd, "packages.R"))
source(file.path(pd, "helpers.R"))
source(file.path(pd, "ui.R"))
source(file.path(pd, "server.R"))
register_plugin(tabs = list(my_ui_tab), server_fn = my_server)
```

**Extra packages:** put installs in `packages.R` and guard them so they only run
when missing:

```r
for (p in c("somePackage")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
```

**Bundled resources** (a lookup table, an image): keep them in your folder and
load by relative path, e.g. `read.csv(file.path(pd, "data/table.csv"))`.

---

## 7. Multiple tabs from one plugin

`register_plugin(tabs = ...)` takes a list, so one plugin can contribute several
tabs:

```r
register_plugin(
  tabs      = list(my_ui_tab, my_settings_tab),
  server_fn = my_server
)
```

All the tabs are driven by the same `server_fn`.

---

## 8. Checklist before you ship

- [ ] Folder copied to `plugins/<name>/` and app shows the new tab.
- [ ] Tab `value` id is unique (not `tab_smlm` / `tab_tracking` / `tab_sim`).
- [ ] Every input/output id is prefixed with the plugin name.
- [ ] Every read of `smlm$...()` is guarded with `validate`/`req`.
- [ ] Plugin behaves sensibly **before** any file is loaded (no red errors).
- [ ] Extra packages are installed from `packages.R`, guarded by `requireNamespace`.
- [ ] Removing the folder cleanly removes the tab with no side effects.

---

## 9. Reference: how the core app wires plugins in

You don't need to edit these, but it helps to know where the hooks are:

- **`global.R`** — defines `register_plugin()`, then discovers and sources every
  `plugins/*/plugin.R` (setting `.plugin_dir` first). Collects tabs into
  `.plugin_tabs` and server functions into `.plugin_servers`.
- **`ui.R`** — appends `.plugin_tabs` to the `navbarPage` after the core tabs.
- **`server.R`** — after starting the core servers, calls each plugin server:
  `.srv(input, output, session, smlm = smlm_exports)`.

That's the entire contract. If your `plugin.R` defines a tab, defines a server
function with signature `function(input, output, session, smlm)`, and calls
`register_plugin()`, it will load.
