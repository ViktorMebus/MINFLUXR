## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## deadcode_smlm.R  —  retired SMLM code, kept for reference only
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##
## WHAT THIS IS
##   Blocks removed from server_smlm.R because they were SHADOWED: the same
##   name was assigned again later in the same function body, so R silently
##   replaced these definitions before any of them could ever run.
##
## THIS FILE IS NOT SOURCED. global.R sources by explicit name and this file is
## deliberately not in that list.
##
## WHY THEY WERE DEAD (this is a shadowing bug, not merely unused code)
##   In R, `x <- 1; x <- 2` leaves x = 2. Shiny server bodies are ordinary
##   function bodies, so assigning a reactive twice keeps only the LAST one.
##   All of these sat at the top level of smlm_server() and were re-assigned
##   further down the same file:
##
##     nn_table                   was 2504-2570, shadowed by the def at 4639
##     secondary_nn_lines         was 2573-2611, shadowed by the def at 4735
##     output$download_nn_table   was 2614-2618, shadowed by the def at 4728
##     output$download_nn_table2  was 2621-2625, shadowed by the def at 4844
##     output$download_cluster2_table was 2737-2762, shadowed by 4828
##
##   Verified with getParseData(): all the duplicate pairs share the same
##   parent block (the body of smlm_server), so the shadowing is unconditional.
##
## WHY THIS ONE WAS ACTIVELY MISLEADING
##   These blocks sat under a banner reading "2D-Only server parts" and were
##   commented "Enhanced NN table (2D)". Reading the file, you would reasonably
##   conclude that 2D used these and 3D used the later ones. It never did: the
##   LATER definitions handled BOTH cases and were the only ones that ran.
##   The surviving versions are strictly more capable — they add the z / 3D
##   median, the multiplex channel scope (nn_group_scope), analysis-group
##   handling and cluster2_color, none of which these older copies had.
##
## Removing them changes NO behaviour: they never executed.
##
## HOW TO RESTORE
##   These live inside smlm_server(input, output, session). The wrapper below is
##   never called and exists only so this file is valid R. To restore, copy a
##   block back into server_smlm.R — but note you would then have to DELETE the
##   later definition, or you would simply recreate the shadowing.
##
## Removed 2026-07-20.
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.deadcode_smlm_server_blocks <- function(input, output, session) {

  # ---- Enhanced NN table (2D) with metadata ----
  nn_table <- reactive({
    req(plot_trigger())
    df <- plottable_data()
    
    if (isTRUE(input$enable_clustering)) {
      # Use cluster medians
      df_geo <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          id        = dplyr::first(cluster),
          tid       = dplyr::first(tid),
          x         = stats::median(x, na.rm = TRUE),
          y         = stats::median(y, na.rm = TRUE),
          cluster2  = dplyr::first(cluster2),
          color = dplyr::first(color),
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups   = "drop"
        )
      id_label <- "cluster"
    } else {
      # Use tid medians
      df_geo <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          id        = dplyr::first(tid),
          x         = stats::median(x, na.rm = TRUE),
          y         = stats::median(y, na.rm = TRUE),
          cluster   = dplyr::first(cluster),
          cluster2  = dplyr::first(cluster2),
          color = dplyr::first(color),
          mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
          .groups   = "drop"
        )
      id_label <- "tid"
    }

    if (nrow(df_geo) < 2) return(df_geo[0, ])

    coords <- as.matrix(df_geo[, c("x","y")])
    pm <- .nn_pairs_scoped(df_geo, coords, input)

    rows <- list()
    for (r in seq_len(nrow(pm))) {
      i <- pm[r, 1]; ni <- pm[r, 2]; dist <- pm[r, 3]
      if (dist >= input$distance_range[1] && dist <= input$distance_range[2]) {
        rows[[length(rows) + 1]] <- data.frame(
          source_id = df_geo$id[i],
          target_id = df_geo$id[ni],
          distance  = round(dist, 2),
          source_x  = df_geo$x[i],
          source_y  = df_geo$y[i],
          target_x  = df_geo$x[ni],
          target_y  = df_geo$y[ni],
          source_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[i] else NA,
          source_cluster2 = df_geo$cluster2[i],
          target_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[ni] else NA,
          target_cluster2 = df_geo$cluster2[ni],
          stringsAsFactors = FALSE
        )
      }
    }
    out <- if (length(rows)) dplyr::bind_rows(rows) else df_geo[0, ]
    names(out)[names(out) == "source_id"] <- paste0("source_", id_label)
    names(out)[names(out) == "target_id"] <- paste0("target_", id_label)
    out
  })

  # ---- Secondary cluster median-to-median NN lines (green, 2D) ----
  secondary_nn_lines <- reactive({
    df <- plottable_data()
    req(isTRUE(input$enable_cluster2))
    
    df_geo <- df %>%
      dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
      dplyr::group_by(cluster2) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        color = dplyr::first(color),    
        .groups = "drop"
      )
    if (nrow(df_geo) < 2) return(list())
    
    coords <- as.matrix(df_geo[, c("x","y")])
    k_eff  <- min(input$nn_k_cluster2 + 1, nrow(coords))
    nn     <- RANN::nn2(coords, k = k_eff)
    
    lines <- list()
    for (i in 1:nrow(coords)) {
      for (j in 2:k_eff) {
        ni   <- nn$nn.idx[i, j]
        dist <- as.numeric(nn$nn.dists[i, j])
        if (dist >= input$distance_range_cluster2[1] && dist <= input$distance_range_cluster2[2]) {
          nn_off <- suppressWarnings(as.numeric(input$c2_label_offset %||% 5))
          if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
          midpoint <- perp_offset_xy(coords[i, ], coords[ni, ], nn_off)
          lines[[length(lines) + 1]] <- list(
            x = c(coords[i, 1], coords[ni, 1]),
            y = c(coords[i, 2], coords[ni, 2]),
            text = round(dist, 2),
            midpoint = midpoint
          )
        }
      }
    }
    lines
  })

  # ---- Downloads ----
  output$download_nn_table <- downloadHandler(
    filename = function() paste0("nn_distances_", Sys.Date(), ".xlsx"),
    content  = function(file) openxlsx::write.xlsx(nn_table(), file, row.names = FALSE),
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  
  output$download_nn_table2 <- downloadHandler(
    filename = function() paste0("nn_distances_cluster2_", Sys.Date(), ".xlsx"),
    content  = function(file) openxlsx::write.xlsx(nn_table_cluster2(), file, row.names = FALSE),
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )

  output$download_cluster2_table <- downloadHandler(
    filename = function() paste0("cluster2_metrics_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      tbl <- cluster2_table_store()
      if (is.null(tbl) || nrow(tbl) == 0) stop("No Median Cluster metrics to download.")
      
      membrane_mode <- isTRUE(input$membrane_protein_options)
      
      if (!membrane_mode) {
        # Drop membrane-specific columns if they exist
        drop_cols <- c(
          "perceived_radius", "bending_radius",
          "distance_to_plane", "angle_to_plane",
          "plane_normal_x", "plane_normal_y", "plane_normal_z"
        )
        tbl <- tbl[, setdiff(names(tbl), drop_cols), drop = FALSE]
      }
      
      openxlsx::write.xlsx(
        x         = tbl,
        file      = file,
        sheetName = "Median Cluster metrics",
        rowNames  = FALSE
      )
    }
  )

}  # end .deadcode_smlm_server_blocks (never called)

