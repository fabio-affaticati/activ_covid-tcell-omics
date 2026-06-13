  res_path_s <- "results/R_scripts_plots/scaled/"

  res_path_notjoined <- "results/R_scripts_plots/notjoined/"

  script_path_from_args <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) > 0) {
      return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
    }
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      active_path <- rstudioapi::getActiveDocumentContext()$path
      if (!is.null(active_path) && nzchar(active_path)) {
        return(normalizePath(active_path, mustWork = FALSE))
      }
    }
    NA_character_
  }

  script_path <- script_path_from_args()
  script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
  repo_root <- if (basename(script_dir) == "R_scripts") dirname(script_dir) else script_dir
  setwd(repo_root)

  # Create directories if they do not exist
  if (!file.exists(res_path_s)) dir.create(res_path_s, recursive = TRUE)
  if (!file.exists(res_path_notjoined)) dir.create(res_path_notjoined, recursive = TRUE)

  library(tidyverse)
  library(pheatmap)
  library(grid)
  library(gridExtra)
  library(circlize)

  source(file.path(repo_root, "R_scripts", "cytometry_metafeatures_functions.R"))


  set.seed(42)
  
  ### ----------------------------
  ### LOAD DATA
  ### ----------------------------
  
  unstim_erc_cd4_base <- read.csv("data/5_per_file_MC_freq_CD4_ERC_neg.csv")
  unstim_erc_cd8_base <- read.csv("data/5_per_file_MC_freq_CD8.csv")
  Meta <- read.csv("data/entero_clusters_cytometrycomparison.csv")
  
  
  
  ### ----------------------------
  ### PROCESS CD4 & CD8
  ### ----------------------------
  
  cd4_notjoined <- process_cytometry(unstim_erc_cd4_base) %>%
    filter(!str_detect(sample_id, "TYF")) %>%
    enforce_55_metaclusters()
  
  cd4 <- inner_join(cd4_notjoined, Meta, by = "sample_id") %>% 
    select(-CMVIgG)
  
  cd8_notjoined <- process_cytometry(unstim_erc_cd8_base) %>%
    filter(!str_detect(sample_id, "TYF")) %>%
    enforce_55_metaclusters()
  
  cd8 <- inner_join(cd8_notjoined, Meta, by = "sample_id") %>% 
    select(-CMVIgG) 
  
  
  ### ----------------------------
  ### COLOR DEFINITIONS
  ### ----------------------------
  
  cluster_colors <- c(
    "~Prevotella"      = "#e31a1c",
    "~Rumino&Christen" = "#FFA500",
    "~Bacteroidaceae"  = "#33a02c",
    "~Lachnospiraceae" = "#1f78b4"
  )
  ann_colors <- list(Cluster = cluster_colors)
  
  
  ### ----------------------------
  ### Loop all timepoints – filtered
  ### ----------------------------
  message("Starting SCALED heatmap generation...")
  
  clustering_results_scaled <- list()
  
  all_TP <- sort(unique(cd4$Timepoint))
  
  for (TP in all_TP) {
    
    message("Processing TP ", TP)
    
    cd4_tp <- cd4 %>% filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
    cd8_tp <- cd8 %>% filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
  
    if (nrow(cd4_tp) == 0 | nrow(cd8_tp) == 0) next
    
    # Build matrices
    mat4 <- make_matrix(cd4_tp)
    mat8 <- make_matrix(cd8_tp)
    
    # Build scaled matrices
    mat4_scaled <- t(scale(t(mat4)))
    mat8_scaled <- t(scale(t(mat8)))
    
    # remove attributes
    attr(mat4_scaled, "scaled:center") <- NULL; attr(mat4_scaled, "scaled:scale") <- NULL
    attr(mat8_scaled, "scaled:center") <- NULL; attr(mat8_scaled, "scaled:scale") <- NULL
    
    # Annotations
    ann4 <- cd4_tp %>% select(sample_id, Cluster) %>% distinct() %>% column_to_rownames("sample_id")
    ann8 <- cd8_tp %>% select(sample_id, Cluster) %>% distinct() %>% column_to_rownames("sample_id")
    
    # ----------------------------
    # Scaled heatmaps
    # ----------------------------
    ph4s <- pheatmap(mat4_scaled, annotation_col = ann4, annotation_colors = ann_colors,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_method = "ward.D2",
                     border_color = NA,
                     legend = TRUE, fontsize = 10, fontsize_row = 8, fontsize_col = 6,
                     main = paste0("CD4 Metaclusters SCALED (", TP, ")"), silent = TRUE)
    ph8s <- pheatmap(mat8_scaled, annotation_col = ann8, annotation_colors = ann_colors,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_method = "ward.D2",
                     border_color = NA,
                     legend = TRUE, fontsize = 10, fontsize_row = 8, fontsize_col = 6,
                     main = paste0("CD8 Metaclusters SCALED (", TP, ")"), silent = TRUE)
    
    
    # ===========================================
    # EXTRACT CLUSTERS FROM HEATMAP (k = 4)
    # ===========================================
    
    # CD4 clusters
    clust4 <- cutree(ph4s$tree_col, k = 4)
    df_clust4 <- data.frame(
      sample_id = names(clust4),
      cluster = clust4,
      Timepoint = TP,
      Lineage = "CD4"
    )
    
    # CD8 clusters
    clust8 <- cutree(ph8s$tree_col, k = 4)
    df_clust8 <- data.frame(
      sample_id = names(clust8),
      cluster = clust8,
      Timepoint = TP,
      Lineage = "CD8"
    )
    
    
    # Store in list (outside loop define: clustering_results <- list())
    clustering_results_scaled[[length(clustering_results_scaled)+1]] <- df_clust4
    clustering_results_scaled[[length(clustering_results_scaled)+1]] <- df_clust8
    # ===========================================
    
    out_file_s <- paste0(res_path_s, "metacluster_heatmap_CD4_CD8_scaled_filtered_", TP, ".png")
    png(out_file_s, width = 3600, height = 2400, res = 300)
    grid.arrange(ph4s$gtable, ph8s$gtable, ncol = 2)
    dev.off()
    
    message("Saved SCALED filtered paired heatmap: ", out_file_s)
  }
  
  message("All timepoints processed (scaled).")
  
  
  ### ----------------------------
  ### Loop all timepoints – filtered
  ### ----------------------------
  message("Starting SCALED heatmap generation...")
  
  clustering_results_scaled <- list()
  
  all_TP <- sort(unique(cd4$Timepoint))
  
  for (TP in all_TP) {
    
    message("Processing TP ", TP)
    
    cd4_tp <- cd4 %>% filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
    cd8_tp <- cd8 %>% filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
    
    if (nrow(cd4_tp) == 0 | nrow(cd8_tp) == 0) next
    
    # Build matrices
    mat4 <- make_matrix(cd4_tp)
    mat8 <- make_matrix(cd8_tp)
    
    # Build scaled matrices
    mat4_scaled <- t(scale(t(mat4)))
    mat8_scaled <- t(scale(t(mat8)))
    
    # remove attributes
    attr(mat4_scaled, "scaled:center") <- NULL; attr(mat4_scaled, "scaled:scale") <- NULL
    attr(mat8_scaled, "scaled:center") <- NULL; attr(mat8_scaled, "scaled:scale") <- NULL
    
    # Annotations
    ann4 <- cd4_tp %>% select(sample_id, Cluster) %>% distinct() %>% column_to_rownames("sample_id")
    ann8 <- cd8_tp %>% select(sample_id, Cluster) %>% distinct() %>% column_to_rownames("sample_id")
    
    # ----------------------------
    # CD4 Heatmap
    # ----------------------------
    col_hc4 <- hclust(dist(t(mat4_scaled)), method = "ward.D2")
    col_clusters4 <- cutree(col_hc4, k = 4)
    
    col_ha4 <- HeatmapAnnotation(
      Enterotype = ann4$Cluster,
      DendroCluster = factor(paste0("C", col_clusters4)),
      col = list(
        Enterotype = cluster_colors,
        DendroCluster = c("C1"="darkblue","C2"="darkgreen","C3"="darkred","C4"="purple")
      )
    )
    
    ht4 <- Heatmap(
      mat4_scaled,
      name = "Scaled value",
      top_annotation = col_ha4,
      cluster_rows = TRUE,
      cluster_columns = col_hc4,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_title = paste0("CD4 Metaclusters SCALED (", TP, ")"),
      row_names_gp = gpar(fontsize = 8),
      heatmap_legend_param = list(title = "Metacluster value")
    )
    
    # Extract clusters from dendrogram
    df_clust4 <- data.frame(
      sample_id = names(col_clusters4),
      cluster = col_clusters4,
      Timepoint = TP,
      Lineage = "CD4"
    )
    
    # ----------------------------
    # CD8 Heatmap
    # ----------------------------
    col_hc8 <- hclust(dist(t(mat8_scaled)), method = "ward.D2")
    col_clusters8 <- cutree(col_hc8, k = 4)
    
    col_ha8 <- HeatmapAnnotation(
      Enterotype = ann8$Cluster,
      DendroCluster = factor(paste0("C", col_clusters8)),
      col = list(
        Enterotype = cluster_colors,
        DendroCluster = c("C1"="darkblue","C2"="darkgreen","C3"="darkred","C4"="purple")
      )
    )
    
    ht8 <- Heatmap(
      mat8_scaled,
      name = "Scaled value",
      top_annotation = col_ha8,
      cluster_rows = TRUE,
      cluster_columns = col_hc8,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_title = paste0("CD8 Metaclusters SCALED (", TP, ")"),
      row_names_gp = gpar(fontsize = 8),
      heatmap_legend_param = list(title = "Metacluster value")
    )
    
    # Extract clusters from dendrogram
    df_clust8 <- data.frame(
      sample_id = names(col_clusters8),
      cluster = col_clusters8,
      Timepoint = TP,
      Lineage = "CD8"
    )
    
    # Save clustering results
    clustering_results_scaled[[length(clustering_results_scaled)+1]] <- df_clust4
    clustering_results_scaled[[length(clustering_results_scaled)+1]] <- df_clust8
    
    # Save heatmap as PNG
    out_file <- paste0(res_path_s, "Heatmap_CD4_CD8_scaled_", TP, ".png")
    png(out_file, width = 3600, height = 2400, res = 300)
    draw(ht4 + ht8, merge_legend = TRUE)
    dev.off()
    
    message("Saved SCALED paired heatmap: ", out_file)
  }
  
  
  message("All timepoints processed (scaled).")
  
  
  ### ============================================================
  ### HEATMAPS FOR NOT-JOINED MATRICES (SCALED ONLY)
  ### ============================================================
  
  
  message("Starting NOT-JOINED heatmap generation (SCALED ONLY)...")
  
  all_TP <- sort(unique(cd4_notjoined$Timepoint))
  
  
  clustering_results_notjoined <- list()
  
  for (TP in all_TP) {
    
    message("Processing NOT-JOINED TP ", TP)
    
    cd4_tp <- cd4_notjoined %>% 
      filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
    cd8_tp <- cd8_notjoined %>% 
      filter(Timepoint == TP) %>%
      remove_samples_percentage()
    
    if (nrow(cd4_tp) == 0 | nrow(cd8_tp) == 0) next
    
    # ---- matrices ----
    mat4 <- make_matrix(cd4_tp)
    mat8 <- make_matrix(cd8_tp)
    
    # ---- SCALED ----
    mat4_scaled <- t(scale(t(mat4)))
    mat8_scaled <- t(scale(t(mat8)))
    attr(mat4_scaled, "scaled:center") <- NULL; attr(mat4_scaled, "scaled:scale") <- NULL
    attr(mat8_scaled, "scaled:center") <- NULL; attr(mat8_scaled, "scaled:scale") <- NULL
    
    # ------------------------------------------------------------
    # SCALED NOT-JOINED HEATMAPS (NO COLOR BAR)
    # ------------------------------------------------------------
    
    ph4s <- pheatmap(mat4_scaled,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_method = "ward.D2",
                     border_color = NA,
                     legend = TRUE, fontsize = 10, fontsize_row = 8, fontsize_col = 6,
                     main = paste0("CD4 Metaclusters SCALED not-joined (", TP, ")"),
                     silent = TRUE)
    
    ph8s <- pheatmap(mat8_scaled,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_method = "ward.D2",
                     border_color = NA,
                     legend = TRUE, fontsize = 10, fontsize_row = 8, fontsize_col = 6,
                     main = paste0("CD8 Metaclusters SCALED not-joined (", TP, ")"),
                     silent = TRUE)
    
    
    
    # ===========================================
    # EXTRACT CLUSTERS FROM HEATMAP (k = 4)
    # ===========================================
    
    # CD4 clusters
    clust4 <- cutree(ph4s$tree_col, k = 4)
    df_clust4 <- data.frame(
      sample_id = names(clust4),
      cluster = clust4,
      Timepoint = TP,
      Lineage = "CD4"
    )
    
    # CD8 clusters
    clust8 <- cutree(ph8s$tree_col, k = 4)
    df_clust8 <- data.frame(
      sample_id = names(clust8),
      cluster = clust8,
      Timepoint = TP,
      Lineage = "CD8"
    )
    
    # Store in list (outside loop define: clustering_results <- list())
    clustering_results_notjoined[[length(clustering_results_notjoined)+1]] <- df_clust4
    clustering_results_notjoined[[length(clustering_results_notjoined)+1]] <- df_clust8
    # ===========================================
    
    
    
    out_file_s <- paste0(
      res_path_notjoined,
      "metacluster_heatmap_CD4_CD8_scaled_notjoined_", TP, ".png"
    )
    
    png(out_file_s, width = 3600, height = 2400, res = 300)
    grid.arrange(ph4s$gtable, ph8s$gtable, ncol = 2)
    dev.off()
    
    message("Saved NOT-JOINED SCALED heatmap: ", out_file_s)
  }
  
  message("All NOT-JOINED SCALED heatmaps saved.")
  
  
  cluster_df <- bind_rows(clustering_results_notjoined)
  
  
  # =====================================================
  # CREATE CHORD PLOTS ACROSS TIMEPOINTS (NOT-JOINED)
  # =====================================================
  
  # Ensure correct ordering
  cluster_df$Timepoint <- factor(cluster_df$Timepoint, levels = sort(unique(cluster_df$Timepoint)))
  
  # Separate CD4 and CD8
  cd4_clusters <- cluster_df %>% filter(Lineage == "CD4")
  cd8_clusters <- cluster_df %>% filter(Lineage == "CD8")
  
  
  # ----------------------------
  # CD4 Chord plot
  # ----------------------------
  
  # Prepare chord table as you did
  cd4_clusters <- cd4_clusters %>%
    filter(Timepoint %in% c("T1", "T2")) %>%
    group_by(sample_id, Timepoint) %>%
    summarise(cluster = first(cluster), .groups = "drop")
  
  long_df <- cd4_clusters %>%
    pivot_wider(names_from = Timepoint, values_from = cluster) %>%
    drop_na() %>%
    mutate(across(all_of(c("T1", "T2")), as.character))
  
  chord <- long_df %>%
    group_by(across(all_of(c("T1", "T2")))) %>%
    summarise(value = n(), .groups = "drop") %>%
    rename(from = !!"T1", to = !!"T2") %>%
    mutate(
      from = paste0("C ", from, " T1"),
      to   = paste0("C ", to, " T2")
    )
  
  png(paste0(res_path_s, "Chord_CD4_notjoined_int_ticks.png"), width = 7, height = 7, units = "in", res = 600)
  circos.clear()
  circos.par(start.degree = 90)
  
  # Pre-allocate space for axis labels
  chordDiagram(chord, big.gap = 30, small.gap = 12,
               order = unique(c(sort(chord$to, decreasing = TRUE), sort(chord$from))))
  
  # Add integer axis for each sector
  sectors <- get.all.sector.index()
  for (s in sectors) {
    xlim <- get.cell.meta.data("xlim", sector.index = s)
    circos.axis(h = "top", sector.index = s, labels.cex = 0.6,
                major.at = seq(0, ceiling(xlim[2]), by = 1),
                labels.niceFacing = TRUE)
  }
  
  dev.off()
  circos.clear()
  
  
  # ----------------------------
  # CD8 Chord plot
  # ----------------------------
  
  # Prepare chord table as you did
  cd8_clusters <- cd8_clusters %>%
    filter(Timepoint %in% c("T1", "T2")) %>%
    group_by(sample_id, Timepoint) %>%
    summarise(cluster = first(cluster), .groups = "drop")
  
  long_df <- cd8_clusters %>%
    pivot_wider(names_from = Timepoint, values_from = cluster) %>%
    drop_na() %>%
    mutate(across(all_of(c("T1", "T2")), as.character))
  
  chord <- long_df %>%
    group_by(across(all_of(c("T1", "T2")))) %>%
    summarise(value = n(), .groups = "drop") %>%
    rename(from = !!"T1", to = !!"T2") %>%
    mutate(
      from = paste0("C ", from, " T1"),
      to   = paste0("C ", to, " T2")
    )
  
  png(paste0(res_path_s, "Chord_CD8_notjoined_int_ticks.png"), width = 7, height = 7, units = "in", res = 600)
  circos.clear()
  circos.par(start.degree = 90)
  
  # Pre-allocate space for axis labels
  chordDiagram(chord, big.gap = 30, small.gap = 12,
               order = unique(c(sort(chord$to, decreasing = TRUE), sort(chord$from))))
  
  # Add integer axis for each sector
  sectors <- get.all.sector.index()
  for (s in sectors) {
    xlim <- get.cell.meta.data("xlim", sector.index = s)
    circos.axis(h = "top", sector.index = s, labels.cex = 0.6,
                major.at = seq(0, ceiling(xlim[2]), by = 1),
                labels.niceFacing = TRUE)
  }
  
  dev.off()
  circos.clear()
  
  
  
  
  
  
  
  
  # =====================================================
  # CREATE CHORD PLOTS ACROSS TIMEPOINTS (SCALED)
  # =====================================================
  
  cluster_df_scaled <- bind_rows(clustering_results_scaled)
  
  # Ensure correct ordering
  cluster_df_scaled$Timepoint <- factor(cluster_df_scaled$Timepoint, levels = sort(unique(cluster_df_scaled$Timepoint)))
  
  # Separate CD4 and CD8
  cd4_clusters <- cluster_df_scaled %>% filter(Lineage == "CD4")
  cd8_clusters <- cluster_df_scaled %>% filter(Lineage == "CD8")
  
  
  # ----------------------------
  # CD4 Chord plot
  # ----------------------------
  
  # Ensure only one cluster per sample per timepoint
  cd4_clusters <- cd4_clusters %>%
    filter(Timepoint %in% c("T1", "T2")) %>%
    group_by(sample_id, Timepoint) %>%
    summarise(cluster = first(cluster), .groups = "drop")
  
  # Pivot wider
  long_df <- cd4_clusters %>%
    pivot_wider(names_from = Timepoint, values_from = cluster) %>%
    drop_na() %>%
    mutate(across(all_of(c("T1", "T2")), as.character))
  
  
  chord <- long_df %>%
    group_by(across(all_of(c("T1", "T2")))) %>%
    summarise(value = n(), .groups = "drop") %>%
    rename(from = !!"T1", to = !!"T2") %>%
    # Add cluster prefix and timepoint
    mutate(
      from = paste0("C ", from, " T1"),
      to   = paste0("C ", to, " T2")
    )
  
  
  
  
  png(paste0(res_path_s, "Chord_CD4_scaled.png"), width = 7, height = 7, units = "in", res = 600)
  circos.clear()
  circos.par(start.degree = 90)
  
  # Pre-allocate space for axis labels
  chordDiagram(chord, big.gap = 30, small.gap = 12,
               order = unique(c(sort(chord$to, decreasing = TRUE), sort(chord$from))))
  
  # Add integer axis for each sector
  sectors <- get.all.sector.index()
  for (s in sectors) {
    xlim <- get.cell.meta.data("xlim", sector.index = s)
    circos.axis(h = "top", sector.index = s, labels.cex = 0.6,
                major.at = seq(0, ceiling(xlim[2]), by = 1),
                labels.niceFacing = TRUE)
  }
  dev.off()
  circos.clear()
  
  
  # ----------------------------
  # CD8 Chord plot
  # ----------------------------
  
  # Ensure only one cluster per sample per timepoint
  cd8_clusters <- cd8_clusters %>%
    filter(Timepoint %in% c("T1", "T2")) %>%
    group_by(sample_id, Timepoint) %>%
    summarise(cluster = first(cluster), .groups = "drop")
  
  # Pivot wider
  long_df_cd8 <- cd8_clusters %>%
    pivot_wider(names_from = Timepoint, values_from = cluster) %>%
    drop_na() %>%
    mutate(across(all_of(c("T1", "T2")), as.character))
  
  # Create chord table
  chord_cd8 <- long_df_cd8 %>%
    group_by(across(all_of(c("T1", "T2")))) %>%
    summarise(value = n(), .groups = "drop") %>%
    rename(from = !!"T1", to = !!"T2") %>%
    mutate(
      from = paste0("C ", from, " T1"),
      to   = paste0("C ", to, " T2")
    )
  
  # Save chord plot
  png(paste0(res_path_s, "Chord_CD8_scaled.png"), width = 7, height = 7, units = "in", res = 600)
  circos.clear()
  circos.par(start.degree = 90)
  
  # Pre-allocate space for axis labels
  chordDiagram(chord_cd8, big.gap = 30, small.gap = 12,
               order = unique(c(sort(chord_cd8$to, decreasing = TRUE), sort(chord_cd8$from))))
  
  # Add integer axis for each sector
  sectors <- get.all.sector.index()
  for (s in sectors) {
    xlim <- get.cell.meta.data("xlim", sector.index = s)
    circos.axis(h = "top", sector.index = s, labels.cex = 0.6,
                major.at = seq(0, ceiling(xlim[2]), by = 1),
                labels.niceFacing = TRUE)
  }
  dev.off()
  circos.clear()




# ------------------------------------------------------------
# Can we reproduce that certain enterotypes are related to
# certain exhausted/senescent MC (new data only)
# ------------------------------------------------------------


