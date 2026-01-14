if (!require(pacman)) {
  install.packages(
    "pacman",
    repos = "http://cran.us.r-project.org"
  )
}
p_load(
  arrow,
  readr,
  dplyr,
  tibble,
  yaml,
  colorspace,
  ggplot2,
  pheatmap,
  proxy,
  RColorBrewer,
  tidyr,
  forcats,
  ggtext,
  reshape2,
  stringr,
  jsonlite
)

config <- read_yaml(file.path("config", "config.yaml"))
models <- unlist(strsplit(config$models, " "))
models <- c(models, "wec")

pr <- list(
  name = "pr",
  x = "Recall",
  y = "Precision"
)

roc <- list(
  name = "roc",
  x = "FPR",
  y = "TPR"
)

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
pan_json <- args[2]
pan_tsv <- args[3]
pan_csv <- args[4]
snv_effects <- args[5]
pathe <- file.path("results", "models", name, "stats")

met <- data.frame() # Optimal F1 metrics for all models, all graphs

for (model in models) {
  metrics <- read_tsv(
    file.path(pathe, paste0(model, "_", "crossval_results.tsv")),
    col_names = TRUE,
    show_col_types = FALSE
  )
  metrics$model <- model
  met <- rbind(met, metrics)
}
met$fold <- factor(met$Fold + 1)

model_labels <- setNames(models, seq_along(models))

for (graphtype in list(pr, roc)) {
  # Accumulate point coordinates for graphtype from separate files
  l <- data.frame()
  for (model in models) {
    model_data <- read_tsv(
      file.path(pathe, paste0(model, "_", graphtype[["name"]], ".tsv")),
      col_names = TRUE,
      show_col_types = FALSE
    )

    model_data$model <- model
    l <- rbind(l, model_data)
  }
  l$fold <- factor(l$fold + 1)
  l$model <- match(l$model, models)

  met_labels <- met %>%
    mutate(
      auc_score = get(paste0(toupper(graphtype[["name"]]), "_AUC")),
      tpr = TP / (TP + FN),
      fpr = FP / (FP + TN),
      label = paste0("AUC: ", sprintf("%.3f", auc_score))
    ) %>%
    select(-auc_score)
  met_labels$model <- match(met_labels$model, models)

  plotvar <- ggplot(
    l, aes(
      x = .data[[tolower(graphtype[["x"]])]],
      y = .data[[tolower(graphtype[["y"]])]]
    )
  ) +
    (if (graphtype[["name"]] == "roc") {
      geom_segment(
        aes(x = 0, y = 0, xend = 1, yend = 1, linetype = "Baseline"),
        color = "gray"
      )
    } else {
      NULL
    }) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_line(aes(color = fold), linewidth = 1.2, alpha = 0.5) +
    geom_point(
      data = met_labels %>%
        select(fold, model, fpr, tpr, `Precision@Thresh`, `Recall@Thresh`) %>%
        mutate(
          !!tolower(graphtype[["x"]]) :=
            if (graphtype[["name"]] == "roc") fpr else `Recall@Thresh`,
          !!tolower(graphtype[["y"]]) :=
            if (graphtype[["name"]] == "roc") tpr else `Precision@Thresh`
        ),
      aes(
        color = fold
      ), size = 3,
      shape = 21,
      fill = "white",
      stroke = 1.5,
    ) +
    geom_smooth(
      aes(linetype = "Average (polynomial)"),
      method = "lm",
      se = FALSE,
      color = "blue",
      formula = y ~ poly(x, 4)
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Baseline" = "dotted", "Average (polynomial)" = "solid"),
    ) +
    geom_text(
      data = met_labels,
      aes(
        x = 1,
        y = 0.35 - (0.07 * (as.numeric(fold) - 1)),
        label = label,
        color = fold
      ),
      hjust = 1,
      size = 2.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    facet_wrap(~model, ncol = 2, labeller = labeller(model = model_labels)) +
    labs(
      x = graphtype[["x"]],
      y = graphtype[["y"]],
      title = paste0(name, " ", toupper(graphtype[["name"]]))
    ) +
    theme_minimal() +
    theme(
      legend.box = "vertical",
      legend.position = c(0.85, 0.05),
      legend.justification = c(1, 0),
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        title = "Fold",
        order = 1
      ),
      linetype = guide_legend(
        nrow = 1,
        title = "Reference",
        order = 2
      )
    )

  ggsave(
    plot = plotvar,
    filename = file.path(
      pathe, "figs", paste0(name, "_", graphtype[["name"]], ".png")
    ),
    width = 10.4,
    height = 8,
    dpi = 300,
    bg = "white",
    create.dir = TRUE
  )
}

l <- data.frame()

featuremodels <- intersect(models, c("logistic", "xgboost"))

for (model in featuremodels) {
  for (fold in 0:4) {
    fold_data <- read_tsv(
      file.path(pathe, paste0(model, "_f", fold, "_features.tsv")),
      col_names = TRUE,
      show_col_types = FALSE
    )

    fold_data$fold <- as.character(fold + 1)
    fold_data$model <- model
    l <- rbind(l, fold_data)
  }
}


l_newscores <- l %>%
  filter(abs(value) > 0) %>%
  group_by(model, feature) %>%
  summarise(
    avg = mean(value),
    n_folds = n(),
    score = avg * n_folds
  ) %>%
  ungroup()


# Bakta file: contig -> panaroo id
pan_seq <- fromJSON(pan_json)$sequences %>%
  select(simple_id, orig_id)

# Bakta file: contig -> annotation
pan_tsv <- read_tsv(
  pan_tsv,
  show_col_types = FALSE, comment = "# "
)

# Old annotations
# Panaroo file:
old_gpa_annot <- read_csv(pan_csv, show_col_types = FALSE) %>%
  select(Gene, `Non-unique Gene name`, Annotation)


# SNV effects
snv_eff <- read_tsv(snv_effects, show_col_types = FALSE) %>%
  filter(!is.na(LOCUS_TAG))

annotate_features <- function(
  l_newscores, pan_seq, pan_tsv, old_gpa_annot, snv_eff
) {
  # Joining: panaroo id -> annotation
  pan_tsv_original <- pan_tsv %>%
    left_join(pan_seq, join_by(`#Sequence Id` == simple_id))

  # Annotating part of gpa
  l_gpa_annot_inter <- l_newscores %>%
    left_join(pan_tsv_original, join_by(feature == orig_id)) %>%
    mutate(
      Gene = if_else(
        (!is.na(Type) & Type != "cds"),
        if_else(
          !is.na(`Locus Tag`),
          `Locus Tag`,
          feature
        ),
        Gene
      )
    )

  # l_gpa_annot_inter[
  #   duplicated(l_gpa_annot_inter[c("model", "feature")]) |
  #     duplicated(l_gpa_annot_inter[c("model", "feature")], fromLast = TRUE),
  # ]

  gpa_new_annot <- l_gpa_annot_inter %>%
    filter(!is.na(`#Sequence Id`)) %>%
    select(-DbXrefs, -`#Sequence Id`, -Type, -Start, -Stop, -Strand)


  l_gpa_annot_inter2 <- l_gpa_annot_inter %>%
    filter(is.na(`#Sequence Id`)) %>%
    left_join(old_gpa_annot, join_by(feature == Gene))


  gpa_old_annot <- l_gpa_annot_inter2 %>%
    filter(!is.na(Annotation)) %>%
    select(where(~ !all(is.na(.)))) %>%
    rename(
      Gene = `Non-unique Gene name`,
      Product = Annotation
    ) %>%
    mutate(
      Gene = gsub("^;+|;+$", "", Gene)
    )

  annotated_gpa <- gpa_new_annot %>%
    bind_rows(gpa_old_annot) %>%
    mutate(
      Effect = "presence",
      Final_name = if_else(
        !is.na(Gene),
        Gene,
        if_else(
          !grepl("group_", feature),
          feature,
          if_else(
            !is.na(`Locus Tag`),
            `Locus Tag`,
            feature
          )
        )
      )
    ) %>%
    select(-`Locus Tag`, -Gene)


  l_snv <- l_gpa_annot_inter2 %>% filter(is.na(Annotation))


  l_snv_split <- l_snv %>%
    select(where(~ !all(is.na(.)))) %>%
    mutate(
      snv_locus = str_split_i(feature, "-", 1),
      snv_pos_alt = str_split_i(feature, "-", 2),
      snv_pos = as.double(str_split_i(snv_pos_alt, "_", 1)),
      snv_alt = str_split_i(snv_pos_alt, "_", 2)
    ) %>%
    select(-snv_pos_alt)

  l_snv_annot <- l_snv_split %>%
    left_join(
      snv_eff, join_by(
        snv_locus == LOCUS_TAG,
        snv_pos == POS,
        snv_alt == ALT
      )
    ) %>%
    rename(
      Product = PRODUCT,
      Gene = GENE,
      Effect = EFFECT
    ) %>%
    select(-REF) %>%
    mutate(
      Effect = str_split_i(str_split_i(Effect, "&", 1), " ", 1)
    )

  l_snv_v_annot <- l_snv_annot %>% filter(!is.na(Effect))

  l_snv_wt_annot <- l_snv_annot %>%
    filter(is.na(Effect)) %>%
    select(where(~ !all(is.na(.)))) %>%
    left_join(
      snv_eff, join_by(
        snv_locus == LOCUS_TAG,
        snv_pos == POS,
        snv_alt == REF
      )
    ) %>%
    select(-ALT, -EFFECT) %>%
    rename(
      Product = PRODUCT,
      Gene = GENE
    ) %>%
    mutate(
      Effect = "wild-variant"
    ) %>%
    distinct(.keep_all = TRUE)

  annotated_snvs <- rbind(l_snv_wt_annot, l_snv_v_annot) %>%
    mutate(
      Final_name = if_else(
        !is.na(Gene),
        Gene,
        snv_locus,
      )
    ) %>%
    select(
      -snv_locus, -snv_pos, -snv_alt, -Gene
    )


  annotated_scored_features <- rbind(annotated_gpa, annotated_snvs)

  return(annotated_scored_features)
}

annotated_scored_features <-
  annotate_features(l_newscores, pan_seq, pan_tsv, old_gpa_annot, snv_eff) %>%
  mutate(abs_score = abs(score))

annotated_ranked_features <- annotated_scored_features %>%
  # filter(!grepl("hypothetical", Final_name)) %>%
  group_by(model) %>%
  mutate(rating = rank(-abs_score, ties.method = "first")) %>%
  ungroup()

plot_data <- annotated_ranked_features %>%
  filter(rating <= 50) %>%
  group_by(Final_name) %>%
  mutate(InBoth = as.integer(n_distinct(model) > 1)) %>%
  ungroup() %>%
  group_by(model) %>%
  arrange(rating) %>%
  mutate(
    Final_name_ordered = paste0(Final_name, "___", model),
    Final_name_ordered = factor(
      Final_name_ordered,
      levels = unique(Final_name_ordered)
    ),
    # Create the border variables in the data
    bar_color = ifelse(InBoth == 1, "black", NA),
    bar_linewidth = ifelse(InBoth == 1, 1, 0)
  ) %>%
  ungroup()

# Create plot
plotvar <- ggplot(
  plot_data,
  aes(x = Final_name_ordered, y = score, fill = Effect)
) +
  geom_bar(
    stat = "identity",
    aes(color = bar_color),
    linewidth = 1
  ) +
  scale_color_identity() +
  facet_wrap(~model, ncol = 2, scales = "free") +
  scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
  labs(
    x = "Feature",
    y = "Value",
    title = paste0(name, " feature importance")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  ) +
  coord_flip()
# export to file
ggsave(
  plot = plotvar,
  filename = file.path(pathe, "figs", paste0(name, "_feature_importance.png")),
  width = 9,
  height = 11,
  dpi = 300,
  bg = "white",
  create.dir = TRUE
)

get_feature_correlations <- function(
  feature_name, corr_df, include_transitive = TRUE
) {
  all_features <- character(0)
  visited <- character(0)

  find_correlations <- function(current_feature) {
    if (current_feature %in% visited) {
      return()
    }
    visited <<- c(visited, current_feature)
    all_features <<- c(all_features, current_feature)

    row_data <- corr_df %>% filter(index == current_feature)
    if (nrow(row_data) > 0) {
      correlations <- row_data %>% select(-index)
      for (col_name in names(correlations)) {
        if (!is.na(correlations[[col_name]])) {
          all_features <<- c(all_features, col_name)
          if (include_transitive) {
            find_correlations(col_name)
          }
        }
      }
    }

    # Find reverse correlations (where current_feature is a column)
    if (current_feature %in% names(corr_df)) {
      col_data <- corr_df %>%
        filter(!is.na(.data[[current_feature]])) %>%
        pull(index)

      all_features <<- c(all_features, col_data)
      if (include_transitive) {
        for (feat in col_data) {
          find_correlations(feat)
        }
      }
    }
  }

  find_correlations(feature_name)
  all_features <- unique(all_features)

  cols_to_keep <- c("index", intersect(names(corr_df), all_features))

  result <- corr_df %>%
    filter(index %in% all_features) %>%
    select(all_of(cols_to_keep))

  return(result)
}

get_and_melt <- function(feature_name, corr_df, include_transitive = FALSE) {
  subset_df <- get_feature_correlations(
    feature_name, corr_df, include_transitive
  )

  if (ncol(subset_df) <= 1) {
    return(tibble(
      input_feature = character(), feature_i = character(),
      feature_j = character(), correlation = numeric()
    ))
  }

  # Melt to long format
  melted <- subset_df %>%
    pivot_longer(
      cols = -index,
      names_to = "feature_j",
      values_to = "correlation"
    ) %>%
    rename(feature_i = index) %>%
    filter(!is.na(correlation)) %>%
    mutate(input_feature = feature_name, .before = 1)
  # Add input feature as first column

  return(melted)
}

get_multiple_features_melted <- function(
  feature_list, corr_df, include_transitive = FALSE
) {
  feature_list %>%
    lapply(function(feat) {
      tryCatch(
        get_and_melt(feat, corr_df, include_transitive),
        error = function(e) {
          warning(paste("Error processing feature:", feat, "-", e$message))
          tibble(
            input_feature = character(), feature_i = character(),
            feature_j = character(), correlation = numeric()
          )
        }
      )
    }) %>%
    bind_rows()
}


# correlation_df <- read_feather(
#   file.path("results", "condensed_data", paste0(name, "_hicorr.feather"))
# )

# melted <- get_multiple_features_melted(
#   unique(annotated_scored_features$feature), correlation_df
# )

# standardized_corr <- melted %>%
#   mutate(
#     feature = if_else(
#       feature_i == input_feature,
#       feature_j,
#       feature_i
#     )
#   ) %>%
#   select(
#     -feature_i, -feature_j
#   )

# scored_corr <- standardized_corr %>%
#   left_join(
#     annotated_scored_features %>%
#       select(model, feature, n_folds, score),
#     join_by(input_feature == feature)
#   ) %>%
#   mutate(
#     score = if_else(
#       model == "logistic",
#       score * correlation,
#       score * abs(correlation)
#     )
#   )

# annot_corr <- annotate_features(
#   scored_corr, pan_seq, pan_tsv, old_gpa_annot, snv_eff
# )
# total_features <- annotated_scored_features %>% bind_rows(annot_corr)

##########################################

# toptops <- correlation[appears_in_both, ] %>%
#   t() %>%
#   as.data.frame()

# filtered <- toptops[
#   toptops %>%
#     mutate(
#       has_high = apply(., 1, function(row) any(abs(row) > 0.9, na.rm = TRUE))
#     ) %>%
#     pull(has_high),
# ] %>%
#   t()

# melted_filtered <- melt(filtered)

# heatmap_plot <- ggplot(melted_filtered, aes(x=Var1, y=Var2, fill=value)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "red",
#     mid = "white",
#     high = "blue",
#     midpoint = 0,
#     limits = c(-1, 1)
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     axis.title = element_blank()
#   ) +
#   labs(fill = "Value")

# ggsave(
#   plot = heatmap_plot,
#   filename = paste0(pathe, "/figs/", name, "_corr.png"),
#   width = 15,
#   height = 15,
#   dpi = 300,
#   bg = "white",
#   create.dir = TRUE
# )
# if (!is.null(filtered) && nrow(filtered) > 0 && ncol(filtered) > 0) {
#   pheatmap(
#     filtered,
#     main = paste0(name, "\ntop feature correlations"),
#     color = colorRampPalette(c("red", "white", "blue"))(100),
#     breaks = seq(-1, 1, length.out = 101),
#     cluster_rows = nrow(filtered) > 2,
#     cluster_cols = ncol(filtered) > 2,
#     display_numbers = TRUE,
#     number_format = "%.2f",
#     fontsize_number = 6,
#     legend = TRUE,
#     fontsize = 9,
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     angle_col = 45,
#     border_color = NA,
#     width = 4 + (ncol(filtered) * 0.5),
#     height = 2.5 + (nrow(filtered) * 0.5),
#     cellwidth = 20,
#     cellheight = 20,
#     filename = file.path(pathe, "figs", paste0(name, "_corr.png"))
#   )
# }
