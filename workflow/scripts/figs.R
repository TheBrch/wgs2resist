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
  reshape2
)

pr <- list(
  name = 'pr',
  x = 'Recall',
  y = 'Precision'
)

roc <- list(
  name = 'roc',
  x = 'FPR',
  y = 'TPR'
)

args <- commandArgs(trailingOnly = TRUE)
name = args[1]
pathe <- file.path("results", "models", name, "stats")

for (graphtype in list(pr, roc)){
  l <- data.frame()
  met <- data.frame()
  for (model in c("logistic", "gaussian", "svm", "xgboost")){
    model_data <- read_tsv(
      file.path(pathe, paste0(model, "_", graphtype[["name"]],".tsv")),
      col_names = TRUE,
      show_col_types = FALSE
    )
      
    model_data$model <- model
    l <- rbind(l, model_data)
    
    metrics <- read_tsv(
      file.path(pathe, paste0(model, "_", "crossval_results.tsv")),
      col_names = TRUE,
      show_col_types = FALSE
    )
    metrics$model <- model
    met <- rbind(met, metrics)
  }
  l$fold <- factor(l$fold + 1)
  met$fold <- factor(met$Fold + 1)
  
  met_labels <- met %>%
    mutate(
      auc_score = get(paste0(toupper(graphtype[["name"]]), "_AUC")),
      label = paste0("AUC: ", sprintf("%.3f", auc_score))
    )
  
  plotvar <- ggplot(l, aes(x = .data[[tolower(graphtype[["x"]])]], y = .data[[tolower(graphtype[["y"]])]])) +
    ( if (graphtype[["name"]] == "roc") {
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                   color = "gray",
                   linetype = "dotted")
      } else {
        NULL
      }) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_line(aes(color = fold)) +
    geom_point(aes(color = fold)) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ poly(x, 3)) +
    geom_text(
      data = met_labels,
      aes(
        x = 0.95,
        y = 0.25 - (0.07 * (as.numeric(fold) - 1) ),
        label = label,
        color = fold
      ),
      hjust = 1,
      size = 3,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ model, ncol = 2) +
    labs(
      x = graphtype[["x"]],
      y = graphtype[["y"]],
      title = paste0(name, " ", toupper(graphtype[["name"]]))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(
    plot = plotvar,
    filename = file.path(pathe, "figs", paste0(name, "_", graphtype[["name"]], ".png")),
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white",
    create.dir = TRUE
  )
}

l <- data.frame()
for (model in c("logistic", "xgboost")){
  for (fold in 0:4) {
    
    fold_data <- read_tsv(
      file.path(pathe, paste0(model, "_f", fold, "_features.tsv")),
      col_names = TRUE,
      show_col_types = FALSE
    )
    
    fold_data$fold <- as.character(fold+1)
    fold_data$model <- model
    l <- rbind(l, fold_data)
  }
}

#filter out features that appear only once
features_in_multiple_folds <- l %>%
  filter(value != 0) %>%
  distinct(feature, model, fold) %>%
  count(feature, model) %>%
  filter(n > 1) %>%
  pull(feature)

# filter by selected feature names
l <- l %>%
  filter(feature %in% features_in_multiple_folds)

# sort the features by normalized values, limit the number to 10%
feature_order <- l %>%
  group_by(model) %>%
  mutate(norm_value = value / max(abs(value), na.rm = TRUE)) %>%
  group_by(feature, model) %>%
  summarise(max_abs = max(abs(norm_value)), .groups = "keep") %>%
  group_by(model) %>%
  slice_max(order_by = max_abs, prop = 0.1, with_ties = FALSE) %>%
  group_by(feature) %>%
  summarise(max_abs = max(max_abs), number=n(), .groups = "keep")

appears_in_both <- feature_order %>%
  filter(number > 1) %>%
  pull(feature)
  
#fill in the data from the main dataset, remove 0 value entries
l <- l %>%
  right_join(feature_order, by = "feature") %>%
  mutate(feature = fct_reorder(feature, max_abs, .desc = TRUE)) %>%
  filter(value != 0)

#plot out
plotvar <- ggplot(l, aes(x = feature, y = value, fill = fold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  labs(
    x = "Feature",
    y = "Value",
    title = paste0(name, " feature importance")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% appears_in_both, 
             paste0("<span style='color:red'>", x, "</span>"), 
             x)
    }
  ) +
  theme(axis.text.x = element_markdown())

#export to file
ggsave(
  plot = plotvar,
  filename = file.path(pathe, "figs", paste0(name, "_feature_importance.png")),
  width = 11,
  height = 6,
  dpi = 300,
  bg = "white",
  create.dir = TRUE
)

correlation <- read_feather(
  file.path("results", "condensed_data", paste0(name, "_hicorr.feather"))
) %>%
  column_to_rownames('index') %>%
  mutate(across(everything(), ~ round(as.numeric(.), 3)))

toptops <- correlation[appears_in_both,] %>%
  t() %>%
  as.data.frame()

filtered <- toptops[
  toptops %>%
    mutate(
      has_high = apply(., 1, function(row) any(abs(row) > 0.9, na.rm = TRUE))
      ) %>%
    pull(has_high),
] %>%
  t()

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
if (!is.null(filtered) && nrow(filtered) > 0 && ncol(filtered) > 0) {
  pheatmap(
    filtered,
    main = paste0(name, "\ntop feature correlations"),
    color = colorRampPalette(c("red", "white", "blue"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = nrow(filtered) > 2,
    cluster_cols = ncol(filtered) > 2,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 6,
    legend = TRUE,
    fontsize = 9,
    show_rownames = TRUE,
    show_colnames = TRUE,
    angle_col = 45,
    border_color = NA,
    width = 4 + (ncol(filtered) * 0.5),
    height = 2.5 + (nrow(filtered) * 0.5),
    cellwidth = 20,
    cellheight = 20,
    filename = file.path(pathe, "figs", paste0(name, "_corr.png"))
  )
}

