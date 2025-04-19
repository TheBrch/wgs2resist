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
  forcats
)

prc <- list(
  name = 'prc',
  x = 'Recall',
  y = 'Precision'
)

roc <- list(
  name = 'roc',
  x = 'FPR',
  y = 'TPR'
)

# args <- commandArgs(trailingOnly = TRUE)
# name = args[1]
name = "Trimethoprim_sulfamethoxazole"
pathe <- paste0("models/", name, "/stats/")

for (graphtype in list(prc, roc)){
  l <- data.frame()
  for (model in c("logistic", "gaussian", "svm", "xgboost")){
    for (fold in 0:4) {
      
      fold_data <- read_tsv(
        paste0(pathe, model, "_f", fold, "_", graphtype[["name"]],".tsv"),
        col_names = TRUE,
        show_col_types = FALSE
      )
      
      fold_data$fold <- as.character(fold+1)
      fold_data$model <- model
      l <- rbind(l, fold_data)
    }
  }
  
  plotvar <- ggplot(l, aes(x = .data[[tolower(graphtype[["x"]])]], y = .data[[tolower(graphtype[["y"]])]])) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                 color = "gray", 
                 linetype = "dotted") +
    geom_line(aes(color = fold)) +
    geom_point(aes(color = fold)) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ poly(x, 3)) +
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
    filename = paste0(pathe, "/figs/", name, "_", graphtype[["name"]], ".png"),
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
      paste0(pathe, model, "_f", fold, "_features.tsv"),
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
  distinct(feature, fold) %>%
  count(feature) %>%
  filter(n > 1) %>%
  pull(feature)


# select all feature names that appear more than once throughout all folds, for each model
selected <- l %>%
  filter(feature %in% features_in_multiple_folds) %>%
  group_by(model, feature) %>%
  summarise(max_value = max(value, na.rm = TRUE), .groups = "drop") %>%
  # group_by(model) %>%
  # slice_max(order_by = max_value, n = 10, with_ties = FALSE) %>%
  # ungroup() %>%
  pull(feature)

# filter by selected feature names
l <- l %>%
  filter(feature %in% selected)

# sort the features by normalized values, limit the number
feature_order <- l %>%
  group_by(model) %>%
  mutate(norm_value = value / max(abs(value), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(feature) %>%
  summarise(max_abs = max(abs(norm_value)), .groups = "drop") %>%
  slice_max(order_by = max_abs, n = 30, with_ties = FALSE)

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
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#export to file
ggsave(
  plot = plotvar,
  filename = paste0(pathe, "/figs/", name, "_feature_importance.png"),
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white",
  create.dir = TRUE
)

# correlation <- read_feather(
#   "condensed_data/Sulfonamide_hicorr.feather"
# ) %>%
#   column_to_rownames('index') %>%
#   mutate(across(everything(), ~ round(as.numeric(.), 3)))
# 
# 
# filtered <- correlation[
#   correlation %>%
#     mutate(
#       has_high = apply(., 1, function(row) any(row > 0.9, na.rm = TRUE))
#       ) %>%
#     pull(has_high),
# ] %>%
#   as.matrix()
# 
# pheatmap(
#   filtered,
#   filename = "dang.png",
#   fontsize = 7,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   show_colnames = TRUE,
#   show_rownames = TRUE,
#   color = colorRampPalette(c("white","blue"))(100)
# )



