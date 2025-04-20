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

args <- commandArgs(trailingOnly = TRUE)
name = args[1]
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
  ungroup() %>%
  group_by(feature, model) %>%
  summarise(max_abs = max(abs(norm_value)), .groups = "keep") %>%
  ungroup() %>%
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
  filename = paste0(pathe, "/figs/", name, "_feature_importance.png"),
  width = 11,
  height = 6,
  dpi = 300,
  bg = "white",
  create.dir = TRUE
)




correlation <- read_feather(
  paste0("condensed_data/", name, "_hicorr.feather")
) %>%
  column_to_rownames('index') %>%
  mutate(across(everything(), ~ round(as.numeric(.), 3)))


filtered <- correlation[
  correlation %>%
    mutate(
      has_high = apply(., 1, function(row) any(row > 0.9, na.rm = TRUE))
      ) %>%
    pull(has_high),
] %>%
  # rownames_to_column("rownames") %>%
  as.matrix()

melted_filtered <- melt(filtered)

# pheatmap(
#   filtered,
#   # filename = "dang.png",
#   # fontsize = 7,
#   # cluster_rows = TRUE,
#   # cluster_cols = TRUE,
#   show_colnames = TRUE,
#   show_rownames = TRUE,
#   color = colorRampPalette(c("white","blue"))(100)
# )



heatmap_plot <- ggplot(melted_filtered, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradientn(colors=colorRampPalette(c("white","blue"))(100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate column labels for better readability
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),  # Remove border
    axis.title = element_blank()    # Remove axis titles
  ) +
  labs(fill = "Value")

ggsave(
  plot = heatmap_plot,
  filename = paste0(pathe, "/figs/", name, "_corr.png"),
  width = 15,
  height = 15,
  dpi = 300,
  bg = "white",
  create.dir = TRUE
)
