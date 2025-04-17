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

name = "Sulfonamide"
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
      
      fold_data$fold <- as.character(fold)
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
    
    fold_data$fold <- as.character(fold)
    fold_data$model <- model
    l <- rbind(l, fold_data)
  }
}
  
l <- l %>%
  filter(value != 0) %>%
  group_by(model) %>%
  mutate(feature = fct_reorder(feature, abs(value), .desc = TRUE)) %>%
  slice_max(abs(value), n = 10) %>%
  ungroup()

plotvar <- ggplot(l, aes(x = fct_reorder(feature, abs(value), .desc = TRUE), y = value, fill = fold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ model, ncol = 2) +
  labs(
    x = "Feature",
    y = "Value",
    title = paste0(name, " feature importance")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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



