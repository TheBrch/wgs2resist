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
  grid,
  gridExtra
)

folders <- basename(list.dirs("results/models/", full.names = TRUE, recursive = FALSE))

l <- data.frame()
for (name in folders) {
  pathe <- paste0("results/models/", name, "/stats/")
  for (model in c("logistic", "gaussian", "svm", "xgboost")){
    metrics <- read_tsv(
      paste0(pathe, model, "_", "crossval_results.tsv"),
      col_names = TRUE,
      show_col_types = FALSE
    )
    metrics$model <- model
    metrics$ab <- name
    
    m <- metrics %>%
      group_by(model, ab) %>%
      summarise(
        accuracy_at_thresh = mean(Score, na.rm = TRUE),
        f1_at_thresh = mean(F1, na.rm = TRUE),
        roc_auc = mean(ROC_AUC, na.rm = TRUE),
        pr_auc = mean(PR_AUC, na.rm = TRUE),
        precision_at_thresh = mean(`Precision@Thresh`, na.rm = TRUE),
        recall_at_thresh = mean(`Recall@Thresh`, na.rm = TRUE),
        best_thresh = mean(`Thresh`, na.rm = TRUE),
        .groups = "keep"
      )
    l <- rbind(l, m)
  }
}

average_metrics <- l %>%
  group_by(model) %>%
  summarise(
    accuracy_at_thresh = mean(accuracy_at_thresh, na.rm = TRUE),
    f1_at_thresh = mean(f1_at_thresh, na.rm = TRUE),
    roc_auc = mean(roc_auc, na.rm = TRUE),
    pr_auc = mean(pr_auc, na.rm = TRUE),
    precision_at_thresh = mean(precision_at_thresh, na.rm = TRUE),
    recall_at_thresh = mean(recall_at_thresh, na.rm = TRUE),
    best_thresh = mean(best_thresh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ab = "Average")

l <- rbind(l, average_metrics)

props <- colnames(l)[-c(1, 2)]


heatmap_plots <- list()

for (prop in props) {

  matrix <- dcast(l, ab ~ model, value.var = prop) %>%
    column_to_rownames("ab")
  
  row_order <- rownames(matrix)
  if ("Average" %in% row_order) {
    other_rows <- row_order[row_order != "Average"]
    row_order <- c(other_rows, "Average")
    matrix <- matrix[row_order, ]
  }
  
  matrix_width <- 2.3 + (ncol(matrix) * 0.5)
  matrix_height <- (nrow(matrix) * 0.4)
  
  plot <- pheatmap(
    as.matrix(matrix),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.3f",
    fontsize_number = 6,
    number_color = "black",
    fontsize = 8,
    width = matrix_width,
    height = matrix_height,
    cellwidth = 20,
    cellheight = 20,
    angle_col = 45,
    main = paste0(toupper(prop), " Score"),
    color = colorRampPalette(c("red", "white", "blue"))(100),
    breaks = seq(min(matrix, na.rm = TRUE), 1, length.out = 101),
    filename = paste0("results/models/", prop, "_summary.png")
  )
  
  if (prop == "roc_auc") {
    heatmap_plots[[1]] <- plot$gtable
  } else if (prop == "pr_auc") {
    heatmap_plots[[2]] <- plot$gtable
  } else if (prop == "f1_at_thresh") {
    heatmap_plots[[3]] <- plot$gtable
  } else if (prop == "accuracy_at_thresh") {
    heatmap_plots[[4]] <- plot$gtable
  }
  
}

png("results/models/combined_metrics_summary.png", width = 2 * matrix_width, height = 2 * matrix_height, units = "in", res = 300)

combined_plot <- grid.arrange(
  heatmap_plots[[1]], heatmap_plots[[2]],
  heatmap_plots[[3]], heatmap_plots[[4]],
  ncol = 2, nrow = 2
)

grid.text("A", x = 0.01, y = 0.99, just = c("left", "top"), 
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("B", x = 0.51, y = 0.99, just = c("left", "top"), 
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("C", x = 0.01, y = 0.49, just = c("left", "top"), 
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("D", x = 0.51, y = 0.49, just = c("left", "top"), 
          gp = gpar(fontsize = 14, fontface = "bold"))

dev.off()

