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
  gridExtra,
  ComplexUpset
)

folders <- basename(
  list.dirs(
    file.path("results", "models"),
    full.names = TRUE, recursive = FALSE
  )
)

config <- read_yaml(file.path("config", "config.yaml"))
models <- unlist(strsplit(config$models, " "))

l <- data.frame()
correctness <- data.frame()
for (name in folders) {
  pathe <- file.path("results", "models", name, "stats")
  for (model in models) {
    metrics <- read_tsv(
      file.path(pathe, paste0(model, "_", "crossval_results.tsv")),
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
  ab_correctness <- read_tsv(
    file.path(pathe, "correctness.tsv")
  ) %>% select(-sample_id)
  correctness <- rbind(correctness, ab_correctness)
  up_set_plots <- list()
  up_set <- upset(
    df,
    intersect = sort(colnames(df)),
    sort_sets = FALSE,
    name = "ML models",
    base_annotations = list(
      "Correct predictions" = (
        intersection_size(
          counts = TRUE
        ) +
          ggtitle(name)
      )
    ),
    set_sizes = (
      upset_set_size(
        mapping = aes(fill = group),
        geom = geom_bar(width = 0.6),
      ) +
        geom_text(
          aes(label = after_stat(count)),
          hjust = -0.2, stat = "count"
        ) +
        scale_fill_hue() +
        guides(fill = "none") +
        ylab("Correct predictions")
    ),
    themes = upset_modify_themes(
      list(
        "overall" = theme_minimal()
      )
    )
  )
  c(my_list, list(up_set))
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
    color = if (prop == "best_thresh") {
      colorRampPalette(c("blue", "white", "red3"))(100)
    } else {
      colorRampPalette(c("red3", "gold", "forestgreen"))(100)
    },
    breaks = seq(0, 1, length.out = 101),
    filename = file.path("results", "models", paste0(prop, "_summary.png"))
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

png(
  file.path("results", "models", "combined_metrics_summary.png"),
  width = 2 * matrix_width, height = 2 * matrix_height, units = "in", res = 300
)

n_cols <- 2
n_rows <- ceiling(length(heatmap_plots) / n_cols)

combined_plot <- do.call(
  grid.arrange,
  c(heatmap_plots, ncol = n_cols, nrow = n_rows)
)

labels <- LETTERS[seq_along(heatmap_plots)]

for (i in seq_along(labels)) {
  row <- ceiling(i / n_cols)
  col <- (i - 1) %% n_cols + 1

  x <- (col - 1) * (1 / n_cols) + 0.01
  y <- 1 - (row - 1) * (1 / n_rows) - 0.01

  grid.text(labels[i],
    x = x, y = y, just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
}

dev.off()

png(
  file.path("results", "models", "combined_up_set.png"),
  width = 3 * matrix_width, height = 4 * matrix_height, units = "in", res = 300
)

n_rows <- ceiling(length(up_set_plots) / n_cols)

combined_plot <- do.call(
  grid.arrange,
  c(up_set_plots, ncol = n_cols, nrow = n_rows)
)

labels <- LETTERS[seq_along(up_set_plots)]

for (i in seq_along(labels)) {
  row <- ceiling(i / n_cols)
  col <- (i - 1) %% n_cols + 1

  x <- (col - 1) * (1 / n_cols) + 0.01
  y <- 1 - (row - 1) * (1 / n_rows) - 0.01

  grid.text(labels[i],
    x = x, y = y, just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
}



dev.off()
