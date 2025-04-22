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

folders <- basename(list.dirs("models/", full.names = TRUE, recursive = FALSE))

l <- data.frame()
for (name in folders) {
  pathe <- paste0("models/", name, "/stats/")
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
        .groups = "keep"
      )
    l <- rbind(l, m)
  }
}

props <- colnames(l)[-c(1, 2)]

for (prop in props) {

  matrix <- dcast(l, ab ~ model, value.var = prop) %>%
    column_to_rownames("ab")
  
  pheatmap(
    as.matrix(matrix),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    fontsize_number = 6,
    number_color = "black",
    fontsize = 8,
    width = 2.5 + (ncol(matrix) * 0.5),
    height = (nrow(matrix) * 0.5),
    cellwidth = 20,
    cellheight = 20,
    angle_col = 45,
    main = paste0(toupper(prop), " Score"),
    color = colorRampPalette(c("red", "white", "blue"))(100),
    breaks = seq(0, 1, length.out = 101),
    filename = paste0("models/", prop, "_summary.png")
  )
}

