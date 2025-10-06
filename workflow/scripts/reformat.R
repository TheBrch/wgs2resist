if (!require(pacman)) {
  install.packages(
    "pacman",
    repos = "http://cran.us.r-project.org"
  )
}
p_load(
  readr,
  dplyr,
  tibble,
  yaml
)

config <- read_yaml(file.path("config", "config.yaml"))

file <- config$source_table
sus <- config$suscept_table

sus_data <- read_tsv(sus, col_names = TRUE, show_col_types = FALSE) %>%
  column_to_rownames("Antibiotic") %>%
  mutate(
    across(
      everything(),
      ~ ifelse(
        grepl("^(R|I|S|R\\(Inh\\))$", .),
        as.numeric(gsub("^(R|I|R\\(Inh\\))$", "0", gsub("^S$", "1", .))),
        .
      )
    )
  ) %>%
  t() %>%
  as.data.frame()

variances <- apply(sus_data, 2, function(x) var(x, na.rm = TRUE))
valid <- variances > 0.05

print("The following sample sets are not diverse enough to be significant:")
print(variances[variances <= 0.05])
sus_data <- sus_data[, valid]

print(colSums(!is.na(sus_data)))

tsv_data <- read_tsv(file, col_names = TRUE, show_col_types = FALSE) %>%
  mutate(combined = paste(CHR, POS, sep = "-")) %>%
  column_to_rownames("combined") %>%
  select(-CHR, -POS, -REF) %>%
  `colnames<-`(gsub("^PAN_", "", colnames(.))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("rownames")

if (!dir.exists("./training_data")) {
  dir.create("./training_data", recursive = TRUE)
}

for (i in colnames(sus_data)) {
  n <- gsub("\\/", "_", i)
  fn <- paste0("./training_data/", n, ".tsv")
  data <- sus_data %>%
    select(all_of(i)) %>%
    na.omit() %>%
    rename("Susceptible" = i) %>%
    rownames_to_column("rownames")
  merged_data <- data %>%
    left_join(tsv_data, by = "rownames")
  assign(n, merged_data)
  write.table(
    merged_data,
    file = fn,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}
