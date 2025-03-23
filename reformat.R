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

config <- read_yaml("config.yaml")

file <- config$source_table
sus <- config$suscept_table

sus_data <- read_tsv(sus, col_names = TRUE, show_col_types = FALSE) %>%
  column_to_rownames("Antibiotic") %>%
  select(contains("_PCD_")) %>%
  mutate(
    across(
      everything(),
      ~ ifelse(
        grepl("^(R|I|S|R\\(Inh\\))$", .),  # Check if value contains R, I, or S
        as.numeric(gsub('^(R|I|R\\(Inh\\))$', '0', gsub('^S$', '1', .))),
        .
      )
    )
  ) %>%
  #filter(complete.cases(.)) %>%
  filter(!apply(., 1, function(x) all(x == x[1]) | all(is.na(x)))) %>%
  t()

tsv_data <- read_tsv(file, col_names = TRUE, show_col_types = FALSE) %>%
  mutate(combined = paste(CHR, POS, sep = "-")) %>%
  column_to_rownames("combined") %>%
  select(-CHR, -POS, -REF) %>%
  `colnames<-`(gsub("^PAN_", "", colnames(.))) %>%
  t()

merged_data <- merge(tsv_data, sus_data, by = "row.names") %>%
  column_to_rownames("Row.names")

write.table(
  merged_data,
  file = "training_matrix.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

