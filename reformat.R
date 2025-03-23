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
  filter(!apply(., 1, function(x) all(x == x[1]) | all(is.na(x)))) %>%
  t() %>%
  as.data.frame()

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

colSums(!is.na(sus_data))

for (i in colnames(sus_data)) {
  n <- gsub("\\/", "_", i)
  fn <- paste0("./training_data/", n, ".tsv")
  data <- sus_data %>%
    select(all_of(i)) %>%
    na.omit() %>%
    rownames_to_column("rownames")
  merged_data <- data %>%
    left_join(tsv_data, by = "rownames") %>%
    column_to_rownames("rownames")
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



