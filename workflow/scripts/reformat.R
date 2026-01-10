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
  yaml,
  arrow
)

config <- read_yaml(file.path("config", "config.yaml"))

args <- commandArgs(trailingOnly = TRUE)
snvs <- args[1]
unitigs <- args[2]
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

variance_thresh <- 0.16

variances <- apply(sus_data, 2, function(x) var(x, na.rm = TRUE))
valid <- variances > variance_thresh

print("The following sample sets are not diverse enough to be significant:")
print(variances[variances <= variance_thresh])
sus_data <- sus_data[, valid]

print(colSums(!is.na(sus_data)))

snv_data <- read_tsv(snvs, col_names = TRUE, show_col_types = FALSE) %>%
  mutate(combined = paste(CHR, POS, sep = "-")) %>%
  column_to_rownames("combined") %>%
  select(-CHR, -POS, -REF) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("rownames")

unitig_data <- read_tsv(unitigs, col_names = TRUE, show_col_types = FALSE) %>%
  column_to_rownames("Unitig_sequence") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("rownames")

if (!dir.exists(file.path("results", "training_data"))) {
  dir.create(file.path("results", "training_data"), recursive = TRUE)
}

if (!dir.exists(file.path("results", "joint_data"))) {
  dir.create(file.path("results", "joint_data"), recursive = TRUE)
}

for (i in colnames(sus_data)) {
  n <- gsub("\\/", "_", i)
  fn <- file.path("results", "training_data", paste0(n, ".tsv"))
  fn_uni <- file.path("results", "training_data", paste0(n, "_unitig.feather"))
  is_susceptible <- sus_data %>%
    select(all_of(i)) %>%
    na.omit() %>%
    rename("Susceptible" = i) %>%
    rownames_to_column("rownames")
  merged_snv_data <- is_susceptible %>%
    left_join(snv_data, by = "rownames")
  write.table(
    merged_snv_data,
    file = fn,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  merged_unitig_data <- is_susceptible %>%
    left_join(unitig_data, by = "rownames") %>%
    select(-Susceptible)
  write_feather(merged_unitig_data, fn_uni)
}
