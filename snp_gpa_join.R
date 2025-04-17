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
  yaml
)

config <- read_yaml("config.yaml")
gpam <- config$gpam_table

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
bin_data <- read_feather(file)

gpam_data <- read_tsv(gpam, col_names = TRUE, show_col_types = FALSE) %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("rownames")

joint <- gpam_data %>%
  right_join(bin_data, by = "rownames") %>%
  column_to_rownames("rownames")

if (!dir.exists("./joint_data")) {
  dir.create("./joint_data", recursive = TRUE)
}

filename <- gsub(".*/(.*)", "\\1", file)
write_feather(joint, paste0("./joint_data/", filename))
