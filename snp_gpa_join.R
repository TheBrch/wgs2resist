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
gpam <- config$gpam_table