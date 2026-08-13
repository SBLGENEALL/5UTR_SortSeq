#!/usr/bin/env Rscript

packages <- c("ggplot2", "scales", "patchwork", "viridisLite", "ragg")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All required R packages are already installed.")
  quit(status = 0)
}

message("Installing missing R packages: ", paste(missing, collapse = ", "))
install.packages(missing, repos = "https://cloud.r-project.org")

still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("R package installation did not complete: ", paste(still_missing, collapse = ", "))
}
