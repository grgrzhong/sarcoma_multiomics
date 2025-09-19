## https://github.com/lima1/PureCN

# vignette("Quick", package = "PureCN")
# vignette("PureCN", package = "PureCN")
library(PureCN)
library(fs)

purecn_dir <- "data/WES/PureCN"

normal_coverage_files <- dir_ls(
    purecn_dir,
    glob = "*-N_recalibrated_coverage_loess.txt.gz",
    recurse = TRUE
)
normal_coverage_files <- as.vector(normal_coverage_files)

purecn_path <- system.file("extdata", package = "PureCN")
normal.coverage.file <- system.file("extdata", "example_normal.txt.gz",
package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt.gz",
package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)

normalDB <- createNormalDatabase(normal_coverage_files)
