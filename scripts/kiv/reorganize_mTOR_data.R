library(tidyverse)
library(fs)

raw_dir <- "/mnt/s/Dataset/Sanger Sequencing/SFT mTOR/SFT mTOR Sanger Data/Original data"
data_dir <- "/mnt/s/Dataset/Sanger Sequencing/SFT mTOR/SFT mTOR Sanger Data"

# Get all files in the raw directory (excluding .zip files)
all_files <- dir_ls(
    raw_dir,
    recurse = TRUE,
    type = "file"
) |>
    str_subset("\\.zip$", negate = TRUE)

# Extract SFT case numbers from filenames and format with hyphen
extract_sft_case <- function(filepath) {
    filename <- basename(filepath)
    sft_match <- str_extract(filename, "SFT-\\d{3}")
    return(sft_match)
}

# Create a data frame with file paths and their corresponding SFT cases
file_df <- tibble(
    filepath = all_files,
    filename = basename(all_files),
    sft_case = map_chr(all_files, extract_sft_case)
) |>
    filter(!is.na(sft_case)) |>
    arrange(sft_case)

# Get unique SFT cases
unique_cases <- file_df |>
    pull(sft_case) |>
    unique() |>
    sort()

cat("Found", length(unique_cases), "unique SFT cases:\n")
cat(paste(unique_cases, collapse = ", "), "\n\n")

# Create directories for each case in the data directory
dir_create(data_dir)

for (case in unique_cases) {
    case_dir <- file.path(data_dir, case)
    dir_create(case_dir)
    cat("Created directory:", case_dir, "\n")
}

# Move files to their respective case directories
cat("\nCopying files to case directories...\n")
for (i in 1:nrow(file_df)) {
    source_file <- file_df$filepath[i]
    case <- file_df$sft_case[i]
    filename <- file_df$filename[i]

    destination_file <- file.path(data_dir, case, filename)

    # Try to copy file to the case directory with error handling
    tryCatch(
        {
            file_copy(source_file, destination_file, overwrite = TRUE)
            cat("Copied:", filename, "to", case, "directory\n")
        },
        error = function(e) {
            cat("Error copying", filename, ":", e$message, "\n")
            cat("Attempting alternative approach...\n")

            # Alternative: try using system cp command
            system_cmd <- paste("cp", shQuote(source_file), shQuote(destination_file))
            result <- system(system_cmd, intern = FALSE)

            if (result == 0) {
                cat("Successfully copied using system command:", filename, "\n")
            } else {
                cat("Failed to copy:", filename, "\n")
            }
        }
    )
}

cat("\nFile organization complete!\n")
cat("Organized files are located in:", data_dir, "\n")