ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

print_banner <- function(section_name) {
  cat("\n")
  cat(strrep("=", 50), "\n")
  cat("SECTION: ", section_name, "\n")
  cat(strrep("=", 50), "\n\n")
}


# Define amplification threshold
amplification_threshold <- 0

# Define ERBB2 Amplification column name
ERBB2_Amp_Col <- "ERBB2_Amp"

# Define the standard significnce threshold
SIGNIFICANCE_LEVEL <- 0.05

# ERBB2 Amplification Status Factor Constants
ERBB2_Factor_Levels <- c(0, 1)
ERBB2_Factor_Labels <- c("ERBB2-", "ERBB2+")
ERBB2_Negative_Colour <- "darkblue"
ERBB2_Positive_Colour <- "red"