# Install the required packages 
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
library(BiocManager)
BiocManager::install("limma")
install.packages(
    "https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", 
    repos = NULL, 
    type = "source"
)
BiocManager::install("edgeR")
install.packages("CIDER", repos = "http://cran.us.r-project.org")

# Stop if CIDER is not installed successfully
if (!"CIDER" %in% installed.packages()[, "Package"]) {
  stop("Package CIDER not installed successfully.")
}

# Output a dummy file to indicate that the script has completed
file.create("../results/install_confirmation.txt")
