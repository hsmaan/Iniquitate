## Imbalanced integration guidelines 

![Example Image](assets/Figure_8_Revised.png)

This folder contains an rmarkdown example for using the guidelines shown in the [Imbalanced integration manuscript](https://www.biorxiv.org/content/10.1101/2022.10.06.511156v3). The tutorial is available as a rendered rmarkdown html document (guidelines.html), which can be downloaded and viewed with a web browser.

The full details of the imbalanced integration guidelines can be found in the manuscript in Results Section: Guidelines for imbalanced single-cell data integration, and the associated Supplementary Table 2.

### Viewing the rendered markdown

#### Through the html

1. Download the `docs/guidelines.html` file

2. Open the file with a web browser (chrome, safari, firefox, etc.)

#### Through the pdf 

1. Download the 'docs/guidelines.pdf' file and view with a PDF browser

Please note that the pdf rendering of the code might not be formatted well as the file was originally rendered as an html vignette

### Running the tutorial

1. Switch to the `docs` folder. Install the conda environment for the guidelines:
```
conda env install -f envs/tutorial.yaml
```

2. Activate the conda environment:
```
conda activate
```

3. Activate R and Install rmarkdown and tinytex
```
R
```
```
install.packages("rmarkdown", dep = TRUE)
install.packages("tinytex")
tinytex::install_tinytex() 
```

4. Run the rmarkdown chunk by chunk in [Rstudio](https://posit.co/downloads/) 
***
