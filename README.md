# recDplot

`recDplot` is an R package for visualizing recombination in viral sequences. This package may be used both for general analysis of recombination along with phylogenetic compatibility matrices and for inspection of specific genome regions. 

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)


## Introduction

`recDplot` provides two approaches for recombination analysis - pairwise distance corresponce plots (PDC plots) and pairwise distance deviation matrices (PDD matrices).  

### PDC plots

PDC plots utilize the pairwise genetic distances in different genome regions. In the absence of recombination, the substitutions are accumulated proportionally, and the distances in two regions will follow a linear relationship.


<img src="https://github.com/v-julia/recDplot/blob/master/images/PDCP_norec_realdata.jpg" align="center" width=500/>


If there are recombinant viruses in the dataset, the pairwise genetic distances between recombinant virus and its parents will diverge from the regression line. 

<img src="https://github.com/v-julia/recDplot/blob/master/images/PDCP_rec_realdata.jpg" align="center" width=500/>

### PDD matrices

To assess the degree of correspondence between genetic distances in different genome regions, PDC plots can be further summarized as pairwise distance deviation matrices (PDD matrices). To build such matrix, the alignment of genomes is divided into windows, then for each window pairwise genetic distances are built. For each pair of regions we estimate the root mean square error of linear regression and visualize it as a heatmap. The higher RMSE, the more pairwise genetic distances deviate from regression line, that might be the consequence of recombination.


<img src="https://github.com/v-julia/recDplot/blob/master/images/PDD_matrix.jpg" align="center" width=300/>

## Installation

`recDplot` can be installed with the `devtools` package:
```R
install.packages("devtools")
library(devtools)
install_github("v-julia/recDplot")
```

## Requirements

`recDplot` was tested with R-4.1.3 version.

## Usage




## Web version

`recDplot` functions are also available in shiny app  http://v-julia.shinyapps.io/recdplot_app.

## Contact
Your comments, bug reports, and suggestions are very welcomed. You can leave your comment via email: vjulia94@gmail.com
