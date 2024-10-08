Shiny App ADNet is an R project that allows users to visualize interactively in a Shiny-based web application the gene co-expression of the TCX AD/Control Consensus Networks 
by visualizing the gene co-expression of all genes within three AD and four control clusters, 
additionally it offers detailed information about the networks in a user-friendly table format, giving you quick access to key metrics and insights for further analysis.

Key features of Shiny App ADNet include:

    1. Written in R and uses the Shiny package, allowing for easy sharing on online platforms e.g. shinyapps.io and Amazon Web Services (AWS), or be hosted via Shiny Server

    2. Information about the AD network including the core genes and the AD KEGG pathway

    3. Information about the gene co-expression clusters of the AD and the control networks including gene descriptions, first neighbor, core, and AD Pathway genes

    4. Visualization of the AD Network 

    5. Interactive plots of the gene co-expression networks for each cluster, with selectable parameters: AD pathway genes, core genes, first neighbor genes, and strength of wTO correlation.


Table of Contents and Additional Tutorials
_____________________________________________________

1. Installation
2. Quick Start Guide
3. Frequently Asked Questions


Installation 

RStudio is preferable to use.
First, users can run the following code to check if the packages required by Shiny App ADNet exist and install them if required:

reqPkg = c("data.table", "Matrix", "igraph", "qgraph", "DT", 
           "dplyr", "shinydashboard", "biomaRt")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}




Quick Start Guide


Shiny App ADNet can then be installed from GitHub as follows:

Download Zip folder from https://github.com/mpstar92/Shiny_App_ADnet
Extract the folder and open the project Shiny_app_ADnet.Rproj
Open file script/shiny_app_ADnet.R
Start the application by clicking on "Run App" in the top right corner. (App will start in a browser after a ~ 30 seconds)

In the App:
