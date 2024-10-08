**ADNet User Guide**

ADNet is a Shiny app designed to facilitate exploration of the temporal cortex Alzheimer's disease (AD) and control gene co-expression networks. It provides a user-friendly interface to examine detailed network information, offering valuable insights for research and further analysis. For more in-depth information about the research and methodology behind ADNet, please refer to the associated paper: [paper title-link].

In ADNet, additional key data is embedded, such as each gene's known function, degree, hub status, differential expression (DEG), gene type (protein-coding or lncRNA), clusters and involvement in Alzheimer's disease pathways according to the KEGG database. This information is presented in a clear, easily navigable table format, allowing quick access to essential metrics and simplifying the process of exploring the gene networks.


Table of Contents and Additional Tutorials
_____________________________________________________

1. Installation
2. Quick Start Guide






1.Installation 


RStudio is preferable to use.
First, users can run the following code to check if the packages required by Shiny App ADNet exist and install them if required:

    reqPkg = c("data.table", "Matrix", "igraph", "qgraph", "DT", 
           "dplyr", "shinydashboard", "biomaRt")
    newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
    if(length(newPkg)){install.packages(newPkg)}






2.Quick Start Guide



Shiny App ADNet can then be installed from GitHub as follows:

https://github.com/mpstar92/Shiny_App_ADnet  
![Alt text](data/images/code.jpg)

Click on the green button <>Code on the top right and download Zip

Extract the folder and open the project Shiny_app_ADnet.Rproj
Open file script/shiny_app_ADnet.R
Start the application by clicking on "Run App" in the top right corner. (App will start in a browser after a ~ 30 seconds)


In the App:

In the top center you can click through the different tabs.

In the AD Network tab you can change between the gene co-expression network, information about the network, the AD KEGG pathway and the core genes.
![Alt text](data/images/tabs_ad.jpg)




In the Clusters tab you can select a cluster, and for each cluster is a tab for 



 cluster information including basic information about the selected gene
 ![Alt text](data/images/cluster_information.jpg)
 
 the first neighbor, core and AD KEGG pathway genes 

 ![Alt text](data/images/core.jpg)

 In each table you can search for specific genes or terms.

 
 And the respective interactive plot.

![Alt text](data/images/network.jpg)




To get a visualization of your cluster of choosing:

Select the Tabs Clusters; the cluster; Network visualization:

![Alt text](data/images/tabs_clusters.jpg)


Select in the left sidebar: 

![Alt text](data/images/sidebar.jpg)

which node attributes to display and the gene of interest in that cluster,
Press PLOT! at the bottom of the sidebar to visualize the cluster (it may take a couple of seconds and if you change parameters press PLOT! again)

the color of the border of a node:

    black = the selected gene
    red = up-regulated
    blue = down regulated
 
the size of the nodes from big to small:

    core genes
    selected gene
    AD pathway genes
    first neighbor genes
    no attribute

the color code for the associated GO terms

![Alt text](data/images/colorcode.jpg)


the links between nodes are colored:

    red = AD-specific
    green = control specific
    grey else

if you select to show the links of selected gene with wTO range the graph will only show the links from the selected genes and not AD/control-specific links.



