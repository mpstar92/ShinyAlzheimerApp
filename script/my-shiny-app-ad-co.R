#library(XML)
#library(rentrez)

library(shiny, warn.conflicts = FALSE)
library(shinydashboard, warn.conflicts = FALSE) 

suppressWarnings({
library(igraph, warn.conflicts = FALSE)
library(DT, warn.conflicts = FALSE)
library("data.table", warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

})



#the shiny app with the preprocessed data , the app needs a minute to start, the fastest cluster is 4 , you can click on the info box and go direct to the gencard site for the gene

# show in table if gene is gatekeeper hub alzheimer gene

# can i publish the code of the app





#########
#cluster1= read.csv(file="../data/cluster/cluster1_APP_Tau_PSEN2.txt", sep = '\t')
#colnames(cluster1)[1:3]<-c("from","to","wTO")
#cluster2= read.csv(file="../data/cluster/cluster2_PSENEN_BACE1_APH1b.txt", sep = '\t')
#colnames(cluster2)[1:3]<-c("from","to","wTO")
#cluster3= read.csv(file="../data/cluster/cluster3_PSEN1_APH1A_ADAM10_ADAM17_NCSTN_BACE2.txt", sep = '\t')
#colnames(cluster3)[1:3]<-c("from","to","wTO")
#tabble= read.csv(file="../data/cluster/cluster4_Apoe4.txt", sep = '\t')

##data input

alzheimergenes_main<-read.csv(file="../data/info/AD_key_genes.txt", sep = '\t')
alzheimergenes<-read.csv(file="../data/info/AD_kegg_genes.txt", sep = '\t')

hubgenes<-read.csv(file="../data/info/hubgenes.txt", sep = ',')
gatekeepergenes<-read.csv(file="../data/info/gatekeepergenes.txt", sep = ',')


cluster1info<-read.csv(file = "../data/info/cluster1info.txt")
gc1<- read_graph("../data/cluster/clustergraphs/cluster1_APP_Tau_PSEN2_graph.txt",format = "gml")
gc1<-delete_vertex_attr(gc1,"id")
gc1_size_links <-gsize(gc1)
gc1_size_nodes <-length(V(gc1))

cluster2info<-read.csv(file = "../data/info/cluster2info.txt")
gc2<- read_graph("../data/cluster/clustergraphs/cluster2_PSENEN_BACE1_APH1b_graph.txt",format = "gml")
gc2<-delete_vertex_attr(gc2,"id")
gc2_size_links <-gsize(gc2)
gc2_size_nodes <-length(V(gc2))

cluster3info<-read.csv(file = "../data/info/cluster3info.txt")
gc3<- read_graph("../data/cluster/clustergraphs/cluster3_PSEN1_APH1A_ADAM10_ADAM17_NCSTN_BACE2_graph.txt",format = "gml")
gc3<-delete_vertex_attr(gc3,"id")
gc3_size_links <-gsize(gc3)
gc3_size_nodes <-length(V(gc3))

cluster4info<-read.csv(file = "../data/info/cluster4info.txt")
g<- read_graph("../data/cluster/clustergraphs/cluster4_Apoe4_graph.txt",format = "gml")
g<-delete_vertex_attr(g,"id")
gc4_size_links <-gsize(g)
gc4_size_nodes <-length(V(g))
##########data input end




### function for hyperlink
createLink <- function(val) {
  sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank" class="btn btn-primary">Info</a>',val)
  
}


###plot function start

plotxfunc<-function(gcx,highlight_selectx,SelectPositiv,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect){   #bool1pos ,2 neg ,x hub ,4 alz , 5 gate, 6 connect..
  ###defintion of special tables
  alztablex<-rep(NA,0)
  dataConnectx<-rep(NA,0)
  hubtablex<-rep(NA,0)
  hubtablex$genes<-""
  gateTablex<-rep(NA,0)
  
  E(gcx)$weight<-edge.attributes(gcx)$wTO
  
  
  
  E(gcx)$color<-adjustcolor("grey",.09)#edgecolor grey and opacity down
  E(gcx)$width<-0.1
  ##vertex shape triangle who are lnc coding
  vertex.attributes(gcx)$shape=ifelse(vertex.attributes(gcx)$genebiotype%in%"lncRNA","triangle","circle")
  
  
  
  E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$color<-adjustcolor("lightblue",2) #all links from selected gene blue
  E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$width<-0.5
  
  
  ## postiv correlated genes darkgreen
  
  some<-as_edgelist(gcx)[,1]  #nodes from
  some2<-as_edgelist(gcx)[,2]#nodes to
  some3<-E(gcx)$wTO       # wto
  
  if(SelectPositiv==TRUE){
    v=1
    while(v<=length(some)){
      if(some3 [v]>=range2){
          if(some[v] == highlight_selectx|some2[v] == highlight_selectx){
        E(gcx)[v]$color<-adjustcolor("darkgreen",.9)
        E(gcx)[v]$width<-abs(E(gcx)$wTO[v])*2
      }
      }
      v=v+1
    }

  }
  ####negativ correlated genes red
  if(SelectNegativ==TRUE){
    w=1
 
    while(w<=length(some)){
      if(some3 [w] <= range1){
        
        if(some[w] == highlight_selectx|some2[w] == highlight_selectx){
          E(gcx)[w]$color<-adjustcolor("red",.9)
          E(gcx)[w]$width<-abs(E(gcx)$wTO[w])*2
        }
      }
      w=w+1
    }
  }
  
  if(checkboxHub==TRUE){#hubgenes
    hubtablex<-completetablex[order(completetablex$degree,decreasing = TRUE ),]
  }
  
  if(checkboxAlzheimer==TRUE){##alzheimer
    w=1
    while(w<=length(alzheimergenes$hgnc_symbol)){
      alztablex<-append(alztablex,alzheimergenes$ensembl_gene_id[w])
      w=w+1
    }
  }
  if(GateKeep==TRUE){   #gatekeepers
    s=1
    while(s<=length(networkx$Gene)){
      gateTablex<-append(gateTablex,networkx$Gene[s])
      s=s+1
    }
  }
  if(checkboxConnect==TRUE){##correlated genes
    dataConnectx<-append(dataConnectx,ends(gcx, E(gcx)[from(highlight_selectx)])[,1])  #"Gene from"
    dataConnectx<-append(dataConnectx,ends(gcx, E(gcx)[from(highlight_selectx)])[,2] ) #"Gene to"
    
    
    #i=1
    #while(i< length(clusterx$from)){
     # if(clusterx$from[i]==highlight_selectx){
    #    dataConnectx<-append(dataConnectx,clusterx$to[i])
     # }
     # if(clusterx$to[i]==highlight_selectx){
     #   dataConnectx<-append(dataConnectx,clusterx$from[i])
    #  }
     # i=i+1
   # }
    
  }##correlated genes
  hub<-ceiling(length(V(gcx))*0.01)
  ##size and color of nodes  
  ## size  selected gene 5,   hub, alzhheimer, gatekeepers genes 4 , correlated genes x
  ## color selected gene orange, hub gene pink, alzheimer genes black , gatekeepers blue
  vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,5,
                                     ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:hub],4,
                                            ifelse(vertex.attributes(gcx)$name%in%alztablex,4,
                                                   ifelse(vertex.attributes(gcx)$name%in%gateTablex,4,
                                                          ifelse(vertex.attributes(gcx)$name%in%dataConnectx,3,1)))))
  
  vertex.attributes(gcx)$color=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,"orange",
                                      ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:hub],"pink",
                                             ifelse(vertex.attributes(gcx)$name%in%alztablex,"black",
                                                    ifelse(vertex.attributes(gcx)$name%in%gateTablex,"blue","white"))))
  
  
  #####
 
  
  #####
  x<-paste("wTO < ",range1,sep = " ")
  z<-paste("wTO > ",range2,sep = "  ")
 
  y2<-paste(completetablex$names[completetablex$genes==highlight_selectx])
 
  plot.igraph(gcx,vertex.label=NA,edge.width=E(gcx)$width)
  #plot with labels
  #plot.igraph(gcx,vertex.label=ifelse(vertex.attributes(gcx)$name%in%c(hubtablex$genes[1:10],alztablex,gateTablex),V(gcx)$hgncsymbol,NA),vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=1,edge.width=E(gcx)$width)
  
  legend("topleft", legend=c(x,z,"Hub Genes","AD Pathway Genes(KEGG)","Gatekeeper Genes",y2,"lnc RNA coding Genes","Protein coding Genes"),lty=c(1,1,0,0,0,0,0,0),pch = c(NA,NA,19,19,19,19,6,19),col = c("red","darkgreen","pink","black","blue","orange","grey","grey"),bg=rgb(1,1,0,alpha=0.05),pt.cex=c(2,2,2,2,2,2,2,2,2))
  #legend("topright", legend=c(x,z,"Hubgenes","AlzheimerGenes","Gatkeepers",highlight_selectx,"LNC Coding","Protein Coding"),lty=c(1,1,1,0,0,0,0,0,0),pch = c(NA,NA,NA,19,19,19,19,6,19),col = c("red","lightblue","darkgreen","pink","black","blue","orange","grey","grey"),bg=rgb(1,1,0,alpha=0.05),pt.cex=c(2,2,2,2,2,2,2,2,2))
  
}  


#####plot function end

###correlated genes Function++++ STart

correlatedGenesFunc<-function(graph,clusterxinfo,highlight_selectx){
  wTO<-E(graph)[from(highlight_selectx)]$wTO  #0                                
  Gene<-ends(graph, E(graph)[from(highlight_selectx)])[,1]  #"Gene from"
  GenesTo<-ends(graph, E(graph)[from(highlight_selectx)])[,2]  #"Gene to"
  dataconnectx<-data.frame(Gene,GenesTo,wTO)
 
  ##order the table so that all selected genes are in row two
  ## highlightselect = ens name
  
  
  for(i in 1:length(dataconnectx$Gene)){
    if (dataconnectx$Gene[i] == highlight_selectx ) dataconnectx$Gene[i]<- dataconnectx$GenesTo [i]
  }
  
  
  dataconnectx<-dataconnectx[order(dataconnectx$wTO,decreasing=TRUE),]### highest wto score first
  dataconnectx$GeneType<-""
  dataconnectx$description<-""
  dataconnectx$Genecard<-""
  e=1
  while(e<length(dataconnectx$Gene)+1){
    dataconnectx$GeneType[e]   =clusterxinfo$gene_biotype[clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    if(highlight_selectx == dataconnectx$Gene[e]){ dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    dataconnectx$Genecard[e]   =createLink(dataconnectx$GenesTo[e])}
    else{                                              dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
    dataconnectx$Genecard[e]=createLink(dataconnectx$Gene[e])}
    dataconnectx$Gene[e]       =clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e]    ]
    dataconnectx$GenesTo[e    ]=clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    
    e=e+1
  }    #add description and gene name 
  
 
  dataconnectx[,-2]
  
  #output dataconnect table with all information
}


###correlated gens function +++ END


####triangle shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

##################################################



## betweenes

betweenFunction<-function(graph){
  
  # use abs values of wto/ won't work with negative values
  E(graph)$weight <-abs(E(graph)$wTO)
  E(graph)$wTO <-E(graph)$weight
  btw <- betweenness(graph, directed=F)
  V(graph)$betweenness<-btw
  #other betweeness type
  btw_edge <- edge_betweenness(graph, directed=F)
  E(graph)$btw_edge <- btw_edge
  
  ###############################################
  ## Find GATEKEEPERS
  #### Prepare dataframe
  
  network <- as.data.frame(btw)
  network<- setDT(network, keep.rownames = TRUE)[]
  colnames(network)<-c("Gene", 'score')
  
  # Filter for top 1%
  Top_betw_network <- network[network$score >= quantile(network$score, 0.99),]
  Top_betw_network<-Top_betw_network[order(Top_betw_network$score , decreasing = TRUE ),]
  
  Top_betw_network
}



###data refining

completetable1<-data.frame(genes=vertex.attributes(gc1)$name,degree=degree(gc1),names=vertex.attributes(gc1)$hgncsymbol,genebiotype=vertex.attributes(gc1)$genebiotype,description="")
completetable1<-completetable1[order(completetable1$genes), ]
completetable1$description<-cluster1info$description
completetable1<-completetable1[order(completetable1$degree), ]
network1<-betweenFunction(gc1)

completetable2<-data.frame(genes=vertex.attributes(gc2)$name,degree=degree(gc2),names=vertex.attributes(gc2)$hgncsymbol,genebiotype=vertex.attributes(gc2)$genebiotype,description="")
completetable2<-completetable2[order(completetable2$genes), ]
completetable2$description<-cluster2info$description
completetable2<-completetable2[order(completetable2$degree), ]
network2<-betweenFunction(gc2)

completetable3<-data.frame(genes=vertex.attributes(gc3)$name,degree=degree(gc3),names=vertex.attributes(gc3)$hgncsymbol,genebiotype=vertex.attributes(gc3)$genebiotype,description="")
completetable3<-completetable3[order(completetable3$genes), ]
completetable3$description<-cluster3info$description
completetable3<-completetable3[order(completetable3$degree), ]
network3<-betweenFunction(gc3)


completetable4<-data.frame(genes=vertex.attributes(g  )$name,degree= degree(g)  ,names=vertex.attributes(g  )$hgncsymbol,genebiotype=vertex.attributes(g  )$genebiotype,description="")
completetable4<-completetable4[order(completetable4$genes), ]
completetable4$description<-cluster4info$description
completetable4<-completetable4[order(completetable4$degree), ]
network4<-betweenFunction(g)


hubbie<-"The Hub Genes are the top 1 % Genes with the most links to other Genes"
gaties<-"Gatekeeper Genes are the top 1% of Genes with the highest betweeness centrality "
gaties2<-"the betweeness centrality is the number of shortest paths that pass throgh the vertex"
num_nodes<-2899
num_edges<-717243


# Define UI ----
ui <- fluidPage(
  
  
  titlePanel("Co-Expression Alzheimer App"),
  
  sidebarLayout(
    sidebarPanel(
      ###
      #tags$a(href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=ECE1", "Click Here!"),
      ####
      
      h5("1. Select your Cluster of interest from the tabs"),
      h5("2. Select Network display options"),
      h5("3. Select your Gene of interst from the cluster Gene list"),
      h5("4. Press   PLOT!   to render the plot with your parameters"),
      h3(""),
      h4("AD Pathway Genes (KEGG)"),
      checkboxInput("checkboxAlzheimer","Yes",value = FALSE),
      h4("Hub Genes"),
      checkboxInput("checkboxHub", "Yes", value = FALSE),
      h4("Gatekeeper Genes"),
      checkboxInput("GateKeep","Yes",value=FALSE),
      
      

      sliderInput("range", 
                  label = "Show wTO above this range:",
                  min = -1, max = 1, value = c(-0.8, 0.8),step = 0.1),
      

      selectInput("select1",
                  h4("Cluster 1 Gene list"), 
                  choices = c(completetable1$names), selected = "DDOST",multiple=FALSE),
      selectInput("select2",
                  h4("Cluster 2 Gene list"), 
                  choices = c(completetable2$names), selected = "ARL2-SNX15",multiple=FALSE),
      selectInput("select3",
                  h4("Cluster 3 Gene list"), 
                  choices = c(completetable3$names), selected = "C1QTNF3-AMACR",multiple=FALSE),
      selectInput("select",
                  h4("Cluster 4 Gene list"), 
                  choices = c(completetable4$names), selected = "ECE1",multiple=FALSE),
      
      h4("Highlight positively links for selected Gene"),
      checkboxInput("checkboxSelectPositiv", "Positiv", value = TRUE),
      h4("Highlight negatively links for selected Gene"),
      checkboxInput("checkboxSelectNegativ", "Negativ", value = TRUE),
      h4("First Neighbor Genes for selected Gene"),
      checkboxInput("checkboxConnect","Yes",value = FALSE),
      
      actionButton("goButton", "PLOT!"),
      
      
   
      ########
      
    ),
    
    mainPanel(
  

               
                tabsetPanel(type = "pills",    #1
                  tabPanel("AD Network",tabsetPanel(
                      tabPanel("Network", img(src = "Ad_network_graph.png") ) ,  
                      tabPanel("Network Information",h4("Number of Genes",num_nodes),
                               h4("Number of links", num_edges),br(),
                               h4("Split in 4 clusters, named after the AD Pathway Genes(KEGG)"),br(),
                               h4("Cluster 1: APP, Tau, PSEN2"),br(),
                               h4("Cluster 2: PSENEN, BACE1, APH1b"),br(),
                               h4("Cluster 3: PSEN1, APH1A, ADAM10, ADAM17, NCSTN, BACE2"),br(),
                               h4("Cluster 4: Apoe4 ")),
                      tabPanel("AD Pathway Genes (KEGG)",br(),uiOutput("kegg_hyper"),br(),img(src = "kegg_pathway.png"),br(),dataTableOutput("alzheimerGenesX")),
                      tabPanel("Hub Genes",h4(hubbie),br(),dataTableOutput("hubgenesX")),
                      tabPanel("Gatekeeper Genes",h4(gaties),h4(gaties2),br(),dataTableOutput("gatekeepersX")),
                       
                      )),
                  tabPanel("Clusters",tabsetPanel(
                    tabPanel("Cluster1",tabsetPanel(
                      tabPanel("Hub Genes",dataTableOutput("hubGeneDesciptions1")),
                      tabPanel("Gatekeeper Genes",dataTableOutput("Gatekeeper1")),
                      tabPanel("First Neighbor Genes",dataTableOutput("connectedGenes1")),
                      tabPanel("Cluster 1 Information",h4("Cluster 1: APP, Tau, PSEN2"),br(),h4("Gene information of selected Gene"),dataTableOutput("gene1desc"),br(),h4("#Genes: ",gc1_size_nodes),h4("#Links: ",gc1_size_links)),
                      tabPanel("Network",plotOutput("graphic1"))
                    )),
                    
                    
                    tabPanel("Cluster2",tabsetPanel(
                      tabPanel("Hub Genes",dataTableOutput("hubGeneDesciptions2")),
                      tabPanel("Gatekeeper Genes",dataTableOutput("Gatekeeper2")),
                      tabPanel("First Neighbor Genes",dataTableOutput("connectedGenes2")),
                      tabPanel("Cluster 2 Information",h4("Cluster 2: PSENEN, BACE1, APH1b"),br(),h4("Gene information of selected Gene"),dataTableOutput("gene2desc"),br(),h4("#Genes: ",gc2_size_nodes),h4("#Links: ",gc2_size_links)),
                      tabPanel("Network",plotOutput("graphic2"))
                    )),
                    tabPanel("Cluster3",tabsetPanel(
                      tabPanel("Hub Genes",dataTableOutput("hubGeneDesciptions3")),
                      tabPanel("Gatekeeper Genes",dataTableOutput("Gatekeeper3")),
                      tabPanel("First Neighbor Genes",dataTableOutput("connectedGenes3")),
                      tabPanel("Cluster 3 Information",h4("Cluster 3: PSEN1, APH1A, ADAM10, ADAM17, NCSTN, BACE2"),br(),h4("Gene information of selected Gene"),dataTableOutput("gene3desc"),br(),h4("#Genes: ",gc3_size_nodes),h4("#Links: ",gc3_size_links)),
                      tabPanel("Network",plotOutput("graphic3"))
                    )),
                    tabPanel("Cluster4",tabsetPanel(
                      tabPanel("Hub Genes",dataTableOutput("hubGeneDesciptions4")),
                      tabPanel("Gatekeeper Genes",dataTableOutput("Gatekeeper4")),
                      tabPanel("First Neighbor Genes",h3("First neighbor Genes"),dataTableOutput("connectedGenes")),
                      tabPanel("Cluster 4 Information",h4("Cluster 4: Apoe4 "),br(),h4("Gene information of selected Gene"),dataTableOutput("gene4desc"),br(),h4("#Genes: ",gc4_size_nodes),h4("#Links: ",gc4_size_links)),
                      tabPanel("Network",plotOutput("graphic4"))
                    ))
                    
                  ))
                  #
                  #tabPanel("Test",dataTableOutput("test")),
                  #tabPanel("Test2",dataTableOutput("test2")),
                  #
                  
                 
                 
                  
               ) #1
            
               
      )
  )
)

# Define server logic ----
options(shiny.maxRequestSize=90*1024^2)
server <- function(input, output) {
  ##all
  
  ###
  output$alzheimerGenesX<-renderDataTable({
    for(x in 1:length(alzheimergenes_main$hgnc_symbol))alzheimergenes_main$genecard[x]<-createLink(alzheimergenes_main$hgnc_symbol[x])
    alzheimergenes_main
    
  }, options=list(pageLength =25,scrollY =TRUE),escape =FALSE )
  
  
  output$hubgenesX<-renderDataTable({
    for(x in 1:length(hubgenes$hgnc_symbol))hubgenes$genecard[x]<-createLink(hubgenes$hgnc_symbol[x])
    hubgenes<-hubgenes[order(hubgenes$hub_score,decreasing = TRUE ),]
    hubgenes
    
  }, options=list(pageLength =50,scrollY =TRUE),escape =FALSE )
  
  output$gatekeepersX<-renderDataTable({
    for(x in 1:length(gatekeepergenes$hgnc_symbol))gatekeepergenes$genecard[x]<-createLink(gatekeepergenes$hgnc_symbol[x])
    gatekeepergenes<-gatekeepergenes[order(gatekeepergenes$between_score,decreasing = TRUE ),]
    gatekeepergenes
    
  }, options=list(pageLength =50,scrollY =TRUE),escape =FALSE )
  ####################
  output$test3<-renderDataTable({
    for(x in 1:length(alzheimergenes$hgnc_symbol))alzheimergenes$genecard[x]<-createLink(alzheimergenes$hgnc_symbol[x])
    alzheimergenes
  }, escape = FALSE)
  ################
  
  
  ##table from input
  output$tbs <-renderDataTable({
    data <-input$file
    if(is.null(data)){return()}
    read.table(data$datapath,sep="")
  })
  
  ###cluster 1 gene desription text
  output$gene_description <- renderText({##cluster4
    highlight_select<-cluster4info$ensembl_gene_id[cluster4info$hgnc_symbol==input$select]
    paste0(cluster4info$hgnc_symbol[cluster4info$ensembl_gene_id==highlight_select],cluster4info$description[cluster4info$ensembl_gene_id==highlight_select],highlight_select,sep="\n")
    #res <- entrez_search(db = "gene", term = highlight_select,api_key =key1) #search gene db for ids
    #recs <- entrez_fetch(db = "gene", id = res$ids[1:25], rettype = "xml", parsed = TRUE) # get fullrecord with ids
    #esums <- entrez_summary(db = "gene", id = res$ids[1])
    #paste(highlight_select,esums$summary,sep="\n")
    #searches for id and than the id in ncbi for the desrcription of the gene
  })
  
  ###### What are Hub , gatekeeper and alzheimer genes?
  
  output$whatis_hub<-renderText(
    paste("What is a Hub Gene","something",sep = "\n")
  )
  output$kegg_hyper <- renderUI({
    url<-a("AD Pathway Genes KEGG", href="https://www.genome.jp/pathway/hsa05010", target="_blank" )
    tagList( url)
  })
  
  
  
  
  ####################
  output$gene4desc<-renderDataTable({
    highlight_select<-cluster4info$ensembl_gene_id[cluster4info$hgnc_symbol==input$select]
    infoTable<-data.frame(Ensamble_ID="",Gene_name="",Gene_description="")
    infoTable$Ensamble_ID<-highlight_select
    infoTable$Gene_name<-cluster4info$hgnc_symbol[cluster4info$ensembl_gene_id==highlight_select]
    infoTable$Gene_description<-cluster4info$description[cluster4info$ensembl_gene_id==highlight_select]
    infoTable$Genecard<- createLink(input$select)
    
    return(infoTable)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  output$gene1desc<-renderDataTable({
    highlight_select1<-cluster1info$ensembl_gene_id[cluster1info$hgnc_symbol==input$select1]
    infoTable1<-data.frame(Ensamble_ID="",Gene_name="",Gene_description="")
    infoTable1$Ensamble_ID<-highlight_select1
    infoTable1$Gene_name<-cluster1info$hgnc_symbol[cluster1info$ensembl_gene_id==highlight_select1]
    infoTable1$Gene_description<-cluster1info$description[cluster1info$ensembl_gene_id==highlight_select1]
    infoTable1$Genecard<- createLink(input$select1)
    infoTable1
    ##put gene description4 in table  cluster1
  }, escape = FALSE)
  output$gene2desc<-renderDataTable({
    
    highlight_select2<-cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2]
    
    infoTable2<-data.frame(Ensamble_ID="",Gene_name="",Gene_description="")
    infoTable2$Ensamble_ID<-highlight_select2
    infoTable2$Gene_name<-cluster2info$hgnc_symbol[cluster2info$ensembl_gene_id==highlight_select2]
    infoTable2$Gene_description<-cluster2info$description[cluster2info$ensembl_gene_id==highlight_select2]
    infoTable2$Genecard<- createLink(input$select2)
    infoTable2
    ##put gene description4 in table  cluster2
  }, escape = FALSE)
  output$gene3desc<-renderDataTable({
    highlight_select3<-cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3]
    infoTable3<-data.frame(Ensamble_ID="",Gene_name="",Gene_description="")
    infoTable3$Ensamble_ID<-highlight_select3
    infoTable3$Gene_name<-cluster3info$hgnc_symbol[cluster3info$ensembl_gene_id==highlight_select3]
    infoTable3$Gene_description<-cluster3info$description[cluster3info$ensembl_gene_id==highlight_select3]
    infoTable3$Genecard<- createLink(input$select3)
    infoTable3
    ##put gene description4 in table  cluster3
  }, escape = FALSE)
 
  
  
  
  
  ###cluster4
  ##top 10 genes with most degree plus information Cluster 4
  output$hubGeneDesciptions4<-renderDataTable({
    hubtable4<-completetable4[order(completetable4$degree,decreasing = TRUE ),]
    for(x in 1:10)hubtable4$genecard[x]<-createLink(hubtable4$genes[x])
    hubtable4<-subset(hubtable4,select = -1)
   
    head(hubtable4,ceiling(length(V(g))*0.01))
   }, escape = FALSE)
  
  ##cluster 1
  output$hubGeneDesciptions1<-renderDataTable({
    hubtable1<-completetable1[order(completetable1$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable1$genecard[x]<-createLink(hubtable1$genes[x])
    hubtable1<-subset(hubtable1,select = -1)
    
    head(hubtable1,ceiling(length(V(gc1))*0.01))
  }, escape = FALSE)
  ##cluster 2
  output$hubGeneDesciptions2<-renderDataTable({
    hubtable2<-completetable2[order(completetable2$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable2$genecard[x]<-createLink(hubtable2$genes[x])
    hubtable2<-subset(hubtable2,select = -1)
    
    head(hubtable2,ceiling(length(V(gc2))*0.01))
  }, escape = FALSE)
  ##cluster 3
  output$hubGeneDesciptions3<-renderDataTable({
    hubtable3<-completetable3[order(completetable3$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable3$genecard[x]<-createLink(hubtable3$genes[x])
    hubtable3<-subset(hubtable3,select = -1)
    
    head(hubtable3,ceiling(length(V(gc3))*0.01))
  }, escape = FALSE)
  #########hub####
  ###GAtekeeper####
  
  output$Gatekeeper1<-renderDataTable({
    gateTable1<-network1
    q=1
    while(q<=length(gateTable1$Gene)){
      gateTable1$name[q]<-cluster1info$hgnc_symbol[cluster1info$ensembl_gene_id==gateTable1$Gene[q]]
      gateTable1$genebiotype[q]<-cluster1info$gene_biotype[cluster1info$ensembl_gene_id==gateTable1$Gene[q]]
      gateTable1$description[q]<-cluster1info$description[cluster1info$ensembl_gene_id==gateTable1$Gene[q]]
      gateTable1$genecard[q]<- createLink(gateTable1$Gene[q])
      q=q+1
    }
    
    gateTable1
  }, escape = FALSE)
  output$Gatekeeper2<-renderDataTable({
    gateTable2<-network2
    q=1
    while(q<=length(gateTable2$Gene)){
      gateTable2$name[q]<-cluster2info$hgnc_symbol[cluster2info$ensembl_gene_id==gateTable2$Gene[q]]
      gateTable2$genebiotype[q]<-cluster2info$gene_biotype[cluster2info$ensembl_gene_id==gateTable2$Gene[q]]
      gateTable2$description[q]<-cluster2info$description[cluster2info$ensembl_gene_id==gateTable2$Gene[q]]
      gateTable2$genecard[q]<- createLink(gateTable2$Gene[q])
      q=q+1
    }

    gateTable2
  }, escape = FALSE)
  output$Gatekeeper3<-renderDataTable({
    gateTable3<-network3
    q=1
    while(q<=length(gateTable3$Gene)){
      gateTable3$name[q]<-cluster3info$hgnc_symbol[cluster3info$ensembl_gene_id==gateTable3$Gene[q]]
      gateTable3$genebiotype[q]<-cluster3info$gene_biotype[cluster3info$ensembl_gene_id==gateTable3$Gene[q]]
      gateTable3$description[q]<-cluster3info$description[cluster3info$ensembl_gene_id==gateTable3$Gene[q]]
      gateTable3$genecard[q]<- createLink(gateTable3$Gene[q])
      q=q+1
    }

    gateTable3
  }, escape = FALSE)
  output$Gatekeeper4<-renderDataTable({
    gateTable4<-network4
    q=1
    while(q<=length(gateTable4$Gene)){
      gateTable4$name[q]<-cluster4info$hgnc_symbol[cluster4info$ensembl_gene_id==gateTable4$Gene[q]]
      gateTable4$genebiotype[q]<-cluster4info$gene_biotype[cluster4info$ensembl_gene_id==gateTable4$Gene[q]]
      gateTable4$description[q]<-cluster4info$description[cluster4info$ensembl_gene_id==gateTable4$Gene[q]]
      gateTable4$genecard[q]<- createLink(gateTable4$Gene[q])
      q=q+1
    }
    
    gateTable4
  }, escape = FALSE)
  
  ###Gatekeeper#########
  
  #### Table with geneinformation For The connected genes and all degree of selected gene
  ##cluster4
  output$connectedGenes <-renderDataTable({
    highlight_select<- cluster4info$ensembl_gene_id[cluster4info$hgnc_symbol==input$select]

    
      dataconnect4<-correlatedGenesFunc(g,cluster4info,highlight_select)
      #######################################################
      
      
      
      
      ############################  g tabble
    
  }, escape = FALSE)

  
  ###cluster 1 correlated genes
  output$connectedGenes1 <-renderDataTable({
    highlight_select1<- cluster1info$ensembl_gene_id[cluster1info$hgnc_symbol==input$select1]
    
    
      dataconnect1<-correlatedGenesFunc(gc1,cluster1info,highlight_select1)
    
  }, escape = FALSE)
  ###cluster 2 correlated genes
  output$connectedGenes2 <-renderDataTable({
    highlight_select2<- cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2]
    
    
      dataconnect2<-correlatedGenesFunc(gc2,cluster2info,highlight_select2)
    
  }, escape = FALSE)
  
  ###cluster 3 correlated genes
  output$connectedGenes3 <-renderDataTable({
    highlight_select3<- cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3]
    
   
      dataconnect3<-correlatedGenesFunc(gc3,cluster3info,highlight_select3)
    
  }, escape = FALSE)
  
  
  ########### information ends##########
  
  
  #####
  
  ################graphs####start#######
  
  ###graph1
  output$graphic1 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select1=isolate(cluster1info$ensembl_gene_id[cluster1info$hgnc_symbol==input$select1])
    
    plot1<-isolate(plotxfunc(gc1,highlight_select1,input$checkboxSelectPositiv,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable1,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network1,input$checkboxConnect))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)

  },height = 1000,width = 1000)
  
  ##graph2 
  output$graphic2 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select2=isolate(cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2])
    
    plot2<-isolate(plotxfunc(gc2,highlight_select2,input$checkboxSelectPositiv,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable2,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network2,input$checkboxConnect))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
  
  },height = 1000,width = 1000)
  
  ####graph3
  output$graphic3 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    
    highlight_select3=isolate(cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3])
    
    plot3<-isolate(plotxfunc(gc3,highlight_select3,input$checkboxSelectPositiv,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable3,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network3,input$checkboxConnect))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
 
  },height =1500,width = 1500)
  ###
 
  ###graph4
  output$graphic4 <-renderPlot({  ###
    set.seed(1277)
    
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select=isolate(cluster4info$ensembl_gene_id[cluster4info$hgnc_symbol==input$select])
   
    
    plot4<-isolate(plotxfunc(g,highlight_select,input$checkboxSelectPositiv,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable4,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network4,input$checkboxConnect) )
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1000,width = 1000)
  
  #############test
  
 
  
  ##########test

}



# Run the app ----
shinyApp(ui = ui, server = server,options =list (launch.browser =TRUE))


