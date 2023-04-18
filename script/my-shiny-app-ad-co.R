library(shiny)
library(shinydashboard)
library(igraph)
library(XML)
library(rentrez)
library("data.table")
library(dplyr)
library(DT)
#the shiny app with the preprocessed data , the app needs a minute to start, the fastest cluster is 4 , you can click on the info box and go direct to the gencard site for the gene




#####HUBGENES reparieren


set_entrez_key("275b491fc56280b6c69e8b23b1073a142108")
key1<-"275b491fc56280b6c69e8b23b1073a142108"
#########
###cluster         in    /data/cluster/
###info+ alzheimer in    /data/info/
###graphs          in    /data/cluster/clustergraphs/


##data input
alzheimergenes<-read.csv(file="../data/info/AD_key_genes.txt", sep = '\t')

cluster1= read.csv(file="../data/cluster/cluster1_APP_Tau_PSEN2.txt", sep = '\t')
colnames(cluster1)[1:3]<-c("from","to","wTO")

cluster1info<-read.csv(file = "../data/info/cluster1info.txt")

gc1<- read_graph("../data/cluster/clustergraphs/cluster1_APP_Tau_PSEN2_graph.txt",format = "gml")
gc1<-delete_vertex_attr(gc1,"id")


cluster2= read.csv(file="../data/cluster/cluster2_PSENEN_BACE1_APH1b.txt", sep = '\t')
colnames(cluster2)[1:3]<-c("from","to","wTO")

cluster2info<-read.csv(file = "../data/info/cluster2info.txt")

gc2<- read_graph("../data/cluster/clustergraphs/cluster2_PSENEN_BACE1_APH1b_graph.txt",format = "gml")
gc2<-delete_vertex_attr(gc2,"id")

cluster3= read.csv(file="../data/cluster/cluster3_PSEN1_APH1A_ADAM10_ADAM17_NCSTN_BACE2.txt", sep = '\t')
colnames(cluster3)[1:3]<-c("from","to","wTO")

cluster3info<-read.csv(file = "../data/info/cluster3info.txt")

gc3<- read_graph("../data/cluster/clustergraphs/cluster3_PSEN1_APH1A_ADAM10_ADAM17_NCSTN_BACE2_graph.txt",format = "gml")
gc3<-delete_vertex_attr(gc3,"id")

tabble= read.csv(file="../data/cluster/cluster4_Apoe4.txt", sep = '\t')

cluster4info<-read.csv(file = "../data/info/cluster4info.txt")

g<- read_graph("../data/cluster/clustergraphs/cluster4_Apoe4_graph.txt",format = "gml")
g<-delete_vertex_attr(g,"id")

##########data input end




### function for hyperlink
createLink <- function(val) {
  sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank" class="btn btn-primary">Info</a>',val)
  
}


###plot function start

plotxfunc<-function(gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect){   #bool1pos ,2 neg ,x hub ,4 alz , 5 gate, 6 connect..
  ###defintion of special tables
  alztablex<-rep(NA,0)
  dataConnectx<-rep(NA,0)
  hubtablex<-rep(NA,0)
  hubtablex$genes<-""
  gateTablex<-rep(NA,0)
  
  
  
  
  E(gcx)$color<-adjustcolor("grey",.09)#edgecolor grey and opacity down
  E(gcx)$width<-0.1
  ##vertex shape triangle who are lnc coding
  vertex.attributes(gcx)$shape=ifelse(vertex.attributes(gcx)$genebiotype%in%"lncRNA","triangle","circle")
  
  
  
  E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$color<-adjustcolor("lightblue",2) #all links from selected gene blue
  E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$width<-0.5
  
  
  ## postiv correlated genes darkgreen
  if(SelectPositiv==TRUE){
    v=1
    while(v<=length(clusterx$from)){
      if(clusterx$wTO[v]>range2){
        if(clusterx$from[v]==highlight_selectx|clusterx$to[v]==highlight_selectx){
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
    while(w<=length(clusterx$from)){
      if(clusterx$wTO[w]<range1){
        if(clusterx$from[w]==highlight_selectx|clusterx$to[w]==highlight_selectx){
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
    i=1
    while(i< length(clusterx$from)){
      if(clusterx$from[i]==highlight_selectx){
        dataConnectx<-append(dataConnectx,clusterx$to[i])
      }
      if(clusterx$to[i]==highlight_selectx){
        dataConnectx<-append(dataConnectx,clusterx$from[i])
      }
      i=i+1
    }
    
  }##correlated genes
  
  ##size and color of nodes  
  ## size  selected gene 5,   hub, alzhheimer, gatekeepers genes 4 , correlated genes x
  ## color selected gene orange, hub gene pink, alzheimer genes black , gatekeepers blue
  vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,5,
                                     ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:10],4,
                                            ifelse(vertex.attributes(gcx)$name%in%alztablex,4,
                                                   ifelse(vertex.attributes(gcx)$name%in%gateTablex,4,
                                                          ifelse(vertex.attributes(gcx)$name%in%dataConnectx,3,1)))))
  
  vertex.attributes(gcx)$color=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,"orange",
                                      ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:10],"pink",
                                             ifelse(vertex.attributes(gcx)$name%in%alztablex,"black",
                                                    ifelse(vertex.attributes(gcx)$name%in%gateTablex,"blue","white"))))
  
  x<-paste("wTO < ",range1,sep = " ")
  y<-paste(range1," < wTO > ", range2,sep = " ")
  z<-paste("wTO > ",range2,sep = " ")
  
  plot.igraph(gcx,vertex.label=NA,edge.width=E(gcx)$width)
  #plot with labels
  #plot.igraph(gcx,vertex.label=ifelse(vertex.attributes(gcx)$name%in%c(hubtablex$genes[1:10],alztablex,gateTablex),V(gcx)$hgncsymbol,NA),vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=1,edge.width=E(gcx)$width)
  
  legend("topleft", legend=c(x,y,z,"Hubgenes","AlzheimerGenes","Gatkeepers",highlight_selectx,"LNC Coding","Protein Coding"),lty=c(1,1,1,0,0,0,0,0,0),pch = c(NA,NA,NA,19,19,19,19,6,19),col = c("red","lightblue","darkgreen","pink","black","blue","orange","grey","grey"),bg=rgb(1,1,0,alpha=0.05),pt.cex=c(2,2,2,2,2,2,2,2,2))
  #legend("topright", legend=c(x,y,z,"Hubgenes","AlzheimerGenes","Gatkeepers",highlight_selectx,"LNC Coding","Protein Coding"),lty=c(1,1,1,0,0,0,0,0,0),pch = c(NA,NA,NA,19,19,19,19,6,19),col = c("red","lightblue","darkgreen","pink","black","blue","orange","grey","grey"),bg=rgb(1,1,0,alpha=0.05),pt.cex=c(2,2,2,2,2,2,2,2,2))
  
}  


#####plot function end

###correlated genes Function++++ STart

correlatedGenesFunc<-function(clusterx,clusterxinfo,highlight_selectx){
  wTO<-0                                
  Gene<-"Gene from"
  GenesTo<-"Gene to"
  dataconnectx<-data.frame(Gene,GenesTo,wTO)
  i=1
  while(i<= length(clusterx$from)){
    if(clusterx$from[i]==highlight_selectx){
      ##fill the columns with the genes and wto score
      dataconnectx<-rbind(dataconnectx,c(clusterx$from[i],clusterx$to[i],clusterx$wTO[i]))
    }
    if(clusterx$to[i]==highlight_selectx){
      ##if selected gene is in the to section add thes to the dataframe
      dataconnectx<-rbind(dataconnectx,c(clusterx$to[i],clusterx$from[i],clusterx$wTO[i]))
    }
    i=i+1
  }
  
  dataconnectx<-dataconnectx[-1,] ##deletes first line that are not genes
  dataconnectx<-dataconnectx[order(dataconnectx$wTO,decreasing=TRUE),]### highest wto score first
  dataconnectx$GeneType<-""
  dataconnectx$description<-""
  dataconnectx$Genecard<-""
  e=1
  while(e<length(dataconnectx$Gene)+1){
    dataconnectx$GeneType[e]   =clusterxinfo$gene_biotype[clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    dataconnectx$Gene[e]       =clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==highlight_selectx]
    dataconnectx$GenesTo[e    ]=clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    dataconnectx$Genecard[e]=createLink(dataconnectx$GenesTo[e])
    e=e+1
  }    #add description and gene name 
  dataconnectx
  
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

###################################################


###degree

#xyz<-aFunction(tabble)
#xyz1<-aFunction(cluster1)
#xyz2<-aFunction(cluster2)
#xyz3<-aFunction(cluster3)




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




# Define UI ----
ui <- fluidPage(
  
  
  titlePanel("The Alzheimer App"),
  
  sidebarLayout(
    sidebarPanel(
      ###
      #tags$a(href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=ECE1", "Click Here!"),
      ####
      
      h4("Select a Gene for each Cluster"),
      h6("Select Tab with Cluster  Network"),
      h6("Choose your Parameters for the Plot"),
      h6("Press   PLOT!   to render the Plot with your Parameters"),
      actionButton("goButton", "PLOT!"),
      
      
      
      selectInput("select1",
                  h4("Select Cluster1 Gene"), 
                  choices = c(completetable1$names), selected = "DDOST",multiple=FALSE),
      selectInput("select2",
                  h4("Select Cluster2 Gene"), 
                  choices = c(completetable2$names), selected = "ARL2-SNX15",multiple=FALSE),
      selectInput("select3",
                  h4("Select Cluster3 Gene"), 
                  choices = c(completetable3$names), selected = "C1QTNF3-AMACR",multiple=FALSE),
      selectInput("select",
                  h4("Select Cluster4 Gene"), 
                  choices = c(completetable4$names), selected = "ECE1",multiple=FALSE),
      
      
      
      sliderInput("range", 
                  label = "Show Correlation in this range:",
                  min = -1, max = 1, value = c(-0.8, 0.8),step = 0.1),
      
      h4("Highlight Positive Related Genes of Selected Gene"),
      checkboxInput("checkboxSelectPositiv", "Positiv", value = TRUE),
      h4("Highlight Negative Related Genes of Selected Gene"),
      checkboxInput("checkboxSelectNegativ", "Negativ", value = TRUE),
      
      h4("Show HUB Genes"),
      checkboxInput("checkboxHub", "Yes", value = FALSE),
      h4("Show Correlated Genes for Selected Gene"),
      checkboxInput("checkboxConnect","Yes",value = FALSE),
      
      h4("Alzheimer Genes"),
      checkboxInput("checkboxAlzheimer","Yes",value = FALSE),
      h4("Gatekeepers"),
      checkboxInput("GateKeep","Yes",value=FALSE)
      
      
   
      ########
      
    ),
    
    mainPanel(
  
      
      tabsetPanel(type = "pills",
                  tabPanel("Gene Information",tabsetPanel(
                      tabPanel("Cluster 1 Gene",dataTableOutput("gene1desc")),
                      tabPanel("Cluster2",dataTableOutput("gene2desc")),
                      tabPanel("Cluster3",dataTableOutput("gene3desc")),
                      tabPanel("Cluster4",dataTableOutput("gene4desc")))),
                  
                  tabPanel("AlzheimerGenes",dataTableOutput("alzheimerGenesX")),
                  
                  #
                  #tabPanel("Test",dataTableOutput("test")),
                  #tabPanel("Test2",dataTableOutput("test2")),
                  #
                  
                  tabPanel("Cluster 1",tabsetPanel(
                                tabPanel("HubGenes",dataTableOutput("hubGeneDesciptions1")),
                                tabPanel("Network",plotOutput("graphic1")),
                                tabPanel("Correlated Genes",dataTableOutput("connectedGenes1")),
                                tabPanel("Gatekeeper",dataTableOutput("Gatekeeper1"))
                  )),
                  tabPanel("Cluster 2",tabsetPanel(
                                tabPanel("HubGenes",dataTableOutput("hubGeneDesciptions2")),
                                tabPanel("Network",plotOutput("graphic2")),
                                tabPanel("Correlated Genes",dataTableOutput("connectedGenes2")),
                                tabPanel("Gatekeeper",dataTableOutput("Gatekeeper2"))
                  )),
                  tabPanel("Cluster 3",tabsetPanel(
                                tabPanel("HubGenes",dataTableOutput("hubGeneDesciptions3")),
                                tabPanel("Network",plotOutput("graphic3")),
                                tabPanel("Correlated Genes",dataTableOutput("connectedGenes3")),
                                tabPanel("Gatekeeper",dataTableOutput("Gatekeeper3"))
                  )),
                  tabPanel("Cluster 4",tabsetPanel(
                                tabPanel("HubGenes",dataTableOutput("hubGeneDesciptions4")),
                                tabPanel("Network",plotOutput("graphic4")),
                                tabPanel("Correlated Genes",dataTableOutput("connectedGenes")),
                                tabPanel("Gatekeeper",dataTableOutput("Gatekeeper4"))
                  ))
                  
               )
      )
  )
)

# Define server logic ----
options(shiny.maxRequestSize=90*1024^2)
server <- function(input, output) {
  ##all
  
  ###
  output$alzheimerGenesX<-renderDataTable({
    for(x in 1:length(alzheimergenes$hgnc_symbol))alzheimergenes$genecard[x]<-createLink(alzheimergenes$hgnc_symbol[x])
    alzheimergenes
  }, escape = FALSE)
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
  
  
  output$gene4desc<-renderDataTable({
    highlight_select<-cluster4info$ensembl_gene_id[cluster4info$hgnc_symbol==input$select]
    infoTable<-data.frame(ens="",name="",desc="")
    infoTable$ens<-highlight_select
    infoTable$name<-cluster4info$hgnc_symbol[cluster4info$ensembl_gene_id==highlight_select]
    infoTable$desc<-cluster4info$description[cluster4info$ensembl_gene_id==highlight_select]
    infoTable$genecard<- createLink(input$select)
    
    return(infoTable)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  output$gene1desc<-renderDataTable({
    highlight_select1<-cluster1info$ensembl_gene_id[cluster1info$hgnc_symbol==input$select1]
    infoTable1<-data.frame(ens="",name="",desc="")
    infoTable1$ens<-highlight_select1
    infoTable1$name<-cluster1info$hgnc_symbol[cluster1info$ensembl_gene_id==highlight_select1]
    infoTable1$desc<-cluster1info$description[cluster1info$ensembl_gene_id==highlight_select1]
    infoTable1$genecard<- createLink(input$select1)
    infoTable1
    ##put gene description4 in table  cluster1
  }, escape = FALSE)
  output$gene2desc<-renderDataTable({
    
    highlight_select2<-cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2]
    
    infoTable2<-data.frame(ens="",name="",desc="")
    infoTable2$ens<-highlight_select2
    infoTable2$name<-cluster2info$hgnc_symbol[cluster2info$ensembl_gene_id==highlight_select2]
    infoTable2$desc<-cluster2info$description[cluster2info$ensembl_gene_id==highlight_select2]
    infoTable2$genecard<- createLink(input$select2)
    infoTable2
    ##put gene description4 in table  cluster2
  }, escape = FALSE)
  output$gene3desc<-renderDataTable({
    highlight_select3<-cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3]
    infoTable3<-data.frame(ens="",name="",desc="")
    infoTable3$ens<-highlight_select3
    infoTable3$name<-cluster3info$hgnc_symbol[cluster3info$ensembl_gene_id==highlight_select3]
    infoTable3$desc<-cluster3info$description[cluster3info$ensembl_gene_id==highlight_select3]
    infoTable3$genecard<- createLink(input$select3)
    infoTable3
    ##put gene description4 in table  cluster3
  }, escape = FALSE)
 
  
  
  
  
  ###cluster4
  ##top 10 genes with most degree plus information Cluster 4
  output$hubGeneDesciptions4<-renderDataTable({
    hubtable4<-completetable4[order(completetable4$degree,decreasing = TRUE ),]
    for(x in 1:10)hubtable4$genecard[x]<-createLink(hubtable4$genes[x])
    hubtable4<-subset(hubtable4,select = -1)
    head(hubtable4,10)
   }, escape = FALSE)
  
  ##cluster 1
  output$hubGeneDesciptions1<-renderDataTable({
    hubtable1<-completetable1[order(completetable1$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable1$genecard[x]<-createLink(hubtable1$genes[x])
    hubtable1<-subset(hubtable1,select = -1)
    head(hubtable1,10)
  }, escape = FALSE)
  ##cluster 2
  output$hubGeneDesciptions2<-renderDataTable({
    hubtable2<-completetable2[order(completetable2$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable2$genecard[x]<-createLink(hubtable2$genes[x])
    hubtable2<-subset(hubtable2,select = -1)
    head(hubtable2,10)
  }, escape = FALSE)
  ##cluster 3
  output$hubGeneDesciptions3<-renderDataTable({
    hubtable3<-completetable3[order(completetable3$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable3$genecard[x]<-createLink(hubtable3$genes[x])
    hubtable3<-subset(hubtable3,select = -1)
    head(hubtable3,10)
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

    
      dataconnect4<-correlatedGenesFunc(tabble,cluster4info,highlight_select)
    
  }, escape = FALSE)

  
  ###cluster 1 correlated genes
  output$connectedGenes1 <-renderDataTable({
    highlight_select1<- cluster1info$ensembl_gene_id[cluster1info$hgnc_symbol==input$select1]
    
    
      dataconnect1<-correlatedGenesFunc(cluster1,cluster1info,highlight_select1)
    
  }, escape = FALSE)
  ###cluster 2 correlated genes
  output$connectedGenes2 <-renderDataTable({
    highlight_select2<- cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2]
    
    
      dataconnect2<-correlatedGenesFunc(cluster2,cluster2info,highlight_select2)
    
  }, escape = FALSE)
  
  ###cluster 3 correlated genes
  output$connectedGenes3 <-renderDataTable({
    highlight_select3<- cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3]
    
   
      dataconnect3<-correlatedGenesFunc(cluster3,cluster3info,highlight_select3)
    
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
    
    plot1<-isolate(plotxfunc(gc1,highlight_select1,input$checkboxSelectPositiv,cluster1,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable1,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network1,input$checkboxConnect))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)

  },height = 1000,width = 1000)
  
  ##graph2 
  output$graphic2 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select2=isolate(cluster2info$ensembl_gene_id[cluster2info$hgnc_symbol==input$select2])
    
    plot2<-isolate(plotxfunc(gc2,highlight_select2,input$checkboxSelectPositiv,cluster2,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable2,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network2,input$checkboxConnect))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
  
  },height = 1000,width = 1000)
  
  ####graph3
  output$graphic3 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    
    highlight_select3=isolate(cluster3info$ensembl_gene_id[cluster3info$hgnc_symbol==input$select3])
    
    plot3<-isolate(plotxfunc(gc3,highlight_select3,input$checkboxSelectPositiv,cluster3,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable3,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network3,input$checkboxConnect))
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
    
    
    plot4<-isolate(plotxfunc(g,highlight_select,input$checkboxSelectPositiv,tabble,input$range[2],input$checkboxSelectNegativ,input$range[1],input$checkboxHub,completetable4,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,network4,input$checkboxConnect) )
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1000,width = 1000)
  
  #############test
  
 
  
  ##########test

}

# Run the app ----
shinyApp(ui = ui, server = server,options =list (launch.browser =TRUE))

