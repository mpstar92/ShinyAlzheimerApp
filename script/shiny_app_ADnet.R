#library(XML)
#library(rentrez)

library(shiny, warn.conflicts = FALSE)
library(shinydashboard, warn.conflicts = FALSE) 

suppressWarnings({
  library(igraph, warn.conflicts = FALSE)
  library(qgraph)
  library(DT, warn.conflicts = FALSE)
  library("data.table", warn.conflicts = FALSE)
  library(dplyr, warn.conflicts = FALSE)
  library(ggraph)
})



##data input

alzheimergenes_main<-read.csv(file="../data/info/AD_key_genes.txt", sep = '\t')
alzheimergenes<-read.csv(file="../data/info/AD_kegg_genes.txt", sep = '\t')

cor_genes <- read.csv(file="../data/shiny_data/data/AD_core_hubs.csv", sep = ',')

AD_all_genes <- read.csv(file="../data/shiny_data/data/AD_all_genes.csv", sep = ',')


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
    
    

    
  }##correlated genes
  hub<-ceiling(length(V(gcx))*0.01)
  ##size and color of nodes  
  ## size  selected gene 5,   hub, alzhheimer, gatekeepers genes 4 , correlated genes x
  ## color selected gene orange, hub gene pink, alzheimer genes black , gatekeepers blue
 # vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,3,
  #                                   ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:hub],2,
   #                                         ifelse(vertex.attributes(gcx)$name%in%alztablex,2,
    #                                               ifelse(vertex.attributes(gcx)$name%in%gateTablex,2,
     #                                                     ifelse(vertex.attributes(gcx)$name%in%dataConnectx,1.5,1)))))
  
  
  
  vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%cor_genes$ensembl_gene_id,cor_genes$Degree_centrality*0.01,1)
  
  
  vertex.attributes(gcx)$color=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,"orange",
                                      ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:hub],"pink",
                                             ifelse(vertex.attributes(gcx)$name%in%alztablex,"#c51b8a",
                                                    ifelse(vertex.attributes(gcx)$name%in%gateTablex,"lightblue",
                                                        ifelse(vertex.attributes(gcx)$name%in%dataConnectx,"lightyellow","white")))))
  
  

  
  #####
  x<-paste("wTO < ",range1,sep = " ")
  z<-paste("wTO > ",range2,sep = "  ")
  
  y2<-paste(completetablex$names[completetablex$genes==highlight_selectx])
  
  
  ##################
  
  #####
  node_degree <- degree(gcx, mode = "all")
  
  # Set the size of the nodes proportional to their degree
 # V(gcx)$size <- node_degree * 0.015

  e <- get.edgelist(gcx,names=FALSE)

  l <- qgraph.layout.fruchtermanreingold(e,weights = edge.attributes(gcx)$wTO,vcount=vcount(gcx),
                                         area=8*(vcount(gcx)^2),repulse.rad=0.01*(vcount(gcx)^3.1))

  plot(gcx,vertex.label=NA,edge.width=E(gcx)$width,layout=l )#,theme(legend.key.size = unit(1.5, "cm")))
 
  #plot(gcx,vertex.label=NA,edge.width=E(gcx)$width,vertex.size = V(gcx)$core,layout=l )#,theme(legend.key.size = unit(1.5, "cm")))
  
  #########
  
  #plot.igraph(gcx,vertex.label=NA,edge.width=E(gcx)$width,layout = coords_jittered )
  #plot with labels
  #plot.igraph(gcx,vertex.label=ifelse(vertex.attributes(gcx)$name%in%c(hubtablex$genes[1:10],alztablex,gateTablex),V(gcx)$hgncsymbol,NA),vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=1,edge.width=E(gcx)$width)
  
  legend("topleft",
         text.width = 0.05,
         bty = "n",
         legend=c(x,z,"Hub Genes","AD Pathway Genes(KEGG)","Gatekeeper Genes",y2,"lnc RNA coding Genes","Protein coding Genes"),
        
         lty=c(1,1,0,0,0,0,0,0),
         pch = c(NA,NA,19,19,19,19,6,19),
         col = c("red","darkgreen","pink","#c51b8a","lightblue","orange","grey","grey"),
         bg=rgb(1,1,0,alpha=0.05),
         pt.cex=c(2,2,2,2,2,2,2,2,2)
         )
  #legend("topright", legend=c(x,z,"Hubgenes","AlzheimerGenes","Gatkeepers",highlight_selectx,"LNC Coding","Protein Coding"),lty=c(1,1,1,0,0,0,0,0,0),pch = c(NA,NA,NA,19,19,19,19,6,19),col = c("red","lightblue","darkgreen","pink","black","blue","orange","grey","grey"),bg=rgb(1,1,0,alpha=0.05),pt.cex=c(2,2,2,2,2,2,2,2,2))
  
}  

plotxfunc2<-function(gcx,highlight_selectx,range2,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect,checkboxwTO){   #bool1pos ,2 neg ,x hub ,4 alz , 5 gate, 6 connect..
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
  
  
  
  
  #E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$color<-adjustcolor("lightblue",2) #all links from selected gene blue
  #E(gcx)[from(highlight_selectx )|to(highlight_selectx )]$width<-0.5
  
  
  ## postiv correlated genes darkgreen
  
  some<-as_edgelist(gcx)[,1]  #nodes from
  some2<-as_edgelist(gcx)[,2]#nodes to
  some3<-E(gcx)$wTO       # wto
  
  
  
  if(checkboxwTO==TRUE){
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
  if(checkboxwTO==TRUE){
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
  
  
  
  
  
  if(checkboxHub==TRUE){#corgenes
    
    #hubtablex<-cor_genes$ensembl_gene_id#completetablex[order(completetablex$degree,decreasing = TRUE ),]
    hubtablex <- completetablex$genes[completetablex$Core == "core"]
    
  }
  
  if(checkboxAlzheimer==TRUE){##alzheimer
    w=1
    while(w<=length(alzheimergenes$hgnc_symbol)){
      alztablex<-append(alztablex,alzheimergenes$ensembl_gene_id[w])
      w=w+1
    }
    
  }
  
  if(checkboxConnect==TRUE){##correlated genes
    dataConnectx<-append(dataConnectx,ends(gcx, E(gcx)[from(highlight_selectx)])[,1])  #"Gene from"
    dataConnectx<-append(dataConnectx,ends(gcx, E(gcx)[from(highlight_selectx)])[,2] ) #"Gene to"
    
   
    
    
  }##correlated genes
  hub<-ceiling(length(V(gcx))*0.01)
  ##size and color of nodes  
  ## size  selected gene 5,   hub, alzhheimer, gatekeepers genes 4 , correlated genes x
  ## color selected gene orange, hub gene pink, alzheimer genes black , gatekeepers blue
  # vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,3,
  #                                   ifelse(vertex.attributes(gcx)$name%in%hubtablex$genes[1:hub],2,
  #                                         ifelse(vertex.attributes(gcx)$name%in%alztablex,2,
  #                                               ifelse(vertex.attributes(gcx)$name%in%gateTablex,2,
  #                                                     ifelse(vertex.attributes(gcx)$name%in%dataConnectx,1.5,1)))))
  
  
  
  #vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%cor_genes$ensembl_gene_id,cor_genes$Degree_centrality*0.01,1)
  #vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%cor_genes$ensembl_gene_id,8,1)
  

  vertex.attributes(gcx)$size=ifelse(vertex.attributes(gcx)$name%in%highlight_selectx,8,
                                     ifelse(vertex.attributes(gcx)$name%in%alztablex,4,
                                     ifelse(vertex.attributes(gcx)$name%in%hubtablex,11,
                                     ifelse(vertex.attributes(gcx)$name%in%dataConnectx,3,1))))
  
  

  
  ##################
  
  #####
  node_degree <- degree(gcx, mode = "all")
  
  # Set the size of the nodes proportional to their degree
  # V(gcx)$size <- node_degree * 0.015
  
  e <- get.edgelist(gcx,names=FALSE)
  
  
  # Define a color palette
  group_colors <- c("Mitochondrial energy production" = "#A5C8E4", 
                    "Apoptosis" = "#C0ECCC",
                    "intracellular trafficking" = "#F4CDA6",
                    "Synapse signaling" = "#F9F0C1",
                    "Others" = "#F2E7DC",
                    "No enriched GO" = "#D9D9D9",
                    "Cytoskeleton organization" = "#78DAFA",
                    "gene expression" = "#D0B5EA",
                    "Neurogenesis" = "#F6A8A6",
                    "cellular homeostasis" = "#F67CFA"
  )
  
  assign_pie_colors <- function(properties, group_colors) {
    # Split properties into a vector
    prop_vector <- strsplit(properties, ",")[[1]]
    
    # Create a color pie based on presence of properties
    pie_colors <- sapply(prop_vector, function(prop) group_colors[prop])
    
    # Return pie chart specification
    return(pie_colors)
  }
  
  
  
  
  vertex_pie_colors <- lapply(V(gcx)$GOcategory, assign_pie_colors, group_colors)
  

  
  # Count the number of properties for each node (for pie chart segments)
  vertex_pie_values <- lapply(V(gcx)$GOcategory, function(props) {
    prop_vector <- strsplit(props, ",")[[1]]
    pie_values <- rep(1, length(prop_vector))  # Equal segments for each property
    return(pie_values)
  })
  
  #####
  
 ##edge color on ad ctrl specificity
  
  group_colors_edges<-c(	"g.AD_TCX" ="#fb8072",
                         "g.Ctrl_TCX" ="#8dd3c7",
                         "a" ="#D9D9D9"
  )
  
   E(gcx)$color <- group_colors_edges[E(gcx)$Phitilde]
  
   
   
   
   
   
   ####frame color DEG
   
   V(gcx)$frame.color[V(gcx)$name%in%alztablex] <- "violet"
   
    frame_colors <- c("Up-regulated" = "red", "Down-regulated" = "blue","non DEG" ="grey")
    V(gcx)$frame.color <- frame_colors[V(gcx)$DEG]
  
  
  
  ###########
    
   # V(gcx)$frame.color[V(gcx)$name == highlight_selectx] <- "black"
    
    V(gcx)$frame.color[V(gcx)$name == highlight_selectx] <- ifelse(
      V(gcx)$frame.color[V(gcx)$name == highlight_selectx] == "red", 
      "darkred", 
      ifelse(
        V(gcx)$frame.color[V(gcx)$name == highlight_selectx] == "blue", 
        "purple", 
        "black"
      )
    )
    
  #####
    l <- qgraph.layout.fruchtermanreingold(e,weights = edge.attributes(gcx)$wTO,vcount=vcount(gcx),
                                           area=18*(vcount(gcx)^3),repulse.rad=0.01*(vcount(gcx)^3.1))
    
  
  ###
  
  plot(gcx,
       vertex.shape = "pie",                # Set vertex shape to pie
       vertex.pie = vertex_pie_values,      # Assign pie segments
       vertex.pie.color = vertex_pie_colors,
       vertex.label=NA,
       vertex.frame.color = V(gcx)$frame.color,
       vertex.frame.width = 7,
       edge.width=E(gcx)$width,
       edge.color = E(gcx)$color,
       layout=l )#,theme(legend.key.size = unit(1.5, "cm")))
  
  #plot(gcx,vertex.label=NA,edge.width=E(gcx)$width,vertex.size = V(gcx)$core,layout=l )#,theme(legend.key.size = unit(1.5, "cm")))
  
  #########
  
  #plot.igraph(gcx,vertex.label=NA,edge.width=E(gcx)$width,layout = coords_jittered )
  #plot with labels
  #plot.igraph(gcx,vertex.label=ifelse(vertex.attributes(gcx)$name%in%c(hubtablex$genes[1:10],alztablex,gateTablex),V(gcx)$hgncsymbol,NA),vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=1,edge.width=E(gcx)$width)
  #####
  x<-paste("Up-regulated genes")
  z<-paste("Down-regulated genes")
  
  y2<-paste(completetablex$names[completetablex$genes==highlight_selectx])
  
  
  legend("topleft",
         #inset=c(-0,0),
         text.width = 0.05,
         bty = "n",
         legend=c(x,z,"common genes",y2,
                  "",
                  "Mitochondrial energy production",
                  "Apoptosis",
                  "Intracellular trafficking",
                  "Synapse signaling",
                  "Others",
                  "No enriched GO",
                  "Cytoskeleton organization",
                  "Gene expression",
                  "Neurogenesis",
                  "Cellular homeostasis",
                  "AD-specific links",
                  "Control-specific links",
                  "Common links",
                  "lncRNA coding Genes","Protein coding Genes"),
         
         lty=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0),
         pch = c(21,   #
                 21,
                 21,
                 21,  #filled circle
                 19,
                 19,19,19,19,19,19,19,19,19,19,NA,NA,NA,6,19),
         col = c("red","blue","grey","black","white",
                 "#A5C8E4",
                 "#C0ECCC",
                 "#F4CDA6",
                 "#F9F0C1",
                 "#F2E7DC",
                 "#D9D9D9",
                 "#78DAFA",
                 "#D0B5EA",
                 "#F6A8A6",
                 "#F67CFA",
                 "#fb8072","#8dd3c7","#D9D9D9",
                 "grey","grey"),
         #bg=rgb(1,1,0,alpha=0.05),
         bg=c("white","white","white","white","white",
              "#A5C8E4",
              "#C0ECCC",
              "#F4CDA6",
              "#F9F0C1",
              "#F2E7DC",
              "#D9D9D9",
              "#78DAFA",
              "#D0B5EA",
              "#F6A8A6",
              "#F67CFA",
              "#fb8072","#8dd3c7","#D9D9D9",
              "grey","grey"),
         lwd= c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
         pt.cex=c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
  )
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
  dataconnectx$kegg<-""
  dataconnectx$GO_category<-""
  dataconnectx$Core<-""
  dataconnectx$DEG<-""
  dataconnectx$description<-""
  dataconnectx$Genecard<-""
  e=1
  while(e<length(dataconnectx$Gene)+1){
    dataconnectx$GeneType[e]   =clusterxinfo$gene_biotype[clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    if(highlight_selectx == dataconnectx$Gene[e]){ 
      dataconnectx$kegg[e] = clusterxinfo$kegg [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$GO_category[e] = clusterxinfo$GO_category [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$Core[e] = clusterxinfo$Core [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$DEG[e] = clusterxinfo$DEG [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      
      dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$Genecard[e]   =createLink(dataconnectx$GenesTo[e])}
    else{  
      dataconnectx$kegg[e]       =clusterxinfo$kegg        [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$GO_category[e] =clusterxinfo$GO_category [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$Core[e] =clusterxinfo$Core [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$DEG[e] =clusterxinfo$DEG [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      
      dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$Genecard[e]=createLink(dataconnectx$Gene[e])}
      dataconnectx$Gene[e]       =clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e]    ]
      dataconnectx$GenesTo[e    ]=clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    
    e=e+1
  }    #add description and gene name 
  
  
  dataconnectx[,-2]
  
  #output dataconnect table with all information
}

#### new function
correlatedGenesFunc2<-function(graph,clusterxinfo,highlight_selectx){
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
  dataconnectx$kegg<-""
  dataconnectx$Parent_Term<-""
  dataconnectx$description<-""
  dataconnectx$Genecard<-""

  e=1
  while(e<length(dataconnectx$Gene)+1){
    dataconnectx$GeneType[e]   =clusterxinfo$gene_biotype[clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    if(highlight_selectx == dataconnectx$Gene[e]){ 
      dataconnectx$kegg[e] = clusterxinfo$kegg [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$Parent_Term[e] = clusterxinfo$Manual_added_parentTerm [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
     # dataconnectx$Term[e] = clusterxinfo$term [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      #dataconnectx$graph_cluster= clusterxinfo$graph_cluster [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
      dataconnectx$Genecard[e]   =createLink(dataconnectx$GenesTo[e])}
    
    else{      

      dataconnectx$kegg[e]       =clusterxinfo$kegg        [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$Parent_Term[e] =clusterxinfo$Manual_added_parentTerm [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      #dataconnectx$Term[e] =clusterxinfo$term [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
     # dataconnectx$graph_cluster[e] =clusterxinfo$graph_cluster [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$description[e]=clusterxinfo$description [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e] ]
      dataconnectx$Genecard[e]=createLink(dataconnectx$Gene[e])}
      dataconnectx$Gene[e]       =clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$Gene[e]    ]
      dataconnectx$GenesTo[e    ]=clusterxinfo$hgnc_symbol [clusterxinfo$ensembl_gene_id==dataconnectx$GenesTo[e] ]
    
      
    e=e+1
  }    #add description and gene name 
  
  
  dataconnectx[,-2]
  
  #output dataconnect table with all information
}

#dataconnect5<-correlatedGenesFunc2(g_AD1,clusterAD1,highlight_select5)
#############

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




###AD cluster 1---- 
clusterAD1<- read.table("../data/shiny_data/data/annotated/AD_cluster1_info.txt", fill = TRUE ,header=T, sep = "," )
g_AD1<- read_graph("../data/shiny_data/data/graphs/AD_cluster1_graph.gml",format = "gml")
g_AD1<-delete_vertex_attr(g_AD1,"id")
g_AD1_size_links <-gsize(g_AD1)
g_AD1_size_nodes <-length(V(g_AD1))
########

completetableAD1<-data.frame(genes=vertex.attributes(g_AD1)$name,degree=degree(g_AD1),names=vertex.attributes(g_AD1)$hgncsymbol,genebiotype=vertex.attributes(g_AD1)$genebiotype,kegg=vertex.attributes(g_AD1)$kegg,Core = vertex.attributes(g_AD1)$core,DEG = vertex.attributes(g_AD1)$DEG,description="")
completetableAD1<-completetableAD1[order(completetableAD1$genes), ]
completetableAD1$description<-clusterAD1$description
completetableAD1<-completetableAD1[order(completetableAD1$degree), ]
networkAD1<-betweenFunction(g_AD1)



###AD cluster 2---- 
clusterAD2<- read.table("../data/shiny_data/data/annotated/AD_cluster2_info.txt", fill = TRUE ,header=T, sep = "," )
g_AD2<- read_graph("../data/shiny_data/data/graphs/AD_cluster2_graph.gml",format = "gml")
g_AD2<-delete_vertex_attr(g_AD2,"id")
g_AD2_size_links <-gsize(g_AD2)
g_AD2_size_nodes <-length(V(g_AD2))
########

completetableAD2<-data.frame(genes=vertex.attributes(g_AD2)$name,degree=degree(g_AD2),names=vertex.attributes(g_AD2)$hgncsymbol,genebiotype=vertex.attributes(g_AD2)$genebiotype,kegg=vertex.attributes(g_AD2)$kegg,Core = vertex.attributes(g_AD2)$core,DEG = vertex.attributes(g_AD2)$DEG,description="")
completetableAD2<-completetableAD2[order(completetableAD2$genes), ]
completetableAD2$description<-clusterAD2$description
completetableAD2<-completetableAD2[order(completetableAD2$degree), ]
networkAD2<-betweenFunction(g_AD2)



###AD cluster 3---- 
clusterAD3<- read.table("../data/shiny_data/data/annotated/AD_cluster3_info.txt", fill = TRUE ,header=T, sep = "," )
g_AD3<- read_graph("../data/shiny_data/data/graphs/AD_cluster3_graph.gml",format = "gml")
g_AD3<-delete_vertex_attr(g_AD3,"id")
g_AD3_size_links <-gsize(g_AD3)
g_AD3_size_nodes <-length(V(g_AD3))
########

completetableAD3<-data.frame(genes=vertex.attributes(g_AD3)$name,degree=degree(g_AD3),names=vertex.attributes(g_AD3)$hgncsymbol,genebiotype=vertex.attributes(g_AD3)$genebiotype,kegg=vertex.attributes(g_AD3)$kegg,Core = vertex.attributes(g_AD3)$core,DEG = vertex.attributes(g_AD3)$DEG,description="")
completetableAD3<-completetableAD3[order(completetableAD3$genes), ]
completetableAD3$description<-clusterAD3$description
completetableAD3<-completetableAD3[order(completetableAD3$degree), ]
networkAD3<-betweenFunction(g_AD3)




#### clust control 1
cluster_con_1<- read.table("../data/shiny_data/data/annotated/ctrl_cluster1_info.txt", fill = TRUE ,header=T, sep = "," )
g_con1<- read_graph("../data/shiny_data/data/graphs/control_cluster1_graph.gml",format = "gml")
g_con1<-delete_vertex_attr(g_con1,"id")
g_con1_size_links <-gsize(g_con1)
g_con1_size_nodes <-length(V(g_con1))
########

completetable_con1<-data.frame(genes=vertex.attributes(g_con1)$name,degree=degree(g_con1),names=vertex.attributes(g_con1)$hgncsymbol,genebiotype=vertex.attributes(g_con1)$genebiotype,kegg=vertex.attributes(g_con1)$kegg,Core = vertex.attributes(g_con1)$core,DEG = vertex.attributes(g_con1)$DEG,description="")
completetable_con1<-completetable_con1[order(completetable_con1$genes), ]
completetable_con1$description<-cluster_con_1$description
completetable_con1<-completetable_con1[order(completetable_con1$degree), ]
networkcon1<-betweenFunction(g_con1)


#### clust control 2 ---cluster 4
cluster_con_2<- read.table("../data/shiny_data/data/annotated/ctrl_cluster2_info.txt", fill = TRUE ,header=T, sep = "," )
g_con2<- read_graph("../data/shiny_data/data/graphs/control_cluster2_graph.gml",format = "gml")
g_con2<-delete_vertex_attr(g_con2,"id")
g_con2_size_links <-gsize(g_con2)
g_con2_size_nodes <-length(V(g_con2))
########

completetable_con2<-data.frame(genes=vertex.attributes(g_con2)$name,degree=degree(g_con2),names=vertex.attributes(g_con2)$hgncsymbol,genebiotype=vertex.attributes(g_con2)$genebiotype,kegg=vertex.attributes(g_con2)$kegg,Core = vertex.attributes(g_con2)$core,DEG = vertex.attributes(g_con2)$DEG,description="")
completetable_con2<-completetable_con2[order(completetable_con2$genes), ]
completetable_con2$description<-cluster_con_2$description
completetable_con2<-completetable_con2[order(completetable_con2$degree), ]
networkcon2<-betweenFunction(g_con2)






###Control cluster 3 ---cluster 1

cluster_con_3<- read.table("../data/shiny_data/data/annotated/ctrl_cluster3_info.txt", fill = TRUE ,header=T, sep = "," )
g_con3<- read_graph("../data/shiny_data/data/graphs/control_cluster3_graph.gml",format = "gml")
g_con3<-delete_vertex_attr(g_con3,"id")
g_con3_size_links <-gsize(g_con3)
g_con3_size_nodes <-length(V(g_con3))
########

completetable_con3<-data.frame(genes=vertex.attributes(g_con3)$name,degree=degree(g_con3),names=vertex.attributes(g_con3)$hgncsymbol,genebiotype=vertex.attributes(g_con3)$genebiotype,kegg=vertex.attributes(g_con3)$kegg,Core = vertex.attributes(g_con3)$core,DEG = vertex.attributes(g_con3)$DEG,description="")
completetable_con3<-completetable_con3[order(completetable_con3$genes), ]
completetable_con3$description<-cluster_con_3$description
completetable_con3<-completetable_con3[order(completetable_con3$degree), ]
networkcon3<-betweenFunction(g_con3)
#########

#control cluster 4

cluster_con_4<- read.table("../data/shiny_data/data/annotated/ctrl_cluster4_info.txt", fill = TRUE ,header=T, sep = "," )
g_con4<- read_graph("../data/shiny_data/data/graphs/control_cluster4_graph.gml",format = "gml")
g_con4<-delete_vertex_attr(g_con4,"id")
g_con4_size_links <-gsize(g_con4)
g_con4_size_nodes <-length(V(g_con4))
########

completetable_con4<-data.frame(genes=vertex.attributes(g_con4)$name,degree=degree(g_con4),names=vertex.attributes(g_con4)$hgncsymbol,genebiotype=vertex.attributes(g_con4)$genebiotype,kegg=vertex.attributes(g_con4)$kegg,Core = vertex.attributes(g_con4)$core,DEG = vertex.attributes(g_con4)$DEG,description="")
completetable_con4<-completetable_con4[order(completetable_con4$genes), ]
completetable_con4$description<-cluster_con_4$description
completetable_con4<-completetable_con4[order(completetable_con4$degree), ]
networkcon4<-betweenFunction(g_con4)
##########data input end

######


gene_list_sort<-completetableAD1[order(completetableAD1$names),]
gene_list_ad_1_names<-gene_list_sort$names

gene_list_sort<-completetableAD2[order(completetableAD2$names),]
gene_list_ad_2_names<-gene_list_sort$names

gene_list_sort<-completetableAD3[order(completetableAD3$names),]
gene_list_ad_3_names<-gene_list_sort$names

gene_list_sort<-completetable_con1[order(completetable_con1$names),]
gene_list_ctrl_1_names<-gene_list_sort$names

gene_list_sort<-completetable_con2[order(completetable_con2$names),]
gene_list_ctrl_2_names<-gene_list_sort$names

gene_list_sort<-completetable_con3[order(completetable_con3$names),]
gene_list_ctrl_3_names<-gene_list_sort$names

gene_list_sort<-completetable_con4[order(completetable_con4$names),]
gene_list_ctrl_4_names<-gene_list_sort$names


#####

hubbie<-"The Hub Genes are the top 1 % Genes with the most links to other Genes"
gaties<-"Gatekeeper Genes are the top 1% of Genes with the highest betweeness centrality "
gaties2<-"the betweeness centrality is the number of shortest paths that pass throgh the vertex"


###################################

# Define UI ----
ui <- fluidPage(
  
  
  titlePanel("ADNet: Interactive Visualization of TCX AD/Control Consensus Networks"),
  
  sidebarLayout(
    sidebarPanel(
      ###
      #tags$a(href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=ECE1", "Click Here!"),
      ####
      h5(HTML("<b> How to use the app:</b>")),
      h5("1- Select a cluster from the top tabs."),
      h5("2- Choose node attributes from the boxes below."),
      h5("3- Pick a gene from the list."),
      h5("4. Click PLOT! to generate the network based on your selections."),
      h3(""),
      h4(HTML("<b>Node attributes:</b>")),
      checkboxInput("checkboxAlzheimer","AD pathway in KEGG",value = FALSE),
      checkboxInput("checkboxHub", "Core genes", value = TRUE),
      
      
      h4(HTML("<b>Select gene of interest from the gene lists below:</b>")),
      
      selectInput("select5",
                  h4("AD cluster 1"), 
                  choices = c(gene_list_ad_1_names), selected = "AAGAB",multiple=FALSE),
      selectInput("select2",
                  h4("AD cluster 2"), 
                  choices = c(gene_list_ad_2_names), selected = "AAK1",multiple=FALSE),
      selectInput("select7",
                  h4("AD cluster 3"), 
                  #change chioces
                  choices = c(gene_list_ad_3_names), selected = "AASDHPPT",multiple=FALSE),
      selectInput("select3",
                  h4("Control cluster 1"), 
                  choices = c(gene_list_ctrl_1_names), selected = "AAGAB",multiple=FALSE),
      selectInput("select",
                  h4("Control cluster 2"), 
                  choices = c(gene_list_ctrl_2_names), selected = "AAAS",multiple=FALSE),
      selectInput("select1",
                  h4("Control cluster 3"), 
                  choices = c(gene_list_ctrl_3_names), selected = "AAMDC",multiple=FALSE),
      selectInput("select8",
                  h4("Control cluster 4"), 
                  #change choices
                  choices = c(gene_list_ctrl_4_names), selected = "AIF1L",multiple=FALSE),
      
      ###########################
      
      ########################
      checkboxInput("checkboxConnect","Highlight the first neighbors of the selected gene",value = FALSE),
      checkboxInput("checkboxwTO","show links of selcted gene with wTO range",value = FALSE),
      
      
      sliderInput("range", 
                  label = "Set wTO to show values above a range",
                  min = -1, max = 1, value = c(-0.8, 0.8),step = 0.1),
      actionButton("goButton", "PLOT!"),
      
      
      
      ########
      width = 2
      ),
    
    mainPanel(
      
      
                             
      tabsetPanel(type = "pills",    #1
                  tabPanel("AD Network",tabsetPanel(
                    tabPanel("Network", img(src = "AD_complet_network.png", height="80%", width="80%", align="left") ) , 
                    tabPanel("Network Information",
                        h4("The AD network has 3 clusters and the control network has 4 clusters"),br(),
                        h4("AD Cluster 1 has "),h4(g_AD1_size_nodes,"  Genes and "),h4(g_AD1_size_links, " links"),br(),
                        h4("AD Cluster 2 has "),h4(g_AD2_size_nodes,"  Genes and "),h4(g_AD2_size_links, " links"),br(),
                        h4("AD Cluster 3 has "),h4(g_AD3_size_nodes,"  Genes and "),h4(g_AD3_size_links, " links"),br(),
                      
                        h4("Control Cluster 1 has "),h4(g_con1_size_nodes,"  Genes and "),h4(g_con1_size_links, " links"),br(),
                        h4("Control Cluster 2 has "),h4(g_con2_size_nodes,"  Genes and "),h4(g_con2_size_links, " links"),br(),
                        h4("Control Cluster 3 has "),h4(g_con3_size_nodes,"  Genes and "),h4(g_con3_size_links, " links"),br(),
                        h4("Control Cluster 4 has "),h4(g_con4_size_nodes,"  Genes and "),h4(g_con4_size_links, " links")
                        ),
                        
                    
                    tabPanel("AD pathway in KEGG",br(),uiOutput("kegg_hyper"),br(),img(src = "kegg_pathway.png"),br(),dataTableOutput("alzheimerGenesX")),
                    tabPanel("AD Network Core Genes",dataTableOutput("corGenesAD") ),
                    
                    
                    
                  )),
                  tabPanel("Clusters",tabsetPanel(   #cluster5
                    tabPanel("AD cluster 1",tabsetPanel(
                      tabPanel("AD cluster 1 Information" ,h4("Information about the selected gene"),dataTableOutput("gene5desc"),br(),h4("Number of Genes: ",g_AD1_size_nodes),h4("Number of Links: ",g_AD1_size_links)),
                      tabPanel("First neighbors of the selected gene",dataTableOutput("connectedGenes5")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ad1")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes5")),
                      tabPanel("Network visualization",plotOutput("graphic5"))
                    )),
                    
                    
                    tabPanel("AD cluster 2",tabsetPanel(   #cluster 2
                      tabPanel("AD cluster 2 Information",h4("Information about the selected gene"),dataTableOutput("gene2desc"),br(),h4("Number of Genes: ",g_AD2_size_nodes),h4("Number of Links: ",g_AD2_size_links)),
                      tabPanel("First neighbors of the selected gene",dataTableOutput("connectedGenes2")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ad2")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes2")),
                      tabPanel("Network visualization",plotOutput("graphic2"))
                    )),
                    #change to cluster 3
                    tabPanel("AD cluster 3",tabsetPanel(   #cluster 3
                      tabPanel("AD cluster 3 Information" ,h4("Information about the selected gene"),dataTableOutput("gene_ad3desc"),br(),h4("Number of Genes: ",g_AD3_size_nodes),h4("Number of Links: ",g_AD3_size_links)),
                      tabPanel("First neighbors of the selected gene" ,dataTableOutput("connectedGenes_ad3")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ad3")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes_ad3")),
                      tabPanel("Network visualization",plotOutput("graphic_ad3"))
                    )),
                    
                    ###
                    tabPanel("Control cluster 1",tabsetPanel(   #cluster 3
                      tabPanel("Control cluster 1 Information",h4("Information about the selected gene"),dataTableOutput("gene3desc"),br(),h4("Number of Genes: ",g_con1_size_nodes),h4("Number of Links: ",g_con1_size_links)),
                      tabPanel("First neighbors of the selected gene",dataTableOutput("connectedGenes3")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ctrl1")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes3")),
                      tabPanel("Network visualization",plotOutput("graphic3"))
                    )),
                    tabPanel("Control cluster 2",tabsetPanel(   #cluster 4
                      tabPanel("Control cluster 2 Information",h4("Information about the selected gene"),dataTableOutput("gene4desc"),br(),h4("Number of Genes: ",g_con2_size_nodes),h4("Number of Links: ",g_con2_size_links)),
                      tabPanel("First neighbors of the selected gene",h3("First neighbor Genes"),dataTableOutput("connectedGenes")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ctrl2")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes4")),
                      tabPanel("Network visualization",plotOutput("graphic4"))
                    )),
                    tabPanel("Control cluster 3",tabsetPanel(    #cluster1
                      tabPanel("Control cluster 3 Information",h4("Information about the selected gene"),dataTableOutput("gene1desc"),br(),h4("Number of Genes: ",g_con3_size_nodes),h4("Number of Links: ",g_con3_size_links)),
                      tabPanel("First neighbors of the selected gene",dataTableOutput("connectedGenes1")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ctrl3")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenes1")),
                      tabPanel("Network visualization",plotOutput("graphic1"))
                    )),
                    tabPanel("Control cluster 4",tabsetPanel(   #cluster 
                      tabPanel("Control cluster 4 Information",h4("Information about the selected gene"),dataTableOutput("genectrl4desc"),br(),h4("Number of Genes: ",g_con4_size_nodes),h4("Number of Links: ",g_con4_size_links)),
                      tabPanel("First neighbors of the selected gene",h3("First neighbor Genes"),dataTableOutput("connectedGenesctrl4")),
                      tabPanel("Core Genes", dataTableOutput("corGenes_ctrl4")),
                      tabPanel("AD Pathway Genes",dataTableOutput("ADGenesctrl4")),
                      tabPanel("Network visualization" ,plotOutput("graphicctrl4"))
                    ))
                    
                  ))
                  
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
  
  ####
  output$corGenesAD<-renderDataTable({
    for(x in 1:length(cor_genes$hgnc_symbol))cor_genes$genecard[x]<-createLink(cor_genes$hgnc_symbol[x])
    cor_genes[order(-cor_genes$Degree_centrality),]
    
    
  }, options=list(pageLength =25,scrollY =TRUE),escape =FALSE )
  #output$hubgenesX<-renderDataTable({
   # for(x in 1:length(hubgenes$hgnc_symbol))hubgenes$genecard[x]<-createLink(hubgenes$hgnc_symbol[x])
  #  hubgenes<-hubgenes[order(hubgenes$hub_score,decreasing = TRUE ),]
  #  hubgenes
    
  #}, options=list(pageLength =50,scrollY =TRUE),escape =FALSE )
  
 # output$gatekeepersX<-renderDataTable({
  #  for(x in 1:length(gatekeepergenes$hgnc_symbol))gatekeepergenes$genecard[x]<-createLink(gatekeepergenes$hgnc_symbol[x])
  #  gatekeepergenes<-gatekeepergenes[order(gatekeepergenes$between_score,decreasing = TRUE ),]
  #  gatekeepergenes
    
  #}, options=list(pageLength =50,scrollY =TRUE),escape =FALSE )
  ####################
 
  
  ###cluster 1 gene desription text
  output$gene_description <- renderText({##cluster4
    highlight_select<-cluster_con_2$ensembl_gene_id[cluster_con_2$hgnc_symbol==input$select]
    paste0(cluster_con_2$hgnc_symbol[cluster_con_2$ensembl_gene_id==highlight_select],cluster_con_2$description[cluster_con_2$ensembl_gene_id==highlight_select],highlight_select,sep="\n")
    #res <- entrez_search(db = "gene", term = highlight_select,api_key =key1) #search gene db for ids
    #recs <- entrez_fetch(db = "gene", id = res$ids[1:25], rettype = "xml", parsed = TRUE) # get fullrecord with ids
    #esums <- entrez_summary(db = "gene", id = res$ids[1])
    #paste(highlight_select,esums$summary,sep="\n")
    #searches for id and than the id in ncbi for the desrcription of the gene
  })
  
  ###### What are Hub , gatekeeper and alzheimer genes?
  
  output$kegg_hyper <- renderUI({
    url<-a("AD Pathway Genes KEGG", href="https://www.genome.jp/pathway/hsa05010", target="_blank" )
    tagList( url)
  })
  
  
  
  
  ####AD cluster 1
  ################# cluster 5 #
  output$gene5desc<-renderDataTable({
    highlight_select5<-clusterAD1$ensembl_gene_id[clusterAD1$hgnc_symbol==input$select5]
    infoTable5<-data.frame(Ensamble_ID="",Gene_name="",kegg="")
    infoTable5$Ensamble_ID<-highlight_select5
    infoTable5$Gene_name<-clusterAD1$hgnc_symbol[clusterAD1$ensembl_gene_id==highlight_select5]

    infoTable5$GO_Term<-clusterAD1$GO_category[clusterAD1$ensembl_gene_id==highlight_select5]
    infoTable5$Core<-clusterAD1$Core[clusterAD1$ensembl_gene_id==highlight_select5]
    infoTable5$DEG<-clusterAD1$DEG[clusterAD1$ensembl_gene_id==highlight_select5]
    
    infoTable5$kegg             <-clusterAD1$kegg[clusterAD1$ensembl_gene_id==highlight_select5]
    infoTable5$Gene_description<-clusterAD1$description[clusterAD1$ensembl_gene_id==highlight_select5]
    infoTable5$Genecard<- createLink(input$select5)
    infoTable5
    return(infoTable5)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  ###########
  
  
  #########
  output$connectedGenes5 <-renderDataTable({
    highlight_select5<- clusterAD1$ensembl_gene_id[clusterAD1$hgnc_symbol==input$select5]
    
    dataconnect5<-correlatedGenesFunc(g_AD1,clusterAD1,highlight_select5)
    
  }, escape = FALSE)
  

  #############
  #ad genes
  output$ADGenes5 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes5<- clusterAD1 %>% filter(kegg == "kegg")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes5$ensembl_gene_id)){
      alzheimergenes5$genecard[q]<- createLink(alzheimergenes5$ensembl_gene_id[q])
      q=q+1
    }
 
    alzheimergenes5
  }, escape = FALSE)

  output$corGenes_ad1 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ad1<- clusterAD1 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ad1$ensembl_gene_id)){
      corgenes_ad1$genecard[q]<- createLink(corgenes_ad1$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ad1
  }, escape = FALSE)
  ############
  
  output$graphic5 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    
    highlight_select5=isolate(clusterAD1$ensembl_gene_id[clusterAD1$hgnc_symbol==input$select5])
    
    plot5<-isolate(plotxfunc2(g_AD1,highlight_select5,input$range[2],input$range[1],input$checkboxHub,completetableAD1,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkAD1,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)

  },height = 1300,width = 1900)
  
  #########not used anymore
  output$hubGeneDesciptions5<-renderDataTable({
    hubtable5<-completetableAD1[order(completetableAD1$degree,decreasing = TRUE ),]
    for(x in 1:10) hubtable5$genecard[x]<-createLink(hubtable5$genes[x])
    hubtable5<-subset(hubtable5,select = -1)
    head(hubtable5,ceiling(length(V(g_AD1))*0.01))
  }, escape = FALSE)
  
  #####
  output$Gatekeeper5<-renderDataTable({
    gateTable5<-networkAD1
    q=1
    while(q<=length(gateTable5$Gene)){
      gateTable5$name[q]<-clusterAD1$hgnc_symbol[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      gateTable5$genebiotype[q]<-clusterAD1$gene_biotype[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      gateTable5$kegg[q]<-clusterAD1$kegg[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      gateTable5$ParentTerm[q]<-clusterAD1$Manual_added_parentTerm[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      #gateTable5$term[q]<-clusterAD1$term[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      #gateTable5$graph_cluster[q]<-clusterAD1$graph_cluster[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      gateTable5$description[q]<-clusterAD1$description[clusterAD1$ensembl_gene_id==gateTable5$Gene[q]]
      gateTable5$genecard[q]<- createLink(gateTable5$Gene[q])
      q=q+1
    }
    
    gateTable5
  }, escape = FALSE)
  
 #################################
 
  
  ##AD cluster 2
  output$gene2desc<-renderDataTable({
    highlight_select2<-clusterAD2$ensembl_gene_id[clusterAD2$hgnc_symbol==input$select2]
    infoTable2<-data.frame(Ensamble_ID="",Gene_name="",kegg="")
    infoTable2$Ensamble_ID<-highlight_select2
    infoTable2$Gene_name<-clusterAD2$hgnc_symbol[clusterAD2$ensembl_gene_id==highlight_select2]
    
    infoTable2$GO_Term<-clusterAD2$GO_category[clusterAD2$ensembl_gene_id==highlight_select2]
    infoTable2$Core<-clusterAD2$Core[clusterAD2$ensembl_gene_id==highlight_select2]
    infoTable2$DEG<-clusterAD2$DEG[clusterAD2$ensembl_gene_id==highlight_select2]
    
    infoTable2$kegg             <-clusterAD2$kegg[clusterAD2$ensembl_gene_id==highlight_select2]
    infoTable2$Gene_description<-clusterAD2$description[clusterAD2$ensembl_gene_id==highlight_select2]
    
    infoTable2$Genecard<- createLink(input$select2)
    infoTable2
    return(infoTable2)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  ##cluster 2
  
  
  ###cluster 2 correlated genes
  output$connectedGenes2 <-renderDataTable({
    highlight_select2<- clusterAD2$ensembl_gene_id[clusterAD2$hgnc_symbol==input$select2]
    
    
    dataconnect2<-correlatedGenesFunc(g_AD2,clusterAD2,highlight_select2)
    
  }, escape = FALSE)
  #############
  #ad genes
  output$ADGenes2 <-renderDataTable({
    
    
    #alzheimergenes2 <- clusterAD2 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes2<- clusterAD2 %>% filter(kegg == "kegg")
    #alzheimergenes2 <- alzheimergenes2 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes2$ensembl_gene_id)){
      alzheimergenes2$genecard[q]<- createLink(alzheimergenes2$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes2
  }, escape = FALSE)
  
  output$corGenes_ad2 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ad2<- clusterAD2 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ad2$ensembl_gene_id)){
      corgenes_ad2$genecard[q]<- createLink(corgenes_ad2$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ad2
  }, escape = FALSE)
  
  ############
  ##graph2 
  output$graphic2 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select2=isolate(clusterAD2$ensembl_gene_id[clusterAD2$hgnc_symbol==input$select2])
    
    plot2<-isolate(plotxfunc2(g_AD2,highlight_select2,input$range[2],input$range[1],input$checkboxHub,completetableAD2,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkAD2,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  
  
  
  #####ad cluster 3
  
  
  output$gene_ad3desc<-renderDataTable({
    highlight_select7<-clusterAD3$ensembl_gene_id[clusterAD3$hgnc_symbol==input$select7]
    infoTable3<-data.frame(Ensamble_ID="",Gene_name="",kegg="")
    infoTable3$Ensamble_ID<-highlight_select7
    infoTable3$Gene_name<-clusterAD3$hgnc_symbol[clusterAD3$ensembl_gene_id==highlight_select7]
    
    
    infoTable3$GO_Term<-clusterAD3$GO_category[clusterAD3$ensembl_gene_id==highlight_select7]
    infoTable3$Core<-clusterAD3$Core[clusterAD3$ensembl_gene_id==highlight_select7]
    infoTable3$DEG<-clusterAD3$DEG[clusterAD3$ensembl_gene_id==highlight_select7]
    infoTable3$kegg             <-clusterAD3$kegg[cluster_con_1$ensembl_gene_id==highlight_select7]
    
    infoTable3$Gene_description<-clusterAD3$description[clusterAD3$ensembl_gene_id==highlight_select7]
    infoTable3$Genecard<- createLink(input$select7)
    infoTable3
    return(infoTable3)
    
  }, escape = FALSE)
  
  output$connectedGenes_ad3 <-renderDataTable({
    highlight_select7<- clusterAD3$ensembl_gene_id[clusterAD3$hgnc_symbol==input$select7]
    
    
    dataconnect3<-correlatedGenesFunc(g_AD3,clusterAD3,highlight_select7)
    
  }, escape = FALSE)
  
  output$ADGenes_ad3 <-renderDataTable({
  
    alzheimergenes3<- clusterAD3 %>% filter(kegg == "kegg")
    q=1
    while(q<=length(alzheimergenes3$ensembl_gene_id)){
      alzheimergenes3$genecard[q]<- createLink(alzheimergenes3$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes3
  }, escape = FALSE)
  
  output$corGenes_ad3 <-renderDataTable({
    
    
    corgenes_ad3<- clusterAD3 %>% filter(Core == "core")
    
    q=1
    while(q<=length(corgenes_ad3$ensembl_gene_id)){
      corgenes_ad3$genecard[q]<- createLink(corgenes_ad3$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ad3
  }, escape = FALSE)
  
  
  
  output$graphic_ad3 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select7=isolate(clusterAD3$ensembl_gene_id[clusterAD3$hgnc_symbol==input$select7])
   
    plot_ad3<-isolate(plotxfunc2(g_AD3,highlight_select7,input$range[2],input$range[1],input$checkboxHub,completetableAD3,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkAD3,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  
  ########
  
  ########Control cluster 1
  output$gene3desc<-renderDataTable({
    highlight_select3<-cluster_con_1$ensembl_gene_id[cluster_con_1$hgnc_symbol==input$select3]
    infoTable3<-data.frame(Ensamble_ID="",Gene_name="",kegg="")
    infoTable3$Ensamble_ID<-highlight_select3
    infoTable3$Gene_name<-cluster_con_1$hgnc_symbol[cluster_con_1$ensembl_gene_id==highlight_select3]
    
    
    infoTable3$GO_Term<-cluster_con_1$GO_category[cluster_con_1$ensembl_gene_id==highlight_select3]
    infoTable3$Core<-cluster_con_1$Core[cluster_con_1$ensembl_gene_id==highlight_select3]
    infoTable3$DEG<-cluster_con_1$DEG[cluster_con_1$ensembl_gene_id==highlight_select3]
    infoTable3$kegg             <-cluster_con_1$kegg[cluster_con_1$ensembl_gene_id==highlight_select3]
    infoTable3$Gene_description<-cluster_con_1$description[cluster_con_1$ensembl_gene_id==highlight_select3]
    
    infoTable3$Genecard<- createLink(input$select3)
    infoTable3
    return(infoTable3)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  
 
  ###cluster 3 correlated genes
  output$connectedGenes3 <-renderDataTable({
    highlight_select3<- cluster_con_1$ensembl_gene_id[cluster_con_1$hgnc_symbol==input$select3]
    
    
    dataconnect3<-correlatedGenesFunc(g_con1,cluster_con_1,highlight_select3)
    
  }, escape = FALSE)
  #############
  #ad genes
  output$ADGenes3 <-renderDataTable({
    
    
    #alzheimergenes3 <- cluster_con_1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes3 <- cluster_con_1 %>% filter(kegg == "kegg")
    #alzheimergenes3 <- alzheimergenes3 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes3$ensembl_gene_id)){
      alzheimergenes3$genecard[q]<- createLink(alzheimergenes3$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes3
  }, escape = FALSE)
  

  output$corGenes_ctrl1 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ctrl1<- cluster_con_1 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ctrl1$ensembl_gene_id)){
      corgenes_ctrl1$genecard[q]<- createLink(corgenes_ctrl1$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ctrl1
  }, escape = FALSE)
  ############
  ####graph3
  output$graphic3 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    
    highlight_select3=isolate(cluster_con_1$ensembl_gene_id[cluster_con_1$hgnc_symbol==input$select3])
    
    plot3<-isolate(plotxfunc2(g_con1,highlight_select3,input$range[2],input$range[1],input$checkboxHub,completetable_con1,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkcon1,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  ###
  
  
  
  #####Control cluster 2----cluster4

  output$gene4desc<-renderDataTable({
    highlight_select<-cluster_con_2$ensembl_gene_id[cluster_con_2$hgnc_symbol==input$select]
    infoTable4<-data.frame(Ensamble_ID="",Gene_name="",kegg="")
    infoTable4$Ensamble_ID<-highlight_select
    infoTable4$Gene_name<-cluster_con_2$hgnc_symbol[cluster_con_2$ensembl_gene_id==highlight_select]
    infoTable4$kegg             <-cluster_con_2$kegg[cluster_con_2$ensembl_gene_id==highlight_select]
    
    infoTable4$GO_Term<-cluster_con_2$GO_category[cluster_con_2$ensembl_gene_id==highlight_select]
    infoTable4$Core<-cluster_con_2$Core[cluster_con_2$ensembl_gene_id==highlight_select]
    infoTable4$DEG<-cluster_con_2$DEG[cluster_con_2$ensembl_gene_id==highlight_select]
    
    infoTable4$Gene_description<-cluster_con_2$description[cluster_con_2$ensembl_gene_id==highlight_select]
    
    infoTable4$Genecard<- createLink(input$select)
    infoTable4
    return(infoTable4)
    ##put gene description4 in table  Cluster4
  }, escape = FALSE)
  
  
 
  #### Table with geneinformation For The connected genes and all degree of selected gene
  ##cluster4
  output$connectedGenes <-renderDataTable({
    highlight_select<- cluster_con_2$ensembl_gene_id[cluster_con_2$hgnc_symbol==input$select]

    dataconnect4<-correlatedGenesFunc(g_con2,cluster_con_2,highlight_select)
    
  }, escape = FALSE)
  
  #############
  #ad genes
  output$ADGenes4 <-renderDataTable({
    
    
    #alzheimergenes4 <- cluster_con_2 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes4 <- cluster_con_2 %>% filter(kegg == "kegg")
    #alzheimergenes4 <- alzheimergenes4 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes4$ensembl_gene_id)){
      alzheimergenes4$genecard[q]<- createLink(alzheimergenes4$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes4
  }, escape = FALSE)
  
  output$corGenes_ctrl2 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ctrl2<- cluster_con_2 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ctrl2$ensembl_gene_id)){
      corgenes_ctrl2$genecard[q]<- createLink(corgenes_ctrl2$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ctrl2
  }, escape = FALSE)
  

  ############

  ###graph4
  output$graphic4 <-renderPlot({  ###
    set.seed(1277)
    
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select=isolate(cluster_con_2$ensembl_gene_id[cluster_con_2$hgnc_symbol==input$select])
    
    
    plot4<-isolate(plotxfunc2(g_con2,highlight_select,input$range[2],input$range[1],input$checkboxHub,completetable_con2,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkcon2,input$checkboxConnect,input$checkboxwTO) )
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  
  ###############Control cluster 3
  
  output$gene1desc<-renderDataTable({
    highlight_select1<-cluster_con_3$ensembl_gene_id[cluster_con_3$hgnc_symbol==input$select1]
    infoTable1<-data.frame(Ensamble_ID="",Gene_name="",KEGG="")
    infoTable1$Ensamble_ID<-highlight_select1
    infoTable1$Gene_name<-cluster_con_3$hgnc_symbol[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$GO_Term<-cluster_con_3$GO_category[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$Core<-cluster_con_3$Core[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$DEG<-cluster_con_3$DEG[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$kegg <-cluster_con_3$kegg[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$Gene_description<-cluster_con_3$description[cluster_con_3$ensembl_gene_id==highlight_select1]
    infoTable1$Genecard<- createLink(input$select1)
    infoTable1
    ##put gene description4 in table  cluster1
  }, escape = FALSE)
  

  
  ###cluster 1 correlated genes
  output$connectedGenes1 <-renderDataTable({
    highlight_select1<- cluster_con_3$ensembl_gene_id[cluster_con_3$hgnc_symbol==input$select1]
    
    
    dataconnect1<-correlatedGenesFunc(g_con3,cluster_con_3,highlight_select1)
    
  }, escape = FALSE)
  #############
  #ad genes
  output$ADGenes1 <-renderDataTable({
    
    #alzheimergenes1 <- cluster_con_3 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes1 <- cluster_con_3 %>% filter(kegg == "kegg")

    #alzheimergenes1<- alzheimergenes1 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes1$ensembl_gene_id)){
      alzheimergenes1$genecard[q]<- createLink(alzheimergenes1$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes1
  }, escape = FALSE)
  

  output$corGenes_ctrl3 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ctrl3<- cluster_con_3 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ctrl3$ensembl_gene_id)){
      corgenes_ctrl3$genecard[q]<- createLink(corgenes_ctrl3$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ctrl3
  }, escape = FALSE)
  ############
  
  ###graph1
  output$graphic1 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select1=isolate(cluster_con_3$ensembl_gene_id[cluster_con_3$hgnc_symbol==input$select1])
    
    plot1<-isolate(plotxfunc2(g_con3,highlight_select1,input$range[2],input$range[1],input$checkboxHub,completetable_con3,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkcon3,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  
 ####control cluster 4
  output$genectrl4desc<-renderDataTable({
    highlight_select8<-cluster_con_4$ensembl_gene_id[cluster_con_4$hgnc_symbol==input$select8]
    
    
    infoTable8<-data.frame(Ensamble_ID="",Gene_name="",KEGG="")
    infoTable8$Ensamble_ID<-highlight_select8
    infoTable8$Gene_name<-cluster_con_4$hgnc_symbol[cluster_con_4$ensembl_gene_id==highlight_select8]
    
    infoTable8$GO_Term<-cluster_con_4$GO_category[cluster_con_4$ensembl_gene_id==highlight_select8]
    infoTable8$Core<-cluster_con_4$Core[cluster_con_4$ensembl_gene_id==highlight_select8]
    infoTable8$DEG<-cluster_con_4$DEG[cluster_con_4$ensembl_gene_id==highlight_select8]
    
    infoTable8$kegg <-cluster_con_4$kegg[cluster_con_4$ensembl_gene_id==highlight_select8]
    infoTable8$Gene_description<-cluster_con_4$description[cluster_con_4$ensembl_gene_id==highlight_select8]
    infoTable8$Genecard<- createLink(input$select8)
    infoTable8
    ##put gene description4 in table  cluster1
  }, escape = FALSE)
  
  
  
  ###cluster 1 correlated genes
  output$connectedGenesctrl4 <-renderDataTable({
    highlight_select8<- cluster_con_4$ensembl_gene_id[cluster_con_4$hgnc_symbol==input$select8]
    
    
    dataconnect8<-correlatedGenesFunc(g_con4,cluster_con_4,highlight_select8)
    
  }, escape = FALSE)
  #############
  #ad genes
  output$ADGenesctrl4 <-renderDataTable({
    
    #alzheimergenes1 <- cluster_con_3 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    alzheimergenes8 <- cluster_con_4 %>% filter(kegg == "kegg")
    
    #alzheimergenes1<- alzheimergenes1 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(alzheimergenes8$ensembl_gene_id)){
      alzheimergenes8$genecard[q]<- createLink(alzheimergenes8$ensembl_gene_id[q])
      q=q+1
    }
    
    alzheimergenes8
  }, escape = FALSE)
  
  
  output$corGenes_ctrl4 <-renderDataTable({
    
    
    #alzheimergenes5 <- clusterAD1 %>% filter(ensembl_gene_id %in% alzheimergenes$ensembl_gene_id)
    corgenes_ctrl4<- cluster_con_4 %>% filter(Core == "core")
    #alzheimergenes5 <- alzheimergenes5 %>% select(-Manual_added_parentTerm, -term)
    q=1
    while(q<=length(corgenes_ctrl4$ensembl_gene_id)){
      corgenes_ctrl4$genecard[q]<- createLink(corgenes_ctrl4$ensembl_gene_id[q])
      q=q+1
    }
    
    corgenes_ctrl4
  }, escape = FALSE)
  ############
  
  ###graph1
  output$graphicctrl4 <-renderPlot({  ###
    set.seed(1277)
    ##########
    if(input$goButton == 0) return()
    ###########
    highlight_select8=isolate(cluster_con_4$ensembl_gene_id[cluster_con_4$hgnc_symbol==input$select8])
    
    plot8<-isolate(plotxfunc2(g_con4,highlight_select8,input$range[2],input$range[1],input$checkboxHub,completetable_con4,input$checkboxAlzheimer,alzheimergenes,input$GateKeep,networkcon4,input$checkboxConnect,input$checkboxwTO))
    #               (gcx,highlight_selectx,SelectPositiv,clusterx,range2,SelectNegativ,range1,checkboxHub,completetablex,checkboxAlzheimer,alzheimergenes,GateKeep,networkx,checkboxConnect)
    
  },height = 1300,width = 1500)
  
  
  
  ############


  
}



# Run the app ----
shinyApp(ui = ui, server = server,options =list (launch.browser =TRUE))


