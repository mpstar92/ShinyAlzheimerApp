### Data Preperation

#biomart
library(biomaRt)
library(igraph)

suppressWarnings({

  library(DT, warn.conflicts = FALSE)
  library("data.table", warn.conflicts = FALSE)
  library(dplyr, warn.conflicts = FALSE)
  
})
#library(gprofiler2)
#library(XML)
#library(rentrez)
#set_entrez_key("275b491fc56280b6c69e8b23b1073a142108")
#key1<-"275b491fc56280b6c69e8b23b1073a142108"


setwd("../ShinyApp_CoExp_Networks")

AD_cluster1 <- read.table("data/shiny_data/Network_table/AD_cluster1.txt",header=T)
AD_cluster1_nodes <- read.csv(file="data/shiny_data/Network_table/AD_cluster1_nodes.csv", sep = ',')
AD_cluster2 <- read.table("data/shiny_data/Network_table/AD_cluster2.txt",header=T)
AD_cluster2_nodes <- read.csv(file="data/shiny_data/Network_table/AD_cluster2_nodes.csv", sep = ',')

control_cluster1 <- read.table("data/shiny_data/Network_table/control_cluster1.txt",header=T)
control_cluster1_nodes <- read.csv(file="data/shiny_data/Network_table/control_cluster1_nodes.csv", sep = ',')
control_cluster2 <- read.table("data/shiny_data/Network_table/control_cluster2.txt",header=T)
control_cluster2_nodes <- read.csv(file="data/shiny_data/Network_table/control_cluster2_nodes.csv", sep = ',')
control_cluster3 <- read.table("data/shiny_data/Network_table/control_cluster3.txt",header=T)
control_cluster3_nodes <- read.csv(file="data/shiny_data/Network_table/control_cluster3_nodes.csv", sep = ',')


AD_cluster1_att <- read.table("data/shiny_data/node_attribute/AD_cluster1_nodes.txt", fill = TRUE ,header=T, sep = "\t" )

AD_cluster2_att <- read.table("data/shiny_data/node_attribute/AD_cluster2_nodes.txt", fill = TRUE ,header=T, sep = "\t" )


control_cluster1_att <- read.table("data/shiny_data/node_attribute/control_cluster1_nodes.txt", fill = TRUE ,header=T, sep = "\t" )

control_cluster2_att <- read.table("data/shiny_data/node_attribute/control_cluster2_nodes.txt", fill = TRUE ,header=T, sep = "\t" )

###############
merged_ad2 <- merge(AD_cluster2_nodes, AD_cluster2_att[, -c(1,2,4)], by = "ensembl_gene_id", all.x = TRUE)

# delete all rows where ens = NA
dt_filtered <- AD_cluster2_att[is.na(AD_cluster2_att$ensembl_gene_id), ]

AD_cluster2_att_comp <- merge(AD_cluster2_nodes, dt_filtered[, -c(2,3,4)], by = "hgnc_symbol", all.x = TRUE)


#######################




###############
#########################
add_description<-function(cluster_info_table,filename){
  
  my_genes<-cluster_info_table$ensembl_gene_id
  
  Sys.setenv("http_proxy" = "http://my.proxy.org:9999")
  #ensembl <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  lookup3 <- getBM(attributes= c('ensembl_gene_id','description'),filters = "ensembl_gene_id", values = my_genes, mart = ensembl, useCache = FALSE)
  
  
  
  lookup <- merge(cluster_info_table, lookup3, by = "ensembl_gene_id", all.x = TRUE)
  
  #writes information without _ underline
  #filename2<-paste("data/",filename,".txt",sep = "")
  #write_graph(g,filename2,format = "gml")
  filename<-paste("data/shiny_data/data/annotated/",filename,"_info",".txt",sep = "")
  write.csv(lookup,filename,row.names = FALSE)
  lookup
}


###############
AD_all_cluster_nodes <- read.csv(file="data/shiny_data/data/AD_all_genes.csv", sep = ',')
ctrl_all_cluster_nodes <- read.csv(file="data/shiny_data/data/control_all_genes.csv", sep = ',')
############


AD_cluster1 <- read.table("data/shiny_data/data/AD_cluster1.txt",header=T)
xyz_cluster1 <- unique(c(AD_cluster1$Node.1, AD_cluster1$Node.2))
AD_cluster1_node <- AD_all_cluster_nodes[AD_all_cluster_nodes$ensembl_gene_id %in% xyz_cluster1, ]
#
AD_cluster2<- read.table("data/shiny_data/data/AD_cluster2.txt",header=T)
xyz_cluster2 <- unique(c(AD_cluster2$Node.1, AD_cluster2$Node.2))
AD_cluster2_node <- AD_all_cluster_nodes[AD_all_cluster_nodes$ensembl_gene_id %in% xyz_cluster2, ]
#
AD_cluster3 <- read.table("data/shiny_data/data/AD_cluster3.txt",header=T)
xyz_cluster3 <- unique(c(AD_cluster3$Node.1, AD_cluster3$Node.2))
AD_cluster3_node <- AD_all_cluster_nodes[AD_all_cluster_nodes$ensembl_gene_id %in% xyz_cluster3, ]

#
#
ctrl_cluster1 <- read.table("data/shiny_data/data/control_cluster1.txt",header=T)
cxyz_cluster1 <- unique(c(ctrl_cluster1$Node.1, ctrl_cluster1$Node.2))
ctrl_cluster1_node <- ctrl_all_cluster_nodes[ctrl_all_cluster_nodes$ensembl_gene_id %in% cxyz_cluster1, ]
#
ctrl_cluster2 <- read.table("data/shiny_data/data/control_cluster2.txt",header=T)
cxyz_cluster2 <- unique(c(ctrl_cluster2$Node.1, ctrl_cluster2$Node.2))
ctrl_cluster2_node <- ctrl_all_cluster_nodes[ctrl_all_cluster_nodes$ensembl_gene_id %in% cxyz_cluster2, ]
#
ctrl_cluster3 <- read.table("data/shiny_data/data/control_cluster3.txt",header=T)
cxyz_cluster3 <- unique(c(ctrl_cluster3$Node.1, ctrl_cluster3$Node.2))
ctrl_cluster3_node <- ctrl_all_cluster_nodes[ctrl_all_cluster_nodes$ensembl_gene_id %in% cxyz_cluster3, ]
#
ctrl_cluster4 <- read.table("data/shiny_data/data/control_cluster4.txt",header=T)
cxyz_cluster4 <- unique(c(ctrl_cluster4$Node.1, ctrl_cluster4$Node.2))
ctrl_cluster4_node <- ctrl_all_cluster_nodes[ctrl_all_cluster_nodes$ensembl_gene_id %in% cxyz_cluster4, ]



AD_cluster1_comp <- add_description(AD_cluster1_node,filename="AD_cluster1") 
AD_cluster2_comp <- add_description(AD_cluster2_node,filename="AD_cluster2") 
AD_cluster3_comp <- add_description(AD_cluster3_node,filename="AD_cluster3") 

ctrl_cluster1_comp <- add_description(ctrl_cluster1_node,filename="ctrl_cluster1") 
ctrl_cluster2_comp <- add_description(ctrl_cluster2_node,filename="ctrl_cluster2") 
ctrl_cluster3_comp <- add_description(ctrl_cluster3_node,filename="ctrl_cluster3") 
ctrl_cluster4_comp <- add_description(ctrl_cluster4_node,filename="ctrl_cluster4") 
####

#####


cor_genes <- read.csv(file="data/shiny_data/data/AD_core_hubs.csv", sep = ',')

cor_genes[order(-cor_genes$Degree_centrality),]

#####
# create igraph graph for the networks and add information as vertex attributes
###with 5 attributes new
make_graph_Func<-function(cluster_table,cluster_info_table){
  
  colnames(cluster_table)[1:5]<-c("Gene1","Gene2","Phi","Phi_tilde","wTO")
  g<-graph_from_data_frame(cluster_table,directed=FALSE)
  
  i=1
  j=1
  while(i<=length(vertex.attributes(g)$name )){
    j=1
    while(j<=length(cluster_info_table$ensembl_gene_id)){
      
      if(vertex.attributes(g)$name[i]==cluster_info_table$ensembl_gene_id[j]){
        V(g)$hgnc_symbol[i]=cluster_info_table$hgnc_symbol[j]
        V(g)$gene_biotype[i]=cluster_info_table$gene_biotype[j]
        V(g)$GO_category[i]=cluster_info_table$GO_category[j]
        V(g)$kegg[i]=cluster_info_table$kegg[j]
        V(g)$core[i]=cluster_info_table$Core[j]
        V(g)$DEG[i]=cluster_info_table$DEG[j]
        V(g)$cluster[i]=cluster_info_table$Cluster[j]
       # V(g)$description[i]=cluster_info_table$description[j]
      }
     
      j=j+1
     }
    i=i+1
  }
  g
}
############
make_graph_Func2<-function(cluster_table,cluster_info_table){
  
  colnames(cluster_table)[1:3]<-c("Gene1","Gene2","wTO")
  g<-graph_from_data_frame(cluster_table,directed=FALSE)
  
  i=1
  j=1
  while(i<=length(vertex.attributes(g)$name )){
    j=1
    while(j<=length(cluster_info_table$ensembl_gene_id)){
      
      if(vertex.attributes(g)$name[i]==cluster_info_table$ensembl_gene_id[j]){
        V(g)$hgnc_symbol[i]=cluster_info_table$hgnc_symbol[j]
        V(g)$gene_biotype[i]=cluster_info_table$gene_biotype[j]
        V(g)$kegg[i]=cluster_info_table$kegg[j]
        #V(g)$Go_category[i]=cluster_info_table$Go_category[j]
        #V(g)$term[i]=cluster_info_table$term[j]
        #V(g)$graph_cluster[i]=cluster_info_table$graph_cluster[j]
        V(g)$description[i]=cluster_info_table$description[j]
      }
      j=j+1
    }
    i=i+1
  }
  g
}############

g_AD1 <- make_graph_Func(AD_cluster1,AD_cluster1_comp )

g_AD2 <- make_graph_Func(AD_cluster2,AD_cluster2_comp )
g_AD3 <- make_graph_Func(AD_cluster3,AD_cluster3_comp )

g_control1 <- make_graph_Func(ctrl_cluster1,ctrl_cluster1_comp)
g_control2 <- make_graph_Func(ctrl_cluster2,ctrl_cluster2_comp )
g_control3 <- make_graph_Func(ctrl_cluster3,ctrl_cluster3_comp )
g_control4 <- make_graph_Func(ctrl_cluster4,ctrl_cluster4_comp )
#########


###############


#writes information without _ underline
write_graph(g_AD1,"data/shiny_data/data/graphs/AD_cluster1_graph.gml",format = "gml")
write_graph(g_AD2,"data/shiny_data/data/graphs/AD_cluster2_graph.gml",format = "gml")
write_graph(g_AD3,"data/shiny_data/data/graphs/AD_cluster3_graph.gml",format = "gml")

write_graph(g_control1,"data/shiny_data/data/graphs/control_cluster1_graph.gml",format = "gml")
write_graph(g_control2,"data/shiny_data/data/graphs/control_cluster2_graph.gml",format = "gml")
write_graph(g_control3,"data/shiny_data/data/graphs/control_cluster3_graph.gml",format = "gml")
write_graph(g_control4,"data/shiny_data/data/graphs/control_cluster4_graph.gml",format = "gml")



gs<- read_graph("data/shiny_data/graphs/control_cluster3_graph.gml",format = "gml")
gs<-delete_vertex_attr(gs,"id")

#################



combined_ad <- rbind(AD_cluster1, AD_cluster2, AD_cluster3)
g_AD_c <- make_graph_Func(combined_ad,AD_all_genes )

write_graph(g_AD_c,"data/shiny_data/data/graphs/AD_all_graph.gml",format = "gml")



combined_ctrl <- rbind(ctrl_cluster1, ctrl_cluster2, ctrl_cluster3,ctrl_cluster4)
g_ctrl_c <- make_graph_Func(combined_ctrl,ctrl_all_cluster_nodes )

write_graph(g_ctrl_c,"data/shiny_data/data/graphs/ctrl_all_graph.gml",format = "gml")


##################


core_hubs <- dplyr::slice_max(order_by = Corness)

Degree_hubs <- df %>%
  arrange(desc(Degree_centrality)) %>%
  slice(1:round(0.05 * n())) 

Betweenness_hubs <- df %>%
  arrange(desc(Betweenness)) %>%
  slice(1:round(0.05 * n()))



#########










################

colorcode <- read.csv(file="data/shiny_data/data/color code.csv", sep = ',')
library(ggraph)

my_colors <- c("black","black","red","green","yellow","blue","black","black")
names(my_colors) <- c("animals","pets","wild animals","rabbit","dog","cat","polar bear","panda bear")



# Create a sample graph
g <- graph.ring(10,directed = FALSE)
E(g)$direction <- sample(c("up", "down"), ecount(g), replace = TRUE)

# Add a node attribute (e.g., group) with values 'A' and 'B'
#V(g)$group <- c("Mitochondrial energy production", "B","B", "B","B", "B","B", "B","B","B" )
V(g)$group <- c(colorcode$Gos[1:10])
V(g)$group2 <- rep(c(colorcode$Gos[13:15]),length.out=10)
V(g)$GO_category <- c("intracellular trafficking,Mitochondrial energy production,Synapse signaling", "Apoptosis,intracellular trafficking,Others", "Apoptosis,intracellular trafficking", "Others", "Neurogenesis,No enriched GO", "Cytoskeleton organization,No enriched GO", "Neurogenesis,No enriched GO", "Neurogenesis", "Apoptosis", "No enriched GO")

V(g)$properties<-0
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

group_colors_edges<-c(	"AD-specific" ="#fb8072",
                  "control specific" ="#8dd3c7",
                  "Conserved" ="#D9D9D9"
)

edge_colors <- c("up" = "red", "down" = "blue")

# Map the group attribute to border colors
V(g)$frame.color <- group_colors_edges[V(g)$group2]

#V(g)$color <- group_colors[V(g)$group]
E(g)$color <- edge_colors[E(g)$direction]
# Plot the graph


assign_pie_colors <- function(properties, group_colors) {
  # Split properties into a vector
  prop_vector <- strsplit(properties, ",")[[1]]
  
  # Create a color pie based on presence of properties
  pie_colors <- sapply(prop_vector, function(prop) group_colors[prop])
  
  # Return pie chart specification
  return(pie_colors)
}

vertex_pie_colors <- lapply(V(g_AD3)$GO_category, assign_pie_colors, group_colors)

# Count the number of properties for each node (for pie chart segments)
vertex_pie_values <- lapply(V(g_AD3)$GO_category, function(props) {
  prop_vector <- strsplit(props, ",")[[1]]
  pie_values <- rep(1, length(prop_vector))  # Equal segments for each property
  return(pie_values)
})


plot(g_AD3,
     vertex.shape = "pie",                # Set vertex shape to pie
     vertex.pie = vertex_pie_values,      # Assign pie segments
     vertex.pie.color = vertex_pie_colors,# Assign pie colors
     vertex.size = 10,                    # Node size
     vertex.label = NA,                   # No labels for simplicity
    # vertex.frame.color = V(g)$frame.color,
    # vertex.frame.width = 2,
    # edge.arrow.size = 0.5,               # Arrow size
     edge.color = E(g)$color,
     layout = layout_with_fr,             # Layout algorithm
     main = "Node Colors as Pie Charts Based on Vertex Attributes"
)







vertex_attribute_names <- vertex_attr_names(g_AD3)
print(vertex_attribute_names)
edge_attribute_names <- edge_attr_names(g_AD3)
print(edge_attribute_names)


#####################



# Create a graph with 10 nodes
g <- make_ring(10)

# Assign a vertex attribute with multiple properties
V(g)$properties <- c("A,B", "B,C", "A,D", "E", "A,B,C", "D,E", "A", "B", "C,D", "A,B,E")

# Define colors for each property
color_map <- c(
  A = "red",
  B = "green",
  C = "blue",
  D = "purple",
  E = "orange"
)

# Function to assign pie chart colors for each node
assign_pie_colors <- function(properties, color_map) {
  # Split properties into a vector
  prop_vector <- strsplit(properties, ",")[[1]]
  
  # Create a color pie based on presence of properties
  pie_colors <- sapply(prop_vector, function(prop) color_map[prop])
  
  # Return pie chart specification
  return(pie_colors)
}

# Create a list of pie chart colors for each vertex
vertex_pie_colors <- lapply(V(g)$properties, assign_pie_colors, color_map)

# Count the number of properties for each node (for pie chart segments)
vertex_pie_values <- lapply(V(g)$properties, function(props) {
  prop_vector <- strsplit(props, ",")[[1]]
  pie_values <- rep(1, length(prop_vector))  # Equal segments for each property
  return(pie_values)
})

# Plot the graph with vertex.shape as pie charts
plot(g,
     vertex.shape = "pie",                # Set vertex shape to pie
     vertex.pie = vertex_pie_values,      # Assign pie segments
     vertex.pie.color = vertex_pie_colors,# Assign pie colors
     vertex.size = 20,                    # Node size
     vertex.label = NA,                   # No labels for simplicity
     edge.arrow.size = 0.5,               # Arrow size
     edge.color = "gray",                 # Edge color
     layout = layout_with_fr,             # Layout algorithm
     main = "Node Colors as Pie Charts Based on Vertex Attributes"
)

###########################




