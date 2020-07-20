# install required packages
if (!require("ggplot2"))
  install.packages("ggplot2")
if (!require("alakazam"))
  install.packages("alakazam")
if (!require("shazam"))
  install.packages("shazam")
if (!require("tidyverse"))
  install.packages("tidyverse")
if (!require("dplyr"))
  install.packages("dplyr")
if (!require("stringr"))
  install.packages("stringr")
if (!require("scales"))
  install.packages("scales")
if (!require("igraph"))
  install.packages("igraph")
if (!require("svglite"))
  install.packages("svglite")
if (!require("Cairo"))
  install.packages("Cairo")

# this script will make plots for the topology analysis for individual lineage trees

# set working directory
# setwd("~/dev/immcantation-kransdorf")

# load BCRE file
# read BCR file
BCRE1 <- read.csv ("./irepertoire_fasta/BCRE1_CATGTC_primers-pass_igblast_db-pass_clone-pass_germ-pass_airr.tsv",sep="\t")

# identify sample_id
BCRE1['sample_id'] = BCRE1

# count the clones in each BCR file to generate seq_count and seq_freq columns
clones <- countClones(BCRE1, group=c("c_call"))

# subset clones for IGHG
clones_IGHG = subset(clones,clones$c_call=="IGHG")

# sort the seq_freq column from highest to lowest
clones_IGHG_freq <- clones_IGHG[order(clones_IGHG$seq_freq, decreasing = TRUE),]

# select the top 10 clones based on the seq_freq column
IGHG_top_10 <- clones_IGHG_freq[9:9,]

# generate all lineage tree plots
nclones = 10
for (clone_id in IGHG_top_10$clone_id) {
  
  # select for IGHG in the c_call column in each BCRE file
  BCRE1_IGHG = subset(BCRE1,BCRE1$c_call=="IGHG")
  
  # subset BCRE into individual clone IDs
  BCRE1_clone <- BCRE1_IGHG[which(BCRE1_IGHG$clone_id == clone_id), ]
  
  nclones = nclones - 1
  if (nclones <= 0) {
    break
  }
  
  # Error in buildPhylipLineage(clone, phylip_exec, rm_temp = TRUE) :
  # The germline and input sequences are not the same length for clone
  
  if (clone_id == "14290") {
    next
  }
  
  print (paste0("Clone_ID: ",toString(clone_id)," Count: ",toString(nrow(BCRE1_clone))))
  
  # make lineage tree from whole tree
  
  # convert c_call column from "factor" to "character" to match sample_id column in BCRE_clone data
  BCRE1_clone$c_call <- as.character(BCRE1_clone$c_call)
  
  # make ChangeoClone
  clone <- makeChangeoClone(BCRE1_clone, text_fields=c("sample_id", "c_call"),pad_end="TRUE",add_count = "TRUE")
  print ("ChangeoClone created")
  
  # run PHYLIP
  phylip_exec <- "/usr/local/bin/dnapars"
  graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
  
  print ("Phylip phylogram generated")
  
  # The graph has shared annotations for the clone
  data.frame(clone_id=graph$clone,
             junction_length=graph$junc_len,
             v_gene=graph$v_gene,
             j_gene=graph$j_gene)
  
  # The vertices have sequence specific annotations
  data.frame(sequence_id=V(graph)$name, 
             c_call=V(graph)$c_call,
             collapse_count=V(graph)$collapse_count)
  
  # Make plot
  # plot(graph)
  
  # Modify graph and plot attributes
  V(graph)$color <- "steelblue"
  V(graph)$color[V(graph)$name == "Germline"] <- "black"
  V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
  V(graph)$label <- V(graph)$c_call
  E(graph)$label <- ""
  
  # Remove large default margins
  par(mar=c(0, 0, 0, 0) + 0.1)
  # Plot graph
  plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
       vertex.label.color="black", vertex.size=10)
  # Add legend
  legend("topleft", c("Germline", "Inferred", "Sample"), 
         fill=c("black", "white", "steelblue"), cex=0.75)
  
  ### this is where the topology vignette code begins 
  
  # plot lineage tree using igraph
  # plot(graph, layout=layout_as_tree)
  
  # And add some annotation complexity to the tree
  V(graph)$sample_id[c(2, 7)] <- "BCRE1"
  V(graph)$c_call[c(2, 7)] <- "IGHG"
  
  # Make a list of example trees excluding multi-isotype trees
  # graph_list <- graph[sapply(graph, function(x) !any(grepl(",", V(x)$c_call)))]
  
  # plotting annotations on a tree
  
  # Set node colors
  V(graph)$color[V(graph)$sample_id == "BCRE1"] <- "seagreen"
  # V(graph)$color[V(graph)$sample_id == "+7d"] <- "steelblue"
  V(graph)$color[V(graph)$name == "Germline"] <- "black"
  V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
  
  # Set node labels
  V(graph)$label <- paste(V(graph)$sample_id, V(graph)$c_call, sep=", ")
  V(graph)$label[V(graph)$name == "Germline"] <- ""
  V(graph)$label[grepl("Inferred", V(graph)$name)] <- ""
  
  # Set node shapes
  V(graph)$shape <- "crectangle"
  V(graph)$shape[V(graph)$name == "Germline"] <- "circle"
  V(graph)$shape[grepl("Inferred", V(graph)$name)] <- "circle"
  
  # Set node sizes
  V(graph)$size <- 60
  V(graph)$size[V(graph)$name == "Germline"] <- 30
  V(graph)$size[grepl("Inferred", V(graph)$name)] <- 15 
  
  # Remove large default margins
  par(mar=c(0, 0, 0, 0) + 0.05)
  
  # Plot the example tree
  plot(graph, layout=layout_as_tree, vertex.frame.color="grey", 
       vertex.label.color="black", edge.label.color="black", 
       edge.arrow.mode=0)
  
  # Add legend
  legend("topleft", c("Germline", "Inferred", "-1h", "+7d"), 
         fill=c("black", "white", "seagreen", "steelblue"), cex=0.75)
  
  # summarizing node properties
  
  # Consider all nodes
  getPathLengths(graph, root="Germline")
  
  # Exclude nodes without an isotype annotation from step count
  getPathLengths(graph, root="Germline", field="c_call", exclude=NA)
  
  # calculating subtree properties
  
  # Summarize tree
  df <- summarizeSubtrees(graph, fields=c("sample_id", "c_call"), root="Germline")
  print(df[1:4])
  print(df[c(1, 5:8)])
  print(df[c(1, 9:12)])
  
  ### EVERYTHING WORKS UNTIL THIS POINT ###
  
  # calculating subtree properties
  ### I don't know if it's correct to use "graph" as the input because it is a list

  # Set sample colors
  sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
  
  # Box plots of node outdegree by sample
  p1 <- plotSubtrees(graph, "sample_id", "outdegree", colors=sample_colors, 
                     main_title="Node outdegree", legend_title="Time", 
                     style="box", silent=TRUE)
  # Box plots of subtree size by sample
  p2 <- plotSubtrees(graph, "sample_id", "size", colors=sample_colors, 
                     main_title="Subtree size", legend_title="Time", 
                     style="box", silent=TRUE)
  # Violin plots of subtree path length by isotype
  p3 <- plotSubtrees(graph, "c_call", "pathlength", colors=IG_COLORS, 
                     main_title="Subtree path length", legend_title="Isotype", 
                     style="violin", silent=TRUE)
  
  ## Warning: Ignoring unknown parameters: fun.y
  
  # Violin plots of subtree depth by isotype
  p4 <- plotSubtrees(graph,  "c_call", "depth", colors=IG_COLORS, 
                     main_title="Subtree depth", legend_title="Isotype", 
                     style="violin", silent=TRUE)
  
  ## Warning: Ignoring unknown parameters: fun.y
  
  # Plot in a 2x2 grid
  gridPlot(p1, p2, p3, p4, ncol=2)
  
  ## No summary function supplied, defaulting to `mean_se()`
  ## No summary function supplied, defaulting to `mean_se()`
  
} # end of PHYLIP loop
