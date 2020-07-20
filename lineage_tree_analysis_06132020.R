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


# set working directory
setwd("~/dev/immcantation-kransdorf")

# read BCR file
BCRE1 <- read.csv ("./irepertoire_fasta/BCRE1_CATGTC_primers-pass_igblast_db-pass_clone-pass_germ-pass_airr.tsv",sep="\t")

# identify sample_id

# count the clones in each BCR file to generate seq_count and seq_freq columns
clones <- countClones(BCRE1, group=c("c_call"))

# make two empty columns for seq_freq and seq_count
# BCRE1["seq_freq"] <- NULL
# BCRE1["seq_count"] <- NULL
BCRE1['sample_id'] = 'BCRE1'

clone_ids = unique(BCRE1$clone_id)

# generate all lineage tree plots
nclones = 10
for (clone_id in clones$clone_id) {

  
  # select for IGHG in the c_call column in each BCR file
  BCRE1_IGHG = subset(BCRE1,BCRE1$c_call=="IGHG")
  
  # subset BCRE1 into individual clone IDs
  BCRE1_clone <- BCRE1_IGHG[which(BCRE1_IGHG$clone_id == clone_id), ]
  
  nclones = nclones - 1
  if (nclones <= 0) {
    break
  }
  
  # skip clones where fewer than threshold number of reads
  if (nrow(BCRE1_clone) < 1000) {
    next
  }
  if (clone_id == "14290") {
    next
  }
  
  print (paste0("Clone_ID: ",toString(clone_id)," Count: ",toString(nrow(BCRE1_clone))))
  
  # make lineage tree from whole tree
  # run PHYLIP
  
  clone <- makeChangeoClone(BCRE1_clone, text_fields=c("sample_id", "c_call"),pad_end="TRUE",add_count = "TRUE")
  
  print ("ChangeoClone created")
  
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
       vertex.label.color="black", vertex.size=40)
  # Add legend
  legend("topleft", c("Germline", "Inferred", "Sample"), 
         fill=c("black", "white", "steelblue"), cex=0.75)
  
  # TODO - Truncate sequence_id
  # TODO - Shrink size of circles
  
  # Make phylogram
  
  # convert to phylo
  phylo <- graphToPhylo(graph)
  
  #plot using ape
  
  # pdf("phylo_9585.pdf") # exports plot to PDF

  plot_title = paste0("Phylogram of Clone ",clone_id)
  plot_filename = paste0("./BCRE1_",clone_id,"_phylo.svg")
  print (paste0("Output Plot: ",plot_filename))
  svg(plot_filename, width=20, height=20)
  plot(phylo, main=plot_title, show.node.label=TRUE,cex=0.1)
  dev.off()

  # ggsave(plot_filename, plot=plot_phylo, device="svg", width=10, height=8)
  # Error in UseMethod("grid.draw") : 
  #   no applicable method for 'grid.draw' applied to an object of class "list"
      
}




# select for IGHG in the c_call column in each BCR file
IGHG_clones = subset(clones,clones$c_call=="IGHG")

# sort the seq_freq column from highest to lowest
IGHG_clones_sorted_by_frequency <- IGHG_clones[order(IGHG_clones$seq_freq, decreasing = TRUE),]

# select the top 10 clones based on the seq_freq column
IGHG_top_10 <- IGHG_clones_sorted_by_frequency[1:10,]

# generate a list of the clone_id numbers
print(IGHG_top_10$clone_id)

# lineage tree loop through 10 rows

# create ChangeoClone object for clone
# clone <- makeChangeoClone(BCRE1, text_fields=c("sample_id", "c_call"),pad_end="TRUE",add_count = "TRUE")

# look at abundance of various sequences
# clone@data[, c("sample_id", "c_call","collapse_count")]

# run PHYLIP
# phylip_exec <- "/Users/mariandribus/Desktop/phylip-3.695/exe/dnapars.app/Contents/MacOS/dnapars"
# graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)

# output the lineage trees to separate files




