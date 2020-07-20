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

# set working directory
# setwd("~/dev/immcantation-kransdorf")

# SECTION 1: generate phylograms for top 10 most frequent clones in a single BCRE dataset

# read BCR file
BCRE1 <- read.csv ("./irepertoire_fasta/BCRE1_CATGTC_primers-pass_igblast_db-pass_clone-pass_germ-pass_airr.tsv",sep="\t")

# count the clones in each BCR file to generate seq_count and seq_freq columns
clones <- countClones(BCRE1, group=c("c_call"))

# subset clones for IGHG
clones_IGHG = subset(clones,clones$c_call=="IGHG")

# sort the seq_freq column from highest to lowest
clones_IGHG_freq <- clones_IGHG[order(clones_IGHG$seq_freq, decreasing = TRUE),]

# select the top 10 clones based on the seq_freq column
IGHG_top_10 <- clones_IGHG_freq[1:10,]

# identify sample ID
BCRE1['sample_id'] = 'BCRE1'

clone_ids = unique(BCRE1$clone_id)

# generate all lineage tree plots
nclones = 10
for (clone_id in IGHG_top_10$clone_id) {
  
  # select for IGHG in the c_call column in each BCR file
  BCRE1_IGHG = subset(BCRE1,BCRE1$c_call=="IGHG")
  
  # subset BCRE1 into individual clone IDs
  BCRE1_clone <- BCRE1_IGHG[which(BCRE1_IGHG$clone_id == clone_id), ]
  
  nclones = nclones - 1
  if (nclones <= 0) {
    break
  }
  
  if (clone_id == "14290") {
    next
  }
  
  print (paste0("Clone_ID: ",toString(clone_id)," Count: ",toString(nrow(BCRE1_clone))))
  
  # make lineage tree from whole tree
  
  # convert c_call column from "factor" to "character" to match sample_id column in BCRE1_clone data
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
  
}


# SECTION 2: create a loop to perform the previous steps on each BCRE dataset (in progress)

# loop for BCRE files

# read BCRE files and sample key
BCRE_files <- dir(path="./irepertoire_fasta/", pattern="*_igblast_db-pass_clone-pass_germ-pass_airr.tsv", full.names=TRUE, recursive=FALSE)
sample_key <- read.csv(file="./kransdorf_sample_key.csv", sep=",")

# create data frame from BCRE files

BCRE_frame = data.frame()
for (file in BCRE_files)
{
  BCRE <- read.csv(file, sep="\t")
  print(file)
  pattern <- "(BCRE)\\d+"
  name <- regmatches(file, regexpr(pattern, file))
  print(name)
  toString(name)
  BCRE['sample_id']=name
  
  sample_row <- filter(sample_key, Sample1ID==name | Sample2ID==name | Bulk_Memory_B_Cell_Seq==name, .preserve=TRUE)
  BCRE['study_id'] = sample_row['STUDY.ID']
  BCRE['sample_type'] = sample_row['Group']
  
  # BCRE_frame <- rbind(BCRE_frame, BCRE) # this might not make sense
  # print(nrow(BCRE_frame)) # this might not make sense
  
}

