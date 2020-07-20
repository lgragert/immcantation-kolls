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

# create data frame

BCRE_files <- dir(path="./irepertoire_fasta/", pattern="*.tsv", full.names=TRUE, recursive=FALSE)
sample_key <- read.csv(file="./kransdorf_sample_key.csv", sep=",")

BCRE_frame = data.frame()
for (file in BCRE_files) {
  BCRE <- read.csv(file, sep="\t")
  print(file)
  pattern <- "(BCRE)\\d+"
  name <- regmatches(file, regexpr(pattern, file))
  print(name)
  toString(name)
  BCRE['sample_id']=name
  
  # count the clones in each BCR file to generate seq_count and seq_freq columns
  clones <- countClones(BCRE, group=c("c_call"))
  
  # subset clones for IGHG
  clones_IGHG = subset(clones,clones$c_call=="IGHG")
  
  # sort the seq_freq column from highest to lowest
  clones_IGHG_freq <- clones_IGHG[order(clones_IGHG$seq_freq, decreasing = TRUE),]
  
  # select the top 10 clones based on the seq_freq column
  IGHG_top_10 <- clones_IGHG_freq[1:10,]
  
  # sample_row <- filter(sample_key, Sample1ID==name | Sample2ID==name | Bulk_Memory_B_Cell_Seq==name)
  # BCRE['study_id'] = sample_row['STUDY.ID']
  # BCRE['sample_type'] = sample_row['Group']

  # create lineage tree plots for top 10 most frequent clones in each BCRE data set
  
  clone_ids = unique(BCRE$clone_id)
  nclones = 10
  for (clone_id in IGHG_top_10$clone_id) {
    
    # select for IGHG in the c_call column in each BCRE file
    BCRE_IGHG = subset(BCRE,BCRE$c_call=="IGHG")
    
    # subset BCRE into individual clone IDs
    BCRE_clone <- BCRE_IGHG[which(BCRE_IGHG$clone_id == clone_id), ]
    
    nclones = nclones - 1
    if (nclones <= 0) {
      break
    }
    
    if (clone_id == "14290") {
      next
    }
    
    if (clone_id == "35898") {
      next
    }
    
    print (paste0("Clone_ID: ",toString(clone_id)," Count: ",toString(nrow(BCRE_clone))))
    
    # make lineage tree from whole tree
    
    # convert c_call column from "factor" to "character" to match sample_id column in BCRE_clone data
    BCRE_clone$c_call <- as.character(BCRE_clone$c_call)
    
    # make ChangeoClone
    clone <- makeChangeoClone(BCRE_clone, text_fields=c("sample_id", "c_call"),pad_end="TRUE",add_count = "TRUE")
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
    # plot_filename = paste0("./BCRE", sample_id, "_", clone_id,"_phylo.svg")
    plot_filename = paste0("./", name, "_", clone_id,"_phylo.svg")
    print (paste0("Output Plot: ",plot_filename))
    svg(plot_filename, width=20, height=20)
    plot(phylo, main=plot_title, show.node.label=TRUE,cex=0.1)
    dev.off()
    
  } # end of lineage tree loop
  
  
} # end of BCRE file loop


