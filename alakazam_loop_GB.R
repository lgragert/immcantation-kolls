#Step1: Set working directory and activate enviroment
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

#Step2: Define Function for vingetes and subsetting
Isotype_GeneUsage_Diversity <-function(BCRE_frame, name, subset_on){
  #Generate plot title variables to include description and subset 
  
  if (subset_on == "study_id"){
    sample_name = paste0("PatientID",name)
  }else if (subset_on == "SampleID"){
    sample_name = name
  }else if (subset_on== "sample_type"){
    if (name == "Sens_undergoing_desens"){
      sample_name <- "Desensitized"
    }else if (name == "Sens_no_desens"){
      sample_name <- "Sensitized"
    }else if (name == "Control"){
      sample_name <- "Control"
    }}
  
  plot_title_ID <- (paste0("Isotype_Diversity_",sample_name))
  plot_title_J <- (paste0("IGHJ_Usage_",sample_name))
  plot_title_length <- (paste0("CDR3_Length_",sample_name))
  plot_title_gravy <- (paste0("CDR3_Hydrophobicity_",sample_name))
  plot_title_basic <- (paste0("CDR3_Basic_Residues_",sample_name))
  plot_title_acidic <- (paste0("CDR3_Acidic_Residues_",sample_name))
  
  #Isotype diversity curves and plots
  # View diversity tests at a fixed diversity order
  # q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
  # A 95% confidence interval will be calculated (ci=0.95)
  # 200 resampling realizations are performed (nboot=200)
  y <- alphaDiversity(BCRE_frame, group=("c_call"),min_q=0, max_q=4, step_q=0.1,ci=0.95, nboot=200)
  # q ranges from 0 (min_q=0) to 2 (max_q=2) in 0.05 increments (step_q=0.05)
  # A 95% confidence interval will be calculated (ci=0.95)
  # 200 resampling realizations are performed (nboot=200)
  z <- alphaDiversity(BCRE_frame, group="c_call", min_q=0, max_q=2, step_q=1, nboot=200)
  plot(y, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_Curve.svg"), path="./IGD_plots", device='svg')
  
  # Test diversity at q=0, q=1 and q=2 (equivalent to species richness, Shannon entropy,
  #Simpson's index) across values in the sample_id column
  # Plot results at q=0 and q=2
  # Plot the mean and standard deviations at q=0 and q=2
  plot(z, 0, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_0.svg"), path="./IGD_plots", device='svg')
  plot(z, 2, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_2.svg"), path="./IGD_plots", device='svg')
  
  # Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
  # Indicate number of sequences resampled from each group in the title
  #sample_main <- paste0("Sample diversity")
  #sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
  #plot(y, colors=sample_colors, main_title=sample_main,legend_title="Sample")
  #plot(z, colors=sample_colors, main_title=sample_main,legend_title="Sample")
  
  
  ##Gene Usage for VDJ
  
  IGHV_list <- list("IGHV1", "IGHV2", "IGHV3", "IGHV4")
  IGHD_list <- list("IGHD1", "IGHD2", "IGHD3", "IGHD4", "IGHD5", "IGHD6")
  
  # Plot IGHV gene usage by sample 1-4
  for (IGHV in IGHV_list){
    plot_title_V <- (paste0(IGHV,'_',sample_name))
    Vgene <- countGenes(BCRE_frame, gene="v_call", mode="gene")
    Vg <- Vgene %>% mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
      filter(getFamily(gene) == IGHV)
    Vplot <- ggplot(Vg, aes(x=gene, y=seq_freq)) +
      theme_bw() +
      ggtitle(plot_title_V) +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
      ylab("Percent_of_Repertoire") +
      xlab("") +
      scale_y_continuous(labels=percent) +
      scale_color_brewer(palette="Dark2") +
      geom_point(color="Blue", size=5, alpha=0.8)
    plot(Vplot)
    ggsave(filename=paste0(sample_name,"_gene_usage_",IGHV,".svg"), path="./IGD_plots", device='svg')
  }
  
  # Plot IGHD gene usage by sample 1-6
  for (IGHD in IGHD_list){
    plot_title_D <- (paste0(IGHD,'_',sample_name))
    Dgene <- countGenes(BCRE_frame, gene="d_call", mode="gene")
    Dg <- Dgene %>% mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
      filter(getFamily(gene) == IGHD)
    Dplot <- ggplot(Dg, aes(x=gene, y=seq_freq)) +
      theme_bw() +
      ggtitle(plot_title_D) +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
      ylab("Percent_of_Repertoire") +
      xlab("") +
      scale_y_continuous(labels=percent) +
      scale_color_brewer(palette="Dark2") +
      geom_point(color="Blue", size=5, alpha=0.8)
    plot(Dplot)
    ggsave(filename=paste0(sample_name, "_gene_usage_",IGHD,".svg"), path="./IGD_plots", device='svg')
  }
  
  # Plot IGHJ gene usage on one graph
  Jgene <- countGenes(BCRE_frame, gene="j_call", mode="gene")
  Jg <- ggplot(Jgene, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle(plot_title_J) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent_of_Repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_brewer(palette="Dark2") +
    geom_point(color="Blue", size=5, alpha=0.8)
  plot(Jg)
  ggsave(filename=paste0(sample_name, "_gene_usage_J.svg"), path="./IGD_plots", device='svg')
  
  
  #CDR3 Diversity 
  # amino acid sequence property analysis
  
  # To reduce the junction sequence to the CDR3 sequence we specify the argument trim=TRUE 
  # which will strip the first and last codon (the conserved residues) prior to analysis. 
  c <- aminoAcidProperties(BCRE_frame, seq="junction_aa", nt=FALSE, trim=TRUE, label="CDR3")
  
  # The full set of properties are calculated by default
  dplyr::select(c[1:3, ], starts_with("CDR3"))
  
  # Define a ggplot theme for all plots
  tmp_theme <- theme_bw() + theme(legend.position="bottom")
  
  # Generate plots for a four of the properties
  # scale_fill_manual(name="Isotype", values=IG_COLORS)
  g1 <- ggplot(c, aes(x=c_call, y=CDR3_aa_length)) + tmp_theme +
    ggtitle(plot_title_length) +
    xlab("Isotype") + ylab("Amino_acids") +
    scale_y_continuous(labels=scales::percent) +
    labs(fill="Isotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(fill=c_call))
  ggsave(filename=paste0(sample_name, "_CDR3_aa_length.svg"), path="./IGD_plots", device='svg')
  
  g2 <- ggplot(c, aes(x=c_call, y=CDR3_aa_gravy)) + tmp_theme +
    ggtitle(plot_title_gravy) +
    xlab("Isotype") + ylab("GRAVY") +
    scale_y_continuous(labels=scales::percent) +
    labs(fill="Isotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(fill=c_call))
  ggsave(filename=paste0(sample_name, "_CDR3_aa_gravy.svg"), path="./IGD_plots", device='svg')
  
  g3 <- ggplot(c, aes(x=c_call, y=CDR3_aa_basic)) + tmp_theme +
    ggtitle(plot_title_basic) +
    xlab("Isotype") + ylab("Basic_residues") +
    scale_y_continuous(labels=scales::percent) +
    labs(fill="Isotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(fill=c_call))
  ggsave(filename=paste0(sample_name, "_CDR3_aa_basic.svg"), path="./IGD_plots", device='svg')
  
  g4 <- ggplot(c, aes(x=c_call, y=CDR3_aa_acidic)) + tmp_theme +
    ggtitle(plot_title_acidic) +
    xlab("Isotype") + ylab("Acidic_residues") +
    scale_y_continuous(labels=scales::percent) +
    labs(fill="Isotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(fill=c_call))
  ggsave(filename=paste0(sample_name, "_CDR3_aa_acidic.svg"), path="./IGD_plots", device='svg')
  gridPlot(g1, g2, g3, g4, ncol=2)
  ggsave(filename=paste0(sample_name, "_CDR3_aa.svg"), path="./IGD_plots", device='svg')
}
# end function Isotype_GeneUsage_Diversity


BCRE_subset_call <- function(BCRE_frame, subset_on=NULL, file_name="whole_frame"){
  if (!is.null(subset_on)){ 
    cols <- BCRE_frame[subset_on]
    cols <- unique(cols[[as.name(subset_on)]])
    for (unicorn in cols) {
      BCRE_subset <- filter(.data=BCRE_frame, !!as.name(subset_on)==unicorn)
      Isotype_GeneUsage_Diversity(BCRE_subset, name=unicorn, subset_on) 
    }
  }else {
    Isotype_GeneUsage_Diversity(BCRE_frame, name=file_name, subset_on) 
  }
} 


#Step3: Input data frame
BCRE_files <- dir(path="./irepertoire_fasta/", pattern="*_igblast_db-pass_clone-pass_germ-pass_airr.tsv", full.names=TRUE, recursive=FALSE)
sample_key <- read.csv(file="./kransdorf_sample_key.csv", sep=",")


#Step4:Set Loop to put all BCRE files into one data frame and to classify based on sample key status
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
  
  
  #!# three separate if...then statements, because the conditions, while mutually exclusive, are for
  #!# different variables
  if (name == sample_row['Sample1ID']) {
    BCRE['SampleID'] = 'PBMC_Pre'
  } else if (name == sample_row['Sample2ID']) {
    BCRE['SampleID'] = 'PBMC_Post'
  } else if (name == sample_row['Bulk_Memory_B_Cell_Seq']) {
    BCRE['SampleID'] = 'MBC'    
  }

  
  
  BCRE_frame <- rbind(BCRE_frame, BCRE)
  print(nrow(BCRE_frame))
}


#Step5:Run selected subset
#BCRE_subset_call(BCRE_frame, subset_on=NULL, file_name="whole_frame")
#seperate by patient id
#BCRE_subset_call(BCRE_frame, subset_on="study_id", file_name="whole_frame")
#seperate by treatment
#BCRE_subset_call(BCRE_frame, subset_on="sample_type", file_name="whole_frame")
#seperate by phase of sample
#BCRE_subset_call(BCRE_frame, subset_on="PBMC_Pre", file_name="whole_frame")
BCRE_subset_call(BCRE_frame, subset_on="SampleID", file_name="whole_frame")
#BCRE_subset_call(BCRE_frame, subset_on="MBC", file_name="whole_frame")