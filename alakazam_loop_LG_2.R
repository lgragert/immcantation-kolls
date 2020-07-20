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
Isotype_GeneUsage_Diversity <-function(BCRE_frame, name, subset_on) {
  
  #Generate plot title variables to include description and subset
  if (subset_on == "study_id") {
    sample_name = paste0("PatientID",name)
  } else if (subset_on == "sample_type") {
    sample_name = paste0(name) # Control, Pre, Desens
  } else if (subset_on == "sampleID") {
    sample_name <- (paste0("SampleID",name))
  }
  
  plot_title_ID <- (paste0("Isotype_Diversity_",sample_name))
  plot_title_J <- (paste0("IGHJ_Usage_",sample_name))
  plot_title_length <- (paste0("CDR3_Length_",sample_name))
  plot_title_gravy <- (paste0("CDR3_Hydrophobicity_",sample_name))
  plot_title_basic <- (paste0("CDR3_Basic_Residues_",sample_name))
  plot_title_acidic <- (paste0("CDR3_Acidic_Residues_",sample_name))
  
  #Isotype diversity curves and plots
  y <- alphaDiversity(BCRE_frame, group=("c_call"),min_q=0, max_q=4, step_q=0.1,ci=0.95, nboot=200)
  z <- alphaDiversity(BCRE_frame, group="c_call", min_q=0, max_q=2, step_q=1, nboot=200)
  plot(y, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_Curve.svg"), path="./IGD_plots", device='svg')
  plot(z, 0, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_0.svg"), path="./IGD_plots", device='svg')
  plot(z, 2, colors=IG_COLORS, main_title=plot_title_ID,legend_title="Isotype")
  ggsave(filename=paste0(sample_name,"_Isotype_Diversity_2.svg"), path="./IGD_plots", device='svg')
  
  
  # #Gene Usage for VDJ
  
  IGHV_list <- list("IGHV1", "IGHV2", "IGHV3", "IGHV4")
  IGHD_list <- list("IGHD1", "IGHD2", "IGHD3", "IGHD4", "IGHD5", "IGHD6")
  
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
  c <- aminoAcidProperties(BCRE_frame, seq="junction_aa", nt=FALSE, trim=TRUE, label="CDR3")
  dplyr::select(c[1:3, ], starts_with("CDR3"))
  tmp_theme <- theme_bw() + theme(legend.position="bottom")
  g1 <- ggplot(c, aes(x=c_call, y=CDR3_aa_length)) + tmp_theme +
    ggtitle(plot_title_length) +
    xlab("Isotype") + ylab("Amino_acids") +
    scale_y_continuous(labels=scales::percent) +
    labs(fill="Isotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(aes(fill=c_call))
  ggsave(filename=paste0(name, "_CDR3_aa_length.svg"), path="./IGD_plots", device='svg')
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
} # end function Isotype_GeneUsage_Diversity

BCRE_subset_call <- function(BCRE_frame, subset_on=NULL, file_name="whole_frame"){
  cols <- BCRE_frame[subset_on]
  cols <- unique(cols[[as.name(subset_on)]])
  
  # make plots for each study_id
  for (unicorn in cols) { # study_ids is subset_on='study_id'
    BCRE_subset <- filter(.data=BCRE_frame, !!as.name(subset_on)==unicorn) # subsets to 1
    Isotype_GeneUsage_Diversity(BCRE_subset, name=unicorn, subset_on) 
  }
  
  ifelse(name=paste0("Patient_id", name), subset_on="study_id", subset_on="_")
      #Error in grid.newpage() : svglite only supports one page
  if (is.variable(unicorn)){
    gsgsub(unicorn= "study_id", name=(paste0("Patient_ID_",name, Isotype_GeneUsage_Diversity)))
  }
      #Error in grid.newpage() : svglite only supports one page
  if (is.variable(unicorn == "study_id")){
    gsgsub(name, name=(paste0("Patient_ID_",name, Isotype_GeneUsage_Diversity)))
  }
      #Error in grid.newpage() : svglite only supports one page
  if (is.function(BCRE_subset_call(subset_on = "study_id"))){
    gsgsub(name, name=(paste0("Patient_ID_",name, Isotype_GeneUsage_Diversity)))
  }
      #gives graphs w/o name
  else if (!is.null(subset_on="study_id")){gsgsub(name, name=paste0("Patient_ID_",name, Isotype_GeneUsage_Diversity))}
      #error of r in grid.newpage() : svglite only supports one page
  while (subset_on= "study_id") { name= paste0("Patient_ID_",unicorn)}
      #error of r in grid.newpage() : svglite only supports one page

  Isotype_GeneUsage_Diversity(BCRE_frame, name=file_name, subset_on) 
  
} 


#Step3: Input data frame
BCRE_files <- dir(path="./irepertoire_fasta/", pattern="*_igblast_db-pass_clone-pass_germ-pass_airr.tsv", full.names=TRUE, recursive=FALSE)
sample_key <- read.csv(file="./kransdorf_sample_key.csv", sep=",")


#Step4:Set Loop
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
  
  BCRE_frame <- rbind(BCRE_frame, BCRE)
  print(nrow(BCRE_frame))
}


#Step5:Run selected subset
#BCRE_subset_call(BCRE_frame, subset_on=NULL, file_name="whole_frame")
BCRE_subset_call(BCRE_frame, subset_on="study_id", file_name="whole_frame")
BCRE_subset_call(BCRE_frame, subset_on="sample_type", file_name="whole_frame")