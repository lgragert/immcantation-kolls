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

#Step2: Input data frame
BCRE_files <- dir(path="./irepertoire_fasta/", pattern="*_igblast_db-pass_clone-pass_germ-pass_airr.tsv", full.names=TRUE, recursive=FALSE)
sample_key <- read.csv(file="./kransdorf_sample_key.csv", sep=",")

#Step3: Define Function for vingetes and subsetting
Isotype_GeneUsage_Diversity <-function(x, name, subset_on){
        if (subset_on == "study_id") {
          sample_name = paste0("PatientID",name)
        }
        else if (subset_on == "sample_type") {
          sample_name = paste0(name) # Control, Pre, Desens
        }
        else { # sampleID
          sample_name <- (paste0("SampleID",name))
        }
        plot_title_V1 <- (paste("IGHV2 Usage",sample_name))
        y <- alphaDiversity(x, group=("c_call"),min_q=0, max_q=4, step_q=0.1,ci=0.95, nboot=200)
        z <- alphaDiversity(x, group="c_call", min_q=0, max_q=2, step_q=1, nboot=200)
        plot(y, colors=IG_COLORS, main_title="Isotype Diversity",legend_title="Isotype")
        ggsave(filename=paste0(name,"_Isotype_Diversity_Curve.svg"), path="./", device='svg')
        plot(z, 0, colors=IG_COLORS, main_title="Isotype Diversity",legend_title="Isotype")
        ggsave(filename=paste0(name,"_Isotype_Diversity_0.svg"), path="./", device='svg')
        plot(z, 2, colors=IG_COLORS, main_title="Isotype Diversity",legend_title="Isotype")
        ggsave(filename=paste0(name,"_Isotype_Diversity_2.svg"), path="./", device='svg')
        Vgene <- countGenes(x, gene="v_call", mode="gene")
         IGHV1 <- Vgene %>%
           mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
           filter(getFamily(gene) == "IGHV1")
         IGHV2 <- Vgene %>%
           mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
           filter(getFamily(gene) == "IGHV2")
         IGHV3 <- Vgene %>%
           mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
           filter(getFamily(gene) == "IGHV3")
         IGHV4 <- Vgene %>%
           mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
           filter(getFamily(gene) == "IGHV4")
         V1 <- ggplot(IGHV1, aes(x=gene, y=seq_freq)) +
           theme_bw() +
           ggtitle(plot_title_V1) +
           theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
           ylab("Percent of Repertoire") +
           xlab("") +
           scale_y_continuous(labels=percent) +
           scale_color_brewer(palette="Dark2") +
           geom_point(color="Blue", size=5, alpha=0.8)
         plot(V1)
         ggsave(filename=paste0(name,"_gene_usage_V1.svg"), path="./", device='svg')
         V2 <- ggplot(IGHV2, aes(x=gene, y=seq_freq)) +
           theme_bw() +
           ggtitle("IGHV2 Usage") +
           theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
           ylab("Percent of Repertoire") +
           xlab("") +
           scale_y_continuous(labels=percent) +
           scale_color_brewer(palette="Dark2") +
           geom_point(color="Orange", size=5, alpha=0.8)
         plot(V2)
         ggsave(filename=paste0(name,"_gene_usage_V2.svg"), path="./", device='svg')
         V3 <- ggplot(IGHV3, aes(x=gene, y=seq_freq)) +
           theme_bw() +
           ggtitle("IGHV3 Usage") +
           theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
           ylab("Percent of Repertoire") +
           xlab("") +
           scale_y_continuous(labels=percent) +
           scale_color_brewer(palette="Dark2") +
           geom_point(color="Red", size=5, alpha=0.8)
         plot(V3)
         ggsave(filename=paste0(name, "_gene_usage_V3.svg"), path="./", device='svg')
         V4 <- ggplot(IGHV4, aes(x=gene, y=seq_freq)) +
           theme_bw() +
           ggtitle("IGHV4 Usage") +
           theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
           ylab("Percent of Repertoire") +
           xlab("") +
           scale_y_continuous(labels=percent) +
           scale_color_brewer(palette="Dark2") +
           geom_point(color="Green", size=5, alpha=0.8)
         plot(V4)
         ggsave(filename=paste0(name, "_gene_usage_V4.svg"), path="./", device='svg')
        Dgene <- countGenes(x, gene="d_call", mode="gene")
        IGHD1 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD1")
        IGHD2 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD2")
        IGHD3 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD3")
        IGHD4 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD4")
        IGHD5 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD5")
        IGHD6 <- Dgene %>%
          mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
          filter(getFamily(gene) == "IGHD6")
        d1 <- ggplot(IGHD1, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD1 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Blue", size=5, alpha=0.8)
        plot(d1)
        ggsave(filename=paste0(name, "_gene_usage_D1.svg"), path="./", device='svg')
        d2 <- ggplot(IGHD2, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD2 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Orange", size=5, alpha=0.8)
        plot(d2)
        ggsave(filename=paste0(name, "_gene_usage_D2.svg"), path="./", device='svg')
        d3 <- ggplot(IGHD3, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD3 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Red", size=5, alpha=0.8)
        plot(d3)
        ggsave(filename=paste0(name, "_gene_usage_D3.svg"), path="./", device='svg')
        d4 <- ggplot(IGHD4, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD4 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Green", size=5, alpha=0.8)
        plot(d4)
        ggsave(filename=paste0(name, "_gene_usage_D4.svg"), path="./", device='svg')
        d5 <- ggplot(IGHD5, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD5 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Blue", size=5, alpha=0.8)
        plot(d5)
        ggsave(filename=paste0(name, "_gene_usage_D5.svg"), path="./", device='svg')
        d6 <- ggplot(IGHD6, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD6 Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Blue", size=5, alpha=0.8)
        plot(d6)
        ggsave(filename=paste0(name, "_gene_usage_D6.svg"), path="./", device='svg')
        dg <- ggplot(Dgene, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHD Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Blue", size=5, alpha=0.8)
        plot(dg)
        ggsave(filename=paste0(name, "_gene_usage_D.svg"), path="./", device='svg')
        Jgene <- countGenes(x, gene="j_call", mode="gene")
        Jg <- ggplot(Jgene, aes(x=gene, y=seq_freq)) +
          theme_bw() +
          ggtitle("IGHJ Usage") +
          theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
          ylab("Percent of Repertoire") +
          xlab("") +
          scale_y_continuous(labels=percent) +
          scale_color_brewer(palette="Dark2") +
          geom_point(color="Blue", size=5, alpha=0.8)
        plot(Jg)
        ggsave(filename=paste0(name, "_gene_usage_J.svg"), path="./", device='svg')
        c <- aminoAcidProperties(x, seq="junction_aa", nt=FALSE, trim=TRUE, label="CDR3")
        dplyr::select(c[1:3, ], starts_with("CDR3"))
        tmp_theme <- theme_bw() + theme(legend.position="bottom")
        g1 <- ggplot(c, aes(x=c_call, y=CDR3_aa_length)) + tmp_theme +
          ggtitle("CDR3 length") +
          xlab("Isotype") + ylab("Amino acids") +
          scale_y_continuous(labels=scales::percent) +
          labs(fill="Isotype") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_boxplot(aes(fill=c_call))
          ggsave(filename=paste0(name, "_CDR3_aa_length.svg"), path="./", device='svg')
        g2 <- ggplot(c, aes(x=c_call, y=CDR3_aa_gravy)) + tmp_theme +
          ggtitle("CDR3 hydrophobicity") +
          xlab("Isotype") + ylab("GRAVY") +
          scale_y_continuous(labels=scales::percent) +
          labs(fill="Isotype") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_boxplot(aes(fill=c_call))
          ggsave(filename=paste0(name, "_CDR3_aa_gravy.svg"), path="./", device='svg')
        g3 <- ggplot(c, aes(x=c_call, y=CDR3_aa_basic)) + tmp_theme +
          ggtitle("CDR3 basic residues") +
          xlab("Isotype") + ylab("Basic residues") +
          scale_y_continuous(labels=scales::percent) +
          labs(fill="Isotype") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_boxplot(aes(fill=c_call))
          ggsave(filename=paste0(name, "_CDR3_aa_basic.svg"), path="./", device='svg')
        g4 <- ggplot(c, aes(x=c_call, y=CDR3_aa_acidic)) + tmp_theme +
          ggtitle("CDR3 acidic residues") +
          xlab("Isotype") + ylab("Acidic residues") +
          scale_y_continuous(labels=scales::percent) +
          labs(fill="Isotype") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_boxplot(aes(fill=c_call))
          ggsave(filename=paste0(name, "_CDR3_aa_acidic.svg"), path="./", device='svg')
        gridPlot(g1, g2, g3, g4, ncol=2)
        ggsave(filename=paste0(name, "_CDR3_aa.svg"), path="./", device='svg')
}

BCRE_subset_call <- function(BCRE_frame, subset_on=NULL, file_name="whole_frame"){
  if (!is.null(subset_on)){ 
    cols <- BCRE_frame[subset_on]
    cols <- unique(cols[[as.name(subset_on)]])
    for (unicorn in cols) {
      BCRE_subset <- filter(.data=BCRE_frame, !!as.name(subset_on)==unicorn)
      Isotype_GeneUsage_Diversity(BCRE_subset, name=unicorn, subset_on)
    }
  } else {
    Isotype_GeneUsage_Diversity(BCRE_frame, name=file_name, subset_on)
  }
}




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


#Step5:Run
#BCRE_subset_call(BCRE_frame, subset_on=NULL, file_name="whole_frame")
BCRE_subset_call(BCRE_frame, subset_on="study_id", file_name="whole_frame")
BCRE_subset_call(BCRE_frame, subset_on="sample_type", file_name="whole_frame")