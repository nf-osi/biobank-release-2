#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(synapser)

# look at RNA samples to WES samples only

##################################################
# get data from Synapse
synLogin()
exome_to_rna <- synGet("syn52578817")$path
colnames(exome_to_rna) <- c("sample_a", "sample_b","relatedness", "ibs2", "hom_concordance", "paired", "source_a", "assay_a", "source_b", "assay_b","comparison_type")

# make relatedness scatter plot
mixed_plot <- ggplot(exome_to_rna, aes(x=relatedness, y=ibs2, color=paired, shape=comparison_type)) + 
  geom_point(alpha=0.5, size=2) + 
  scale_color_manual(values = c("yes" = "#c9182c", "no" = "#0D3B66"), name = "Same Individual") +
  scale_shape_manual(values = c("cell_line" = 9, "tumor" = 16, "xenograft" = 3, "normal"=15, "tumor/normal"=17), labels = c("cell line/tumor or normal", "tumor/tumor", "tumor/normal", "xenograft/tumor or normal"), name = "Compared Tissues") +
  theme_bw() +
  xlim(c(-2.5,1)) +
  xlab("Relatedness Coefficient") +
  ylab("Number of sites with same genotype") +
  ggtitle("Relatedness between WES and RNA-seq samples") +
  geom_vline(xintercept = 0.9, linetype = "dashed")

# save plot as PDF file
pdf("rna_to_wes.relatedness_scatter_plots.cutoff_0.9.pdf", height = 8, width = 14)
plot_grid(rna_plot, exome_plot, mixed_plot, ncol = 2, nrow = 2)
dev.off()


##################################################
# get relatedness matrix and row/column annotation data
matrix_rna_to_wes <- synGet("syn52578817")$path
rownames(matrix_rna_to_wes) <- matrix_rna_to_wes[,1]
matrix_rna_to_wes <- matrix_rna_to_wes[,-1]

row_annot <- synGet("syn52578818")$path
rownames(row_annot) <- row_annot[,1]
row_annot <-row_annot[,-1]
colnames(row_annot) <- c("SampleType", "Assay")

col_annot <- synGet("syn52578822")$path
rownames(col_annot) <- col_annot[,1]
col_annot <-col_annot[,-1]
colnames(col_annot) <- c("SampleType", "Assay")


# make relatedness heatmap
ann_colors = list(SampleType = c(cell_line="#61C9A8",normal_tissue="#0D3B66", primary_tumor="#7E52A0", recurrent_tumor="#F4D35E", blood="#ab3428", xenograft_passage="#f49e4c"),
                  Assay = c(Whole_Exome_Sequencing="#2d728f", Bulk_RNA_sequencing="#ff5d8f"))

rna_to_wes_heatmap <- pheatmap(matrix_rna_to_wes, 
                               cluster_rows = F, 
                               cluster_cols = F, 
                               annotation_row = row_annot, 
                               annotation_col = col_annot, 
                               annotation_colors = ann_colors)


# save plot as PDF file
pdf("rna_to_wes.relatedness_heatmap.pdf", height = 25, width = 28)
rna_to_wes_heatmap
dev.off()

