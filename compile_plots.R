## compile figures from plot panels
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)

## Empty place holder
empty <- ggplot() +
          theme_bw() +
          theme(panel.border = element_blank())

## Load cohort plot
#cohort_plot <- readRDS("cohort_plot.RDS")
#cohort_plot <- cohort_plot + theme(legend.position = "right")

## Load copy number plots
substraction_plot <- readRDS("copy_number_analysis/broad_analysis/plots/genome_paired_substractionPlot.RDS")
focal_rates_clinic <- readRDS("copy_number_analysis/focal_analysis/plots/amp_del_gene_plot_CNA.RDS")
focal_rates_clinic <- focal_rates_clinic + theme(legend.position = "bottom")

## Load signature plots 
signature_stacked_all <- readRDS("copy_number_signatures/plots/signatures_allSamples_stacked_bar.RDS")
signature_stacked_all <- signature_stacked_all + theme(legend.position = "none")

sig_box_wilcox <- readRDS("copy_number_signatures/plots/sig_box_wilcox.RDS")
sig_box_wilcox <- sig_box_wilcox + labs(title = "",caption = "") + theme(legend.position = "none")
                        
sig_box_paired_wilcox <- readRDS("copy_number_signatures/plots/sig_box_paired_wilcox.RDS")
sig_box_paired_wilcox <- sig_box_paired_wilcox +
    labs(title = "",caption = "") + theme(legend.position = "none")

signature_stacked_alt <- readRDS("copy_number_signatures/plots/signatures_stacked_bar_alt.RDS")
signature_stacked_alt <- signature_stacked_alt + theme(legend.position = "none")

signature_paired_line <- readRDS("copy_number_signatures/plots/paired_line_plot.RDS")
signature_paired_line <- signature_paired_line + theme(legend.position = "none")

sig_legend_a <- get_legend(
    signature_stacked_all + 
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)

sig_legend_b <- get_legend(
    sig_box_wilcox + 
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)

## Load tissue plots
cna_tissue_rates <- readRDS("copy_number_analysis/sites_of_relapse/plots/cna_rates_tissue.RDS")
cna_tissue_rates <- cna_tissue_rates + theme(legend.position = "none")
cna_tissue_counts <- readRDS("copy_number_analysis/sites_of_relapse/plots/cna_count_tissue.RDS")
cna_tissue_counts <- cna_tissue_counts + theme(legend.position = "none",axis.text = element_text(angle = 90))
cna_tissue_ith <- readRDS("copy_number_analysis/sites_of_relapse/plots/ith_tissue.RDS")
cna_tissue_ith <- cna_tissue_ith + theme(legend.position = "none")

tissue_legend_a <- get_legend(
  cna_tissue_counts + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
tissue_legend_b <- get_legend(
  cna_tissue_rates + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

## Figure 1
# Thomas

## Figure 2
# Thomas

## Figure 3
fig3 <- plot_grid(substraction_plot,focal_rates_clinic,labels = c("A","B"),nrow = 2)
ggsave2(filename = "plots/figure_3_render.png",plot = fig3,width = 11,height = 8,units = "in",dpi = 300)
pdf(file = "plots/figure_3_pdf_render.pdf",width = 11,height = 8)
fig3
dev.off()

## Figure 4

fig4 <- plot_grid(signature_stacked_all,sig_box_wilcox,plot_grid(sig_legend_a,sig_legend_b,ncol = 2,labels = c("","")),
                  nrow = 3,
                  rel_heights = c(1,1,0.1),
                  labels = c("A","B"))
ggsave2(filename = "plots/figure_4_render.png",plot = fig4,width = 10,height = 11,units = "in",dpi = 300)  

## Figure 5
fig5 <- plot_grid(plot_grid(signature_paired_line,
                  plot_grid(signature_stacked_alt,empty,ncol = 2,labels = c("B","C")),
                  nrow = 2,
                  rel_heights = c(0.4,0.6),
                  labels = c("A","")),nrow = 2,
                  plot_grid(sig_legend_a,sig_legend_b,ncol = 2,labels = c("","")),rel_heights = c(1,0.1))
ggsave2(filename = "plots/figure_5_render.png",plot = fig5,width = 10,height = 9,units = "in",dpi = 300)  

## Figure 6
fig6 <- plot_grid(
          plot_grid(cna_tissue_rates,
                    cna_tissue_counts,align = "hv",
                    labels = c("A","B")),
          plot_grid(tissue_legend_b,tissue_legend_a),
          nrow = 2,rel_heights = c(1,0.1))
ggsave2(filename = "plots/figure_6_render.png",plot = fig6,width = 10,height = 10,units = "in",dpi = 300)

## Supplemental - TCGA-Britroc focal rate comparison 
britroc_tcga_comp_raw <- readRDS("copy_number_analysis/focal_analysis/plots/britroc_tcga_comp_raw.RDS")
britroc_tcga_comp_raw <- britroc_tcga_comp_raw +
                          labs(title = "") + theme(legend.position = "none")
britroc_tcga_comp_norm <- readRDS("copy_number_analysis/focal_analysis/plots/britroc_tcga_comp_norm.RDS")
britroc_tcga_comp_norm <- britroc_tcga_comp_norm +
                           labs(title = "") + theme(legend.position = "bottom")
sup_plot_britroc_tcga_rates <- plot_grid(britroc_tcga_comp_raw,britroc_tcga_comp_norm,nrow = 2,labels = c("A","B"))
ggsave2(filename = "plots/supplemental_figure_tcga_britroc_rates.png",plot = sup_plot_britroc_tcga_rates,width = 8,height = 8,units = "in",dpi = 300)  
