
## **********************************************************
## Load libraries and source files
## **********************************************************

## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html


rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd('skyIslands/')

library(pals)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)

source('dataPrep/src/trees_init.R')
load("indiv.comm16sR0.Rdata")
load('presAbsTable.Rdata')
load('spec_RBCL_16s.Rdata')
load('physeq16s.Rdata')
load('networks/microNets.Rdata')
load('spec_RBCL_16s.Rdata')


setwd("../../skyIslands_saved")

## **********************************************************
## Create microbe phylo trees by genus
## **********************************************************

## Bombus tree
bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus", genus.or.spp='Genus', finalASV, bombus_sites, do_collapse = TRUE)
panelA <- bombus_tree[[1]] + labs(tag="A. Bombus (n=444)")
bombus_meta <- bombus_tree[[2]]
panelA

## Melissodes tree
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', finalASV, melissodes_sites, do_collapse = TRUE)
panelB <- melissodes_tree[[1]] + labs(tag="B. Melissodes (n=51)")
melissodes_meta <- melissodes_tree[[2]]
panelB

## **********************************************************
## Create custom legend
## **********************************************************

## Obligate symbionts
these_obligates <- c("Acetobacteraceae",
                     "Bartonellaceae",
                     "Bifidobacteriaceae",
                     "Lactobacillaceae",
                     "Neisseriaceae",
                     "Orbaceae")


leg_col_order <- c("#F6A600", #Acetobacteriaceae
                   "#604E97", #Bartonellaceae
                   "#B3446C", #Bifidobacteriaceae
                   "#882D17", #Lactobacillaceae
                   "#DCD300", #Neisseriaceae
                   "#8DB600") #Orbaceae

# Create a DataFrame to store legend data
data.leg <- data.frame( 
  Xdata = rnorm(6),
  Ydata = rnorm(6), 
  Family = these_obligates,
  leg_color = leg_col_order)

# Create a Scatter Plot to generate proper legend
gplot <- ggplot(data.leg, aes(Xdata, Ydata, color = Family)) +    
  geom_point(size = 7) +
  scale_color_manual(values=data.leg$leg_color) +
  theme(legend.position='bottom') +
  labs(color='Bacteria Family') +
  guides(colour = guide_legend(nrow = 1)) + theme(legend.key=element_blank(),
                                                   legend.text=element_text(size=12),
                                                  legend.title=element_text(size=12))

## Draw only legend without plot 
panelD  <- get_legend(gplot)                     
plot(get_legend(gplot) )

## **********************************************************
## Create full paneled figure with legend
## **********************************************************

#Set up the layout matrix so that the bottom legend spans both columns
layout <- rbind(c(1, 2),
                c(3, 3)) # The legend will span both columns

# Create the final layout
final_plot <- arrangeGrob(panelA, panelB, panelD,
                          layout_matrix = layout,
                          heights = c(9, 1)) # Adjust the heights as needed

# Center the legend panelD properly within its grid
panelD_centered <- arrangeGrob(panelD, 
                               ncol = 1, 
                               padding = unit(1, "lines"))  # Add padding if necessary

# Open a PDF device to save the plot
pdf("../skyIslands/analysis/microbiome/figures/final/grid_trees.pdf",
    height=8, width=11)  # Adjust width and height as needed

# Plot the final combined figure
grid.arrange(final_plot)

# Close the PDF device to complete saving
dev.off()
