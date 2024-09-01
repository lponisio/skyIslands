############### relative abundance at each site 1/9/23

## tutorial: https://github.com/surh/scip_barplot/blob/a775e0d97b47e8dd40498b6713016bbb7e045003/extended_example.Rmd#L285-L301

##what we want: relative abundance bar charts with:
## a facet for each bee genus at each site
## y axis is rel abund
## x axis is individual bee
## so 7 sites and four bee genera (rows will be site and cols will be genus)
## colored by ASV (family probs)
## maybe will have to subet to each genus then facet only by site but we will see
#rm(list=ls())
wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands'
setwd(wdpath)

library(stringr)
library(tidyverse)
library(ggVennDiagram)
library(gt)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(viridis)
library(UpSetR)
library(ComplexUpset)


source("analysis/microbiome/src/vennFunctions.R")

load("data/spec_RBCL_16s.Rdata")

spec16s <- spec.net %>%
  filter(Apidae == 1) %>%
  select(UniqueID, Site, Genus, GenusSpecies, starts_with('16s')) %>%
  na.omit() %>%
  pivot_longer(-c(UniqueID, Site, Genus, GenusSpecies), names_to = 'Bacteria', values_to = 'Abundance')
     
gensp_summary_table <- spec.net %>%
  filter(Apidae == 1) %>%
  select(GenusSpecies) %>%
  group_by(GenusSpecies) %>%
  summarize(n = n()) %>%
  arrange(desc(n))



apis_subset <- genus_filter(spec16s, 'Apis')
bombus_subset <- genus_filter(spec16s, 'Bombus')
anthophora_subset <- genus_filter(spec16s, 'Anthophora')
megachile_subset <- genus_filter(spec16s, 'Megachile')
melissodes_subset <- genus_filter(spec16s, 'Melissodes')


## adding weights transient and weights obligate
spec.net$WeightsObligate <- if_else(spec.net$PD.obligate!=0 & !is.na(spec.net$PD.obligate),
                                    1, 0)
spec.net$WeightsTransient <- if_else(spec.net$PD.transient!=0 & !is.na(spec.net$PD.transient),
                                     1, 0)

## relative abundance plots
relabund.dat <- spec.net %>%
  filter(Apidae == 1) %>%
  select(UniqueID, Site, Genus, WeightsObligate,
         WeightsTransient, Year, starts_with('16s')) %>%
  na.omit() %>%
  pivot_longer(-c(UniqueID, Site, Genus, WeightsObligate,
                  WeightsTransient, Year), names_to = 'Bacteria', values_to = 'Abundance') %>%
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))




# #change the column names to just include the bacteria family
#relabund.dat$Bacteria <- gsub("; s__.+$","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
#relabund.dat$Bacteria <- gsub('; g__.+$',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

# #change the column names to just include the bacteria order
# relabund.dat$Bacteria <- gsub(".*D_3__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
# relabund.dat$Bacteria <- gsub('\\.D_4__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

## some issues with resolution: for now to get around this i am filtering out samples that included anything that was not resolved to family

library(randomcoloR)

relabund.dat.clean <- relabund.dat %>%
  filter(Genus != 'Agapostemon') %>%
  filter(Abundance > 0.01) ## too many groups -- decide what is the cutoff to show on relabund bars


## commenting out because found a palette i kinda liked
n <- length(unique(relabund.dat.clean$Bacteria))
palette <- distinctColorPalette(n)
strains <- unique(relabund.dat.clean$Bacteria)


#now need to export palette
color_dict <- data.frame(palette, strains)

#write.csv(color_dict, "CustomPalette.csv", row.names=FALSE)

## import color dict

#color_dict <- read.csv('CustomPalette.csv')


## good enough for now probs

## now to make prelim plot

# relabund.dat.clean %>%
#   ggplot(aes(x=UniqueID, y=Abundance)) +
#   geom_bar(aes(fill=Bacteria), stat='identity', position='fill') +
#   facet_grid(rows=vars(relabund.dat.clean$Site), cols=vars(relabund.dat.clean$Genus))

## yuck :( okay needs to do each genus separately

ob_data <- relabund.dat.clean %>%
  filter(WeightsObligate == 1)

trans_data <- relabund.dat.clean %>%
  filter(WeightsTransient == 1)



# plot_genus_by_site_relabund <- function(data, genus){
# 
#   genus_subset <- data %>%
#     filter(Genus == genus) %>%
#     arrange(Bacteria) %>%
#     group_by(Bacteria) %>%
#     mutate(rank = rank(Abundance)) %>%  # rank samples at the level of each peak subtype
#     mutate(UniqueID = reorder(UniqueID, -rank)) %>%   # this reorders samples
#     ungroup()
#     
#   these_colors <- color_dict %>%
#     filter(strains %in% genus_subset$Bacteria) %>%
#     arrange(strains)
# 
#   ggplot(data=genus_subset, aes(x=UniqueID, y=Abundance)) +
#     geom_bar(aes(fill=Bacteria), stat='identity', position='fill', width=1) +
#     facet_grid(~Site, scales='free_x', space='free_x') +
#     theme_classic() +
#     theme(axis.text.x=element_blank(),
#           #axis.ticks.x=element_blank(),
#           axis.text.y = element_text(color = "black"),
#           legend.position='None') +
#     scale_fill_manual(values=as.vector(these_colors$palette), name = "Bacteria Family") +
#     labs(x='Individuals', y='Relative Abundance')
# }

## new funct
plot_genus_by_site_relabund <- function(data, genus) {
  
  # Filter the data for the specified genus
  genus_subset <- data %>%
    filter(Genus == genus)
  
  # Group samples by the most abundant bacterial family (peak class)
  sample_grouping <- genus_subset %>%
    group_by(UniqueID) %>%
    slice_max(order_by = Abundance) %>%
    select(Bacteria, UniqueID) %>%
    rename(peak_class = Bacteria)
  
  # Reorder bars by ranking the samples within each peak class
  genus_reordered <- genus_subset %>%
    inner_join(sample_grouping, by = "UniqueID") %>%
    group_by(peak_class) %>%
    mutate(rank = rank(Abundance)) %>%
    mutate(UniqueID = reorder(UniqueID, -rank)) %>%
    ungroup() %>%
    arrange(Year, rank) %>%
    mutate(UniqueID = factor(UniqueID, levels = unique(UniqueID)))
  
  # Define colors for the plot
  these_colors <- color_dict %>%
    filter(strains %in% genus_reordered$Bacteria) %>%
    arrange(strains)
  
  # Create the plot with reordered bars
  ggplot(data = genus_reordered, aes(x = UniqueID, y = Abundance, fill = Bacteria)) +
    geom_bar(stat = 'identity', position = 'fill', width = 1) +
    facet_grid(Year ~ Site, scales = 'free_x', space = 'free_x') + # Facet by both Year and Site
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          legend.position = 'None',
          strip.background = element_blank(),
          strip.text = element_text(face = "bold") # Make facet labels bold
    ) +
    scale_fill_manual(values = these_colors$palette, name = "Bacteria Family") +
    labs(x = 'Individuals', y = 'Relative Abundance')
}

## TODO beautify plot

# Example usage with your data
# plot_genus_by_site_relabund(my_data, "SomeGenus")


apis_ob_plot <- plot_genus_by_site_relabund(ob_data, 'Apis')
plot(apis_ob_plot)

apis_trans_plot <- plot_genus_by_site_relabund(trans_data, 'Apis')
apis_trans_plot

bombus_ob_plot <- plot_genus_by_site_relabund(ob_data, 'Bombus')
bombus_ob_plot

bombus_trans_plot <- plot_genus_by_site_relabund(trans_data, 'Bombus')
bombus_trans_plot

mell_ob_plot <- plot_genus_by_site_relabund(ob_data, 'Melissodes')
mell_ob_plot

mell_trans_plot <- plot_genus_by_site_relabund(trans_data, 'Melissodes')
mell_trans_plot

# apis_plot <- plot_genus_by_site_relabund(ob_data, 'Apis')
# apis_plot
# 
# bombus_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Bombus')
# bombus_plot
# 
# anthophora_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Anthophora')
# anthophora_plot
# 
# megachile_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Megachile')
# megachile_plot

## lets see how they look combined
library(ggpubr)

ggarrange(
  apis_plot,
  bombus_plot,
  anthophora_plot,
  megachile_plot,
  labels = c("Ap", "B", "An", "M"),
  common.legend = FALSE, legend = "bottom"
)


## venn diagram
library(ggVennDiagram)

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(UniqueID, Site, Genus, starts_with('X16s')) %>%
  na.omit() %>%
  pivot_longer(-c(UniqueID, Site, Genus), names_to = 'Bacteria', values_to = 'Abundance') %>%
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))

##summary count of indiv in each taxa
wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands'
setwd(wdpath)
load("data/spec_RBCL_16s.Rdata")

microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 

level5 <- spec.net %>% select(UniqueID, all_of(microbes))%>% mutate_if(is.numeric, ~1 * (. != 0)) %>%
  select(-UniqueID) %>%summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  t()