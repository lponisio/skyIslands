############### relative abundance at each site 1/9/23

## tutorial: https://github.com/surh/scip_barplot/blob/a775e0d97b47e8dd40498b6713016bbb7e045003/extended_example.Rmd#L285-L301

##what we want: relative abundance bar charts with:
## a facet for each bee genus at each site
## y axis is rel abund
## x axis is individual bee
## so 7 sites and four bee genera (rows will be site and cols will be genus)
## colored by ASV (family probs)
## maybe will have to subet to each genus then facet only by site but we will see
rm(list=ls())
wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands'
setwd(wdpath)

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


venn_data <- list(Anthophora = anthophora_subset,
                  Apis = apis_subset,
                  Bombus = bombus_subset,
                  Megachile = megachile_subset,
                  Melissodes = melissodes_subset)


ggVennDiagram(venn_data,
              label_alpha = 0.5) + 
  scale_color_manual(values = plasma(5)) +
  scale_fill_gradient(low="white",high = "black") +
  scale_x_continuous(expand = expansion(mult = .2))

upset_genus <- upset(fromList(venn_data),
                  colnames(fromList(venn_data)),
                  queries=query_by_degree(
                    fromList(venn_data),
                    colnames(fromList(venn_data)),
                    params_by_degree=list(
                      '0'=list(color='black', fill='black'),
                      '1'=list(color='black', fill='black'),
                      '2'=list(color='black', fill='black'),
                      '3'=list(color='black', fill='black'),
                      '4'=list(color='black', fill='black'),
                      '5'=list(color='red', fill='red')
                    )))
upset_genus

# intersection_vals <- process_region_data(Venn(venn_data))
# 
# shared_otus_all <- unlist(intersection_vals$item[intersection_vals$id == '12345'])
# 
# bombus_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Bombus'])
# 
# apis_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Apis'])
# 
# megachile_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Megachile'])
# 
# melissodes_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Melissodes'])
# 
# anthophora_only_otus <- unlist(intersection_vals$item[intersection_vals$name == 'Anthophora'])
# 
# 
# #now making a nice table to show the shared OTUs between all genera
# shared_table_data <- spec16s %>%
#   filter(Bacteria %in% shared_otus_all) %>%
#   select(UniqueID, Bacteria, Abundance) %>%
#   filter(Abundance > 0) %>%
#   filter(!duplicated(Bacteria)) %>%
#   select(-Abundance) %>%
#   mutate(Feature.ID = row_number()) %>%
#   mutate(Taxon = Bacteria) %>%
#   select(-c(Bacteria, UniqueID)) %>%
#   mutate(Taxon = str_replace_all(Taxon, '16s:', '')) %>%
#   parse_taxonomy() %>%
#   mutate(Kingdom = str_replace_all(Kingdom, 'd__', '')) %>%
#   arrange(Phylum, Class, Order, Family, Genus, Species) %>%
#   mutate(Genus = str_replace_all(Genus, '_', ' ')) %>%
#   mutate(Species = str_replace_all(Species, '_', ' ')) %>%
#   #mutate(which_genera = 'Shared Among All') %>%
#   as_tibble()
# 
# gt_shared <- gt(shared_table_data) %>%
#   sub_missing(columns=1:7, missing_text = '') %>%
#   tab_row_group(label = 'Present in all genera',
#                 rows=1:31)
# 
# 
# gt_shared 



### subsetting different intersections based on how many individuals we have of each
### will do one for ones >50, >20, >10, and > 4


## >= 50 indiv
Amellifera_subset <- species_filter(spec16s, 'Apis mellifera') 
Bcentralis_subset <- species_filter(spec16s, 'Bombus centralis')
Bhuntii_subset <- species_filter(spec16s, 'Bombus huntii')
Bbifarius_subset <- species_filter(spec16s, 'Bombus bifarius')
Mconfusus_subset <- species_filter(spec16s, 'Melissodes confusus')
## >= 20 indiv
Amontana_subset <- species_filter(spec16s, 'Anthophora montana')
Brufocinctus_subset <- species_filter(spec16s, 'Bombus rufocinctus')
Bflavifrons_subset <- species_filter(spec16s, 'Bombus flavifrons')
Mfrigida_subset <- species_filter(spec16s, 'Megachile frigida')
## >= 10 indiv
Bnevadensis_subset <- species_filter(spec16s, "Bombus nevadensis")
Mcomata_subset <- species_filter(spec16s, "Megachile comata")
## >= 4 indiv
Bfervidus_subset <- species_filter(spec16s, "Bombus fervidus")
Aurbana_subset <- species_filter(spec16s, "Anthophora urbana")
Bmixtus_subset <- species_filter(spec16s, "Bombus mixtus")
Bmorrisoni_subset <- species_filter(spec16s, "Bombus morrisoni")
Bsonorus_subset <- species_filter(spec16s, "Bombus sonorus")
Dmaura_subset <- species_filter(spec16s, "Dufourea maura")
Mrelativa_subset <- species_filter(spec16s, "Megachile relativa")



venn_data_50 <- list(Apis_mellifera = Amellifera_subset,
                  Bombus_centralis = Bcentralis_subset,
                  Bombus_huntii = Bhuntii_subset,
                  Bombus_bifarius = Bbifarius_subset,
                  Melissodes_confusus = Mconfusus_subset)

venn_50 <- ggVennDiagram(venn_data_50,
              label_alpha = 0.5) + 
  scale_color_manual(values = plasma(5)) +
  scale_fill_gradient(low="white",high = "black") +
  scale_x_continuous(expand = expansion(mult = .2))
venn_50

upset_50 <- upset(fromList(venn_data_50),
                  colnames(fromList(venn_data_50)),
                  queries=query_by_degree(
                    fromList(venn_data_50),
                    colnames(fromList(venn_data_50)),
                    params_by_degree=list(
                      '0'=list(color='black', fill='black'),
                      '1'=list(color='black', fill='black'),
                      '2'=list(color='black', fill='black'),
                      '3'=list(color='black', fill='black'),
                      '4'=list(color='black', fill='black'),
                      '5'=list(color='red', fill='red')
                      )))
upset_50

### now do overlap between just the species that are in > 20

venn_data_20 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset,
                     Anthophora_montana = Amontana_subset,
                     Bombus_rufocinctus = Brufocinctus_subset,
                     Bombus_flavifrons = Bflavifrons_subset,
                     Megachile_frigida = Mfrigida_subset)

upset_20 <- upset(fromList(venn_data_20),
                  colnames(fromList(venn_data_20)),
                  queries=query_by_degree(
                    fromList(venn_data_20),
                    colnames(fromList(venn_data_20)),
                    params_by_degree=list(
                      '0'=list(color='black', fill='black'),
                      '1'=list(color='black', fill='black'),
                      '2'=list(color='black', fill='black'),
                      '3'=list(color='black', fill='black'),
                      '4'=list(color='black', fill='black'),
                      '5'=list(color='black', fill='black'),
                      '6'=list(color='black', fill='black'),
                      '7'=list(color='black', fill='black'),
                      '8'=list(color='black', fill='black'),
                      '9'=list(color='red', fill='red')
                    )))
upset_20

## now do 10
venn_data_10 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset,
                     Anthophora_montana = Amontana_subset,
                     Bombus_rufocinctus = Brufocinctus_subset,
                     Bombus_flavifrons = Bflavifrons_subset,
                     Megachile_frigida = Mfrigida_subset,
                     Bombus_nevadensis = Bnevadensis_subset,
                     Megachile_comata = Mcomata_subset)

upset_10 <- upset(fromList(venn_data_10), nsets = 11, order.by='freq')
upset_10
## now do 4
venn_data_4 <- list(Apis_mellifera = Amellifera_subset,
                     Bombus_centralis = Bcentralis_subset,
                     Bombus_huntii = Bhuntii_subset,
                     Bombus_bifarius = Bbifarius_subset,
                     Melissodes_confusus = Mconfusus_subset,
                     Anthophora_montana = Amontana_subset,
                     Bombus_rufocinctus = Brufocinctus_subset,
                     Bombus_flavifrons = Bflavifrons_subset,
                     Megachile_frigida = Mfrigida_subset,
                     Bombus_nevadensis = Bnevadensis_subset,
                     Megachile_comata = Mcomata_subset,
                    Bombus_fervidus = Bfervidus_subset,
                    Anthophora_urbana = Aurbana_subset,
                    Bombus_mixtus = Bmixtus_subset,
                    Bombus_morrisoni = Bmorrisoni_subset,
                    Bombus_sonorus = Bsonorus_subset,
                    Dufourea_maura = Dmaura_subset)

upset_4 <- upset(fromList(venn_data_4), nsets = 17, order.by='freq')
upset_4


# relabund.dat <- read.csv('spec_RBCL_16s.csv') %>%
#   filter(Apidae == 1) %>%
#   select(UniqueID, Site, Genus, starts_with('X16s')) %>%
#   na.omit() %>%
#   pivot_longer(-c(UniqueID, Site, Genus), names_to = 'Bacteria', values_to = 'Abundance') %>% 
#   mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
#                                       "SM",
#                                       "SC",
#                                       "MM",
#                                       "HM",
#                                       "PL",
#                                       "CH")))
# 
# # #change the column names to just include the bacteria family
# relabund.dat$Bacteria <- gsub(".*D_4__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
# relabund.dat$Bacteria <- gsub('\\.D_5__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family
# 
# # #change the column names to just include the bacteria order
# # relabund.dat$Bacteria <- gsub(".*D_3__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
# # relabund.dat$Bacteria <- gsub('\\.D_4__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family
# 
# ## some issues with resolution: for now to get around this i am filtering out samples that included anything that was not resolved to family
# 
# library(randomcoloR)
# 
# relabund.dat.clean <- relabund.dat %>%
#   filter(!grepl('X16s', Bacteria)) %>%
#   filter(Genus != 'Agapostemon') %>%
#   filter(Abundance > 0.01) ## too many groups -- decide what is the cutoff to show on relabund bars
# 
# 
# ## commenting out because found a palette i kinda liked 
# # n <- length(unique(relabund.dat.clean$Bacteria))
# # palette <- distinctColorPalette(n)
# # strains <- unique(relabund.dat.clean$Bacteria)
# 
# 
# #now need to export palette 
# #color_dict <- data.frame(palette, strains)
# 
# #write.csv(color_dict, "CustomPalette.csv", row.names=FALSE)
# 
# ## import color dict
# 
# color_dict <- read.csv('CustomPalette.csv')
# 
# 
# ## good enough for now probs
# 
# ## now to make prelim plot
# 
# # relabund.dat.clean %>%
# #   ggplot(aes(x=UniqueID, y=Abundance)) +
# #   geom_bar(aes(fill=Bacteria), stat='identity', position='fill') +
# #   facet_grid(rows=vars(relabund.dat.clean$Site), cols=vars(relabund.dat.clean$Genus))
# 
# ## yuck :( okay needs to do each genus separately
# 
# plot_genus_by_site_relabund <- function(data, genus){
#   
#   genus_subset <- data %>%
#     filter(Genus == genus) %>%
#     arrange(Bacteria)
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
#           axis.text.y = element_text(color = "black")) +
#     scale_fill_manual(values=as.vector(these_colors$palette), name = "Bacteria Family") + 
#     labs(x='Individuals', y='Relative Abundance')
# }
# 
# apis_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Apis')
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
# 
# ## lets see how they look combined 
# library(ggpubr)
# 
# ggarrange(
#   apis_plot, 
#   bombus_plot,
#   anthophora_plot,
#   megachile_plot,
#   labels = c("Ap", "B", "An", "M"),
#   common.legend = FALSE, legend = "bottom"
# )
# 
# 
# ## venn diagram
# library(ggVennDiagram)
# 
# spec16s <- read.csv('spec_RBCL_16s.csv') %>%
#   filter(Apidae == 1) %>%
#   select(UniqueID, Site, Genus, starts_with('X16s')) %>%
#   na.omit() %>%
#   pivot_longer(-c(UniqueID, Site, Genus), names_to = 'Bacteria', values_to = 'Abundance') %>% 
#   mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
#                                       "SM",
#                                       "SC",
#                                       "MM",
#                                       "HM",
#                                       "PL",
#                                       "CH")))

