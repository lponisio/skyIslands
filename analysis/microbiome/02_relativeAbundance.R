############### relative abundance at each site 1/9/23

## tutorial: https://github.com/surh/scip_barplot/blob/a775e0d97b47e8dd40498b6713016bbb7e045003/extended_example.Rmd#L285-L301

##what we want: relative abundance bar charts with:
## a facet for each bee genus at each site
## y axis is rel abund
## x axis is individual bee
## so 7 sites and four bee genera (rows will be site and cols will be genus)
## colored by ASV (family probs)
## maybe will have to subet to each genus then facet only by site but we will see

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

library(tidyverse)


relabund.dat <- read.csv('spec_RBCL_16s.csv') %>%
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

# #change the column names to just include the bacteria family
relabund.dat$Bacteria <- gsub(".*D_4__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
relabund.dat$Bacteria <- gsub('\\.D_5__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

# #change the column names to just include the bacteria order
# relabund.dat$Bacteria <- gsub(".*D_3__","",relabund.dat$Bacteria)#remove all of the column name before and up to D_4__
# relabund.dat$Bacteria <- gsub('\\.D_4__.*',"",relabund.dat$Bacteria) #remove everything after D_5__ to isolate just the bacteria family

## some issues with resolution: for now to get around this i am filtering out samples that included anything that was not resolved to family

library(randomcoloR)

relabund.dat.clean <- relabund.dat %>%
  filter(!grepl('X16s', Bacteria)) %>%
  filter(Genus != 'Agapostemon') %>%
  filter(Abundance > 0.01) ## too many groups -- decide what is the cutoff to show on relabund bars


## commenting out because found a palette i kinda liked 
# n <- length(unique(relabund.dat.clean$Bacteria))
# palette <- distinctColorPalette(n)
# strains <- unique(relabund.dat.clean$Bacteria)


#now need to export palette 
#color_dict <- data.frame(palette, strains)

#write.csv(color_dict, "CustomPalette.csv", row.names=FALSE)

## import color dict

color_dict <- read.csv('CustomPalette.csv')


## good enough for now probs

## now to make prelim plot

# relabund.dat.clean %>%
#   ggplot(aes(x=UniqueID, y=Abundance)) +
#   geom_bar(aes(fill=Bacteria), stat='identity', position='fill') +
#   facet_grid(rows=vars(relabund.dat.clean$Site), cols=vars(relabund.dat.clean$Genus))

## yuck :( okay needs to do each genus separately

plot_genus_by_site_relabund <- function(data, genus){
  
  genus_subset <- data %>%
    filter(Genus == genus) %>%
    arrange(Bacteria)
  
  these_colors <- color_dict %>%
    filter(strains %in% genus_subset$Bacteria) %>%
    arrange(strains)
  
  ggplot(data=genus_subset, aes(x=UniqueID, y=Abundance)) +
    geom_bar(aes(fill=Bacteria), stat='identity', position='fill', width=1) +
    facet_grid(~Site, scales='free_x', space='free_x') +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.text.y = element_text(color = "black")) +
    scale_fill_manual(values=as.vector(these_colors$palette), name = "Bacteria Family") + 
    labs(x='Individuals', y='Relative Abundance')
}

apis_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Apis')
apis_plot

bombus_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Bombus')
bombus_plot

anthophora_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Anthophora')
anthophora_plot

megachile_plot <- plot_genus_by_site_relabund(relabund.dat.clean, 'Megachile')
megachile_plot

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

genus_filter <- function(dataframe, genus){
  new_df <- dataframe %>%
    filter(Genus == genus) %>%
    filter(Abundance > 0) %>%
    select(Bacteria)
  unique(new_df$Bacteria)
}

apis_subset <- genus_filter(relabund.dat, 'Apis')
bombus_subset <- genus_filter(relabund.dat, 'Bombus')
anthophora_subset <- genus_filter(relabund.dat, 'Anthophora')
megachile_subset <- genus_filter(relabund.dat, 'Megachile')


venn_data <- list(Apis = apis_subset,
                  Bombus = bombus_subset,
                  Anthophora = anthophora_subset,
                  Megachile = megachile_subset)

ggVennDiagram(venn_data,
              edge_size = 2,
              label_alpha = 0,
              set_color = c("black","grey","blue","green")) + 
  scale_color_manual(values = c("black","grey","blue","green")) +
  ggplot2::scale_fill_gradient(low="white",high = "red")
