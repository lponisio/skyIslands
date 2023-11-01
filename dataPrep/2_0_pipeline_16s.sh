#!/usr/bin/env bash

## *****************************************************************************
##                                   16s 
## *****************************************************************************
# 1. ASSIGN TAXONOMY 16s
## *****************************************************************************
#navigate to Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw

# We need taxonomic classifications 

# Use a information-rich database that is clustered at 99% sequence similarity at least 
# In our case, using Silva for 16s, and NCBI AND RDP for rbcl
# we have to "train" the classifier dataset just once.

## IF YOU ALREADY HAVE THE CLASSIFIERS FOLDER, SKIP TO 1c

# 1a. Download the newest silva 132 database into a new working directory from https://www.arb-silva.de/download/archive/qiime

# use 99_otus_16S.fasta and  consensus_taxonomy_7_levels.txt to create the training set.
# mkdir 16s-trainingclassifier
# cd 16s-trainingclassifier
# wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
# unzip Silva_132_release.zip
# rm SILVA_132_release.zip


#downloaded from: https://docs.qiime2.org/2023.9/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier
#silva 136 SSURef N99 seq and tax files

# 1b. Train feature classifier
# import reference sequences from silva data as a qiime2 artifact

## go back into qiime
docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1

cd ../mnt/SI_pipeline/R2018/2023_sequence_results_raw/16s-trainingclassifier/

## 99 is 99% match between our seq and the database
# qiime tools import \
# --type 'FeatureData[Sequence]' \
# --input-path silva-138-99-seqs.qza \
# --output-path 99_otus_16S.qza
# 
# # import taxonomy strings. Check and see if your taxonomy file is a tab-seperated file without a header.
# # if it doesnt have a header, specify "headerlessTSVTaxonomyFormat" since the default source formats usually have headers
# 
# qiime tools import \
# --type 'FeatureData[Taxonomy]' \
# --input-format TSVTaxonomyFormat \
# --input-path taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
# --output-path 99_otus_16S_taxonomy.qza

# We now have two Qiime2 artifacts, 99_otus_16s.qza (reference sequences) and 99_otus_16s_taxonomy.qza (taxonomic names). 
# trim silva to my region using my sequencing primers. We tell the algorithm our genomic primer forward and reverse sequences
# we do this because taxonomic classification is more accurate when a naive bayes classifier is trained on the region
# of the 16s sequence that we sequenced (Werner et al. 2012).
qiime feature-classifier extract-reads --i-sequences silva-138-99-seqs.qza --p-f-primer CMGGATTAGATACCCKGG --p-r-primer AGGGTTGCGCTCGTTG --o-reads ref-seqs16s.qza

# visualize:
qiime feature-table tabulate-seqs --i-data ref-seqs16s.qza --o-visualization ref-seqs16s.qzv
#qiime tools view ref-seqs16s.qzv # (but this is big and took forever to load...).

## may need to clean up docker memory usage
#docker system prune

#Train the classifier:
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs16s.qza --i-reference-taxonomy silva-138-99-tax.qza  --o-classifier classifier16s.qza

# 1c. classify rep seqs and get the resulting taxonomic ids 

# may need to install  scikit learn
# cd ../ until you get into your main folder of computer/wherever you install stuff
# pip install -U scikit-learn==0.19.1 #OR WHATEVER VERSION IS COMPATIBLE WITH THE CONTAINER! Might be an old version

#CD to the main SI_pipeline folder where the classifiers are
#cd ../../

qiime feature-classifier classify-sklearn --i-classifier classifier16s.qza --i-reads  ../merged/16s/rep-seqs-16s.qza --o-classification  ../merged/16s/taxonomy16s.qza


## switch to the newest version of qiime
exit
docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

# visualize. navigate back to where you have your taxonomy files
cd ../mnt/SI_pipeline/R2018/2023_sequence_results_raw/merged/16s

qiime metadata tabulate --m-input-file taxonomy16s.qza --o-visualization taxonomy16s.qzv

#cant do this step without a master map, which we don't have for the merged files. 
#skip for now. we still wanna do the other steps on the merged files
# cd ../../
# qiime taxa barplot --i-table merged/16s/table16s.qza --i-taxonomy merged/16s/taxonomy16s.qza --m-metadata-file merged/16s/maps/sky2018map16s.txt --o-visualization merged/16s/taxa-bar-plots.qzv


### ************************************************************************
# 2. FILTERING STEPS (PART A)
### ************************************************************************

# Go through and filter out unwanted reads on the merged tables and repseqs.
# We can't do all the filtering steps on the merged files, 
# but we can do some, then we make a phylogenetic tree, then go back and filter round-specific issues

# 2a. Prefilter step: The silva database does not have Lactobacillus michnerii, so our reads for Lactobacillus kunkeei are likely wrong. 

# first visualize taxonomy16s.qzv. from the visualizer, download taxonomy16s.tsv. 
# make a copy of taxonomy16s.tsv and rename it taxonomy16sfixed.tsv - we will make corrections here
# then visualize repseqs16sfiltered.qzv in qiime2view

## unmerged
# qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv

cd merged/16s/

qiime feature-table tabulate-seqs --i-data rep-seqs-16s.qza --o-visualization rep-seqs-16s.qzv


# using these two documents, select the sequence corresponding with the feature ids for L. kunkeei

# Make corrections to taxonomy16sfixed.tsv. convert it to txt manually on your computer.
# Convert this text file back into a qza object using qiime

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy16s.qza/taxonomy16sfixed.txt \
--output-path taxonomy16s.qza/taxonomy16sfixed.qza

# we didnt have to end making these modifications in sky islands or SI 2019 because we didnt have L. kunkeei or michinerii. but left this for future projects
# if you have to make these changes,  make sure to change the paths below appropraitely. we kept everything labeled taxonomy16s.qza, NOT fixed

#2b. filter 1: out the chloroplast and mitochondria reads 

qiime taxa filter-table --i-table table16s.qza --i-taxonomy taxonomy16s.qza --p-exclude mitochondria,chloroplast --o-filtered-table tablefilt1.qza

qiime feature-table summarize --i-table tablefilt1.qza --o-visualization tablefilt1.qzv 

#2c. filter 2: remove sequences only found in one sample 

qiime feature-table filter-features --i-table tablefilt1.qza --p-min-samples 2 --o-filtered-table tablefilt2.qza

qiime feature-table summarize --i-table tablefilt2.qza --o-visualization tablefilt2.qzv


#2d: Filter our rep seqs file so that you can see what samples you have left after filtering and subsampling
qiime feature-table filter-seqs --i-data rep-seqs-16s.qza --i-table tablefilt2.qza --o-filtered-data rep-seqs-16s-filtered.qza

qiime feature-table tabulate-seqs --i-data rep-seqs-16s-filtered.qza --o-visualization rep-seqs-16s-filtered.qzv

qiime tools view rep-seqs-16s-filtered.qzv

### ************************************************************************
## 3. #GENERATE TREE FOR PHYLOGENETIC DIVERSITY ANALYSES
### ************************************************************************

# generate a tree with the merged data

#Now alignment of the reads using MAFFT.
qiime alignment mafft --i-sequences rep-seqs-16s-filtered.qza --o-alignment aligned_repseqs16s.qza

#Then filter out the unconserved, highly gapped columns from alignment
qiime alignment mask --i-alignment aligned_repseqs16s.qza --o-masked-alignment masked_aligned_repseqs16s.qza

#And make a tree with Fasttree which creates a maximum likelihood phylogenetic trees from aligned sequences
#for more information on fastree, http://www.microbesonline.org/fasttree/
qiime phylogeny fasttree --i-alignment masked_aligned_repseqs16s.qza --o-tree unrooted_tree16s.qza

#now root it
qiime phylogeny midpoint-root --i-tree unrooted_tree16s.qza --o-rooted-tree rooted-tree16s.qza

### *************************************************************************
## 4. SPLIT MERGED DATA TO CONTINUE CLEANING
### ************************************************************************

#########made it here! need to make txt file


# We have 5 rounds of data. We need make a txt file for each round with the list of IDs for samples associated with each one

mkdir split
cd split

# the files are made inside split and called "R0samples.txt", "R1samples.txt", "R2samples.txt", "R2samples.txt", "R4samples.txt" 

# R0 = plate 1 2021
# R1 = plate 2 2021
# R2 = plate 3 2021
# R3 = plate 4 2021
# R4 = plate 5 2021
# R5 = plate 6 2021



# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R0samples.txt --o-filtered-table tableR0.qza

# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R1samples.txt --o-filtered-table tableR1.qza

# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R2samples.txt --o-filtered-table tableR2.qza

# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R3samples.txt --o-filtered-table tableR3.qza

# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R4samples.txt --o-filtered-table tableR4.qza

# qiime feature-table filter-samples --i-table tablefilt2.qza --m-metadata-file R5samples.txt --o-filtered-table tableR5.qza


# Now that you have split the tables you can continue to filter out samples. 
# You don't need to filter rep-seqs because we only used that for making the tree, which is done already

### *************************************************************************
# 5. FILTERING STEPS (PART B)
### *************************************************************************


# 5a. Filter out large taxonomic groups of bacteria that are known to be contaminants

# there are known bacteria are salt-loving and often found in buffers. 
# they include  Halomonas, Shewanella, Oceanospirillales, and the acne bacteria Propionibacterium. 


qiime taxa filter-table --i-table tableR0.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR0_f1.qza

# qiime taxa filter-table --i-table tableR1.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR1_f1.qza

# qiime taxa filter-table --i-table tableR2.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR2_f1.qza

# qiime taxa filter-table --i-table tableR3.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR3_f1.qza

# qiime taxa filter-table --i-table tableR4.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR4_f1.qza

# qiime taxa filter-table --i-table tableR5.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tableR5_f1.qza


# 5b. Let's specifically look at sequences in our controls and remove the bacteria that are in them

# cd ../

qiime taxa barplot --i-table tablefilt3.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2018map16s.txt --o-visualization split/taxa-bar-plots-R0_filt3.qzv

# qiime taxa barplot --i-table split/tableR1_f1.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2019_R1map16s.txt --o-visualization split/taxa-bar-plotsR1_f1.qzv

# qiime taxa barplot --i-table split/tableR2_f1.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2019_R2map16s.txt --o-visualization split/taxa-bar-plotsR2_f1.qzv

# qiime taxa barplot --i-table split/tableR3_f1.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2019_R3map16s.txt --o-visualization split/taxa-bar-plotsR3_f1.qzv

# qiime taxa barplot --i-table split/tableR4_f1.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2019_R4map16s.txt --o-visualization split/taxa-bar-plotsR4_f1.qzv


## *************************************************************
## Needed manual step!!!!!!!! Need to do in qiime 2 view
# for each taxa bar plot, download the .csv (level 7).
## *************************************************************


# There are multiple ways to filter from here. 
# we used code from "2_1_pipeline16s.R" to generate a list of sequences in the controls

# use the following rationale for selecting bacterial sequences to remove:
# look at what bacteria show up our control samples. we then, for each one of these
# bacteria, look up if it is in all the samples or just that control and made a call about whether or not to remove it
# if bacteria is present in just the control sample, you want to remove it
# however, some are also present in a lot of other samples! we don't want to lose data. if bacteria is in the control AND in more than 30 samples just leave it and don't filter it out (10% of samples)
# you may want to make an exception if the bacterial contaminant is obviously present in one just one plate, because even though its in a lot of samples its likely a contaminant
# another exception is if you have the contaminant in a lot of samples BUT also in a lot of the controls, get rid of it

#this is one option, making a metadata file of samples to keep and filter based on that
#qiime feature-table filter-features --i-table dada2-16s/tablefilt3.qza --m-metadata-file dada2-16s/featurestokeep_16s.txt --p-where id --o-filtered-table dada2-16s/tablefilt4.qza

# this is one option, making a metadata file of samples to remove and filter based on that
#qiime feature-table filter-features --i-table dada2-16s/tablefilt3.qza --m-metadata-file dada2-16s/contaminants_16s.txt --p-exclude-ids TRUE --o-filtered-table dada2-16s/tablefilt4.qza

#this works the best with the least troubleshooting

## REMOVE ALL THE ENDING ;___
## sample :"Unassigned"

# R0 (2018) dna controls: "DNActrl1", "DNActrl2"
# R0 illumina controls:"16SIlluminaCtrl1", "16SIlluminaCtrl2", "16SIlluminaCtrl3

qiime taxa filter-table --i-table tablefilt3.qza --i-taxonomy taxonomy16s.qza --p-mode exact --p-exclude "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Micrococcales;D_4__Microbacteriaceae;D_5__Galbitalea;__",\
"D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae;__;__",\
"D_0__Bacteria;D_1__Tenericutes;D_2__Mollicutes;D_3__Entomoplasmatales;D_4__Spiroplasmataceae;D_5__Spiroplasma;__",\
"D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Kineosporiales;D_4__Kineosporiaceae;D_5__Kineosporia;D_6__uncultured;bacterium",\
"D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae;D_5__Segetibacter;__",\
"D_0__Bacteria;__;__;__;__;__;__" --o-filtered-table tableR0filt4.qza


# # R1 illumina controls: 5843, 5844, 5845, 5846, 5847, 5848, 5849, 5850
# # R1 dna controls: 5834, 5835, 
# qiime taxa filter-table --i-table split/tableR1_f1.qza --i-taxonomy taxonomy16s.qza --p-mode exact --p-exclude "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Micrococcales;D_4__Microbacteriaceae;D_5__Galbitalea;__",\
# "D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae;__;__"
# ,\
# "D_0__Bacteria;D_1__Tenericutes;D_2__Mollicutes;D_3__Entomoplasmatales;D_4__Spiroplasmataceae;D_5__Spiroplasma;__"
#  --o-filtered-table split/tableR1_f2.qza

#5c. remove control samples. make a txt file with the list of IDs for samples you want to keep, where you remove control IDs. 
# the files are made inside maps and called "R0samplesNoCtrl.txt", "R1samplesNoCtrl.txt", "R2samplesNoCtrl.txt", "R2samplesNoCtrl.txt", "R4samplesNoCtrl.txt" 

cd split

qiime feature-table filter-samples --i-table tableR0filt4.qza --m-metadata-file maps/R0samplesNoCtrl.txt --o-filtered-table tableR0filt5.qza

# qiime feature-table filter-samples --i-table tableR1_f2.qza --m-metadata-file R1samplesNoCtrl.txt --o-filtered-table tableR1_f3.qza

#5d. check it worked by making new taxabarplots without contaminants. use regular maps to RsamplesNoCtrl?

cd ../

qiime taxa barplot --i-table tableR0filt5.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2018map16s.txt --o-visualization taxa-bar-plots-R0_filt5.qzv
# qiime taxa barplot --i-table split/tableR1_f3.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/SI2019_R1map16s.txt --o-visualization split/taxa-bar-plotsR1_f3.qzv

## make sure to visualize and check to see if all were removed!!!!

### *************************************************************************
#6. DETERMINE SUBSAMPLE DEPTH
### *************************************************************************

#We want to make a rarefaction curve to see how subsampling depths influence our alpha diversity metrics.
#Open visualization in Qiime2 View and look at visualized table to see impact of different depths

#For subsampling R0:  1182 reads.  257 (88.01%) samples at the
#specifed sampling depth.

qiime diversity alpha-rarefaction --i-table tableR0filt5.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/SI2018map16s.txt --o-visualization alphararefact16sR0.qzv

qiime feature-table summarize --i-table tableR0filt5.qza --o-visualization tableR0filt5.qzv

# #For subsampling R1:  1130 depth. kept 368 (98.66%) 

# qiime diversity alpha-rarefaction --i-table split/tableR1_f3.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/SI2019_R1map16s.txt --o-visualization split/alphararefact16sR1.qzv

# qiime feature-table summarize --i-table split/tableR1_f3.qza --o-visualization split/tableR1_f3.qzv

# #For subsampling R2: 1090 depth, kept 370 (99.20%)

# qiime diversity alpha-rarefaction --i-table split/tableR2_f3.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/SI2019_R2map16s.txt --o-visualization split/alphararefact16sR2.qzv

# qiime feature-table summarize --i-table split/tableR2_f3.qza --o-visualization split/tableR2_f3.qzv

# #For subsampling R3: 1729 depth kept 367 (98.13%) 

# qiime diversity alpha-rarefaction --i-table split/tableR3_f3.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/SI2019_R3map16s.txt --o-visualization split/alphararefact16sR3.qzv

# qiime feature-table summarize --i-table split/tableR3_f3.qza --o-visualization split/tableR3_f3.qzv

# #For subsampling R4: 1752. kept 355 (97.80%) 

# qiime diversity alpha-rarefaction --i-table split/tableR4_f3.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/SI2019_R4map16s.txt --o-visualization split/alphararefact16sR4.qzv

# qiime feature-table summarize --i-table split/tableR4_f3.qza --o-visualization split/tableR4_f3.qzv



### *************************************************************************
# 7. ALPHA AND BETA DIVERSITY
### *************************************************************************

#Generate core metrics folder. This command makes the new directory and outputs a bunch of files in it
# make sured youre in "merged/16s"

mkdir final

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table tableR0filt5.qza --p-sampling-depth 1182 --m-metadata-file maps/SI2018map16s.txt --output-dir final/core_metrics16sR0 --verbose

# qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table split/tableR1_f3.qza --p-sampling-depth 1130 --m-metadata-file maps/SI2019_R1map16s.txt --output-dir final/core_metrics16sR1 --verbose

# qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table split/tableR2_f3.qza --p-sampling-depth 1090 --m-metadata-file maps/SI2019_R2map16s.txt --output-dir final/core_metrics16sR2 --verbose

# qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table split/tableR3_f3.qza --p-sampling-depth 1729 --m-metadata-file maps/SI2019_R3map16s.txt --output-dir final/core_metrics16sR3 --verbose

# qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table split/tableR4_f3.qza --p-sampling-depth 1752 --m-metadata-file maps/SI2019_R4map16s.txt --output-dir final/core_metrics16sR4 --verbose

#if you want to specify the number of parallel cores running, add this to the code
--p-n-jobs 8

#now we have lots of files in that core_metrics directory. gotta export rarefied_table.qza so we can get a qzv and convert to a csv
#NOT SURE WHAT MAP TO USE HERE, the map or the RnoCtrl file.....
qiime taxa barplot --i-table final/core_metrics16sR0/rarefied_table.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/R0samplesNoCtrl.txt --o-visualization final/core_metrics16sR0/rarefiedtable16s.qzv

#drag the new qzv into qiime 2 view, and download a csv! put this into "final asv tables" folder

