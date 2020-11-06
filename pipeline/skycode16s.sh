#!/usr/bin/env bash

# ---- We skipped the first chunk of code because  the artifacts described below were already created when we ran through the RBCL code


# qiime2 analysis 

#tutorial for reference: https://docs.qiime2.org/2019.1/tutorials/moving-pictures/
#Another useful tutorial is https://rachaellappan.github.io/VL-QIIME2-analysis/pre-processing-of-sequence-reads.html

#PIPELINE USING DOCKER
#*note: this requires docker to already be installed on the machine that you'll be using for the analyses. 

#*This is much easier to do in person rather than remotely.

#I will start with the original sequence files from the illumina run. To get into Qiime 2 I need to first extract barcodes using qiime1 extract_barcodes.py. This script formats your sequence and barcode data so its compatible for demultiplexing.

#This argument works on the illumina output files, which were called flowcelllane1pair1 and flowcelllane1pair2. usually, pair 1 corresponds to forward and pair 2 corresponds to reverse reads. However, because of how our primers were designed, we need to switch around the order of these. I renamed pair 2 as "rawforward" and pair 1 as "rawreverse". 

#First, log in to the remote computer / cluster where you'll be running the analyses. e.g.:
ssh hc@osmia.dyn.ucr.edu
ssh gsmith@osmia.dyn.ucr.edu
#Build and run a docker container that has qiime1 installed. We used sglim2/qiime-1.9.1. 

#For this command, the -it tells docker to make the root of the container interactive in the terminal. 

 #the v (of -itv) tells docker to link the directories you specify ([local path]:[
container path]) to each other. 

 #note: these paths have to be absolute from root, rather than relative from your current directory
 
#finally, you specify the docker image. 
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline sglim2/qiime-1.9.1


#you are now inside the container! You can navigate around and input commands normally. 
 #*to get back to your local directory in terminal, type "exit". 

#Navigate to the new folder you linked with your local directory to confirm that everything is synched. 
cd mnt/SI_pipeline/16sRBCL
ls
#Now you can use qiime1 commands to extract the barcodes from the fastq files

#### note!, don't extract 

 #*mac archive utility (i.e., double clicking the files) seems to work, though you can also extract them with qiime.
extract_barcodes.py -f forward.fastq.gz -r reverse.fastq.gz  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#cd into parsed_barcodes, then gzip the output using the gzip argument. the asterisk looks for everything that has fastq. this will turn everything into a .gz file, which qiime 2 needs in order to work
cd parsed_barcodes
gzip *.fastq

#in your native working directory, manually rename the reads 1 file as forward.fastq.gz, and read 2 as reverse.fastq.gz and barcodes.fastq.gz

#OPEN ENVIRONMENT for qiime 2 

#exit your qiime1 container (sglim2/qiime-1.9.1) by typing: 
exit

#now run a new container for qiime2. the majority of the remaining steps will be completed in this container.
#as above, we'll associate a volume within the container with the SI_pipeline folder on our local machine.

docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

#Once again you'll be in an interactive terminal within the container. Test qiime2 to make sure it all worked:
qiime --help

## make sure you are in the root directory
#Set working directory to the mounted volume with cd. 
cd mnt/SI_pipeline/16sRBCL

#Check mapping file: 
qiime metadata tabulate --m-input-file sky2018map16s.txt --o-visualization sky2018map16s.qzv

#For this project we have two mapping files. One for 16s and one for RBCL. We are going to
#start this project w/the RBCL map, run through the pipeline, then demultiplex data again
#with the 16s metadata file read in. 

#Use chrome browser and the website  https://view.qiime2.org to view the .qvz files or
qiime tools view sky2018map16s.qzv

#to move back a directory, add ../ before the output object. e.g. you could write "../seq.qza" or you could, before running the code, write cd ../

#you want to specify a file path (instead of just the working directory) for qiime tool imports because the argument works on everything in a directory, 

## rename the files
mv reads1.fastq.gz  forward.fastq.gz
mv reads2.fastq.gz  reverse.fastq.gz

#so just setting the working directory will make it work on everything inside of it
qiime tools import --type EMPPairedEndSequences --input-path /mnt/SI_pipeline/16sRBCL/parsed_barcodes/ --output-path seqs.qza

---- We started this pipeline here because the artifacts described above were already created when we ran through the RBCL code

#go to "DEMULTIPLEX" step next! If didn't use Docker, use the alternative below. 

-----
#DOCKER ALTERNATIVE: PREPARE WORKSPACE WITHOUT DOCKER (DO THIS STEP JUST ONCE)

#Download miniconda and install qiime 2: https://docs.qiime2.org/2019.10/install/
#Download qiime 1, which requires macqiime if you are on a Mac

#OPEN ENVIRONMENT to enter QIIME1: type and enter "macqiime" into terminal. 

#I will start with the original sequence files from the illumina run. To get into Qiime 2 I need to first extract barcodes using qiime1 extract_barcodes.py. This script formats your sequence and barcode data so its compatible for demultiplexing.

#This argument works on the illumina output files, which were called flowcelllane1pair1 and flowcelllane1pair2. usually, pair 1 corresponds to forward and pair 2 corresponds to reverse reads. However, because of how our primers were designed, we need to switch around the order of these. I renamed pair 2 as "rawforward" and pair 1 as "rawreverse". 

extract_barcodes.py -f SI_forward.fastq -r SI_reverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#cd into parsed_barcodes, then gzip the output using the gzip argument. the asterisk looks for everything that has fastq. this will turn everything into a .gz file, which qiime 2 needs in order to work

cd parsed_barcodes
gzip *.fastq

#in your native working directory, manually rename the reads 1 file as forward.fastq.gz, and read 2 as reverse.fastq.gz and barcodes.fastq.gz


#OPEN ENVIRONMENT for qiime 2 

#exit macqiime first
exit

#Activate conda environment 
source activate qiime2-2017.12
source tab-qiime
#Test installation
qiime --help

#Set working directory with cd

#Check mapping file: 
qiime metadata tabulate --m-input-file sky2018map16s.txt --o-visualization sky2018map16s.qzv

#For this project we have two mapping files. One for 16s and one for RBCL. We are going to
#start this project w/the RBCL map, run through the pipeline, then demultiplex data again
#with the 16s metadata file read in. 

#Use chrome browser and the website  https://view.qiime2.org to view the .qvz files or

qiime tools view sky2018map16s.qzv


#to move back a directory, add ../ before the output object. e.g. you could write "../seq.qza" or you could, before running the code, write cd ../
#you want to specify a file path (instead of just the working directory) for qiime tool imports because the argument works on everything in a directory, 
#so just setting the working directory will make it work on everything inside of it

qiime tools import --type EMPPairedEndSequences --input-path /Users/melittology/Documents/ffarsunflower/parsed_barcodes --output-path seqs.qza

---- We started this pipeline here because the artifacts described above were already created when we ran through the RBCL code



-------

#DEMULTIPLEX

# This step splits data into individual, per-sample fastq files.


qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file ../sky2018map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza --p-no-golay-error-correction --output-dir error


#to visualize results:
qiime demux summarize --i-data demux16s.qza --o-visualization demux16s.qzv

#to see in html window:  https://view.qiime2.org/ 
# OR run the command

qiime tools view demux16s.qzv


#SEQUENCE QUALITY CONTROL AND MAKE FEATURE TABLE

#The following steps are part of the DADA2 pipeline for detecting and correcting amplicon sequence data.

#this step denoises paired-end sequences, dereplicates them, and filters chimeras. at this point, phix is already removed from demultiplexing
#Examine the "Interactive Quality Plot" tab in the demux16s.qzv file to determine which parameters to pass in
#--p-trunc-len-f and --p-trunc-len-r will truncate each sequence at position f or r and remove low quality regions of the sequences
#--p-trim-left-m will remove the firm m bases of each sequence

#In our specific case, truncate F reads at 266 and R reads at 288 for the first pass. 
#our first set of bases are high quality so dont need to use trim-left

#outputs: you will get a zipped files of fasta files that have one sequence for each asv. this is your rep-seqdada2.qza file
#you would use this repseq file for performing taxonomy. you could also use it to measure unifrac distances if you decide to perform phylogenetics
#the other output is an asv feature table, with asvs as rows and samples as columns with counts
#note: in the newer versions of qiime2, need to add an output for denoising stats "--o-denoising-stats stats-dada2.qza" 
#FOR RBCL, expect quality drop off around 110ish, for 16s, drop off around 260ish

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demux16s.qza  \
--p-trunc-len-f 266  \
--p-trunc-len-r 288  \
--p-trim-left-f 0  \
--p-n-threads 10  \
--output-dir dada2-16s.output  \
 --o-representative-sequences rep-seqs-dada2-16s.qza  \
 --o-table table16s.qza 
 
#VISUALIZE FEATURE TABLE
#then make a “feature table” visualization to explore the data
#the "feature-table summarize" commands tells you how many sequences are associated with each sample and feature
# the "feature-table tabulate-seqs" command gives a map of feature IDs to sequences and provides links to BLAST each sequence

qiime feature-table summarize --i-table dada2-16s.output/table16s.qza --o-visualization dada2-16s.output/table16s.qzv
qiime tools view dada2-16s.output/table16s.qzv

#And repseqs visualization

qiime feature-table tabulate-seqs --i-data dada2-16s.output/rep-seqs-dada2-16s.qza --o-visualization dada2-16sL.output/rep-seqs-dada2-16s.qzv
qiime tools view dada2-16sL.output/rep-seqs-dada2-16s.qzv


#ASSIGNING TAXONOMOY

#We need reference sequences and their taxonomic classifications. 
#Use a information-rich database that is
clustered at 99% sequence similarity at least (In our case, using Silva)

#I downloaded the newest silva 132 database into a new working directory from 
https://www.arb-silva.de/download/archive/qiime
# and used 99_otus_16S.fasta and 
consensus_taxonomy_7_levels.txt to create the training set. 

mkdir 16s-trainingclassifier
cd 16s-trainingclassifier


wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
unzip Silva_132_release.zip
rm Silva_132_release.zip


cd SILVA_132_QIIME_release

###Train feature classifier: 

#import reference sequences from silva data as a qiime2 artifact
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path 99_otus_16S.qza

#import taxonomy strings. Check and see if your taxonomy file is a tab-seperated file without a header. 
#if it doesnt have a header, specify "headerlessTSVTaxonomyFormat" since the default source formats usually have headers

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
--output-path 99_otus_16S_taxonomy.qza

#We now have two Qiime2 artifacts, 99_otus_16s.qza (reference sequences) and 99_otus_16s_taxonomy.qza (taxonomic names)

#trim silva to my region using my sequencing primers. We tell the algorithm our genomic primer forward and reverse sequences
#we do this because taxonomic classification is more accurate when a naive bayes classifier is trained on the region 
#of the 16s sequence that we sequenced (Werner et al. 2012). 

qiime feature-classifier extract-reads --i-sequences 99_otus_16S.qza --p-f-primer CMGGATTAGATACCCKGG --p-r-primer AGGGTTGCGCTCGTTG --o-reads ref-seqs16s.qza


#visualize: 
qiime feature-table tabulate-seqs --i-data ref-seqs16s.qza --o-visualization ref-seqs16s.qzv

qiime tools view ref-seqs.qzv (but this is big and took forever to load...).

#Train the classifier:
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs16s.qza --i-reference-taxonomy 99_otus_16S_taxonomy.qza  --o-classifier classifier16s.qza

#classify rep seqs and put the resulting taxonomic ids into the training classifier folder
cd ../ until you get into your main folder 

qiime feature-classifier classify-sklearn --i-classifier 16s-trainingclassifier/SILVA_132_QIIME_release/classifier16s.qza --i-reads  dada2-16s.output/rep-seqs-dada2-16s.qza --o-classification  16s-trainingclassifier/taxonomy16s.qza

#cd into 16s-training classifier

#And make a file that I can view…
qiime metadata tabulate --m-input-file taxonomy16s.qza --o-visualization taxonomy16s.qzv

# Visualize taxonomy:
cd ../ until you get to the right directory

qiime taxa barplot --i-table dada2-16s/table16s.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file ffar2018map16s.txt --o-visualization dada2-16s.output/taxa-bar-plots.qzv


#FILTERING STEPS. Go through and filter and THEN subsample and THEN remake trees. Go back and fix the order of things. 

#filter 1: out the chloroplast and mitochondria reads for now, I will then look at the blanks and decide whether to filter out Halomonas and Shewanella reads as well. 

qiime taxa filter-table --i-table dada2-16s.output/table16s.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --p-exclude mitochondria,chloroplast --o-filtered-table dada2-16s.output/tablefilt1.qza

qiime feature-table summarize --i-table dada2-16s.output/tablefilt1.qza --o-visualization dada2-16s.output/tablefilt1.qzv --m-sample-metadata-file ffar2018map16s.txt 


#filter 2: remove sequences only found in one sample 

qiime feature-table filter-features --i-table dada2-16s.output/tablefilt1.qza --p-min-samples 2 --o-filtered-table dada2-16s.output/tablefilt2.qza

qiime feature-table summarize --i-table dada2-16s.output/tablefilt2.qza --o-visualization dada2-16s.output/tablefilt2.qzv 

#filter 3/4: look at controls and remove the bacteria that are in them. These bacteria are often salt-loving and thats why theyre found in buffers. they include 
#Halomonas, Shewanella, Oceanospirillales, and the acne bacteria Propionibacterium

#filter out whole taxonimc groups that don't need to be in there (e.g. Propionibacterium)
#filter out specific strains that came up in our controls. for example, there is an uncultured species of lactobacillus that is likely showing up
#in a plate of our samples and in the dna control. remove it. (need to check)
#filter out samples that are control samples so that we don't include them in analyses


#To look at contaminants: we made a copy of our metadata map and added in a column that indicated if the sample was not a control, a dna control, or a field control
#then we visualized our taxonomic IDs with this file as the metadata file
#We also looked at the filtered taxonomy file as a csv and looked at what bacteria showed up our control samples. we then, for each one of these
#bacteria, looked up if it shows up in all the samples or just that control and made a call about whether or not to remove it
#if bacteria is present in just the control sample, you want to remove it
#however, some are also present in a lot of other samples! we don't want to lose data. if bacteria is in the control AND in more than 30 samples just leave it and don't filter it out (10% of samples)
#you may want to make an exception if the bacterial contaminant is obviously present in one just one plate, because even though its in a lot of samples its likely a contaminant
#another exception is if you have the contaminant in a lot of samples BUT also in a lot of the controls, get rid of it
#We made a removal list called "contaminants.txt" which has the bacteria we want to remove

#run this code and visualize the taxa barplot results. drop down the "taxonomic level" menu to level 7 and then "download csv" 
#use excel to check out and summarize the bacterial reads that are your control samples 
qiime taxa barplot --i-table dada2-16s.output/table16s.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file ffar2018ctrlmap16s.txt --o-visualization dada2-16s.output/taxa-bar-plotsCTRl.qzv

#filter steps

qiime taxa filter-table --i-table dada2-16s.output/tablefilt2.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --p-exclude lautropia,"Curtobacterium;Ambiguous_taxa",rothia,haemophilus,ruminococcaceae,"Corynebacterium 1;Ambiguous_taxa",desulfuromonadales,"Betaproteobacteriales;D_4__Burkholderiaceae;__;__",devosia,"Hymenobacter;Ambiguous_taxa","uncultured Lactobacillus",sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table dada2-16s.output/tablefilt3.qza

#not sure this next step works lol but try it
qiime feature-table filter-samples --i-table dada2-16s.output/tablefilt3.qza --m-metadata-file dada2-16s.output/bacterialcontaminants.txt --p-where "contaminant" --p-exclude-ids --o-filtered-table dada2-16s.output/tablefilt4.qza

#check to make sure your contaminants actually filtered out

qiime taxa barplot --i-table dada2-16s.output/tablefilt4.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file ffar2018ctrlmap16s.txt --o-visualization dada2-16s.output/taxa-bar-plotsnocontam.qzv


#filter 5: remove control samples. make a txt file with the list of IDs for samples you want to keep, where you remove control IDs. called "samplestokeep.txt"

qiime feature-table filter-samples --i-table dada2-16s.output/tablefilt4.qza --m-metadata-file dada2-16s.output/samplestokeep.txt --o-filtered-table dada2-16s.output/tablefilt5.qza

#check it worked
qiime taxa barplot --i-table dada2-16s.output/tablefilt5.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file ffar2018ctrlmap16s.txt --o-visualization dada2-16s.output/taxa-bar-plotsNoCtrl.qzv


#We decided to filter our rep seqs file so that you can see what samples you have left after filtering and subsampling

qiime feature-table filter-seqs --i-data dada2-16s.output/rep-seqs-dada2-16s.qza --i-table dada2-16s.output/tablefilt5.qza --o-filtered-data dada2-16s.output/rep-seqs-dada2-16s-filtered.qza

qiime feature-table tabulate-seqs --i-data dada2-16s.output/rep-seqs-dada2-16s-filtered.qza --o-visualization dada2-16s.output/rep-seqs-dada2-16s-filtered.qzv
qiime tools view dada2-16s.output/rep-seqs-dada2-16s-filtered.qzv


#DETERMINE SUBSAMPLE DEPTH - CHECK THIS TO MAKE SURE WE STILL WANT THESE DEPTHS

I want to make a rarefaction curve to see how subsampling depths influence my alpha diversity metrics. 

qiime diversity alpha-rarefaction \
--i-table dada2-16s.output/tablefilt5.qza \
--i-phylogeny dada2-16s.output/rooted-tree16s.qza \
--p-max-depth 14000 \
--m-metadata-file ffar2018map16s.txt \
--o-visualization dada2-16s.output/alphararefaction16s.qzv

#Open visualization in Qiime2 View

#For subsampling I will use 3,100 reads based on where the rarefaction curve peters off

qiime feature-table summarize --i-table dada2-16s.output/tablefilt5.qza --o-visualization dada2-16s.output/tablefilt5.qzv


When I visualize the table (go to "interactive sample detail", I can see this allows me to keep 250 samples (85.62%)

- STOPPED HERE -------

#GENERATE TREE FOR PHYLOGENETIC DIVERSITY ANALYSES

#Now alignment of the reads using MAFFT.

qiime alignment mafft --i-sequences dada2-16s.output/rep-seqs-dada2-16s.qza --o-alignment dada2-16s.output/aligned_repseqs16s.qza

#Then filter out the unconserved, highly gapped columns from alignment
qiime alignment mask --i-alignment dada2-16s.output/aligned_repseqs16s.qza --o-masked-alignment dada2-16s.output/masked_aligned_repseqs16s.qza

#And make a tree with Fasttree which creates a maximum likelihood phylogenetic trees from aligned sequences 
#for more information on fastree, http://www.microbesonline.org/fasttree/
qiime phylogeny fasttree --i-alignment dada2-16s.output/masked_aligned_repseqs16s.qza --o-tree dada2-16s.output/unrooted_tree16s.qza

#now root it
qiime phylogeny midpoint-root --i-tree dada2-16s.output/unrooted_tree16s.qza --o-rooted-tree dada2-16s.output/rooted-tree16s.qza







#STOPPED HERE. Make sure to use the subsampling depth that we've determined in any future analyses. 
#Before continuing, blast Lactobacillus kunkeii against NCBI database. The silva database does not have Lactobacillus michnerii, but that is likely what this bacteria is. 


-------------QUINNS ANALYSES NOTES BELOW



qiime diversity core-metrics-phylogenetic --i-phylogeny dada2.output/rooted_tree.qza --i-table dada2.output/tablefilt3.qza --p-sampling-depth 5212 --m-metadata-file nosemamapflipbarcode.txt --output-dir core_div_metrics_filt3
qiime diversity beta-group-significance --i-distance-matrix core_div_metrics_filt3/weighted_unifrac_distance_matrix.qza --m-metadata-file nosemamapflipbarcode.txt --m-metadata-category IndResist --p-pairwise --o-visualization core_div_metrics_filt3/wuniIndResistSig.qzv




qiime diversity core-metrics-phylogenetic --i-phylogeny dada2.output/rooted_tree.qza --i-table dada2.output/tablefilt1.qza --p-sampling-depth 10072 --m-metadata-file nosemamapflipbarcode.txt --output-dir core_div_metrics

qiime diversity beta-group-significance --i-distance-matrix core_div_metrics/weighted_unifrac_distance_matrix.qza --m-metadata-file nosemamapflipbarcode.txt --m-metadata-category Treatment --p-pairwise --o-visualization core_div_metrics/weighteduniTrt_group_significance.qzv

Testing between individual bees resistant to Nosema infection: 
qiime diversity beta-group-significance --i-distance-matrix core_div_metrics_2/weighted_unifrac_distance_matrix.qza --m-metadata-file nosemamapflipbarcode.txt --m-metadata-category IndResist --p-pairwise --o-visualization core_div_metrics_2/wuniIndResistSig.qzv
This was significant: PERMANOVA test statistic	pseudo-F sample size	159 number of groups	2 test statistic	3.12899 p-value	0.008 number of permutations	999

Now filter Shewanella, Halomonas, and Propionbacterium:
qiime taxa filter-table --i-table dada2.output/tablefilt1.qza --i-taxonomy dada2.output/taxonomy.qza --p-exclude Shewanella,Oceanospirillales,Propionibacterium --o-filtered-table dada2.output/tablefilt2.qza
qiime feature-table summarize --i-table dada2.output/tablefilt2.qza --o-visualization dada2.output/tablefilt2.qzv
qiime tools view dada2.output/tablefilt2.qzv 
sample to 5,212 for diversity will keep 196 samples (86%).
qiime diversity core-metrics-phylogenetic --i-phylogeny dada2.output/rooted_tree.qza --i-table dada2.output/tablefilt2.qza --p-sampling-depth 5212 --m-metadata-file nosemamapflipbarcode.txt --output-dir core_div_metrics_filt2
qiime diversity beta-group-significance --i-distance-matrix core_div_metrics_filt2/weighted_unifrac_distance_matrix.qza --m-metadata-file nosemamapflipbarcode.txt --m-metadata-category IndResist --p-pairwise --o-visualization core_div_metrics_filt2/wuniIndResistSig.qzv


One last filter, to be sure that including samples without qPCR data does not mess up the analysis: 
qiime feature-table filter-samples --i-table dada2.output/tablefilt3.qza --m-metadata-file nosemamapflipbarcode.txt  --p-where "IndResist IN ('0','1')" --o-filtered-table dada2.output/tablefilt4.qza
qiime feature-table summarize --i-table dada2.output/tablefilt4.qza --o-visualization dada2.output/tablefilt4.qzv
subsample at 5212 to get 192 (92.31%)
Results are exactly the same, no need to worry...
qiime diversity beta-group-significance --i-distance-matrix weighted_unifrac_distance_matrix.qza --m-metadata-file ../nosemamapflipbarcode.txt --m-metadata-category IndResist --p-pairwise --o-visualization wuniIndResistSig.qzv
qiime diversity beta-group-significance --i-distance-matrix weighted_unifrac_distance_matrix.qza --m-metadata-file ../nosemamapflipbarcode.txt --m-metadata-category ColResist --p-pairwise --o-visualization wuniColResistSig.qzv
qiime diversity beta-group-significance --i-distance-matrix we

