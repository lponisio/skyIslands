
#1: On whatever computer will be running the tasks, innitiate a docker container for qiime1

docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline sglim2/qiime-1.9.1


#2: cd into the correct folder within the container, verify that the files you need are there

cd mnt/SI_pipeline/16sRBCL
ls


#3: If your reads are in *fasta.gz format, unzip them with qiime, then
#rename them "forward.fastq" and "reverse.fastq". Quinn designed the
#barcodes in reverse so flowcell1 is the reverse and flowcell2 is the
#forward

gunzip *.gz

mv flowcellpair2.fastq rawforward.fastq
mv flowcellpair1.fastq rawreverse.fastq


#4: parse the barcodes in the files, putting our data into a format qiime2 will be able to use

extract_barcodes.py -f rawforward.fastq -r rawreverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#5: re-zip the output files

cd parsed_barcodes
gzip *.fastq

mv reads1.fastq.gz forward.fastq.gz
mv reads2.fastq.gz reverse.fastq.gz

#5: Exit Qiime1 and use docker to open the environment for Qiime 2. 

exit 
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

#6: Test that the container for Qiime 2 is properly associated, then make sure you are in the root directory using ls, then set working directory to the mounted volume
 
qiime --help
ls

#7: Import your parser barcodes into an object you can demultiplex in Qiime 2. Make sure working directory is set correctly.
qiime tools import --type EMPPairedEndSequences --input-path /mnt/SI_pipeline/16sRBCL/parsed_barcodes/ --output-path /mnt/SI_pipeline/16sRBCL/seqs.qza

#8: Examine your mapping file, which is metadata you create associated with the project. Use qiime to make it into a qzv object.
#If you are examining multiple amplicon types, pick a map associated with one to start with (e.g. 16s)

qiime metadata tabulate --m-input-file sky2018map16s.txt --o-visualization sky2018map16s.qzv
qiime tools view sky2018map16s.qzv

#9: Demultiplex 16s reads. Only works in version Qiime2 2019.1
exit 
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1

cd ../mnt/SI_pipeline/16sRBCL

# qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file sky2018map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza --p-golay-error-correction FALSE --output-dir error


qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file sky2018map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza 


#9a: Visualize Results

qiime demux summarize --i-data demux16s.qza --o-visualization demux16s.qzv
qiime tools view demux16s.qzv

## when interpretating the quality boxes, you can use the bottom of
## the black box as a conservative measure for the phred score (not the
## whiskers and not the middleo f the box)

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demux16s.qza  \
--p-trunc-len-f 180  \
--p-trunc-len-r 220  \
--p-trim-left-f 0  \
--p-n-threads 2  \
--output-dir mnt/SI_pipeline/16sRBCL/dada2-16s  \
 --o-representative-sequences rep-seqs-dada2-16s.qza  \
 --o-table table16s.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv

qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## repeate the steps for RBCL

qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file sky2018mapRBCL.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demuxRBCL.qza 

#9a: Visualize Results

qiime demux summarize --i-data demuxRBCL.qza --o-visualization demuxRBCL.qzv
qiime tools view demuxRBCL.qzv

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demuxRBCL.qza  \
--p-trunc-len-f 180  \
--p-trunc-len-r 218  \
--p-trim-left-f 0  \
--p-n-threads 2  \
--output-dir /mnt/SI_pipeline/16sRBCLdada2-RBCLs  \
 --o-representative-sequences rep-seqs-dada2-RBCL.qza  \
 --o-table tableRBCL.qza

qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCL.qzv

qiime feature-table summarize --i-table dada2-RBCL/tableRBCL.qza --o-visualization dada2-RBCL/tableRBCL.qzv

#check outputs to make sure you didn't lose too many samples. If you did, you may want to retruncate.

## *****************************************************************************
## 16s
## *****************************************************************************

#10. ASSIGN TAXONOMY 16s

#We need reference sequences and their taxonomic classifications.
#Use a information-rich database that is clustered at 99% sequence similarity at least (In our case, using Silva for 16s) 
#we have to "train" the classifier dataset just once.

#10a. Download the newest silva 132 database into a new working directory from https://www.arb-silva.de/download/archive/qiime

# use 99_otus_16S.fasta and  consensus_taxonomy_7_levels.txt to create the training set.
mkdir 16s-trainingclassifier
cd 16s-trainingclassifier
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
unzip Silva_132_release.zip
rm SILVA_132_release.zip
#10b. Train feature classifier
#import reference sequences from silva data as a qiime2 artifact

## go back into qiime
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1
cd ../mnt/SI_pipeline/16sRBCL/16s-trainingclassifier/SILVA_132_QIIME_release


## 99 is 99% match between our seq and the database
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path rep_set/rep_set_16S_only/99/SILVA_132_99_16S.fna \
--output-path 99_otus_16S.qza

#import taxonomy strings. Check and see if your taxonomy file is a tab-seperated file without a header.
#if it doesnt have a header, specify "headerlessTSVTaxonomyFormat" since the default source formats usually have headers

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
--output-path 99_otus_16S_taxonomy.qza

#We now have two Qiime2 artifacts, 99_otus_16s.qza (reference sequences) and 99_otus_16s_taxonomy.qza (taxonomic names). 
#trim silva to my region using my sequencing primers. We tell the algorithm our genomic primer forward and reverse sequences
#we do this because taxonomic classification is more accurate when a naive bayes classifier is trained on the region
#of the 16s sequence that we sequenced (Werner et al. 2012).
qiime feature-classifier extract-reads --i-sequences 99_otus_16S.qza --p-f-primer CMGGATTAGATACCCKGG --p-r-primer AGGGTTGCGCTCGTTG --o-reads ref-seqs16s.qza
#visualize:
qiime feature-table tabulate-seqs --i-data ref-seqs16s.qza --o-visualization ref-seqs16s.qzv
qiime tools view ref-seqs16s.qzv # (but this is big and took forever to load...).

## may need to clean up docker memory usage
docker system prune

#Train the classifier:
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs16s.qza --i-reference-taxonomy 99_otus_16S_taxonomy.qza  --o-classifier classifier16s.qza

#10c. classify rep seqs and put the resulting taxonomic ids into the training classifier folder
# cd ../ until you get into your main folder

## may need to install older version of scikit learn
# pip install -U scikit-learn==0.19.1

cd ../../

qiime feature-classifier classify-sklearn --i-classifier 16s-trainingclassifier/SILVA_132_QIIME_release/classifier16s.qza --i-reads  dada2-16s/rep-seqs-dada2-16s.qza --o-classification  16s-trainingclassifier/taxonomy16s.qza

## switch to the newest version of qiime
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

# 10c visualize
cd ../mnt/SI_pipeline/16sRBCL/16s-trainingclassifier/

qiime metadata tabulate --m-input-file taxonomy16s.qza --o-visualization taxonomy16s.qzv

cd ../ 
qiime taxa barplot --i-table dada2-16s/table16s.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file sky2018map16s.txt --o-visualization dada2-16s/taxa-bar-plots.qzv


# 11 FILTERING STEPS. Go through and filter out unwanted reads. THEN subsample and THEN remake trees.  

# HAVENT DONT THIS YET, BUT NEED TO DISCUSS WHY WE DON'T HAVE A TON OF SPECIFICITY

# 11a. Prefilter step: The silva database does not have Lactobacillus michnerii, so our reads for Lactobacillus kunkeei are likely wrong. 

# first visualize taxonomy16s.qzv. from the visualizer, download taxonomy16s.tsv. 
# make a copy of taxonomy16s.tsv and rename it taxonomy16sfixed.tsv - we will make corrections here
# then visualize repseqs16sfiltered.qzv in qiime2view 
qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv

# using these two documents, select the sequence corresponding with the feature ids for L. kunkeei


# Make correctionss to taxonomy16sfixed.tsv. convert it to txt manually on your computer.
# Convert this text file back into a qza object using qiime

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy16s.qza/taxonomy16sfixed.txt \
--output-path taxonomy16s.qza/taxonomy16sfixed.qza

#DID THIS
#11b. filter 1: out the chloroplast and mitochondria reads 

qiime taxa filter-table --i-table dada2-16s/table16s.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --p-exclude mitochondria,chloroplast --o-filtered-table dada2-16s/tablefilt1.qza

qiime feature-table summarize --i-table dada2-16s/tablefilt1.qza --o-visualization dada2-16s/tablefilt1.qzv --m-sample-metadata-file sky2018map16s.txt

#DID THIS
#11c. filter 2: remove sequences only found in one sample 

qiime feature-table filter-features --i-table dada2-16s/tablefilt1.qza --p-min-samples 2 --o-filtered-table dada2-16s/tablefilt2.qza

qiime feature-table summarize --i-table dada2-16s/tablefilt2.qza --o-visualization dada2-16s/tablefilt2.qzv  --m-sample-metadata-file sky2018map16s.txt 

#HAVENT DONT THIS YET BUT FOUND WHAT WE NEED TO FILTER
#11d. filter 3: look at DNA and Illumina PCR controls and remove the bacteria that are in them. 
#We don't want to remove any contaminant bacteria that are real bee bacterial, so make decisions based on what you know are contaminants 
#there are known bacteria are salt-loving and often found in buffers. they include 
#Halomonas, Shewanella, Oceanospirillales, and the acne bacteria Propionibacterium. 
#visualize taxabarplot.qzv (level 7) and download the .csv from the visualization to see whats in your controls. 

#filtering rationale
# look at what bacteria show up our control samples. we then, for each one of these
#bacteria, look up if it is in all the samples or just that control and made a call about whether or not to remove it
#if bacteria is present in just the control sample, you want to remove it
#h owever, some are also present in a lot of other samples! we don't want to lose data. if bacteria is in the control AND in more than 30 samples just leave it and don't filter it out (10% of samples)
# you may want to make an exception if the bacterial contaminant is obviously present in one just one plate, because even though its in a lot of samples its likely a contaminant
# another exception is if you have the contaminant in a lot of samples BUT also in a lot of the controls, get rid of it

## *****************************************************************************
## move into R to choose bacteria in controls to filter out
## *****************************************************************************
## load in bacteria csv
bact <- read.csv("~/Desktop/taxa-bar-plots.csv")
control.names <- c("DNActrl1", "DNActrl2", "16SIlluminaCtrl1", "16SIlluminaCtrl2", "16SIlluminaCtrl3")

primer.names <- c("forwardbarcode", "revbarcode", "barcodesequence", "forwardgenomicprimer", "revgenomicprimer")

control.bact <- bact[bact$index %in% control.names, !colnames(bact) %in% primer.names]
rownames(control.bact) <- control.bact$index
control.bact$index <- NULL
control.bact <- as.matrix(control.bact)
bact.count <- colSums(control.bact)
bacteria.in.controls <- names(bact.count[bact.count > 0])

sample.bact <- bact[!bact$index %in% control.names, !colnames(bact) %in% primer.names]
rownames(sample.bact) <- sample.bact$index
sample.bact$index <- NULL
sample.bact <- as.matrix(sample.bact)
bact.count.samples <- apply(sample.bact, 2, function(x) sum(x > 0))
hist(bact.count.samples, breaks=100)
bacteria.in.samples <- names(bact.count.samples[bact.count.samples > 0])

bacteria.in.controls[bacteria.in.controls %in% bacteria.in.samples]

#we want to remove the following:
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Frankiales;D_4__Nakamurellaceae;D_5__Nakamurella;D_6__uncultured Nakamurellaceae bacterium
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Kineosporiales;D_4__Kineosporiaceae;D_5__Kineococcus;__
REVISIT D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Kineosporiales;D_4__Kineosporiaceae;D_5__Kineosporia;D_6__uncultured bacterium
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Micrococcales;__;__;__
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Propionibacteriales;D_4__Propionibacteriaceae;D_5__Cutibacterium;__
REVISIT D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;__;__;__;__
REVISIT D_0__Bacteria;D_1__Actinobacteria;D_2__Thermoleophilia;D_3__Gaiellales;D_4__uncultured;__;__
REVISIT D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Chitinophagales;D_4__Chitinophagaceae;D_5__Segetibacter;__
D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Cytophagales;D_4__Hymenobacteraceae;D_5__Hymenobacter;__
D_0__Bacteria;D_1__Chloroflexi;D_2__Ktedonobacteria;D_3__Ktedonobacterales;D_4__Ktedonobacteraceae;D_5__1959-1;D_6__uncultured bacterium
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Staphylococcaceae;D_5__Staphylococcus;__
REVISIT D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Enterococcaceae;D_5__Enterococcus;__
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Leuconostocaceae;D_5__Leuconostoc;__
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Lactococcus;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Acetobacterales;D_4__Acetobacteraceae;__;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Beijerinckiaceae;D_5__Methylobacterium;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Beijerinckiaceae;D_5__Microvirga;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Rhizobiaceae;D_5__Aureimonas;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhodobacterales;D_4__Rhodobacteraceae;D_5__Paracoccus;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Sphingomonadales;D_4__Sphingomonadaceae;D_5__Sphingomonas;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Sphingomonadales;D_4__Sphingomonadaceae;D_5__uncultured;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Pseudomonadales;D_4__Pseudomonadaceae;D_5__Pseudomonas;__
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Micrococcales;D_4__Microbacteriaceae;D_5__Galbitalea;__
D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Propionibacteriales;D_4__Nocardioidaceae;D_5__Nocardioides;D_6__Nocardioides sp. kmb27
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;__;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhodobacterales;D_4__Rhodobacteraceae;D_5__Rubellimicrobium;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Sphingomonadales;D_4__Sphingomonadaceae;D_5__Novosphingobium;__
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Sphingomonadales;D_4__Sphingomonadaceae;__;__
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;D_5__Oceanobacillus;__

# Get out everything under "D_3__Oceanospirillales"

# Make sure the other filter steps got these out
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rickettsiales;D_4__Mitochondria;__;__
D_0__Bacteria;D_1__Cyanobacteria;D_2__Oxyphotobacteria;D_3__Chloroplast;__;__;__

# NOT DONE YET: make a removal list called "contaminants.txt" which has the bacteria we want to remove



# NOT DONE YET
# 11d. fitler 4: remove samples that are control samples (DNA controls and Illumina PCR controls) so that we don't include them in analyses
# make a txt file with the list of IDs for samples you want to keep, where you remove control IDs. called "samplestokeep.txt"

qiime feature-table filter-samples --i-table dada2-16s/tablefilt3.qza --m-metadata-file dada2-16s/samplestokeep.txt --o-filtered-table dada2-16s/tablefilt4.qza

#check it worked
qiime taxa barplot --i-table dada2-16s/tablefilt4.qza --i-taxonomy 16s-trainingclassifier/taxonomy16s.qza --m-metadata-file ffar2018ctrlmap16s.txt --o-visualization dada2-16s/taxa-bar-plotsNoCtrl.qzv

#Filter rep seqs file so that you can see what samples you have left after filtering and subsampling

qiime feature-table filter-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --i-table dada2-16s/tablefilt4.qza --o-filtered-data dada2-16s/rep-seqs-dada2-16s-filtered.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s-filtered.qza --o-visualization dada2-16s/rep-seqs-dada2-16s-filtered.qzv
qiime tools view dada2-16s/rep-seqs-dada2-16s-filtered.qzv



NOT DONE YET, CODE ISNT FIXED. DO WE SUBSAMPLE OR MAKE PHYLOGENY FIRST

# 12. DETERMINE SUBSAMPLE DEPTH - CHECK THIS TO MAKE SURE WE STILL WANT THESE DEPTHS

# make a rarefaction curve to see how subsampling depths influence my alpha diversity metrics. 

qiime diversity alpha-rarefaction \
--i-table dada2-16s/tablefilt4.qza \
--i-phylogeny dada2-16s/rooted-tree16s.qza \
--p-max-depth 14000 \
--m-metadata-file ffar2018map16s.txt \
--o-visualization dada2-16s/alphararefaction16s.qzv

#Open visualization in Qiime2 View

#For subsampling I will use 3,100 reads based on where the rarefaction curve peters off

qiime feature-table summarize --i-table dada2-16s.output/tablefilt5.qza --o-visualization dada2-16s.output/tablefilt5.qzv

When I visualize the table (go to "interactive sample detail", I can see this allows me to keep X NUMBER OF samples (X%)

#13 GENERATE TREE FOR PHYLOGENETIC DIVERSITY ANALYSES

# align of the reads using MAFFT.

qiime alignment mafft --i-sequences dada2-16s/rep-seqs-dada2-16s.qza --o-alignment dada2-16s/aligned_repseqs16s.qza

#Then filter out the unconserved, highly gapped columns from alignment
qiime alignment mask --i-alignment dada2-16/aligned_repseqs16s.qza --o-masked-alignment dada2-16s/masked_aligned_repseqs16s.qza

#And make a tree with Fasttree which creates a maximum likelihood phylogenetic trees from aligned sequences 
#for more information on fastree, http://www.microbesonline.org/fasttree/
qiime phylogeny fasttree --i-alignment dada2-16s/masked_aligned_repseqs16s.qza --o-tree dada2-16s/unrooted_tree16s.qza

#now root it
qiime phylogeny midpoint-root --i-tree dada2-16s/unrooted_tree16s.qza --o-rooted-tree dada2-16s/rooted-tree16s.qza










## *****************************************************************************
## RBCL
## *****************************************************************************

## 12 classify RBCL


#The next steps are to classify taxonomy with a trained RBCL classifier. But Qiime2 does not have a rbcl classifier for download. 

#We will use the Bell et al. RBCL RDP classifier, described here:
#https://github.com/KarenBell/rbcL-dual-index-metabarcoding

#This tutorial describes how to train and apply their RBCL classifier. However, we won't use their entire pipeline because their script processes raw files, then joins them, demultiplexes them, etc, and only assigns taxonomy at the end. All we need to do is use the classifier they developed. Follow the script below, which we modified. 

#we downloaded the classifier, which is already trained, from https://gitlab.umiacs.umd.edu/derek/qiime/tree/master/rdp_classifier_2.2

#go to Bell and Brosi's rbcl reference library: https://figshare.com/collections/rbcL_reference_library/3466311/1

#and download "rbcl_utax_trained.zip" and "rbcl_rdp_trained_reference_database" into a folder called RBCLclassifierRDP in your local working directory. unzip the files.

mkdir RBCLclassifierRDP
cd ../RBCLclassifierRDP

unzip qiime-master-rdp_classifier_2.2.zip

unzip rbcL_rdp_trained_reference_database.zip

unzip rbcL_utax_trained.zip

#we will work on the rep-seqs object which is a qiime artifact currently, so export it to a fasta file. first cd back to the main ffarsunflower folder. NOTE: use syntax below; online docs for qiime tools export have different syntax that doesnt work. 

docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

cd ../mnt/SI_pipeline/16sRBCL/

qiime tools export --input-path dada2-RBCL/rep-seqs-dada2-RBCL.qza --output-path dada2-RBCL
  
#Take the file that was generated "dna sequences.fasta" and rename it "rep-seqs-dada2RBCL.fasta" and move it where you want in your directory

mv dna-sequences.fasta rep-seqs-dada2RBCL.fasta

exit


#CLASSIFY TAXONOMY WITH RDP

#First make sure you have java installed
#exit qiime2 with "exit" command

# We ran two different scripts to generate two different taxonomic files. The first is a fixrank file. the second is an allrank file (the default output). See this link for description: https://github.com/rdpstaff/classifier

#NOTE: you may need to download a java developer kit (JDK) from Oracle for the java command line stuff to work

cd ../ 

#this one gives taxonomy to species (we ended up only using this one after all)
java -Xmx2g -jar RBCLclassifierRDP/qiime-master-rdp_classifier_2.2/rdp_classifier_2.2/rdp_classifier-2.2.jar classify -t RBCLclassifierRDP/rbcL_rdp_trained/rbcL.properties -o RBCLclassifierRDP/rbcl_classified_rdp.txt -q dada2-RBCL/rep-seqs-dada2RBCL.fasta


#Our resulting taxonomic IDs are found in two sheets in the RBCL classifier folder

#CLASSIFY TAXONOMY WITH BLAST / NCBI

#Let's compare sequences we get when using the RDP Classifier and BLAST.  We are going to use BioPerl

#to blast against NCBI, you need to get a package installer that has python arguments and stuff from

## install package installer for your computer from 
## ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


mkdir RBCLclassifierNCBI

#We used the .dmg file (second to last one) for a Mac
#We put it in "RBCLclassifierNCBI"

#When running BLAST against a local database downloaded to your computer, it is storage intensive to use the entire NCBI database. Instead, download just the nucleotide files you need.

#Go to NCBI and manually download RBCL files into your working directory. Go to https://www.ncbi.nlm.nih.gov/ and press download, search rbcl, press nucleotide under "Genomes". select what you want, use "send to" bottom at upper right corner and download as fasta. We downloaded plant nucleotide data and called it "rbcL_only_NCBI.fasta". Put it in "RBCL classifier NCBI" folder.

#then download gi_taxid_nucl.dmp from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
#expand it, put the dmp into the appropriate directory

#NOTE: you have to add the makeblastdb executable from the ncbi-blast-2.10.0+ bin to your $PATH! Just drag and drop it in Users/Ponisio/perl5, for example

#Use the fasta file to create a database that we'll then use in BioPerl
makeblastdb -in -taxid_map gi_taxid_nucl.dmp RBCLclassifierNCBI/rbcL_only_NCBI.fasta -out RBCLclassifierNCBI/rbcL_only_DB -parse_seqids -dbtype nucl


#if this all worked, you should get a series of files called "rbcL_only_DB" with some extension or numerical at the end. These are ALL your database, and one of the files (the .nal) is an alias file that lets everything talk to each other


## open NCBI docker containter, we ran cpan commands within the
## container and they installed correctly. We used find -name to find
## them and then updated the .pl run file with their location on by setting lib (line 2).


## the container includes blastn in usr

docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline ncbi/blast

cpan App::cpanminus
cpanm Bio::SearchIO
cpanm Bio::Tools::Run::StandAloneBlast --force 


cd mnt/SI_pipeline/16sRBCL/

perl RBCLclassifierNCBI/BlastToTableexistingdb.pl dada2-RBCL/rep-seqs-dada2RBCL.fasta RBCLclassifierNCBI/rbcL_only_DB > RBCLclassifierNCBI/rbcl_classified_NCBI.txt


#We looked at the RBCL taxonomic determinations generated by RDP and BLAST and compared the. We generated a set of rules to determine which classifier to use for each query. The rules we used were

#1	If genus matches on NCI and RDP, use that genus, but still check if its the region. if not, done use
#2	If genus does not match, look up samples in cal flora. If one sample is present in calflora and other is abset, use that genus
#3	If both samples are present in calflora, use the genus that occurs in Yolo County (based off Calflora maps and our veg data from our sites)
#4	if both are present in Yolo county, use NCBI unless the RDP classifier is over .75 confidence or if either is already present previously and confirmed
#5	If neither are present in Yolo county, leave blank, use family level ID if it's the same
# 8	There are a couple species that bees were caught on, but their pollen rbcL matched with something else as the top hit. What they had been caught on did match in NCBI though, with only a single snp difference causing the known species to not be the first hit (e.g., 176/179 bp matches for Ipomoea wrightii vs 175/179 for Convolvulus arvensis). We manually changed I. wrightii to Con arv, and Arctotheca calendula to H. annuus.   


#We saved our final determinations into our ffarunsflower folder as a .txt called "taxonomyRBCL.txt" inside "RBCLclassifierRDP" folder
#We need to make sure this file is in the right format for the next steps in qiime to process. Label the first column "Feature ID" and the second column "Taxon". Feature ID is the same as the query
#We removed all samples that didn't have an NCBI match and removed queries occurring in less than two samples


## ************************************************************************************
## move into R for pairing the databases
## ***********************************************************************************
rm(list=ls())
ncbi <- read.table("~/Dropbox/skyIslands_saved/data/raw/RBCL_16s/rbcl_classified_NCBI.txt",
		   sep="\t", header=TRUE)

rdp <- read.table("~/Dropbox/skyIslands_saved/data/raw/RBCL_16s/rbcl_classified_rdp.txt",
		 sep="\t", header=TRUE)

## rdp cleaning
rdp$clean_family <- gsub("f__", "", rdp$family)
rdp$clean_family <- sapply(strsplit(rdp$clean_family, "_"),
			  function(x) x[1])

rdp$clean_genus <- gsub("g__", "", rdp$genus)
rdp$clean_genus <- sapply(strsplit(rdp$clean_genus, "_"),
			  function(x) x[1])

rdp$clean_species <- gsub("s__", "", rdp$species)
rdp$GenusSpecies <- sapply(strsplit(rdp$clean_species, "_"),
			  function(x) x[1])

## ncbi cleaning

ncbi$genus <- sapply(strsplit(ncbi$name, "[ ]"),
		     function(x) x[1])

ncbi$species <- sapply(strsplit(ncbi$name, "[ ]"),
		       function(x) x[2])
ncbi$GenusSpecies <- paste(ncbi$genus, ncbi$species)


## plants IDed by hand at sites
plants <- read.csv("~/Dropbox/skyIslands_saved/data/raw/plants.csv")

veg.genera <- unique(plants$PlantGenus)
veg.genera <- veg.genera[!is.na(veg.genera)]

veg.sp <- unique(plants$GenusSpecies)
veg.sp <- veg.sp[!is.na(veg.sp)]

nrow(ncbi)
nrow(rdp)

ids <- data.frame(Sample=rdp$ID,
		  GenusSpeciesRDP=rdp$GenusSpecies,
		  GenusSpeciesRDPMatch=rdp$species_match,
		  GenusRDP=rdp$clean_genus,
		  GenusRDPMatch=rdp$genus_match)	       

ids$GenusSpeciesNCBI <- ncbi$GenusSpecies[match(ids$Sample,
						ncbi$Query)]
ids$GenusNCBI <- ncbi$genus[match(ids$Sample,
						ncbi$Query)]

ids$GenusRDP_NCBI_match <- ids$GenusNCBI == ids$GenusRDP
ids$SpeciesRDP_NCBI_match <- ids$GenusSpeciesNCBI == ids$GenusSpeciesRDP

ids$FinalGenusSpecies <- NA
ids$FinalGenus <- NA


## genera that match between RDP and NCBI
ids$FinalGenus[ids$GenusRDP_NCBI_match & !is.na(ids$GenusRDP_NCBI_match)]  <- ids$GenusRDP[ids$GenusRDP_NCBI_match & !is.na(ids$GenusRDP_NCBI_match)]

## and are also in the veg data IDed by hand
unique(ids$FinalGenus)
ids$FinalGenusInVeg <- NA
ids$FinalGenusInVeg[ids$FinalGenus %in% veg.genera] <- 1
ids$FinalGenusInVeg[!ids$FinalGenus %in% veg.genera & !is.na(ids$FinalGenus)] <- 0
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0])

## IDs that have a solid match but are not in the hand collected veg dats 
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0][ids$GenusRDPMatch > 0.8])
## "Swertia"      "Scrophularia" "Rubus"        "Ziziphus"   "Glycine"      "Garrya"    

## Garrya, Swertia, Rubus (blackberry?), Scrophularia are entierly possible
## https://www.fs.fed.us/database/feis/plants/shrub/garwri/all.html
## https://www.wildflower.org/plants/result.php?id_plant=SWPE
## https://plants.usda.gov/core/profile?symbol=SCLA

## Ziziphus is not impossible, it is a cultivar planted in NM and AZ
## https://aces.nmsu.edu/pubs/_h/H330/welcome.html

## Glycine, the genus of soybeans not native to the US is very unlikely.

## IDs that have an okay match but are not in the hand collected veg dats 
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0][ids$GenusRDPMatch >0.6 & ids$GenusRDPMatch < 0.8])

# "Cichorium"      "Astragalus"     "Tanacetum"     
# [5] "Alnus"          "Symphoricarpos" "Fragaria"       "Mentzelia"     
# [9] "Ziziphus"       "Veronica"       "Swertia"        "Glycine"       
#[13] "Quercus"        "Carex"          "Galega"         "Stenorrhynchos"
#[17] "Rubus"          "Gleditsia"      "Garrya"         "Strychnos"    # 
#[21] "Ascolepis"      "Cercis"         "Lotus"

## Cichorium is an invaisve species that likes road sides, present in AZ and NM so possible
## https://www.invasive.org/browse/subinfo.cfm?sub=5332

## Astragalus is very possible, "milkvetch/locoweed" https://en.wikipedia.org/wiki/Astragalus

## Tanacetum possible
## https://www.npsnm.org/pdfs/Asterhandout.pdf

## Alnus, alder possible
## https://plants.usda.gov/core/profile?symbol=ALNUS

## Symphoricarpos, snowberry possible
## https://www.wildflower.org/plants/result.php?id_plant=syal

## Fragaria, wild strawberry possible
## https://wnmu.edu/academic/nspages/gilaflora/fragaria_vesca.html

## Mentzelia possible
## https://www.swcoloradowildflowers.com/Yellow%20Enlarged%20Photo%20Pages/mentzelia.htm

## Veronica possible
## http://www.coloradowildbuds.com/wildflowers-by-color/blue-purple/veronica-species/

## Quercus, oak possible

## carex possible
## https://wnmu.edu/academic/nspages/gilaflora/carex_hystericina.html

## Stenorrhynchos, rare orchid found in alligator juniper
## https://www.fs.fed.us/rm/pubs/rmrs_p023/rmrs_p023_095_098.pdf

## Gleditsia, locust likely
## https://www.honey-plants.com/calendar/new-mexico/gleditsia-triacanthos/

## Lotus, likely

## Cercis, redbud ornamental possible.. 

## Galega UNLIKLEY, though invasive in NA not found in NM/AZ
## https://swbiodiversity.org/seinet/collections/map/googlemap.php?usethes=1&taxa=50789

## Strychnos UNLIKELY https://plants.usda.gov/java/nameSearch?mode=sciname&keywordquery=Strychnos

## Ascolepis NOT POSSIBLE, not found in NA


## species that match between RDP and NCBI
ids$FinalGenusSpecies[ids$SpeciesRDP_NCBI_match & !is.na(ids$SpeciesRDP_NCBI_match)]  <- ids$GenusSpeciesRDP[ids$SpeciesRDP_NCBI_match & !is.na(ids$SpeciesRDP_NCBI_match)]

## and are also in the veg data IDed by hand
unique(ids$FinalGenusSpecies)
ids$FinalGenusSpeciesInVeg <- NA
ids$FinalGenusSpeciesInVeg[ids$FinalGenusSpecies %in% veg.sp] <- 1
ids$FinalGenusSpeciesInVeg[!ids$FinalGenusSpecies %in% veg.sp & !is.na(ids$FinalGenusSpecies)] <- 0

unique(ids$FinalGenusSpecies[ids$FinalGenusSpeciesInVeg == 0])
unique(ids$FinalGenusSpecies[ids$FinalGenusSpeciesInVeg == 0][ids$GenusRDPMatch >0.8])

## "Vicia cirrhosa"        "Chlorophytum nimmonii"  "Rubus idaeus" ## none of these are at all likely

## where NCBI didn't have a match, if RDP is > 50%
ids$FinalGenus[is.na(ids$GenusRDP_NCBI_match) & ids$GenusRDPMatch > 0.5] <- ids$GenusRDP[is.na(ids$GenusRDP_NCBI_match) & ids$GenusRDPMatch > 0.5]
## this never happens
