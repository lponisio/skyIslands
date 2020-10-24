
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

#10. ASSIGN TAXONOMY

#We need reference sequences and their taxonomic classifications.
#Use a information-rich database that is clustered at 99% sequence similarity at least (In our case, using Silva for 16s) 
#we have to "train" the classifier dataset just once.

#10a. Download the newest silva 132 database into a new working directory from https://www.arb-silva.de/download/archive/qiime

# EXIT QIIME
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
