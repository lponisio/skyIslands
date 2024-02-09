#!/usr/bin/env bash 

# Pipeline Code for Processing 16s and RBCL reads from fastq files. 

# 1: #1: Download docker on your computer, allows you to run qiime and
# other software without downloading them.

# On whatever computer will be running the tasks, innitiate a docker
# container for qiime1

#only need to do this docker pull step once per user
docker pull mbari/qiime1

#Helpful notes:
# 1. if running on the shared lab computer, you may get errors if your zipped sequence results are not synched 
#     to your dropbox account. Double clicking the file will usually make it appear. If you look at the filesize and 
# it is 0, this is likely because the file is not synched to the local.

# 2. If you get an error about memory or space issues, check that the other users do not have their dropbox files
#     on the local computer, you can fix this by logging in and making all dropbox files online only. This should
#     clear up enough space to fix this error.

## m1 mac specific flag: --platform linux/amd64 ???

docker run -itv /Volumes/bombus/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline mbari/qiime1

source activate qiime1 

# 2 for 2018 cd into the correct folder within the container, verify that the files you need are there

cd ../../mnt/SI_pipeline/
ls

#3: If your reads are in *fasta.gz format, unzip them with qiime, then
#rename them "forward.fastq" and "reverse.fastq". Quinn designed the
#barcodes in reverse so flowcell1 is the reverse and flowcell2 is the
#forward

## Step 3 only need to be run once ever.

cd R2018/
gunzip *.gz

mv 4376_S0_L001_R1_001.fastq rawreverse.fastq
mv 4376_S0_L001_R2_001.fastq rawforward.fastq

## Run 1 2018
gunzip GC3F-JZ-7102---6632_S1_L001_R1_001.fastq.gz
gunzip GC3F-JZ-7102---6632_S1_L001_R2_001.fastq.gz

mv GC3F-JZ-7102---6632_S1_L001_R1_001.fastq rawreverse.fastq
mv GC3F-JZ-7102---6632_S1_L001_R2_001.fastq rawforward.fastq

## Run 2 2020
cd R2023/lane1

gunzip GC3F-JZ-7102---6632_S1_L001_R1_001.fastq.gz
gunzip GC3F-JZ-7102---6632_S1_L001_R2_001.fastq.gz

mv GC3F-JZ-7102---6632_S1_L001_R1_001.fastq rawreverse.fastq
mv GC3F-JZ-7102---6632_S1_L001_R2_001.fastq rawforward.fastq

## Run 3 2020
cd R2023/lane2

gunzip GC3F-JZ-7632---7075_S1_L001_R1_001.fastq.gz
gunzip GC3F-JZ-7632---7075_S1_L001_R2_001.fastq.gz

mv GC3F-JZ-7632---7075_S1_L001_R1_001.fastq rawreverse.fastq
mv GC3F-JZ-7632---7075_S1_L001_R2_001.fastq rawforward.fastq

#4a: parse the barcodes in the files, putting our data into a format
#qiime2 will be able to use. -- NOTE this step takes ~40 minutes to run!
cd R2018/
extract_barcodes.py -f rawforward.fastq -r rawreverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#5a: re-zip the output files
cd parsed_barcodes
gzip *.fastq

mv reads1.fastq.gz forward.fastq.gz
mv reads2.fastq.gz reverse.fastq.gz

#4b: parse the barcodes in the files, putting our data into a format
#qiime2 will be able to use. -- NOTE this step takes ~40 minutes to run!
cd R2023/lane1
extract_barcodes.py -f rawforward.fastq -r rawreverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#5b: re-zip the output files
cd parsed_barcodes
gzip *.fastq

mv reads1.fastq.gz forward.fastq.gz
mv reads2.fastq.gz reverse.fastq.gz

#4c: parse the barcodes in the files, putting our data into a format
#qiime2 will be able to use. -- NOTE this step takes ~40 minutes to run!
cd R2023/lane2
extract_barcodes.py -f rawforward.fastq -r rawreverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#5c: re-zip the output files
cd parsed_barcodes
gzip *.fastq

mv reads1.fastq.gz forward.fastq.gz
mv reads2.fastq.gz reverse.fastq.gz

#6: Exit Qiime1 and use docker to open the environment for Qiime 2. 
exit

docker run -itv /Volumes/bombus/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core
source activate qiime2

#6: Test that the container for Qiime 2 is properly associated, then
#make sure you are in the root directory using ls and/or pwd, then set working
#directory to the mounted volume

#7: Import your parser barcodes into an object you can demultiplex in
#Qiime 2. Make sure working directory is set correctly.
cd ../../
cd mnt/SI_pipeline/

## Run 1 2018
cd R2018/
qiime tools import --type EMPPairedEndSequences --input-path parsed_barcodes/ --output-path seqs.qza

## Run 2 2023
cd R2023/lane1
qiime tools import --type EMPPairedEndSequences --input-path parsed_barcodes/ --output-path seqs.qza

#Run 2
cd R2023/lane2
qiime tools import --type EMPPairedEndSequences --input-path parsed_barcodes/ --output-path seqs.qza


#8: Examine your mapping file, which is metadata you create associated
#with the project. Use qiime to make it into a qzv object.

#If you are examining multiple amplicon types, pick a map associated
#with one to start with (e.g. 16s)

## 2018 run 1
qiime metadata tabulate --m-input-file maps/sky2018map16s.txt --o-visualization sky2018map16s.qzv
qiime tools view sky2018map16s.qzv

## 2020 run 2
qiime metadata tabulate --m-input-file maps/sky2020map16s_1.txt --o-visualization sky2020map16s_1.qzv
qiime tools view sky2020map16s_1.qzv

## 2020 run 3
qiime metadata tabulate --m-input-file maps/sky2020map16s_2.txt --o-visualization sky2020map16s_2.qzv
qiime tools view sky2020map16s_2.qzv


#9: Demultiplex 16s reads first. Only works in version Qiime2 2019.1

exit

docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1

cd ../../mnt/SI_pipeline
#cd R2018/2023_sequence_results_raw/lane1 #or whichever run you're working on
#cd R2018/2023_sequence_results_raw/lane1
cd R2018/2023_sequence_results_raw/lane2
#note: this step takes ~2 hours!
#qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/sky2020map16s_1.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza 

qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/sky2020map16s_2.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza 

#9a: Visualize Results

qiime demux summarize --i-data demux16s.qza --o-visualization demux16s.qzv

qiime tools view demux16s.qzv

## when interpretating the quality boxes, you can use the bottom of
## the black box as a conservative measure for the phred score (not the
## whiskers and not the middleo f the box)

## watch out for mnt

## The truncation length will vary between each run! Make sure to adjust the numbers pasted below. 
## R2018 16s: f = 180 , r = 220
## R2023 16s: f = 180 , r = 220

#note this step takes hours!

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demux16s.qza  \
--p-trunc-len-f 180  \
--p-trunc-len-r 220  \
--p-trim-left-f 0  \
--p-n-threads 2  \
--output-dir dada2-16s  \
 --o-representative-sequences dada2-16s/rep-seqs-dada2-16s.qza  \
 --o-table dada2-16s/table16s.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv

qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## repeat the steps for RBCL
exit

#docker run -itv /Volumes/bombus/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1
docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core:2019.1


cd ../../mnt/SI_pipeline

#cd R2018/2023_sequence_results_raw/lane1
cd R2018/2023_sequence_results_raw/lane2

## 2020 run 2
#cd ../lane2

#qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/sky2020mapRBCL_1.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demuxRBCL.qza 
qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/sky2020mapRBCL_2.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demuxRBCL.qza 

#9a: Visualize Results

qiime demux summarize --i-data demuxRBCL.qza --o-visualization demuxRBCL.qzv
qiime tools view demuxRBCL.qzv

## The truncation length will vary between each run! Make sure to
## adjust the numbers pasted below.
#this step takes like ~1 hr

## R0 RBCL:  f = 180, r = 218

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demuxRBCL.qza  \
--p-trunc-len-f 180  \
--p-trunc-len-r 218  \
--p-trim-left-f 0  \
--p-n-threads 2  \
--output-dir dada2-RBCL  \
 --o-representative-sequences dada2-RBCL/rep-seqs-dada2-RBCL.qza  \
 --o-table dada2-RBCL/tableRBCL.qza

qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCL.qzv

qiime feature-table summarize --i-table dada2-RBCL/tableRBCL.qza --o-visualization dada2-RBCL/tableRBCL.qzv

# check outputs to make sure you didn't lose too many samples. 
# We found being more conservative and doing shorter truncations gives you the same number of sequences, 
# but sorted into fewer features, likely cause trimmed off poor quality reads


## *****************************************************************************
##       MERGE files from runs ONCE THEY EXIST
## *****************************************************************************
# time to merge the files from your different runs.
# NOTE: much of this comes from https://john-quensen.com/tutorials/merging-dada2-results-in-qiime2/ 

# cd back to your main folder (in this case SI_pipeline) that has the separate run folders
# make a couple new directories

cd ../ 
mkdir merged
cd merged
mkdir 16s
mkdir RBCL
cd ../

#### 1. 16s ###
# 1a. first merge the table files. do this from your main SI_pipeline folder
qiime feature-table merge \
      --i-tables lane1/dada2-16s/table16s.qza \
      --i-tables lane2/dada2-16s/table16s.qza \
      --o-merged-table merged/16s/table16s.qza

 # --i-tables R1/dada2-16s/table16s.qza \
 # --i-tables R2/dada2-16s/table16s.qza \
 # --i-tables R3/dada2-16s/table16s.qza \
 # --i-tables R4/dada2-16s/table16s.qza \
 # --i-tables R2018/dada2-16soutput/table16s.qza \
 # --o-merged-table merge/16s/table16s.qza

# 1b: next merge the rep-seqs
qiime feature-table merge-seqs \
       --i-data lane1/dada2-16s/rep-seqs-dada2-16s.qza \
       --i-data lane2/dada2-16s/rep-seqs-dada2-16s.qza \
	--o-merged-data merged/16s/rep-seqs-16s.qza

 # --i-data R1/dada2-16s/rep-seqs-dada2-16s.qza \
 # --i-data R2/dada2-16s/rep-seqs-dada2-16s.qza \
 # --i-data R3/dada2-16s/rep-seqs-dada2-16s.qza \
 # --i-data R4/dada2-16s/rep-seqs-dada2-16s.qza \
 # --i-data R0_2018/dada2-16soutput/rep-seqs-dada2-16s.qza \
 # --o-merged-data merge/16s/rep-seqs-16s.qza
 
### 2. RBCL ### *see note below before proceeding
# 2a: first merge the table files
qiime feature-table merge \
 --i-tables lane1/dada2-RBCL/tableRBCL.qza \
 --i-tables lane2/dada2-RBCL/tableRBCL.qza \
 --o-merged-table merged/RBCL/tableRBCL.qza
      
 # --i-tables R1/dada2-RBCL/tableRBCL.qza \
 # --i-tables R2/dada2-RBCL/tableRBCL.qza \
 # --i-tables R3/dada2-RBCL/tableRBCL.qza \
 # --i-tables R4/dada2-RBCL/tableRBCL.qza \
 # --i-tables R0_2018/dada2RBCLoutput/tableRBCL.qza \
 # --o-merged-table merged/RBCL/tableRBCL.qza
 
# 2b: next merge the rep-seqs
qiime feature-table merge-seqs \
 --i-data lane1/dada2-RBCL/rep-seqs-dada2-RBCL.qza \
 --i-data lane2/dada2-RBCL/rep-seqs-dada2-RBCL.qza \
 --o-merged-data merged/RBCL/rep-seqs-RBCL.qza


 # --i-data R1/dada2-RBCL/rep-seqs-dada2-RBCL.qza \
 # --i-data R2/dada2-RBCL/rep-seqs-dada2-RBCL.qza \
 # --i-data R3/dada2-RBCL/rep-seqs-dada2-RBCL.qza \
 # --i-data R4/dada2-RBCL/rep-seqs-dada2-RBCL.qza \

# * note: When we ran this code, got an error. leaving the code below, where we fixed it, annotated out, cause not relevant for future orojecs
# we have sample 16s 4285 in both in 16s R1 and R4. decided to remove from R4. kept this in the plate maps, but took out by filtering.
# i called the new filtered table and repseqs outputs in R4 "clean", but then manually renamed them to replace the old artifats
# need to filter out sample 4285 from RBCL of round R4 table16s and rep seqs and then merge

# make a txt file with the list of IDs for samples you want to keep, where you remove control IDs. called "samplestokeep.txt"
# qiime feature-table filter-samples --i-table dada2-RBCL/tableRBCL.qza --m-metadata-file dada2-RBCL/samplestokeep.txt --o-filtered-table dada2-RBCL/tableRBCLclean.qza

# filter out repseqs
# qiime feature-table filter-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --i-table dada2-RBCL/tableRBCLclean.qza --o-filtered-data dada2-RBCL/rep-seqs-dada2-RBCLclean.qza
# qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCLclean.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCLclean.qzv
# qiime tools view dada2-RBCL/rep-seqs-dada2-RBCLclean.qzv


