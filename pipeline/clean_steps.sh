
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

