
#1: On whatever computer will be running the tasks, innitiate a docker container for qiime1

docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline sglim2/qiime-1.9.1


#2: cd into the correct folder within the container, verify that the files you need are there

cd mnt/SI_pipeline/16sRBCL
ls


#3: If your reads are in *fasta.gz format, unzip them with qiime, then rename them "forward.fastq" and "reverse.fastq"

gunzip *.gz

mv [file name] forward.fastq
mv [file name] reverse.fastq


#4: parse the barcodes in the files, putting our data into a format qiime2 will be able to use

extract_barcodes.py -f forward.fastq -r reverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes

#5: re-zip the output files

cd parsed_barcodes
gzip *.fastq

#5: Exit Qiime1 and use docker to open the environment for Qiime 2. 

exit 
docker run -itv ~/Dropbox/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline qiime2/core

#6: Test that the container for Qiime 2 is properly associated, then make sure you are in the root directory using ls, then set working directory to the mounted volume
 
qiime --help
ls
cd mnt/SI_pipeline/16sRBCL

#7: Import your parser barcodes into an object you can demultiplex in Qiime 2. Make sure working directory is set correctly.
qiime tools import --type EMPPairedEndSequences --input-path /mnt/SI_pipeline/16sRBCL/parsed_barcodes/ --output-path seqs.qza

#8: Examine your mapping file, which is metadata you create associated with the project. Use qiime to make it into a qzv object.
#If you are examining multiple amplicon types, pick a map associated with one to start with (e.g. 16s)

qiime metadata tabulate --m-input-file sky2018map16s.txt --o-visualization sky2018map16s.qzv
qiime tools view sky2018map16s.qzv

#9: Demultiplex 16s reads

qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file ../sky2018map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza --p-no-golay-error-correction --output-dir error

#9a: Visualize Results

qiime demux summarize --i-data demux16s.qza --o-visualization demux16s.qzv
qiime tools view demux16s.qzv

