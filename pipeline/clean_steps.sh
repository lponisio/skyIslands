
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
