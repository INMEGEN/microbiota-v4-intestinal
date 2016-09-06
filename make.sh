#Pipeline for Illumina data with qiime (Sep 2016)
#Para poder extraer los barcodes de manera sencilla manejare el rev como fwd y vice versa.
#En el proceso de split libraries hay una opcion de escribir el
#reverse complement antes de mandarlo al output file, para corregir
#todo

## Mas info en http://qiime.org/tutorials/processing_illumina_data.html

# usage:
# ./split_library.sh PREFIX Mapping_File.txt


PREFIX=$1
MF=$2


# Validate mapping file
validate_mapping_file.py -m $MF -o mapping_output

# Unzip forward and Reverse fastq files
gunzip ${PREFIX}_R1.fastq.gz
gunzip ${PREFIX}_R2.fastq.gz

# We have to extract the barcode files with: Note that the
# barcode is extracted only from the rev primer

extract_barcodes.py -f ${PREFIX}_R2.fastq \
                    -c barcode_single_end \
                    --bc1_len 6 \
                    -o ${PREFIX}_processed_revseqs


# Once extracted the barcodes then do join the pair end reads with
# this instruction Using the files obtained from extracting the barcode
# (reads and barcode) and the forwards. Used them inverted. This will
# allow the input into the split libraries command. The default
# parameter is around 10 overlap, however here I will use 50 bp. The
# max score diff is 8%

join_paired_ends.py -f ${PREFIX}_processed_revseqs/reads.fastq \
                    -r ${PREFIX}_R1.fastq \
                    -b ${PREFIX}_processed_revseqs/barcodes.fastq \
                    -o ${PREFIX}_fastqjoin_50bp \
                    -j 50

# Once you have the joined paired ends, demultiplex, filter quality
# data and extract quality scores. The Phred is 30 and the max barcode
# error is 0, specify the barcode length to 6

split_libraries_fastq.py -i ${PREFIX}_fastqjoin_50bp/fastqjoin.join.fastq \
                         -b ${PREFIX}_fastqjoin_50bp/fastqjoin.join_barcodes.fastq \
                         --barcode_type 6 \
                         --max_barcode_errors 0 \
                         --rev_comp \
                         -o ${PREFIX}_slout_q30/ \
                         -m $MF \
                         --store_qual_scores -q30

# NOTE: Usually it is fine if you loose up to 10% of the reads. It
# will be fin to have only around 50,000 seq pers sample the American
# Gut project only has around thousends of reads not hundred thousend.
# So try and play with the phred score up to 30

# Filter out chimeric sequences with usearch.
# note: usearch61 must be in the $PATH

GG=`print_qiime_config.py | grep pick_otus_reference_seqs_fp | awk '{print $2}'`

# identify chimeras
identify_chimeric_seqs.py -m usearch61 \
                          -i ${PREFIX}_slout_q30/seqs.fna \
                          -r $GG \
                          -o ${PREFIX}_chimera_checked_seqs_q30/

# remove identified chimeras
filter_fasta.py -f ${PREFIX}_slout_q30/seqs.fna \
                -o ${PREFIX}_seqs_chimeras_filtered.fna \
                -s ${PREFIX}_chimera_checked_seqs_q30/chimeras.txt \
                -n
