#Pipeline for Illumina data with qiime (Sep 2016)
#Para poder extraer los barcodes de manera sencilla manejare el rev como fwd y vice versa.
#En el proceso de split libraries hay una opcion de escribir el reverse complement antes
de mandarlo al output file, para corregir todo
## Mas info en http://qiime.org/tutorials/processing_illumina_data.html

### PARA CORRER QIIME EN MUUK
source ~/.virtualenvs/qiime/bin/activate

##### 1. Create a mapping file

##### 2. Validate mapping file
validate_mapping_file.py -m Fasting_Map.txt -o mapping_output

##### 3. Unzip Reverse (R2) fastq files
gunzip E_S1_L001_R2_001.fastq.gz

##### 3. We have to extract the barcode files with: Note that the barcode is extracted only from the rev primer

extract_barcodes.py -f R_S1_L001_R2_001.fastq -c barcode_single_end --bc1_len 6 -o A1_processed_revseqs
extract_barcodes.py -f S_S2_L001_R2_001.fastq -c barcode_single_end --bc1_len 6 -o B1_processed_revseqs
extract_barcodes.py -f T_S3_L001_R2_001.fastq -c barcode_single_end --bc1_len 6 -o C1_processed_revseqs


#### 4. Once extracted the barcodes then do join the pair end reads with this instruction 
###Using the files obtained from extracting the barcode (reads and barcode) and the forwards. 
###Used them inverted. This will allow the input into the split libraries command. The default parameter is around 10 overlap,
###however here I will use 50 bp. The max score diff is 8%

join_paired_ends.py -f J1_processed_revseqs/reads.fastq -r A_S1_L001_R1_001.fastq -b J1_processed_revseqs/barcodes.fastq -o J1_fastqjoin_50bp -j 50

#### NOTE: this will generate de 16S sequences reversed I still have to figure out if this does not matter for the alignment. Either as clustering or with a database.

#### 5. Once you have the joined paired ends demultiplex, filter quality data and extract quality scores, the Phred is 30 and the max barcode error is 0, specify the barcode length to 6
split_libraries_fastq.py -i AG_fastqjoin_50bp/fastqjoin.join.fastq -b AG_fastqjoin_50bp/fastqjoin.join_barcodes.fastq --barcode_type 6 --max_barcode_errors 0 --rev_comp -o AG_slout_q30/ -m Mapping_file_AG.txt --store_qual_scores -q30

#### 6. To check all the statistics, sequences per sample, discarded sequences etc we only need to check the log file
NOTE: Usually it is fine if you loose up to 10% of the reads. It will be fin to have only around 50,000 seq pers sample
the American Gut project only has around thousends of reads not hundred thousend.
So try and play with the phred score up to 30

### SI YA SE TIENE INSTALADO USEARCH SOLO SE HACE A PARTIR DE ESTE PASO
##Desde muuk:
bash-4.1$  qiime > export PATH=$PATH:/export/home/smoran/ # Como ya se creo el link a usearch61, esto es para ponerlo en mi path                     


# 7. Ahora si a checar quimeras en cada bloque por separado porque si no no alcanza la version gratuita de usearch 32bit
## Esto lo hare en las seqs de new merge, q30
identify_chimeric_seqs.py -m usearch61 -i AH_slout_q30/seqs.fna -r /export/home/smoran/.virtualenvs/qiime/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -o AH_chimera_checked_seqs_q30/

#8. Este comando me va a dar varios log y el archivo de chimeras.txt para poder filtrar mis secuencias con el siguiente comando:
filter_fasta.py -f O_slout_q30/seqs.fna -o O_seqs_chimeras_filtered.fna -s O_chimera_checked_seqs_q30/chimeras.txt -n

# 9. Do this on all 4 sets of reads and concatenate:
cat AK_29JUL16-38462681/AK_seqs_chimeras_filtered.fna BK_29JUL16-38462682/BK_seqs_chimeras_filtered.fna CK_29JUL16-38462683/CK_seqs_chimeras_filtered.fna DK_29JUL16-38462684/DK_seqs_chimeras_filtered.fna > PlacaK_seqs_chimera_free_q30.fna

#### Merge mapping files as well
merge_mapping_files.py -m Mapping_Ricardo1y2.txt,RICARDO3/Mapping_Ricardo3.txt -o Mapping_Ricardo_1_2y3.txt -n 'Data not collected'

# 10. Run pick otus in closed reference and open reference (on the combined sequences (With human gut as it is unlikely to find very new things)
####Check if you loose many sequences and if you do, you can run open reference
pick_closed_reference_otus.py -i PlacaF_seqs_chimera_free_q30.fna -o Closedref_PlacaF_chimera_free_q30 -a -O 5
pick_open_reference_otus.py -i PlacaH_seqs_chimera_free_q30.fna -o Openref_PlacaH_chimera-free_q30 -r /export/home/smoran/.virtualenvs/qiime/local/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -a -O 5

#11. Para saber cuantas lecturas quedaron por muestra:
biom summarize-table -i otu_table_mc2_w_tax_no_pynast_failures.biom -o summary_otu_table_open_PlacaM.txt
biom summarize-table -i otu_table.biom -o summary_otu_table_closed_PlacaM.txt

# Para saber cuantas OTUS hay por muestra al comando anteior le pasamos el paramtero --qualitative
biom summarize-table -i otu_table_SinC25yF4_ABD_E27_no_otus_less10percy01perc_sinG4.biom -o summary_otu_table_SInC25yF4_ABD_E27_no_otus_10perc_y_01perc_OTUS.txt --qualitative

##CORRI OPEN Y CLOSED REF PARA LAS PLACAS F-J Y SE PIERDEN ENTRE EL 9 Y EL 11% DE LAS SECUENCIAS POR LO QUE UTILIZARE LOS RESULTADOS DE CLOSED REFERENCE

### COMO EN LOS MAPPING FILE PUSE LA VE14.009 SIN NINGUNA LETRA REFERIDA A LA PLACA SE LA QUITE AL ARCHIVO BIOM, PARA PODER HACER EL MERGE ENTRE ELLOS
filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_sinControl.biom --sample_id_fp /export/home/smoran/Exclude_control.txt --negate_sample_id_fp

merge_otu_tables.py -i Placa_F/Closedref_PlacaF_chimera_free_q30/otu_table_sinControl.biom,Placa_G/Closedref_PlacaG_chimera_free_q30/otu_table_sinControl.biom,Placa_H/Closedref_PlacaH_chimera_free_q30/otu_table_sinControl.biom,Placa_J/Closedref_PlacaJ_chimera_free_q30/otu_table_sinControl.biom -o Otu_table_closed_Placa_FGHJ_sin_Control.biom

