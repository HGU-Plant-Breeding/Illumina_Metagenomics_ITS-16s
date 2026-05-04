# Illumina_Metagenomics_ITS-16s
Started from the demultiplexed, adapter-trimmed, PCR-primer-trimmed, merged fastq files.
This means the raw data from teh Illumina sequencer has already been basecalled and reads assigned to samples, for these samples the sequencing barcodes have already been removed. Also the 16s (called rDNA in the following) and ITS primers have been removed. 
I assume for 16s 341F-785R and for ITS ITS1f-ITS2 were used as primers but it does not matter for the following

## 1. Quality filtering of input reads
probably not necessary since Illumina sequencers usually provide reads with a very low error rate anyways but it can't hurt. 

So we are using vsearch for this

```--fastq_maxee 1.0``` means we allow for at most 1 error per read 
```--fastq_minlen 100``` means reads need to be at least 100 bp
```--fastq_maxns 0``` means we remove all reads that have unspecified "N" bases
The reads were compressed with bgzip so we ran ```bzcat``` to decompress them and then stream/send them to vsearch with ```|```

```
bzcat input/341F-785R-P01-A01-T8_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/341F-785R-P01-A01-T8_joined-SR.fasta
bzcat input/341F-785R-P01-B01-T25_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/341F-785R-P01-B01-T25_joined-SR.fasta
bzcat input/341F-785R-P01-C01-T34_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/341F-785R-P01-C01-T34_joined-SR.fasta
bzcat input/ITS1f-ITS2-P01-A01-T8_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/ITS1f-ITS2-P01-A01-T8_joined-SR.fasta
bzcat input/ITS1f-ITS2-P01-B01-T25_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/ITS1f-ITS2-P01-B01-T25_joined-SR.fasta
bzcat input/ITS1f-ITS2-P01-C01-T34_joined-SR.fastq.bz2 | vsearch --fastq_filter - --fastq_maxee 1.0 --fastq_minlen 100 --fastq_maxns 0  --fastaout filter/ITS1f-ITS2-P01-C01-T34_joined-SR.fasta
```
Runs only a few seconds per sample. 

Important notice for next time: vsearch converts fastq (with quality scoring) to fasta (just sequences)

## 2. Dereplication of filtered reads
Since PCR was used to amplifiy the sequencing library before sequencing, we now remove the PCR-duplicates generated in that process. 

We also use vsearch for that

```
for i in 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do vsearch \
    --derep_fulllength filter/${i}.fasta \
    --output derep/${i}.fasta \
    --sizeout --minuniquesize 2
 done
```
For simplicity I ran this all in a bash for loop so the ```for i in``` says for every element ```i``` from the filenames following afterwards run this command and the ```${i}``` is then the variable name in the loop. This way we get the same input and output names. We just ake it from the folder `filter` and write it to `derep`

## 3. Removing Chimeric Reads
It can happen during PCR and/or sequencing that independent amplicons get fused together and show up as one read in the sequencing output. To fix this we run `vsearch --uchime`, which splits these chimeras

```
for i in 341F-785R-P01-A01-T8_joined-SR 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do vsearch \
    --uchime_denovo derep/${i}.fasta \
    --nonchimeras nochim/${i}.fasta
  done
```

## 4. Merge Reads belonging to the same Species rDNA/ITS sequence into clusters
Reads that are from the same rDNA gene or ITS spacer are now merged into operational taxonomic units (OTUs)

```
for i in 341F-785R-P01-A01-T8_joined-SR 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do vsearch \
    --cluster_size nochim/${i}.fasta \
    --id 0.99 --centroids otu/${i}.fasta \
    --relabel OTU_
  done
```

We merge them if they are at least 99% identical and we rename the resulting consensus sequences to ```OTU_``` and an increasing number as a unique identifier

## 5. Counting Reads per OTU-Cluster
Now we use the filtered reads from step 1 and map them back to the OTUs
```
for i in 341F-785R-P01-A01-T8_joined-SR 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do vsearch \
    --usearch_global filter/${i}.fasta \
    --db otu/${i}.fasta --id 0.97 otu_map/${i}.txt
  done
```

We accept mappings if they are at least a 97% match.

## 6. Assigning a species to the OTUs

### 6.1 preparing databases
I used the SILVA database for the 16s and the UNITE database for the ITS species assignment. Both can be downloaded online 
silva from here [https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/ARB_files/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz](https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/ARB_files/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz)
unite from here [https://doi.plutof.ut.ee/doi/10.15156/BIO/3301229](https://doi.plutof.ut.ee/doi/10.15156/BIO/3301229)

Since they have very sequence long names that include the complete taxonomy starting from kingdom or domain I used two python scripts to split them in a way that I have only the sequnce ID and sequence in the fasta file and the taxonomy for all IDs in a separate file
```
python reformat_silva_fasta.py -i SILVA_138.2_SSURef_NR99_tax_silva.fasta -o SILVA_reformat.fasta -t SILVA_taxonomy.tsv
python reformat_unite_fasta.py -i sh_general_release_dynamic_19.02.2025.fasta -o unite_reformatted.fasta -t unite_taxonomy.tsv
```

Then I used NCBI-BLAST to make a blast compatible database out of both fasta-files

```
makeblastdb -in SILVA_reformat.fasta -dbtype nucl -out SILVA_reformat -parse_seqids
makeblastdb -in unite_reformatted.fasta -dbtype nucl -out unite_reformatted -parse_seqids
```
### 6.2 Checking OTUs against databases
Since we have three ITS and three 16s samples we check three against the unite and three against the silva databases we build

```
for i in 341F-785R-P01-A01-T8_joined-SR 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR
  do
    blastn \
      -query otu/${i}.fasta \
      -db database/SILVA/SILVA_reformat \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -out blast/${i}.blast.out \
      -max_target_seqs 5 \
      -num_threads 8 \
      -max_hsps 1
  done

for i in ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do
    blastn \
      -query otu/${i}.fasta \
      -db database/unite/unite_reformatted \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -out blast/${i}.blast.out \
      -max_target_seqs 5 \
      -num_threads 8 \
      -max_hsps 1
  done
```
`query` is the file with the OTUs
`db` the database we created
`outfmt` specifies the columns we want in the output file
`out` is the output file name
`max_target_seqs` the maximum number of results we allow per OTU
`num_threads` the number of CPUs the programm is allowed to use
`may-hsps` specifies that a referenc ehit can show up only once per OTU in the results (else sometimes the beginning of a sequnce and its end match the same reference sequence and it shows up twice in the results)

### 6.3 Assigning Species to each OTU
With the count files from step 5, the annotations from 6.2 and the taxonomy from 6.1 we can now compile these results and produce a table for each samples that tells us how many reads belong to each OTU and to which species this OTU belongs. 
We again use a self-made python-script for this

```
for i in 341F-785R-P01-A01-T8_joined-SR 341F-785R-P01-B01-T25_joined-SR 341F-785R-P01-C01-T34_joined-SR
  do python blast_to_count_table.py \
    -d silva -b blast/${i}.blast.out \
    -t database/SILVA/SILVA_taxonomy.tsv \
    -c otu_map/${i}.txt \
    -o count_table_rDNA/${i}.cluster_count_taxonomy.tsv
  done

for i in ITS1f-ITS2-P01-A01-T8_joined-SR ITS1f-ITS2-P01-B01-T25_joined-SR ITS1f-ITS2-P01-C01-T34_joined-SR
  do python blast_to_count_table.py \
    -d unite -b blast/${i}.blast.out \
    -t database/unite/unite_taxonomy.tsv \
    -c otu_map/${i}.txt \
    -o count_table_ITS/${i}.cluster_count_taxonomy.tsv
  done
```

Since the taxonomy files for unite and silva look slightly different we have to specify the input type with `-d` as unite or silva (Unite starts at kindom, Silva at domain). `-b` is our blast results file from 6.2. `-t` the taxonomy from 6.1. `-c` are the counts from step 5. And `-o` our output-files (we want them written to the same folder and ending with `.cluster_count_taxonomy.tsv` for step 7.)

## 7. Build Species Matrix
Finally we build a matrix/table where each line is a species (with full taxonomy) and the columns are the numbers of sequences found belonging to this species in each respective samples. Again we make one for ITS and one for 16s.

```
python combine_sample_taxonomy_tables_from_folder.py -i count_table_rDNA/ -o species_matrix_rDNA.tsv
python combine_sample_taxonomy_tables_from_folder.py -i count_table_ITS/ -o species_matrix_ITS.tsv
```

## 8. Make Bar-Plots
Loading the `species_matrix_rDNA.tsv` and `species_matrix_ITS.tsv` into R and building bar-charts using ab R-Script called `bar-chart-rDNA.R`





