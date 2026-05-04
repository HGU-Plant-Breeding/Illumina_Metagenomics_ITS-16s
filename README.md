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
I used the SILVA database for the 16s and the UNITE database for the ITS species assignment. Both can be downloaded 
