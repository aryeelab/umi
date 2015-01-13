# umi
Preprocessing tools for unique molecular index (UMI) sequencing reads

## Example
The example directory contains undemultipexed data from a subset of an Illumina MiSeq run.

#### Demultiplex reads
    cd test
    python ../demultiplex.py --min_reads 1000 --read1 undemux.r1.fastq.gz --read2 undemux.r2.fastq.gz --index1 undemux.i1.fastq.gz --index2 undemux.i2.fastq.gz --sample_barcodes samplekey.txt
    
#### Add a molecular index (UMI) tag to the header of the R1 and R2 reads
    SAMPLE=sample1
    python ../umitag.py --read1_in $SAMPLE.r1.fastq --read2_in $SAMPLE.r2.fastq --read1_out $SAMPLE.r1.umitagged.fastq --read2_out $SAMPLE.r2.umitagged.fastq --index1 $SAMPLE.i1.fastq --index2 $SAMPLE.i2.fastq
    
#### Consolidate reads with the same molecular index
    python ../consolidate.py $SAMPLE.r1.umitagged.fastq $SAMPLE.r1.consolidated.fastq 15 0.9
    python ../consolidate.py $SAMPLE.r2.umitagged.fastq $SAMPLE.r2.consolidated.fastq 15 0.9



    
