# umi
Preprocessing tools for unique molecular index (UMI) sequencing reads

## Depedencies

 - argparse
 - HTSeq

## Example
The example directory contains undemultipexed data from a subset of an Illumina MiSeq run.

#### Demultiplex reads
    cd example
    python ../demultiplex.py --min_reads 1000 --read1 undemux.r1.fastq.gz --read2 undemux.r2.fastq.gz --index1 undemux.i1.fastq.gz --index2 undemux.i2.fastq.gz --sample_barcodes samplekey.txt
    
#### Add a molecular index (UMI) tag to the header of the R1 and R2 reads
    python ../umitag.py --read1_in mysample.r1.fastq --read2_in mysample.r2.fastq --read1_out mysample.r1.umitagged.fastq --read2_out mysample.r2.umitagged.fastq --index1 mysample.i1.fastq --index2 mysample.i2.fastq
    
#### Consolidate reads with the same molecular index
    python ../consolidate.py mysample.r1.umitagged.fastq mysample.r1.consolidated.fastq 15 0.9
    python ../consolidate.py mysample.r2.umitagged.fastq mysample.r2.consolidated.fastq 15 0.9


    
