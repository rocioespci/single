# SINGLe 
SNIPs In Nanopore reads of Gene Libraries

## Installation
The easiest way to install this package is to run in an R console:
require(devtools)
install_github("rocioespci/single")

You can also download it as a /tar and run, inside R, install.package("path_to_file/single.tar.gz", repos=NULL)

## Run

Pre-processing: 

Run the following separately for the reference and for the library reads 

```
minimap2 -ax map-ont --sam-hit-only  REFERENCE.fasta FILE.fastq >FILE.sam    # ALIGN USING MINIMAP2 TO REFERENCE.fast
samtools view -S -b FILE.sam > FILE.bam                                      # TRANSFORM TO BAM
samtools sort FILE.bam -o FILE.bam                                           # SORT BAM FILE
samtools mpileup -Q 0 FILE.bam > FILE.counts.txt                             # COUNT BASES PER POSITION
```

Recommended: Also do it independently for forward and reverse reads, adding --for-only and --rev-only in minimap2 options

Use SINGLe in an R console

Train model using single_train_model. 

Evaluate the model in a data set using single_evaluate_model. 

Compute consensus sequences for barcoded sequences using single_consensus_byBarcode
   


