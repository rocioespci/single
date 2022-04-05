# SINGLe 
SNIPs In Nanopore reads of Gene Libraries

## Installation

From bioconductor, run in an R console:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("single")


Or from git-hub, run in an R console:
require(devtools)
install_github("rocioespci/single")


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
   


