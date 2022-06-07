# SINGLe 
SNIPs In Nanopore reads of Gene Libraries

   
# Introduction

SINGLe computes consensus sequence of DNA reads by (noisy) nanopore sequencing. It is focused on long amplicons sequencing, and it aims to the reads of gene libraries, typically used in directed evolution experiments.

SINGLe takes advantage that gene libraries are created from an original wild type or reference sequence, and it characterizes the systematic errors made by nanopore sequencing. Then, uses that information to correct the confidence values (QUAL) assigned to each nucleotide read in the mutants library.

Finally, given that you can identify which variant was read in each case (for example by the use of unique molecular identifiers or DNA barcodes), SINGLe groups them and computes the consensus sequence by weighting the frequencies with the corrected confidence values.

For more details, please refer to our pre-print "Accurate gene consensus at low nanopore coverage" doi: https://doi.org/10.1101/2020.03.25.007146 for more information.


# Installation

Using bioconductor:

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("single")
```

Or from git-hub:

```{r}
 require(devtools)
 install_github("rocioespci/single")
```

# How to use

## Before starting

To use SINGLe you must have the following data:

- A fasta file with the reference sequence (ex. a wild type from which you generated the library) (__REF.fasta__).
- Nanopore reads of a reference sequence (__REF_READS.fastq__).
- Nanopore reads of a gene libraries (__LIB_READS.fastq__).
- Identification of each read in the gene library with a gene variant (ex via the use of DNA barcodes on your experiment) (__BC_TABLE.txt__).

## Preprocessing data

Align nanopore reads to the reference sequence and create a sorted bam file.

For the reads of the reference:
```{bash, eval=F}
minimap2 -ax map-ont --sam-hit-only REF.fasta REF_READS.fastq > REF_READS.sam
samtools view -S -b REF_READS.sam > REF_READS.bam
samtools sort REF_READS.bam -o REF_READS.sorted.bam
```

And for the reads of the library:
```{bash, eval=F}
minimap2 -ax map-ont --sam-hit-only  REF.fasta LIB_READS.fastq >LIB_READS.sam
samtools view -S -b LIB_READS.sam > LIB_READS.bam
samtools sort LIB_READS.bam -o LIB_READS.sorted.bams
```

Recommendation: SINGLe works better if you work separately with nanopore reads of the forward and the reverse strand separately. To do that, you can add --for-only and --rev-only in minimap2 options, and follow the downstream analysis independently for each set of aligned reads.


## Run SINGLe in R

SINGLe consists on three steps: train model, evaluate model, compute consensus. As it can be time consuming, I will only analyze a subset of positions
```{r}
library(single)
pos_start <- 1
pos_end <- 10

refseq_fasta <- system.file("extdata", "ref_seq.fasta", package = "single")
ref_seq <- Biostrings::subseq(Biostrings::readDNAStringSet(refseq_fasta), pos_start,pos_end)

```

First, train the model using nanopore reads of the reference (wild type).
```{r}
REF_READS <- system.file("extdata", "train_seqs_500.sorted.bam",package = "single")
train <- single_train(bamfile=REF_READS,
                      output="train",
                      refseq_fasta=refseq_fasta,
                      rates.matrix=mutation_rate,
                      mean.n.mutations=5.4,
                      pos_start=pos_start,
                      pos_end=pos_end,
                      verbose=FALSE,
                      save_partial=FALSE,
                      save_final= FALSE)
print(head(train))
```


Second, evaluate model: use the fitted model to evaluate the reads of your library, and re-weight the QUAL (quality scores).
```{r}
LIB_READS <- system.file("extdata","test_sequences.sorted.bam",package ="single")
corrected_reads <- single_evaluate(bamfile = LIB_READS,
                single_fits = train,
                ref_seq = ref_seq,
                pos_start=pos_start,pos_end=pos_end,
                verbose=FALSE,
                gaps_weights = "minimum",
                save = FALSE)
corrected_reads
```

Finally, use the reads of the library with the corrected QUAL scores to compute a weighted consensus sequences in subsets of reads. The sets of reads corresponding to each variant are indicated in a table (here BC_TABLE) of two columns: SeqID (name of the read) and BCid (barcode or group identity).
```{r}
BC_TABLE = system.file("extdata", "Barcodes_table.txt",package = "single")
consensus <- single_consensus_byBarcode(barcodes_table = BC_TABLE,
                           sequences = corrected_reads,
                           verbose = FALSE)
consensus
```


### Other functions included in the package

Use pileup to create a data.frame with counts by position nucleotide and quality score
```{r}
counts_pnq <- pileup_by_QUAL(bam_file=REF_READS,
                    pos_start=pos_start,
                    pos_end=pos_end)
head(counts_pnq)
```

Compute a priori probability of making errors
```{r}
p_prior_errors <- p_prior_errors(counts_pnq=counts_pnq,
                                  save=FALSE)
p_prior_errors
```

Compute a priori probability of having a mutation
```{r}
p_prior_mutations <- p_prior_mutations(rates.matrix = mutation_rate,
                        mean.n.mut = 5,ref_seq = ref_seq,save = FALSE)
head(p_prior_mutations)
```

Fit SINGLe logistic regression using the prior probabilities and the counts
```{r}
fits <- fit_logregr(counts_pnq = counts_pnq,ref_seq=ref_seq,
                    p_prior_errors = p_prior_errors,
                    p_prior_mutations = p_prior_mutations,
                    save=FALSE)
head(fits)
```

Use the fits to obtain the replacement Qscores after SINGLe fit, for all possible QUAL, nucleotide and position values
```{r}
evaluated_fits <- evaluate_fits(pos_range = c(1,5),q_range = c(0,10),
                                data_fits = fits,ref_seq = ref_seq,
                                save=FALSE,verbose = FALSE)
head(evaluated_fits)
```

Compute one consensus sequence weighted by QUAL values.
```{r}
data_barcode = data.frame(
    nucleotide=unlist(sapply(as.character(corrected_reads),strsplit, split="")),
    p_SINGLe=unlist(1-as(Biostrings::quality(corrected_reads),"NumericList")),
    pos=rep(1:Biostrings::width(corrected_reads[1]),length(corrected_reads)))
consensus_seq <- weighted_consensus(df = data_barcode, cutoff_prob = 0.9)
consensus_seq
another_consensus_seq <- weighted_consensus(df = data_barcode, cutoff_prob = 0.999)
another_consensus_seq
list_mismatches(ref_seq[[1]],another_consensus_seq)
```


