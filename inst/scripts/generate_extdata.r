require(single)
refseq_fasta<- system.file("extdata", "ref_seq.fasta", package = "single")
ref_seq = Biostrings::readDNAStringSet(refseq_fasta)
train_reads_example <- system.file("extdata", "train_seqs_500.sorted.bam",
                                package = "single")


## Generate fits_example.txt
counts_pnq <- pileup_by_QUAL(bam_file=train_reads_example,
                             pos_start=1,
                             pos_end=10)

p_prior_errors <- p_prior_errors(counts_pnq=counts_pnq,
                                 output_file="none",
                                 save=F)

p_prior_mutations   <- p_prior_mutations(rates.matrix = mutation_rate,
                                         mean.n.mut = 5,
                                         ref_seq = ref_seq,
                                         save = FALSE,
                                         output_file="none")

fits_example <- fit_logregr(counts_pnq = counts_pnq,
                    ref_seq=ref_seq,
                    p_prior_errors = p_prior_errors,
                    p_prior_mutations = p_prior_mutations,
                    save=F)
write.table(fits_example, file="single/inst/extdata/fits_example.txt")

## Generate train_example.txt
train_example <- single_train(bamfile=train_reads_example,
                refseq_fasta=refseq_fasta,
                rates.matrix=mutation_rate,mean.n.mutations=5,
                pos_start=1,pos_end=10,)
write.table(train_example, file="single/inst/extdata/train_example.txt")
