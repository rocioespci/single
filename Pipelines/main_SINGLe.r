################### SINGLe ##################
# SINGLe is a package for nanopore sequencing, to be applied after basecalling.
# It reweights the Qscores returned by basecallers, and improves clustering and consensus computations.
# Please refer to our manuscript
# “SINGLe: Accurate detection of single nucleotide polymorphisms using nanopore sequencing in gene libraries”,
# by Espada, Zarevski, Dramé-Maigné, Rondelez, 2020.

#################### INPUTS ####################
require(single)

working_path         =  "~/Documents"                  # path of working directory (where data is located)
filename_wildtype    =  "barcode01.maf"                # .maf file with reads of the reference sequence ( “REFERENCE.MAF”)
filename_to_evaluate =  "barcode05.maf"                # .maf file with reads to analyse ( “SEQUENCES.MAF”)
LAST_reference_label =  "ORF_KlenTaq"                  # Label  used in LAST  (in .maf files) to name the reference sequence.
length_sequencing    =  1829                           # Length of the strands sequenced, integer
ORF_ini              =  82                             # In which position of the sequenced strands the open reading frame of the gene starts
ORF_end              =  1746                           # In which position of the sequenced strands the open reading frame of the gene end*
reference_sequence   =  KT_wildtype                    # string with the reference sequence
mean.n.mut           =  5.4                            #integer, mean number of mutations per sequence
mutationratePCR      =  matrix(c(   NA,  14.1,  25.5,  28.5,    # Numerix matrix, 4x5 with the mutation rates from A/C/G/T to A/C/G/T/deletion
                                   4.7,    NA,   4.1,  17.5,    # It changes depending on the method used for creating a library
                                   17.5,   4.1,    NA,   4.7,   # This one represents mutaenzymeII kit from agilent
                                   28.5,  25.5,  14.1,    NA,
                                   1.2, 1.2, 1.2, 1.2), ncol=5,nrow=4)/100

#*This is in case you sequence a region longer than the one you are interested in.
# In this way you can restrict the analysis to ORF_ini to ORF_end to reduce computational times.
# If you are interested in the full region, use ORF_ini=1 and ORF_end=length_sequencing. ORF_ini and
# ORF_end numbering are in reference to 1...length_sequencing.

################ MAIN #########################
setwd(working_path)
reference_sequence <- strsplit(reference_sequence, "")[[1]]
colnames(mutationratePCR) <-bases; rownames(mutationratePCR) <- bases[-c(5)]

LRminion_train_model   (filename_wildtype = filename_wildtype,LAST_reference_label = LAST_reference_label,length_sequencing =  length_sequencing,ORF_ini =  ORF_ini,ORF_end = ORF_end, gaps_weights = 'mean',homopolymer_positions = NULL, auxiliary_files=NULL)
LRminion_evaluate_model(filename_to_evaluate = filename_to_evaluate, filename_wildtype = filename_wildtype,LAST_reference_label = LAST_reference_label,length_sequencing =  length_sequencing, ORF_ini = ORF_ini,ORF_end = ORF_end,reference_sequence = reference_sequence, gaps_weights = 'mean',homopolymer_positions = NULL, auxiliary_files=NULL)

