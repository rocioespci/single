################### DESCRIPTION v 3 marzo 2020 ##################
# This code performs the correction on minION reads of amplicons
# described in Espada et al. 2020.
# Please include your inputs in the "INPUTS" section before running it
# Take into account that several intermediate files will be stored during the process.

# Input should be sequences aligned with LAST to the wild type reference
# Output will be table of Q scores observed in each library (corrected and not). Columns are:
# "position" "nucleotide" "counts" "qualityQ" "p_error_minION" "wt.base" "isWT" "p_right_priors_model" "p_right_null_model"

# INPUTS
# working_path:            It should have the output of LAST alignments. Several auxiliary files and final results will be saved in it.
# filename_wildtype:       .maf file that will be used for training the model
# filename_to_evaluate:    It should be a vector with the names of all the files that you want to be evaluated
# LAST_reference_label:    When you align using LAST, the rows with the reference sequence (wild type here) is tagged with some name you've chosen. Indicate this tag in LAST_reference_label, that will be used for parsing LAST files
# length_sequencing:       Length of the reference sequence you used in LAST
# ORF_ini and ORF_end:     If the reference sequence is longer than the ORF or the interesting fragment you want to analyse, indicate inital and last positions of the ORF in the reference sequence
# reference_sequence:     Wildtype or reference sequence. Use on string as in the example.
# For priors of mutations you should provide:
#    I. mutationratePCR (a matrix 4x4 with the mutation rate of the PCR, in the order a c g t, provided by the kit you used), and
#    mean.n.mut the mean number of mutations indtroduced. This depends on your experimental characteristics.

#################### INPUTS ####################
require(single)

working_path         =  "/home/robotspm/Documents/Rocio/LRminION2020/RawData" #system.file(package = "minIONLogRegFull") #filename/character
filename_wildtype    =  "barcode01.maf" #filename/character
filename_to_evaluate =  c("barcode05.maf")#,"barcode06.maf","barcode07.maf","barcode08.maf","barcode09.maf","barcode10.maf","barcode11.maf","barcode12.maf") #filename/character
LAST_reference_label =  "ORF_KlenTaq" #character
length_sequencing    =  1829 #integer
ORF_ini              =  82   #integer
ORF_end              =  1746 #integer
reference_sequence   =  KT_wildtype #character: sequence DNA quote
mutationratePCR      =  matrix(c(   NA,  14.1,  25.5,  28.5,
                                   4.7,    NA,   4.1,  17.5,
                                   17.5,   4.1,    NA,   4.7,
                                   28.5,  25.5,  14.1,    NA,
                                   1.2, 1.2, 1.2, 1.2), ncol=5,nrow=4)/100
mean.n.mut           =  5.4 #integer

################ MAIN #########################


setwd(working_path)
reference_sequence <- strsplit(reference_sequence, "")[[1]]
colnames(mutationratePCR) <-bases; rownames(mutationratePCR) <- bases[-c(5)]

LRminion_train_model   (filename_wildtype = filename_wildtype,LAST_reference_label = LAST_reference_label,length_sequencing =  length_sequencing,ORF_ini =  ORF_ini,ORF_end = ORF_end, gaps_weights = 'mean',homopolymer_positions = NULL, auxiliary_files=NULL)
LRminion_evaluate_model(filename_to_evaluate = filename_to_evaluate, filename_wildtype = filename_wildtype,LAST_reference_label = LAST_reference_label,length_sequencing =  length_sequencing, ORF_ini = ORF_ini,ORF_end = ORF_end,reference_sequence = reference_sequence, gaps_weights = 'mean',homopolymer_positions = NULL, auxiliary_files=NULL)

