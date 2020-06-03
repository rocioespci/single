#################################################
# Clustering (Distance + hierarchical)
#################################################
require(single);require(dendextend);require(ggplot2);

#---------------------------------------------------------------------------------------#
#Inputs
working_directory       <- "~/Documents/"                 # path of working directory (where data is located)
files_inputs_fastq      <- "barcode11_guppy.fastq"        # Names fastq files, outputs of main_SINGLe. It can be one file, or a vector of files.
file_output_dist_matrix <- "DistanceMatrix.Rdata"         # Name of R.data file in which store the distance matrix
file_output_plot        <- "Dendrogram.png"               # Name of .png file to save the hierarchical clustering plot  
file_output_clusters    <- "Clusters.txt"                 # Name of txt file to save classification of sequences into clusters
file_output_consensus   <- "ConsensusSequences.fasta"     # Name of .fasta file to save consensus calling
reference_sequence      <- KT_wildtype                    # string with the reference sequence
k_clusters              <- 7                              # Number of desired clusters
#---------------------------------------------------------------------------------------#

setwd(working_directory)

#----------# Load Data
data_all     <- read_combine_fastq(files_inputs_fastq)
nseq_all     <- nrow(data_all$sequences)

data_all$weights             <- plyr::mapvalues(data_all$Qscores, from=ascii$Symbol, to=ascii$P, warn_missing = F)
data_all$weights             <- 1-apply(data_all$weights, 2,as.numeric)
rownames(data_all$sequences) <- rownames(data_all$Qscores) <- rownames(data_all$weights)<- data_all$sequence_id

ref_seq                      <- strsplit(reference_sequence, "")[[1]]
homopolymer_positions_all    <- detect_homopolymer_positions(ref_seq)
homopolymer_positions        <- homopolymer_positions_all[sapply(homopolymer_positions_all, length)>1]

#----------# Compute distance matrix
distance_matrix <- matrix(ncol=nseq_all, nrow=nseq_all)
for(s1 in 1:(nseq_all-1)){for(s2 in (s1+1):nseq_all){
  distance_matrix[s1,s2] <- dist_modified_hamming(seq_1 = data_all$sequences[s1,],seq_2 = data_all$sequences[s2,],
                                                  ref = ref_seq,
                                                  w_1 = data_all$weights[s1,],w_2 = data_all$weights[s2,])} }
distance_matrix <- accomodate_dist(x = distance_matrix,data_all$sequence_id)
colnames(distance_matrix) <- rownames(distance_matrix)
save(distance_matrix, file=file_output_dist_matrix)

#----------# Hierarchical clustering
#load(file=file_output_dist_matrix)
dendrogram_original  <- hclust(as.dist(distance_matrix),  method = "ward.D2")
dendrogram           <- as.dendrogram(dendrogram_original)%>% set("branches_k_color", k=k_clusters)%>% set("branches_lwd", 1) %>% set("labels",NULL)
dendrogram           <- as.ggdend(dendrogram)
plot_dendrogram      <- ggplot(dendrogram) #+ scale_y_continuous(limits = c(-12,72))

ggsave(plot_dendrogram,width=200, height = 80, units = "mm",filename = file_output)

clusters_table <- cutree(dendrogram_original, k_clusters)
clusters_table <- data.frame(IDsequence = names(clusters_table), cluster=clusters_table)
write.table(clusters_table, file=file_output_clusters, row.names = F)

#----------# Consensus calling
# clusters_table <- read.table(file=file_output_clusters, header=T)

for(i in 1:k_clusters){
  sequences_id <- clusters_table$IDsequence[which(clusters_table$cluster==i)]
  consensus_hc_single       <- consensus_sequence(data_all$sequences[sequences_id,], probabilities = data_all$weights[sequences_id,])
  cat(">Consensus ",i,"\n",paste0(consensus_hc_single, collapse = ""), file=file_output_consensus)
}

