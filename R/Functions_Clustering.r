#DISTANCE FUNCTIONS
dist_modified_hamming      <- function(seq_1,seq_2,ref, w_1, w_2){
  #index where sequences are different to reference
  ind_1 <- which(seq_1!=ref & seq_1!= "-");
  ind_2 <- which(seq_2!=ref & seq_2!= "-");
  ind <- unique(sort(c(ind_1, ind_2)))

  #verify if they are the same or different between the two sequences
  ind_equal <- ind[seq_1[ind]==seq_2[ind]]
  ind_diff  <- ind[seq_1[ind]!=seq_2[ind]]

  #weight by weights
  aux1 <- if(length(ind_diff)>0) {sum(w_1[ind_diff]  * w_2[ind_diff], na.rm=T)}else{0}
  aux2 <- if(length(ind_equal)>0){sum(w_1[ind_equal] * w_2[ind_equal], na.rm=T)}else{0}
  distance_value <- aux1-aux2
  return(distance_value)
}
dist_different_minus_equal <- function(seq_1,seq_2,ref){
  #index where sequences are different to reference
  ind_1 <- which(seq_1!=ref & seq_1!= "-");
  ind_2 <- which(seq_2!=ref & seq_2!= "-");
  ind <- unique(sort(c(ind_1, ind_2)))

  #verify if they are the same or different between the two sequences
  ind_equal <- ind[seq_1[ind]==seq_2[ind]]
  ind_diff  <- ind[seq_1[ind]!=seq_2[ind]]

  distance_value <- length(ind_diff)-length(ind_equal)
  return(distance_value)

}
dist_modified_hamming_soft <- function(seq_1,seq_2,ref, w_1, w_2, cut_off=0.9){
  #index where sequences are different to reference
  ind_1 <- which(seq_1!=ref & seq_1!= "-");
  ind_2 <- which(seq_2!=ref & seq_2!= "-");
  ind_unique   <- intersect(ind_1, ind_2)
  #and they are different between the two sequences
  ind_diff  <- ind_unique[seq_1[ind_unique]!=seq_2[ind_unique]]
  #filter with cut_off on weights
  ind <- ind_diff[which(w_1[ind_diff]>cut_off | w_2[ind_diff]>cut_off)]
  #distance is product of weights
  distance_value  <- if(length(ind)>0){ sum(w_1[ind] * w_2[ind])}else{0}
  return(distance_value)
}
dist_equal                 <- function(seq_1,seq_2,ref, w_1, w_2, cut_off=0.9){
  #index where sequences are different to reference
  ind_1 <- which(seq_1!=ref & seq_1!= "-");
  ind_2 <- which(seq_2!=ref & seq_2!= "-");
  ind_both   <- intersect(ind_1, ind_2)
  #and they are different between the two sequences
  ind_diff  <- ind_both[seq_1[ind_both]!=seq_2[ind_both]]
  #filter with cut_off on weights
  ind <- ind_diff[which(w_1[ind_diff]>cut_off | w_2[ind_diff]>cut_off)]

  distance_value <- length(ind)
  return(distance_value)
}

modified_hamming_distance_matrix <- function(filename_to_evaluate, reference_sequence, weighting="priors",
                                             dist_function="dist_modified_hamming", cut_off=0.9,
                                             filename_out=NULL,filename_ids=NULL){
  if(!is.null(filename_ids) & length(filename_to_evaluate)!=length(filename_ids)){stop('Length filename_to_evaluate must be the same as filename_ids')}

  data <- lapply(filename_to_evaluate, load_sequences_from_path)
  if(!is.null(filename_ids)){
    ids_repetitions <- sapply(data, length)
    ids <- rep(filename_ids, ids_repetitions)
  }

  data      <- unlist(data, recursive=FALSE)
  sequences <- t(sapply(data, function(x){x$nucleotide}))

  if(weighting=="raw"){       weights   <- t(sapply(data, function(x){1-x$p_error_minION}))      }
  if(weighting=="naive"){     weights   <- t(sapply(data, function(x){x$p_right_null_model}))    }
  if(weighting=="priors"){    weights   <- t(sapply(data, function(x){x$p_right_priors_model}))  }

  n_sequences     <- length(data)
  distance_matrix <- matrix(ncol=n_sequences, nrow=n_sequences)

  if(dist_function=="dist_modified_hamming"){
    for(s1 in 1:(n_sequences-1)){    for(s2 in (s1+1):n_sequences){
        distance_matrix[s1,s2] <- dist_modified_hamming(sequences[s1,],sequences[s2,], reference_sequence, weights[s1,],weights[s2,])
    }  }
  }
  if(dist_function=="dist_modified_hamming_soft"){
    for(s1 in 1:(n_sequences-1)){    for(s2 in (s1+1):n_sequences){
      distance_matrix[s1,s2] <- dist_modified_hamming_soft(sequences[s1,],sequences[s2,], reference_sequence, weights[s1,],weights[s2,],cut_off = cut_off)
    }  }
  }
  if(dist_function=="dist_different_minus_equal"){
    for(s1 in 1:(n_sequences-1)){    for(s2 in (s1+1):n_sequences){
      distance_matrix[s1,s2] <- dist_different_minus_equal(sequences[s1,],sequences[s2,], reference_sequence)
    }  }
  }
  if(dist_function=="dist_equal"){
    for(s1 in 1:(n_sequences-1)){    for(s2 in (s1+1):n_sequences){
      distance_matrix[s1,s2] <- dist_equal(sequences[s1,],sequences[s2,], reference_sequence, weights[s1,],weights[s2,], cut_off = cut_off)
    }  }
  }

  matrix_names <- unlist(sapply(filename_to_evaluate, dir))
  distance_matrix <- accomodate_dist(x = distance_matrix, matrix_names)
  if(!is.null(filename_out)){  write.table(distance_matrix, file = filename_out)  }
  return(distance_matrix)
}

#GENERAL FUNCTIONS
accomodate_dist <- function(x, names=NULL){
  if(ncol(x)!=nrow(x)){stop('x must be a squared matrix')}
  if(!is.null(names) & length(names)!= ncol(x)){'length(names) must be equal to ncol(x) or NULL'}
  x[is.na(x)] <- 0
  x[x<0] <- 0
  y = x+t(x)
  colnames(y) <- rownames(y) <- names
  return(y)
}
make_labels     <- function(x){
  n <- length(x)
  y <- c()
  if(n>6){l <- 3}else{l <- 2}
  for(i in 1:n){
    y <- paste(y, x[i])
    if(i%%l==0){y <- paste0(y, "\n")}
  }
  return(y)
}
parse_labels    <- function(strings, split, element,name=F,fixed=FALSE){
  y <- strsplit(strings, split=split, fixed=fixed)
  y <- sapply(y, function(x){x[element]})
  if(name==FALSE){names(y) <- NULL}
  return(y)
}
plot_tree       <- function(hierarchical_clustering,ids,k_clusters, colors=rainbow(k_clusters),main,consensus_labels=NULL){
  hc_order  <- hierarchical_clustering$order
  hc        <- as.dendrogram(hierarchical_clustering)

  if(!is.null(consensus_labels)){
    plotting_order <- unique(ids[hc_order])
    plotting_x     <- sapply(1:k_clusters, function(x){mean(which(ids[hc_order]==x))})
    plotting_labels <- consensus_labels[plotting_order]
   # plotting_colors <- colors[plotting_order]
  }

  ids <- as.factor(ids)
  ids <- as.numeric(ids)
  labels_colors(hc) <- colors[ids[hc_order]]
  labels(hc)        <- rep(".", length(ids))

  plot(hc, xlab="", main=main)
  legend(x="topright", inset=0.02, col=colors,pch=19,
         legend =  1:k_clusters, ncol=2,title = "Sequence")

  if(!is.null(consensus_labels)){
    for (i in 1:k_clusters){
      label <- make_labels(unlist(strsplit(plotting_labels[i], " ")))
      mtext(at=plotting_x[i], side=1, line=0.1,text =label , col=colors[i],
            cex=1, padj=1)
      mtext(at=plotting_x[i], side=1, line=-0.1,text = paste0("Cluster ",plotting_order[i]) ,
            col=colors[i], cex=0.5, padj = 0)
    }
  }

  return(0)

}

#CONSENSUS FUNCTIONS
clusters_consensus    <- function(sequences,clusters_vector,reference_sequence,weights_vector=NULL,cutoff_prob=0.2){
  if(is.null(weights_vector)){weights_vector <- matrix(1,ncol=ncol(sequences), nrow=nrow(sequences))}
  k_clusters=length(unique(clusters_vector))
  consensus_for_clusters <- data.frame(cluster = 1:k_clusters, consensus_mutations = rep(NA,k_clusters))
  for(i in 1:k_clusters){
    index_seq <- which(clusters_vector==i)
    consenso <-  consensus_sequence(sequences = sequences[index_seq,],
                                    probabilities = weights_vector[index_seq,],cutoff_prob = cutoff_prob)
    index    <- which(consenso!=reference_sequence)
    if(length(index)>0){
      labels   <- paste0(reference_sequence[index],index,consenso[index])
      labels <- paste(labels, collapse = " ")
    }else{labels <- NA}

    consensus_for_clusters[i,1] <- i
    consensus_for_clusters[i,2] <- labels
  }
  return(consensus_for_clusters)
}
consensus_for_subsets <- function(path_name,nseq,nrep,weighting,cutoff_prob=0,
                                  reference_sequence,reference_mutant,
                                  fragment_homopolymers=F,sort_homopolymers=F,
                                  auxFile=NULL){
  # Modify auxFile name
  if(!is.null(auxFile)){
    if(fragment_homopolymers){auxFile <- paste0(auxFile,"_FragHP")}
    if(sort_homopolymers){auxFile <- paste0(auxFile,"_SortHP")}
    if(!file.exists(auxFile)){cat("Sample Id","Subset size","Mutations", "\n", sep="\t", file=auxFile)}
  }
  if(fragment_homopolymers | sort_homopolymers){
    homopolymer_positions_all <- detect_homopolymer_positions(reference_sequence)
    homopolymer_positions     <- homopolymer_positions_all[sapply(homopolymer_positions_all, length)>1]
  }
  if(fragment_homopolymers){
    reference_sequence <- sapply(homopolymer_positions_all, function(x){paste0(reference_sequence[x], collapse = "")})
    reference_mutant   <- sapply(homopolymer_positions_all, function(x){paste0(reference_mutant[x], collapse = "")})
  }

  #load sequences
  data       <- load_sequences_from_path(path_name)
  seqs       <- t(sapply(data, function(x){x$nucleotide}))
  if(weighting=="raw"){       weights   <- t(sapply(data, function(x){-10*log10(x$p_error_minION)}))      }
  if(weighting=="naive"){     weights   <- t(sapply(data, function(x){-10*log10(1-x$p_right_null_model)}))    }
  if(weighting=="priors"){    weights   <- t(sapply(data, function(x){-10*log10(1-x$p_right_priors_model)}))  }
  weights[is.infinite(weights)]=93

  #Order sequences and probabilities in homopolyers, putting gaps in the end
  if(fragment_homopolymers){
    for(j in length(homopolymer_positions):1){
      pos <- homopolymer_positions[[j]]

      seqs[,pos[1]] <- apply(seqs[,pos], 1, paste0, collapse="")
      seqs          <- seqs[,-c(pos[-1])]

      weights[,pos[1]] <- apply(weights[,pos], 1, prod)
      weights          <- weights[,-c(pos[-1])]
    }
  }
  if(sort_homopolymers){
    for(k in 1:nrow(seqs)){
      for(j in seq_along(homopolymer_positions)){
        hp_pos    <- homopolymer_positions[[j]]
        hp_seq    <- seqs[k,hp_pos]
        anygaps <- any(hp_seq=="-")
        if(!anygaps){next()}
        horder <- c(which(hp_seq!="-"), which(hp_seq=="-"))
        seqs[k,hp_pos]     <- hp_seq[horder]
        weights[k,hp_pos]  <- weights [k,hp_pos[horder]]
      }
    }
  }

  #make subsets of sequences and compute consensus
  subset    <- rep(nseq, each=nrep)
  l         <- length(subset)
  empty_vec <- rep(NA,l)
  results_subset <- data.frame(path=rep(path_name,l), weights=rep(weighting,l),nseq=empty_vec,n_diferences=empty_vec)
  n<-0
  for(p in seq_along(subset)){
    ind_seq <- sample(1:nrow(seqs), subset[p])
    s_min   <- consensus_sequence(sequences = seqs[ind_seq,],probabilities = weights[ind_seq,],cutoff_prob = cutoff_prob)

    #store
    n<-n+1
    if(!is.null(auxFile)){
      aux_indM <- which(s_min!=reference_mutant)
      if(length(aux_indM)>0){
        cat(n,subset[p],paste(reference_mutant[aux_indM],aux_indM, s_min[aux_indM],sep=""),'\n',sep="\t", file=auxFile, append = T)
      }else{
        cat(n,subset[p],"-",'\n',sep="\t", file=auxFile, append = T)
      }
    }
    results_subset$nseq[n]             <- subset[p]
    results_subset$n_diferences[n]     <- sum(s_min!=reference_mutant, na.rm=T)
  }
  results_subset$equal_ref    <- results_subset$n_diferences==0

  #Average
  res_means <- results_subset %>% group_by(path, nseq)%>%
    dplyr::summarise(ntest=n(),
                     mean_differences = (mean(n_diferences, na.rm=T)),
                     sd_differences=sd(n_diferences, na.rm=T),
                     right_sequences = sum(equal_ref, na.rm = T))


  return(res_means)
}

#COLORS
red      <- rgb(218,0,0,maxColorValue = 255)
green    <- rgb(0,122,0,maxColorValue = 255)
red_tr   <- rgb(230,0,0,alpha = 200, maxColorValue = 255)
blue_tr  <- rgb(0,0,218,alpha = 250, maxColorValue = 255)
green_tr <- rgb(0,122,0,alpha = 200, maxColorValue = 255)
blue     <- rgb(0,0,250, maxColorValue = 255)


# ######## FUNCTIONS QUE QUIZAS HALLA QUE TIRAR ###########
# plot_nmut_vs_nseq           <- function(data, bc_id, id, known_muts){
#   g <- ggplot(data %>% filter(bc==bc_id), aes(x=nseq))+
#     ggtitle(name_title[id])+
#     scale_x_continuous(breaks=seq(0,100,by=10))+
#     scale_y_continuous(breaks=seq(10,100,by=10))+
#     theme_bw()+ylab("Number of mutations")+xlab("Number of sequences")+
#     geom_hline(aes(yintercept=0), linetype=3) +
#     geom_point(aes(y=mean_minION)) +  geom_line(aes(y=mean_minION)) +
#     geom_errorbar(aes(ymin=mean_minION-sd_minION/2, ymax=mean_minION+sd_minION/2),width=0.5)+
#     geom_point(aes(y=mean_corrected), col=2) +  geom_line(aes(y=mean_corrected), col=2) +
#     geom_errorbar(aes(ymin=mean_corrected-sd_corrected/2, ymax=mean_corrected+sd_corrected/2), col=2,width=0.5)
#   return(g)
# }
# plot_exactconsensus_vs_nseq <- function(data, bc_id, id){
#   g <-   ggplot(data %>% filter(bc==bc_id), aes(x=nseq)) +
#     theme_bw()+ylab("Number of correct sequences")+xlab("Number of sequences")+
#     ggtitle(name_title[id])+
#     geom_point(aes(y=right_minION)) +  geom_line(aes(y=right_minION)) +
#     geom_point(aes(y=right_corrected), col=2) +  geom_line(aes(y=right_corrected), col=2) #+
#   #ylim(c(0,100))
#   return(g)
# }

