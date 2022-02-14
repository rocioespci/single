#' Save evaluated SINGLE model
#'
#' This is an auxiliary function in single package. It saves the results of the evaluation in a txt file.
#' @param input_sam_file File containing the counts per position returned by samtools mpileup
#' @param single_fits Results of the SINGLE model as returned by single_train(). It can be either the output data.frame or the saved file.
#' @param output_txt_file String. Prefix for output files
#' @param ref_seq Reference sequence: vector of characters, as returned by load_ref_seq
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param comment.char Character used to indicate comments in the input_sam_file
#' @param gaps_weights One of "minimum","none","mean". How to assign qscores to deletions.
#' @param remove.first.reverse Logical. Do you want to remove first position in the reverse reads?
#' @import dplyr
#' @importFrom utils read.table
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export save_txt_file
#' @return Writes txt file
save_txt_file <- function(input_sam_file, single_fits,output_txt_file,ref_seq,
                    pos_start=NULL, pos_end=NULL,comment.char = "@",gaps_weights,remove.first.reverse=FALSE){
    #Verify inputs
    if(is.null(pos_start)){pos_start=1;warning("save_txt_file: pos_start set to 1")}
    if(is.null(pos_end)){pos_end=length(ref_seq);
    warning("save_txt_file: pos_end set to legth(referece_sequence)+pos_start (",length(ref_seq)+pos_start,")\n")}
    if(!file.exists(input_sam_file)){stop("save_txt_file: file ", input_sam_file," does not exist")}
    if(file.exists(output_txt_file)){warning("save_txt_file: overwritting ", output_txt_file,"\n"); file.remove(output_txt_file)}

    # Load data and define general parameters
    if(is.character(single_fits)){
        if(!file.exists(single_fits)){
            stop("save_txt_file: file ", single_fits," does not exist")
        }
        qtable = utils::read.table(single_fits, header=TRUE)
    }else{
        qtable = single_fits
    }
    cigar_reference_counts   <- c("M","D","N","=","X")
    cigar_query_counts <- c("M","I","S","=","X")

    cat(c("SeqID","position","nucleotide","isWT","original_quality","p_right_priors_model\n"), file=output_txt_file, append=FALSE)

    connection <- file(input_sam_file)
    open(connection)
    good_line <- TRUE
    while(good_line){
        #read line and check it
        line <- readLines(connection, n=1)
        if(length(line)==0){good_line <- FALSE;next()}  #file is not over
        if(substr(line,1,1)==comment.char){next()}  #line it's not header

        #parse line into something friendly
        line_data <- scan(text = line,what=character(1), quiet = TRUE,quote = "")

        pos_ini_aliref  <-  as.numeric(line_data[4])
        cigar    <-  line_data[6]
        sequence <-  strsplit(line_data[10],split="")[[1]]
        qscores  <-  strsplit(line_data[11],split="")[[1]]

        #Obtain CIGAR vector
        cigar.table      <- data.frame(n=as.numeric(strsplit(cigar, split="[MIDNSHP=X]")[[1]]),Op= strsplit(cigar, split="[0-9]+")[[1]][-1])
        cigar_vec        <- unlist(apply(cigar.table,1, function(x){rep(x[2],x[1])}))
        names(cigar_vec) <- NULL

        #Data frame for the original sequence, qscores, and reference for the model (full)
        data_seq <- data.frame(cigar=cigar_vec, nucleotide=NA,qscore=NA,quality=NA,pos_in_samalig=NA,position=NA,wt_seq=NA)
        if(remove.first.reverse){
            if(line_data[2]=="16"){
                data_seq <- data_seq[-1,]
            }
        }
        rows_with_query     <- data_seq$cigar%in%cigar_query_counts
        rows_with_reference <- data_seq$cigar%in%cigar_reference_counts

        if(sum(rows_with_query)!= length(sequence)){warning("OJO", line_data[1])}

        data_seq$nucleotide[rows_with_query]  <- sequence     #Nucleotide
        data_seq$qscore[rows_with_query]      <- qscores      #Qscore
        data_seq$quality = dplyr::recode(data_seq$qscore, !!! stats::setNames(ascii$Q,  as.character(ascii$Symbol)))
        data_seq$pos_in_samalig[rows_with_reference] <- pos_ini_aliref + c(0:(sum(rows_with_reference)-1))  #position - according to aligning reference

        data_seq <- data_seq[data_seq$cigar!="H",]
        ## Complete or cut sequence to relevant nucleotides

        #1. Check if there is any information in positions of interest. If not discard.
        # Ref        ***
        # Read  ****
        pos_range <- range(data_seq$pos_in_samalig, na.rm=TRUE)
        if(pos_range[2]<pos_end | pos_range[1]>pos_end){next()}

        #2. Add to reference the soft clipping sequences at the begining and at the end
        # Ref    ********
        # Read   --****--    (--- not matching sequences)
        if(is.na(data_seq$pos_in_samalig[1])){
            pos_not_na <- which(!is.na(data_seq$pos_in_samalig))[1]
            data_seq$pos_in_samalig[seq_len(pos_not_na-1)] <- data_seq$pos_in_samalig[pos_not_na]-rev(seq_len(pos_not_na-1))
        }
        if(is.na(data_seq$pos_in_samalig[nrow(data_seq)])){
            pos_not_na <- which(!is.na(data_seq$pos_in_samalig))
            pos_not_na <- pos_not_na[length(pos_not_na)]

            data_seq$pos_in_samalig[(pos_not_na+1):nrow(data_seq)] <- data_seq$pos_in_samalig[pos_not_na]+seq_along((pos_not_na+1):nrow(data_seq))
        }

        #3. Read is longer than reference. Keep positions between pos_start and pos_end
        # Ref     *****
        # Read  *********
        row_initial <- which(data_seq$pos_in_samalig == pos_start)
        if(length(row_initial)==0){row_initial <- 1}
        row_final <- which(data_seq$pos_in_samalig == pos_end)
        if(length(row_final)==0){row_final <- nrow(data_seq)}
        data_seq <- data_seq[row_initial:row_final,]

        #Positions according to reference and reference seq
        data_seq$position <- data_seq$pos_in_samalig- pos_start + 1  #position - according to single reference
        data_seq$wt_seq[!is.na(data_seq$position)]   <- ref_seq[data_seq$position[!is.na(data_seq$position)]]

        #4. Read is shorter than reference. Add rows.
        # Ref   *********
        # Read    ****
        pos_one <- data_seq$position[1]-1
        if(pos_one>0){
            data_aux <-  data.frame(cigar="D",nucleotide="-",
                                    qscore=NA,quality=NA,pos_in_samalig=NA,
                                    position=seq_len(pos_one),
                                    wt_seq=ref_seq[seq_len(pos_one)])
            data_seq <- rbind(data_aux,data_seq)
        }
        pos_last <- data_seq$position[nrow(data_seq)]
        total_pos <- length(ref_seq)
        if(pos_last < total_pos){
            data_aux <-  data.frame(cigar="D",nucleotide=NA,
                                    qscore=NA,quality=NA,pos_in_samalig=NA,
                                    position=(pos_last+1):total_pos,
                                    wt_seq=ref_seq[(pos_last+1):total_pos])
            data_seq <- rbind(data_seq,data_aux)
        }


        ## Remove insertions
        data_seq <- data_seq[data_seq$cigar != "I",]

        ## Replace deletions by - and provide a Qscore
        del_positions <- which(data_seq$cigar=="D")
        data_seq$nucleotide[del_positions] <- "-"
        nr <- nrow(data_seq)

        #If they are at the beginning, I replace by the first Qscore
        if(length(del_positions)>0){
            if(del_positions[1]==1){
                n = length(del_positions)
                if(n==1){
                    data_seq$quality[1] <- data_seq$quality[2]
                    del_positions <- del_positions[-1]
                }else{
                    ind_dif <- which(del_positions[seq(2,n)] != del_positions[seq_len(n-1)]+1)[1]
                    if(is.na(ind_dif)) {ind_dif <- n}
                    data_seq$quality[seq_len(ind_dif)] <- data_seq$quality[ind_dif+1]
                    del_positions <- del_positions[-seq_len(ind_dif)]
                }
            }
        }

        #If they are at the end, I replace by the last Qscore
        if(length(del_positions)>0 ){
            n = length(del_positions)
            if(del_positions[n]==nr){
                #if all the remaining deletions are the end of the string, replace all
                if(all(del_positions== nr-seq(n-1,0))){
                    data_seq$quality[del_positions] <- data_seq$quality[min(del_positions)-1]
                    del_positions <- NULL
                }else{
                    #detect which are the values on the end of the string
                    del_positions_rev <- rev(del_positions)
                    ind_dif           <- which(del_positions_rev[seq(2,n)] != del_positions_rev[seq_len(n-1)]-1)[1]
                    index             <- rev(del_positions_rev[seq_len(ind_dif)])
                    #replace those
                    data_seq$quality[index] <- data_seq$quality[min(index)-1]
                    del_positions <- del_positions[- seq(n-ind_dif+1,n)]
                }
            }
        }


        if(length(del_positions)>0){
            if(gaps_weights!='none'){
                for(j in del_positions){
                    #skip if there are only gaps before or after
                    start <- data_seq$nucleotide[seq_len(j)]
                    end   <- data_seq$nucleotide[seq(j,nr)]

                    #Replace ~ by the mean of the next non gaps qualities:
                    #-if all positions before are gaps, replace by the first Qscore no gap after
                    if(all(start == "-")){
                        nend   <- which(end[-1]!="-")[1]
                        data_seq$quality[j] <- data_seq$quality[j+nend]
                        next()
                    }
                    #-if all positions after are gaps, replace by the first Qscore no gap before
                    if(all(end == "-")){
                        nstart <- which(rev(start)[-1]!="-")[1]
                        data_seq$quality[j] <- data_seq$quality[j-nstart]
                        next()
                    }
                    #-if it has nucleotides before and after, replace by smthg of the Qscores of nearest neighbours
                    nstart <- which(rev(start)[-1]!="-")[1]
                    nend   <- which(end[-1]       !="-")[1]

                    #Mean of probabilities
                    aux.qual       <- data_seq$quality[c(j-nstart, j+nend)]
                    if(gaps_weights=='mean'){
                        mean.qscore    <- round(mean(aux.qual))
                    }
                    if(gaps_weights=='minimum'){
                        mean.qscore    <- min(aux.qual)
                    }
                    data_seq$quality[j] <- mean.qscore
                }
            }else{
                #* All gaps get ~ score
                data_seq$quality[del_positions] <- 0
            }
        }

        #Use fitted data to re-evaluate the sequence
        df <- data_seq %>%
            left_join(qtable %>% select(c("nucleotide","position","quality","isWT","p_right_priors_model")),
                    by = c("nucleotide","position","quality"))                              # merge with p_right from correction
        df$p_right_priors_model[which(df$isWT==TRUE)] <- 1-10^(-df$quality[which(df$isWT==TRUE)]/10)     # For wildtype keep original Qscore (or values that did not occurred)

        ## Save
        df <- df %>%
            dplyr::mutate(SeqID=line_data[1])%>%
            dplyr::select(.data$SeqID,.data$position,.data$nucleotide,.data$isWT,.data$quality,.data$p_right_priors_model) %>%
            dplyr::rename(original_quality = .data$quality) %>%
            dplyr::filter(.data$position <= (pos_end - pos_start+1))

        utils::write.table(df, append = TRUE, file=output_txt_file, row.names = FALSE,col.names = FALSE)
    }
    close(connection)
    return(0)
}
