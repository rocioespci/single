#' Parse nucleotide' counts file
#'
#' This is an auxiliary function in single package. Reshapes file returned by samtools mpileup into a data.frame.
#' @param input_file File containing the counts per position returned by samtools mpileup
#' @param output_file Name of file to save parsed data.
#' @param pos_start Numeric. Position to start analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param pos_end  Numeric. Position to stop analyzing, counting starts from 1 and it refers to reference used for minimap2 alignment.
#' @param save Logical. Should data be saved in a output_file?
#' @return data frame with columns position nucleotide quality counts
#' @importFrom stringr str_locate_all
#' @importFrom rlang .data
#' @importFrom utils write.table
#' @export parse_nucleotides_per_qscore
#' @examples
#' input_file = system.file("extdata", "example_train_countsPNQ.txt", package = "single")
#' parse_nucleotides_per_qscore(input_file=input_file,output_file=NA, pos_start=1,pos_end=10)
parse_nucleotides_per_qscore   <- function(input_file,output_file, pos_start=NULL,pos_end=NULL, save=FALSE){
    # Load data
    data_input <- readLines(input_file)
    data_input <- strsplit(data_input, split="\t")
    data_input <- do.call(rbind,data_input)
    data_input <- as.data.frame(data_input, stringsAsFactors = FALSE)
    if(ncol(data_input)!=6){stop('count_nucleotides_per_score: input_file does not have the correct size: it should have 6 columns.')}
    # Make a nice dataframe
    data_input <- data_input %>%
        dplyr::rename(ref_seq_id = .data$V1, position=.data$V2,ref_base = .data$V3,n_reads=.data$V4,bases=.data$V5,qscores=.data$V6) %>%
        dplyr::mutate(position = as.numeric(.data$position), n_reads=as.numeric(.data$n_reads))
    # Filter to positions between pos_start and pos_end
    if(!is.null(pos_start) | !is.null(pos_end)){
        if(is.null(pos_start)){pos_start <- 1}
        if(is.null(pos_end)){pos_end <- max(data_input$position)}
        data_input <- data_input %>%
            dplyr::filter(.data$position >=pos_start & .data$position<= pos_end) %>%
            dplyr::mutate(position = .data$position - pos_start+1)
    }
    #All to upper case
    data_input$bases <- toupper(data_input$bases)
    # Reshape data to counts table
    if(save){cat(x=c("position", "nucleotide","quality","counts\n"),file=output_file,append = FALSE)}
    for(i in seq_len(nrow(data_input))){
        n.bases <- data_input$n_reads[i]
        length.bases.col   <- nchar(data_input$bases[i])
        length.qscores.col <- nchar(data_input$qscores[i])
        if(n.bases==0){next()}
        if(length.qscores.col!=n.bases){stop('Different Qscores than length')}
        # Remove characters that are not bases (such as starting or finishing characters)
        if(length.bases.col!=n.bases){
            aux <- toupper(data_input$bases[i])
            #Detect ^. and ^, , character that indicates initial position of a read, and remove them
            match_start <-  gsub(replacement = "",x = aux, pattern = "\\^.",fixed = FALSE) ## hat and the next one
            aux <- match_start
            match_start <-  gsub(replacement = "",x = aux, pattern = "\\^,",fixed = FALSE) ## hat and the next one
            aux <- match_start
            #Detect $, character that indicates final position of a read, and remove them
            match_end <-     gsub(replacement = "",x = aux, pattern = "$",fixed = TRUE) ## hat and the next one
            aux <- match_end

            #Detect -, character indicating a deletion
            match_minus <-  stringr::str_locate_all(string = aux, pattern = "-[0-9]+")
            if(nrow(match_minus[[1]])>0){
                length_match_minus <- apply(match_minus[[1]], 1, function(x){as.numeric(substr(aux,x[1]+1,x[2]))})
                match_minus[[1]][,2] <- match_minus[[1]][,2] + length_match_minus
            }
            #Detect +, character indicating an insertion
            match_plus <-  stringr::str_locate_all(string = aux, pattern = "\\+[0-9]+")
            if(nrow(match_plus[[1]])>0){
                length_match_plus <- apply(match_plus[[1]], 1, function(x){as.numeric(substr(aux,x[1]+1,x[2]))})
                match_plus[[1]][,2] <- match_plus[[1]][,2] + length_match_plus
            }
            #Remove detected - and + (insertion and deletion flags)
            match_all <- rbind(match_minus[[1]], match_plus[[1]])
            if(nrow(match_all)>0){
                remove_characters <- apply(match_all,1,function(x){seq(x[1],x[2],by=1)})
                remove_characters <- sort(unlist(remove_characters))
                aux_out <- strsplit(aux, split="")[[1]]
                aux_out <- aux_out[-c(remove_characters)]
                aux <- paste0(aux_out, collapse = "")
            }
            data_input$bases[i] <- aux
        }
        if(nchar(data_input$bases[i])!=n.bases){stop("Line",i)}

        bases.vec <- strsplit(data_input$bases[i], split="")[[1]]
        qscores <- strsplit(data_input$qscores[i], split="")[[1]]
        df <- data.frame(nucleotide=bases.vec,qualitySymbol=qscores, stringsAsFactors = FALSE)
        df.counts <- df %>%
            dplyr::group_by(.data$nucleotide,.data$qualitySymbol) %>%
            dplyr::tally(name = "counts")%>%
            dplyr::mutate(position=data_input$position[i])
        df.counts$quality = recode(df.counts$qualitySymbol, !!! setNames(ascii$Q,  as.character(ascii$Symbol)))
        df.counts$nucleotide[df.counts$nucleotide=="*"] <- "-"

        if(i==1){
            df_final <- df.counts[,c("position", "nucleotide","quality","counts")]
        }else{
            df_final <-  df_final %>% rbind(df.counts[,c("position", "nucleotide","quality","counts")])
        }
        if(save){utils::write.table(x=df.counts[,c("position", "nucleotide","quality","counts")],
                             file=output_file,
                             append = TRUE, col.names = FALSE,row.names = FALSE)
        }
    }
    return(df_final)
}
