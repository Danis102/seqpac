#' Checks PAC format
#'
#' \code{PAC_check} Quickly checks rownames in PAC.
#'
#' Given a PAC object the function will cross check the rownames and colnames
#' across Pheno, Anno and Counts.
#' 
#' @family PAC generation
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object.
#' 
#' @param sample_check Can be provided either as a data frame or the name of a
#'   data frame in PAC. Will check both the row and column names of the object
#'   for identical sample names in PAC. A simple character vector can also be
#'   provided.
#' 
#' @param seq_check Can be provided either as a data frame or the name of an
#'   data frame in PAC. Will check only the row names of the object for
#'   identical sequence names in PAC. A simple character vector can also be
#'   provided.
#' 
#' @return If PAC format is correct will return logical TRUE, else it will throw an error. 
#'
#' @examples
#' 
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' 
#' PAC_check(pac)  # TRUE
#' 
#' # Remove first row in Anno
#' new_Anno <- pac$Anno[-1,]               
#' 
#' 
#' #PAC_check(pac, seq_check=new_Anno) # Error   
#' 
#' # Add to pac
#' pac$Anno <- new_Anno                    
#' #PAC_check(pac) # Error                      
#' 
#' ## Reload data
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' PAC_check(pac)  # TRUE
#' 
#' # Remove a sample column in normalized table 
#' new_norm <- pac$norm$cpm[,-1]           
#' PAC_check(pac, seq_check=new_norm)      # All sequences are good
#' #PAC_check(pac, sample_check=new_norm)  # Error
#' 
#' pac$norm$cpm <- new_norm                # Add to pac
#' #PAC_check(pac)                         # Error 
#' 
#' ## Reload data
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' PAC_check(pac)  # TRUE
#'
#' seq_names <- rownames(pac$Counts)
#' sample_names <- colnames(pac$Counts) 
#' 
#' PAC_check(pac, seq_check=seq_names)
#' PAC_check(pac, sample_check=sample_names) 
#' 
#' #PAC_check(pac, seq_check=rev(seq_names))
#' #PAC_check(pac, sample_check=rev(sample_names)) 
#' 
#' 
#' @export

PAC_check <- function(PAC, sample_check=NULL, seq_check=NULL){
  sampl <- paste0("Please make sure that all samples are represented, ",
           "\nnamed and ordered the correct way in both tables.") 
  seqs <-  paste0("Please make sure that all sequences are represented, ",
                  "\nnamed and ordered the correct way in both tables.")
  if(!identical(rownames(PAC$Pheno), colnames(PAC$Counts))){
    stop(paste0("\nSample names in Pheno (row names) and Counts",
                "(column names) are not identical.\n", sampl))
  }
  if(!identical(rownames(PAC$Anno), rownames(PAC$Counts))){
    stop(paste0("\nSequence names in Anno (row names) and Counts",
                "(row names) are not identical.\n", seqs))
  }
  if(length(PAC$norm) > 0){ 
    norm_logi_row  <- any(unlist(lapply(PAC$norm, function(x){
      !identical(rownames(PAC$Counts), rownames(x))
    })))
    norm_logi_col  <- any(unlist(lapply(PAC$norm, function(x){
      !identical(colnames(PAC$Counts), colnames(x))
    })))
    if(norm_logi_row){
      stop("\nRow names in normalized tables of PAC$norm are not identical",
           "\nwith row names in Counts.")
    }
    if(norm_logi_col){
      stop("\nColumn names in normalized tables of PAC$norm are not identical",
           "\nwith column names in Counts.")
    }
  }
  if(length(PAC$summary) > 0){ 
    sum_logi_row  <- any(!unlist(
      lapply(lapply(PAC$summary, rownames), function(x){
        identical(rownames(PAC$Counts), x)
    })))
    if(sum_logi_row){
      stop("\nRow names in summary tables of PAC$summary are not identical",
           "\nwith row names in Counts.")
    }
  } 
  res <- NULL
  ### sample_check #############################################################          
  if(!is.null(sample_check)){
    ## Provided table
    if(class(sample_check) %in% c("data.frame", "matrix")){
      if(any(rownames(sample_check)[1] %in% rownames(PAC$Pheno))){
        if(!identical(rownames(PAC$Pheno), rownames(sample_check))){
          stop(paste0("\nSample (row) names in new table are not identical",
                      "\nwith sample names in the main PAC.\n", sampl))
          }
      } 
      if(any(colnames(sample_check)[1] %in% rownames(PAC$Pheno))){
        if(!identical(rownames(PAC$Pheno), colnames(sample_check))){
          stop(paste0("\nSample (column) names in new table are not identical",
                      "\nwith sample names in the main PAC.\n", sampl))
          }
      }
    }else{
      ## Provided name of table in PAC
      if(is.character(sample_check)){
        if(length(sample_check)==1){
          whr  <- which(unlist(lapply(PAC, function(x){  
            if(!class(x)=="list"){return(FALSE)}
            if(class(x)=="list"){any(grepl(sample_check, names(x)))}})))
          if(length(whr)==1){
            if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), 
                         rownames(PAC[[whr]][[sample_check]])))==2){
              if(!identical(rownames(PAC$Pheno), 
                            rownames(PAC[[whr]][[sample_check]]))){
                stop(paste0(
                  "\nSample (row) names in new table are not identical",
                  "\nwith sample names in the main PAC.\n", sampl))
                }
            }else{
              if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), 
                           colnames(PAC[[whr]][[sample_check]])))==2){
                if(!identical(rownames(PAC$Pheno), 
                              colnames(PAC[[whr]][[sample_check]]))){
                  stop(paste0(
                    "\nSample (column) names in new table are not identical",
                    "\nwith sample names in the main PAC.\n", sampl))
                  }
              }else{   
                stop(
                  "\nPAC_check did find a matching object name, but could not",
                  "\nfind the samples.\nIf you meant to check sequences and not",
                  "\nsamples, provide a seq_check object instead.")
              }
             }
            }
          
        }else{
          
          # Provided names only       
          if(length(sample_check)>1){ 
            if(!identical(rownames(PAC$Pheno), sample_check)){
              stop(
                paste0(
                  "\nSample names are not identical with sample names in",
                  "\nthe main PAC.\n", sampl))
            }
          }
        }
        
      }else{
        stop(
          "\nPAC_check did not recognize the sample_check input.",
          "\nDid you provide a data.frame, character vector or the",
          "\nname of a PAC object?")
        }  
    }}
  ### seq_check #############################################################          
  if(!is.null(seq_check)){
    ## Provided table
    if(class(seq_check) %in% c("data.frame", "matrix")){
      if(!identical(rownames(PAC$Anno), rownames(seq_check))){
        stop(
          paste0(
            "\nSequence (row) names in new table are not identical with",
            "\nsequence names in the main PAC.\n", sampl))
        }
    }else{
      ## Provided name of table in PAC
      if(is.character(seq_check)){
        if(length(seq_check)==1){
          whr  <- which(unlist(lapply(PAC, function(x){  
            if(!class(x)=="list"){return(FALSE)}
            if(class(x)=="list"){any(grepl(seq_check, names(x)))}})))
          if(length(whr)==1){
            if(!identical(rownames(PAC$Anno), rownames(PAC[[whr]][[seq_check]]))){
              stop(
                paste0(
                  "\nSequence (row) names in new table are not identical",
                  "\nwith sequence names in the main PAC.\n", sampl))
              }
          }
        }else{
          ## Provided names only        
          if(length(seq_check)>1){ 
            if(!identical(rownames(PAC$Anno), seq_check)){
              stop(paste0(
                "\nSequence names are not identical with sequence names",
                "\nin the main PAC.\n", sampl))
              }
          }}                   
      }else{
        stop(
          "\nThe PAC_check did not recognize the sample_check input.",
          "\nDid you provide a data.frame, character vector or the",
          "\nname of a PAC object?")
        }  
    }}
  return(TRUE)
}
