#' Checks PAC format
#'
#' \code{PAC_check} Quickly checks rownames in PAC.
#'
#' Given a PAC object the function will cross check the rownames and colnames
#' across Pheno, Anno and Counts.
#' 
#' @family
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
#' load("~/OneDrive/Programmering/Programmering/Pipelines/Drosophila/Pipeline_3.1/seqpac/dm_test_PAC.Rdata")
#' 
#' PAC_filt$summary$pheno <- PAC_filt$Pheno
#' PAC_filt$summary$pheno_trans <- t(PAC_filt$Pheno)
#' PAC_filt$summary$pheno_trans2 <- t(PAC_filt$Pheno)[,-1]
#' PAC_filt$summary$pheno_rev <- PAC_filt$Pheno[rev(seq(PAC_filt$Pheno)),]
#' 
#' 
#' PAC_check(PAC_filt, sample_check="pheno") # TRUE
#' PAC_check(PAC_filt, sample_check="pheno_trans") # TRUE (Automatically looking for samples in columns instead)
#' PAC_check(PAC_filt, sample_check="pheno_trans2") #Error with hint
#' PAC_check(PAC_filt, sample_check="pheno_rev") #Error not identical
#' 
#' PAC_filt$summary$anno <- PAC_filt$Anno
#' PAC_filt$summary$anno_trans <- t(PAC_filt$Anno)
#' PAC_filt$summary$anno_mis <- PAC_filt$Anno[-567,]
#' 
#' PAC_check(PAC_filt, seq_check="anno") # TRUE
#' PAC_check(PAC_filt, seq_check="anno_trans") # Error not identical
#' PAC_check(PAC_filt, seq_check="anno_mis") # Error not identical
#' 
#' PAC_check(PAC_filt, seq_check=PAC_filt$summary$anno) #TRUE
#' PAC_check(PAC_filt, seq_check=PAC_filt$summary$anno_trans)  # Error not identical
#' PAC_check(PAC_filt, seq_check=PAC_filt$summary$pheno) #Error not identical
#' PAC_check(PAC_filt, sample_check=PAC_filt$summary$pheno_rev) # Error not identical
#' PAC_check(PAC_filt, sample_check=PAC_filt$summary$pheno_trans) # TRUE (Automatically looking for samples in columns instead)
#' PAC_check(PAC_filt, sample_check=rownames(PAC_filt$summary$pheno)) # TRUE
#' PAC_check(PAC_filt, sample_check=rownames(PAC_filt$summary$pheno_trans)) # Error not identical
#' 
#' PAC_check(PAC_filt, sample_check=as.integer(13))
#' PAC_check(PAC_filt, sample_check=list(13,45))
#' 
#' @export

PAC_check <- function(PAC, sample_check=NULL, seq_check=NULL){
                    sampl <- "Please make sure that all samples are represented, named and ordered the correct way in both tables." 
                    seqs <-  "Please make sure that all sequences are represented, named and ordered the correct way in both tables."
                    if(!identical(rownames(PAC$Pheno), colnames(PAC$Counts))){stop(paste0("\n  Sample names in Pheno (row names) and Counts (column names) are not identical.\n  ", sampl))}
                    if(!identical(rownames(PAC$Anno), rownames(PAC$Counts))){stop(paste0("\n  Sequence names in Anno (row names) and Counts (row names) are not identical.\n  ", seqs))}
                    if(length(PAC$norm) > 0){ 
                                            norm_logi_row  <- any(!unlist(lapply(lapply(PAC$norm, rownames), function(x){identical(rownames(PAC$Counts), x)})))
                                            norm_logi_col  <- any(!unlist(lapply(lapply(PAC$norm, colnames), function(x){identical(colnames(PAC$Counts), x)})))
                                            if(norm_logi_row){stop("Row names in normalized tables of PAC$norm are not identical with row names in Counts.")}
                                            if(norm_logi_col){stop("Column names in normalized tables of PAC$norm are not identical with column names in Counts.")}
                                            }
                    if(length(PAC$summary) > 0){ 
                                            sum_logi_row  <- any(!unlist(lapply(lapply(PAC$summary, rownames), function(x){identical(rownames(PAC$Counts), x)})))
                                            if(norm_logi_row){stop("Row names in summary tables of PAC$summary are not identical with row names in Counts.")}
                                            } 
                    res <- NULL
                    ### sample_check #############################################################          
                    if(!is.null(sample_check)){
                            ## Provided table
                            if(class(sample_check) %in% c("data.frame", "matrix")){
                                                          if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), rownames(sample_check)))==2){
                                                                              if(!identical(rownames(PAC$Pheno), rownames(sample_check))){stop(paste0("\n  Sample (row) names in new table are not identical with sample names in the main PAC.\n  ", sampl))}
                                                          }else{
                                                          if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), colnames(sample_check)))==2){
                                                                              if(!identical(rownames(PAC$Pheno), colnames(sample_check))){stop(paste0("\n  Sample (column) names in new table are not identical with sample names in the main PAC.\n  ", sampl))}
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
                                                       if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), rownames(PAC[[whr]][[sample_check]])))==2){
                                                                                      if(!identical(rownames(PAC$Pheno), rownames(PAC[[whr]][[sample_check]]))){stop(paste0("\n  Sample (row) names in new table are not identical with sample names in the main PAC.\n  ", sampl))}
                                              }else{
                                                       if(sum(grepl(paste(rownames(PAC$Pheno)[1:2], collapse="|"), colnames(PAC[[whr]][[sample_check]])))==2){
                                                                                      if(!identical(rownames(PAC$Pheno), colnames(PAC[[whr]][[sample_check]]))){stop(paste0("\n  Sample (column) names in new table are not identical with sample names in the main PAC.\n  ", sampl))}
                                                       }else{   stop("\nPAC_check did find a matching object name, but could not find the samples.\nIf you ment to check sequences and not samples, provide a seq_check object instead.")}}}
                                            
                                    }else{
                                    
                                    # Provided names only       
                                    if(length(sample_check)>1){ 
                                                       if(!identical(rownames(PAC$Pheno), sample_check)){stop(paste0("\n  Sample names are not identical with sample names in the main PAC.\n  ", sampl))}
                                    }
                                }
                            
                            }else{stop("\nPAC_check did not recognize the sample_check input.\nDid you provide a data.frame, character vector or the name of a PAC object?")}  
                    }}
                    ### seq_check #############################################################          
                    if(!is.null(seq_check)){
                            ## Provided table
                            if(class(seq_check) %in% c("data.frame", "matrix")){
                                                              if(!identical(rownames(PAC$Anno), rownames(seq_check))){stop(paste0("\n  Sequence (row) names in new table are not identical with sequence names in the main PAC.\n  ", sampl))}
                            }else{
                            ## Provided name of table in PAC
                            if(is.character(sample_check)){
                              if(length(seq_check)==1){
                                    whr  <- which(unlist(lapply(PAC, function(x){  
                                                              if(!class(x)=="list"){return(FALSE)}
                                                              if(class(x)=="list"){any(grepl(seq_check, names(x)))}})))
                                    if(length(whr)==1){
                                                              if(!identical(rownames(PAC$Anno), rownames(PAC[[whr]][[seq_check]]))){stop(paste0("\n  Sequence (row) names in new table are not identical with sequence names in the main PAC.\n  ", sampl))}
                                    }
                            }else{
                            ## Provided names only        
                            if(length(seq_check)>1){ 
                                               if(!identical(rownames(PAC$Anno), seq_check)){stop(paste0("\n  Sequence names are not identical with sequence names in the main PAC.\n  ", sampl))}
                              }}                   
                           }else{stop("\nThe PAC_check did not recognize the sample_check input.\nDid you provide a data.frame, character vector or the name of a PAC object?")}  
                    }}
                    return(TRUE)
}
                                      