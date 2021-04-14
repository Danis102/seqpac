#' Binary converter for PAC_filtsep  
#'
#' \code{filtsep_bin} Converts PAC_filtsep data.frame output into a binary
#' (hit-no-hit) data.frame
#' 
#' Given a PAC_filtsep output data.frame, where each column contains the
#' sequences that passed the filter for a specific group specified in
#' pheno_target, filtsep_bin converts this into a data.frame where sequences are
#' reported as hit (=1) or no hit (=0). Such binary coded group occurance can
#' for example be used by UpSetR::upset to generate visualization of overlaps
#' using UpSet plots.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param filtsep_out PAC_filtsep output data.frame, where each column contains the
#' names of sequences that passed the filter for a specific group specified by a
#' pheno_target object.
#'   
#' @return data.frame where each uniques sequence (row names) in filtsep_out are
#'   reported as hit (=1) or no hit (=0) across samples or sample groups (column
#'   names).
#'   
#' @examples
#' 
#' 

filtsep_bin <- function(filtsep_out){
     filtsep_lst  <- as.list(filtsep_out)
     seqs <- unique(unlist(filtsep_lst, use.names = FALSE))
     seqs <- seqs[!is.na(seqs)]
     df <- as.data.frame(matrix(NA, nrow=length(seqs), ncol=length(filtsep_lst)))
     colnames(df) <- names(filtsep_lst)
     rownames(df) <- seqs
     rm(seqs)
     for(i in 1:length(filtsep_lst)){
         df[,i]  <- as.numeric(rownames(df) %in% filtsep_lst[[i]])
     }
     return(df)
}