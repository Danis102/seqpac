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
#' @param filtsep_out PAC_filtsep output data.frame, where each column contains
#'   the names of sequences that passed the filter for a specific group
#'   specified by a pheno_target object.
#'
#' @return data.frame where each uniques sequence (row names) in filtsep_out are
#'   reported as hit (=1) or no hit (=0) across samples or sample groups (column
#'   names).
#'   
#' @examples
#' 
#'  \donttest{
#' 
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' 
#' ## Keep sequences with 5 counts (threshold) in 100% (coverage) of 
#' ## samples in a group:
#'  # Use PAC_filtsep to find sequences 
#'  filtsep <- PAC_filtsep(pac_master, norm="counts", threshold=5, 
#'                         coverage=100, pheno_target= list("stage"))
#'  
#'  # Venn-diagram of the overlap
#'  olap <- reshape2::melt(filtsep, 
#'                         measure.vars = c("Stage1", "Stage3", "Stage5"), 
#'                         na.rm=TRUE)
#'  plot(venneuler::venneuler(data.frame(olap[,2], olap[,1]))) 
#'  
#'  # Filter by unique sequences passing filtsep  
#'  filtsep <- unique(do.call("c", as.list(filtsep)))
#'  pac_filt <- PAC_filter(pac_master, subset_only = TRUE, anno_target= filtsep)
#'  
#' ## Upset plot using the UpSetR package
#'  # (when output="binary" PAC_filtsep uses filtsep_bin for binary conversion
#'  # Use PAC_filtsep with binary output
#'  filtsep_bin <- PAC_filtsep(pac_master, norm="counts", threshold=5, 
#'                             coverage=100, pheno_target= list("stage"), 
#'                             output="binary")
#'  
#'  # Plot with UpSetR 
#'  UpSetR::upset(filtsep_bin, sets = colnames(filtsep_bin), 
#'                mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order=TRUE)
#'                
#' }            
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