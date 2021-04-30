#' Generates normalized values from a PAC object
#'
#' \code{PAC_norm} generates normalized values from a PAC object
#' 
#' Using the counts in a PAC object to generate normalized values in a dataframe
#' with the same rownames as in the original PAC object
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a Count table with raw counts.
#'
#' @param norm Character indicating what type of normalization method that
#'   should be applied to the PAC$Counts. If norm="cpm", counts per million
#'   reads) is returned. Each sequence is then divided against the total number
#'   of counts in a given sample. If norm="vst", PAC$Counts will be imported
#'   into the \code{\link[DESeq2]{varianceStabilizingTransformation}} function
#'   in the DESeq2-package with the options blind=TRUE and fitType="mean". If
#'   norm="rlog", PAC$Counts will instead be imported into the
#'   \code{\link[DESeq2]{rlog}} function also available in
#'   DESeq2-package (options blind=TRUE and fitType="mean") for a log2
#'   transformed version of vst, and is more robust to varying library sizes.
#'   
#' @param merge_pac logical whether the normalized table should be returned and
#'   stored in the PAC$norm 'folder' of the provided PAC object (TRUE) or be
#'   returned as a data frame.
#'   
#' @return A normalized count table, or if pac_merge=TRUE, a PAC object with
#'   normalized counts table added to the norm folder (PAC$norm).
#'   
#' @examples
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' pac_norm  <- PAC_norm(pac, norm="cpm") 
#' df_norm <- PAC_norm(pac, norm = "vst", merge_pac = FALSE)
#' 
#' 
#' @export
#' 
PAC_norm <- function(PAC, norm="cpm", merge_pac=TRUE){
  if(norm %in% c("cpm", "rpm")){
    lib_sizes <- colSums(PAC$Counts)
    counts_cpm <- data.frame(
      matrix(NA,
             nrow=nrow(PAC$Counts),
             ncol=ncol(PAC$Counts)))
    colnames(counts_cpm) <- colnames(PAC$Counts)
    rownames(counts_cpm) <- rownames(PAC$Counts)
    for (i in 1:length(lib_sizes)){
      counts_cpm[,i] <- (PAC$Counts[,i]/(lib_sizes[i]/1000000))
    }
    fin <- counts_cpm
    PAC$norm$cpm <- fin
    
  }
  if(norm=="vst") { 
    
    fin <- DESeq2::varianceStabilizingTransformation(as.matrix(PAC$Counts), 
                                                     blind=TRUE, 
                                                     fitType="parametric")
    PAC$norm$vst <- fin 
  }
  if(norm=="rlog") { 
    fin <- DESeq2::rlogTransformation(as.matrix(PAC$Counts), 
                                      blind=TRUE, fitType="parametric")
    row.names(fin) <- row.names(PAC$Counts)
    PAC$norm$rlog <- fin 
    
  }
  if(merge_pac==TRUE){return(PAC)} else {return(fin)} 
}
