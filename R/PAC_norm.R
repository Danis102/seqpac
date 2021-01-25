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
#' @param norm Character describing what type of normalization method to be
#'   applied to the PAC$Counts. If norm="cpm" will return reads/counts per
#'   million (or sometimes refered to as counts per million reads). Each
#'   sequence is then divided against the total number of counts in a given
#'   sample. If norm="vst" PAC$Counts will instead be imported into the
#'   \code{varianceStabilizingTransformation} function in the DESeq2-package
#'   with the options blind=TRUE and fitType="mean". If norm="rlog" PAC$Counts
#'   will instead be imported into the \code{rlogTransformation} function also
#'   available in DESeq2-package (options blind=TRUE and fitType="mean") for a
#'   log2 transformed version of vst that are more robust to varying library
#'   sizes.
#'   
#' @param PAC_merge logical whether the normalized table should be returned and
#'   stored in the PAC$norm 'folder' of the provided PAC object (TRUE) or be
#'   returned as a data frame.
#'   
#' @return A normalized count table 
#'
#' @examples
#'   df  <- PAC_norm(PAC, norm="cpm") 
#' @export
#' 
PAC_norm <- function(PAC, norm="cpm", PAC_merge=TRUE){
                  if(norm %in% c("cpm", "rpm")){
                          lib_sizes <- colSums(PAC$Counts)
                          counts_cpm <- data.frame(matrix(NA, nrow=nrow(PAC$Counts), ncol=ncol(PAC$Counts)))
                          colnames(counts_cpm) <- colnames(PAC$Counts)
                          rownames(counts_cpm) <- rownames(PAC$Counts)
                          for (i in 1:length(lib_sizes)){
                                            counts_cpm[,i] <- (PAC$Counts[,i]/(lib_sizes[i]/1000000))
                          }
                          fin <- counts_cpm
                          PAC$norm$cpm <- fin
                          
                          }
                  if(norm=="vst") { 

                          fin <- DESeq2::varianceStabilizingTransformation(as.matrix(PAC$Counts), blind=TRUE, fitType="parametric")
                          PAC$norm$vst <- fin 
                  }
                  if(norm=="rlog") { 
                          fin <- DESeq2::rlogTransformation(as.matrix(PAC$Counts), blind=TRUE, fitType="parametric")
                          PAC$norm$rlog <- fin 
                  }
          if(PAC_merge==TRUE){return(PAC)} else {return(fin)} 
           }
