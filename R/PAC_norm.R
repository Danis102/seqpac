#' Generates normalized values from a PAC object
#'
#' \code{PAC_norm} generates RPM values from a PAC object
#' 
#' Using the counts in a PAC object to generate RPM values in a dataframe with
#' the same rownames as in the original PAC object
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a Count table with raw counts.
#'
#' @param type Character describing what type of normalization method to be
#'   applied to the PAC$Counts. If type="rpm" will return reads per million (or
#'   sometimes refered to as counts per million reads). Each sequence is then
#'   divided against the total number of sequences in a given sample. If
#'   type="vst" PAC$Counts will instead be imported into the
#'   \code{varianceStabilizingTransformation} function in the DESeq2-package
#'   with the options blind=TRUE and fitType="mean".
#'   If type="rlog" PAC$Counts will instead be imported into the
#'   \code{rlogTransformation} function in the DESeq2-package
#'   with the options blind=TRUE and fitType="mean".
#'
#' @param PAC_merge logical whether the normalized table should be returned as
#'   stored in the PAC$norm folder of the provided PAC object (TRUE) or be
#'   retruned as a data frame.
#'   
#' @return A normalized count table 
#'
#' @examples
#'   df  <- PAC_norm(PAC, type="rpm") 
#' @export
#' 
PAC_norm <- function(PAC, type="rpm", PAC_merge=TRUE){
                  if(type=="rpm"){
                          lib_sizes <- colSums(PAC$Counts)
                          counts_RPM <- data.frame(matrix(NA, nrow=nrow(PAC$Counts), ncol=ncol(PAC$Counts)))
                          colnames(counts_RPM) <- colnames(PAC$Counts)
                          rownames(counts_RPM) <- rownames(PAC$Counts)
                          for (i in 1:length(lib_sizes)){
                                            counts_RPM[,i] <- (PAC$Counts[,i]/(lib_sizes[i]/1000000))
                          }
                          fin <- counts_RPM
                          PAC$norm$rpm <- fin
                          
                          }
                  if(type=="vst") { 

                          fin <- DESeq2::varianceStabilizingTransformation(as.matrix(pac_master$Counts), blind=TRUE, fitType="parametric")
                          PAC$norm$vst <- fin 
                  }
                  if(type=="rlog") { 
                          fin <- DESeq2::rlogTransformation(as.matrix(pac_master$Counts), blind=TRUE, fitType="parametric")
                          PAC$norm$rlog <- fin 
                  }
          if(PAC_merge==TRUE){return(PAC)} else {return(fin)} 
           }
