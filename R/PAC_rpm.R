#' Generates RPM values from a PAC object
#'
#' \code{generate_rpm} Generates RPM values from a PAC object
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
#' @return dataframe
#'
#' @examples
#'   df  <- PAC_rpm(PAC) 
#' @export
#' 
PAC_rpm <- function(PAC){
                                                      lib_sizes <- colSums(PAC$Counts)
                                                      counts_RPM <- data.frame(matrix(NA, nrow=nrow(PAC$Counts), ncol=ncol(PAC$Counts)))
                                                      colnames(counts_RPM) <- colnames(PAC$Counts)
                                                      rownames(counts_RPM) <- rownames(PAC$Counts)
                                                      for (i in 1:length(lib_sizes)){
                                                                        counts_RPM[,i] <- (PAC$Counts[,i]/(lib_sizes[i]/1000000))
                                                      }
                                                      PAC$norm$rpm <- counts_RPM
                                                      return(PAC)
                                }