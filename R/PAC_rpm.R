#' Generates normalized values from a PAC object (deprecated)
#'
#' \code{PAC_rpm} Generates RPM values from a PAC object (deprecated)
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
#' @return A normalized count table (
#'
#' @examples
#'   df  <- PAC_norm(PAC, type="rpm") 
#' @export
#' 
PAC_rpm<- function(PAC){
                  stop("PAC_rpm has been depricated.\nPlease use PAC_norm instead.")  
                  }
