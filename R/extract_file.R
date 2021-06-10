#' Extract character vectors from text files
#'
#' Simple function that uses the tidyverse package to extracts character vectors
#' from any text file, but with special implementation on quickly retrieving
#' sequence names from large fasta files.
#' 
#' Given the path to a file, this function will search each row in that file for
#' a identifier. Rows containing that identifier will be saved, where the
#' resulting strings can be trimmed shorter at yet another identifier. 
#' 
#'
#' @seealso 
#'   \url{https://github.com/Danis102} for updates on the current package.
#'
#' @param file_path Character string defining the file path.
#'
#' @param mode Character string. If mode="fasta", then the function will
#'   automatically search for rows containing fasta sequence names. If
#'   mode="fasta_first", in addition each fasta name will be trimmed at the
#'   first white space. If mode="custom" (default) users must provide \code{search}
#'   and \code{cut} manually.
#'
#' @param search Character string written as a regular expression. This is the
#'   first identifier that identifies which row to be saved.
#'
#' @param cut Character string written as a regular expression. This is the
#'   second identifier that identifies where to trim/cut each saved row string.
#'
#' @param remove_search Logical whether the \code{search} string should be
#'   removed.
#'   
#' @examples
#' 
#' ## Download ss object from GtRNAdb ##
#' dest_path <- file.path("/home/danis31/Desktop/Temp/trna.fa")
#' download.file(url="http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa", destfile=dest_path)
#' file_path <- "/home/danis31/Desktop/Temp/trna.fa"
#' 
#' ## Extracts fasta names upto first white space ##
#' fst_names <- extract_file(file_path, mode="fasta_first")
#' 
#' ## Extracts full fasta names  ##
#' fst_names <- extract_file(file_path, mode="fasta")
#' 
#' ## Extracts custom rows  ##
#' fst_names <- extract_file(file_path, search="tRNA-Val-TAC")
#' 
#' ## Extracts custom rows but with a nicer finish  ##
#' fst_names <- extract_file(file_path, search="tRNA-Val-TAC", cut="Sc:")
#' 
#' @export
extract_file <- function(file_path, mode="custom", search=NULL, cut= NULL, 
                          remove_search=FALSE){
  if(mode=="fasta"){
    search <- "^>"
    remove_search=TRUE
  }
  if(mode=="fasta_first"){
    search <- "^>"
    cut= " "
    remove_search=TRUE
  }
  `%>%` <- dplyr::`%>%`
  text <- readr::read_delim(file_path, col_names = FALSE,  delim="¤¤")
  
  if(is.null(cut)){
    text <- text %>% 
      dplyr::filter(grepl(search, X1)) %>%
         dplyr::pull(X1)
  }
  if(!is.null(cut)){     
    text <- text %>% 
      dplyr::filter(grepl(search, X1)) %>%
        tidyr::separate(col=X1, sep=cut, into=c("X1"), 
                        remove=TRUE, extra="drop") %>%
          dplyr::pull(X1)   
  }
  if(remove_search==TRUE){
        text <-  gsub(pattern = search, replacement = "", text)
  }
  return(text)
}


  