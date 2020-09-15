#' piRNA analysis of PAC object
#'
#' Analyzing and plotting piRNAs between groups.
#'
#' Piwi interacting RNA are often classified according to what type of genomic
#' element they silence. Sileincing targets may differ from cell-type to
#' cell-type, which makes piRNA analysis challenging. Given a PAC object, this
#' function automates the mapping and annotation of piRNA.
#' 
#' @family PAC analysis
#'   
#' @seealso \url{https://github.com/Danis102} for updates.
#'
#' @param PAC PAC-list object. The Anno object needs to contain a piRNA column
#'   generated during the \family{PAC reannotation} workflow
#'   (\code{map_reanno},\code{import_reanno},\code{make_reanno},
#'   \code{add_reanno}.
#'
#' @param data_pirna
#' 
#' @param data_custom
#' 
#' @param graphs Logical whether graphs should be reported or not.
#'   
#' @return List of ggplot2 plots and the data used for generating the plots. Use
#'   ls.str() to explore each level.
#'   
#' @examples
#' 
#' library(seqpac)
#' 
#' ### Prepare pirbase data file  ###
#' # - Download pirbase data from pirbase. Make sure you get a table where dataset identity is specified in a seperate column.
#' # - You may also use your own file, with custom   
#' # - Resulting data.table should have 4 columns named: "name", "seq", "dataset", "n_evidence"
#' #   
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Humans/pirBase/hg38/piR_hsa.txt.gz"), delim="\t", col_names = TRUE)
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.txt.gz"), delim="\t", col_names = TRUE)
#' 
#' datset <- strsplit(pirbase_dat$dataset, " ")
#' datset_var <- lapply(datset, function(x){return(unlist(lapply(strsplit(x, ":"), function(y){y[1]})))}) # Takes a few minutes
#' n_evidence  <- unlist(lapply(datset_var, length))
#' dats  <- unlist(lapply(datset_var, function(x){paste0(x, collapse=";")}))
#' 
#' dat_pirna <- data.table::data.table(name=pirbase_dat$name, seq=pirbase_dat$sequence, dataset=dats, n_evidence=n_evidence)
#' dat_pirna <- dat_pirna[!is.na(ref_pirna$seq),]
#' 
#' dat_pirna <- dat_pirna
#' 
#' ### Prepare repeatMasker data file ###
#' # 
#' 
#' rm_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Humans/pirBase/hg38/piR_hsa.txt.gz"), delim="\t", col_names = TRUE)
#' 
#' 
#' 
#' ### Prepare repeatMasker data file ###
#' # 
#' 
#' load("/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/PAC_merge_10counts100.Rdata")
#' PAC <- PAC_merge_10counts100
#' slct <- as.character(unique(PAC$Pheno$Groups2))
#' PAC_sub <- PAC_filter(PAC, pheno_target=list("Groups2", slct[c(5,6,1,3)]))
#' v_list <- PAC_filtsep(PAC_sub, type = "rpm", threshold = 10, coverage = 100, pheno_target=list("Groups2"))
#' PAC_sub <- PAC_filter(PAC_sub, anno_target= unique(do.call("c", as.list(v_list))))
#' 
#' PAC <- PAC_sub
#' 
#' 
#' 
#' 
#' @export


PAC_pirna <- function(PAC, ref_pirna=NULL, ref_custom=NULL, )

                      
                      


pirbase_dat


test3

185230  175900