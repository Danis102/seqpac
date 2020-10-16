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
#' @param mismatches 
#' 
#' @param pirna
#' 
#' @param genome
#' 
#' @param fasta_list
#' 
#' @param gtf_list
#' 
#' @param nuc_bias Integer vector indicating which nucleotide to perform
#'   nucluotide bias analysis on. As default, nuc_bias = c(1,10). Nucleotide
#'   bias analysis is then performed on the 1st and 10th nucleotide.
#'   
#' @param graphs Logical whether graphs should be reported or not.
#'  
#' @param threads Integer indiciating the number of parallel processes.   
#'     
#' @return List of ggplot2 plots and the data used for generating the plots. Use
#'   ls.str() to explore each level.
#'   
#' @examples
#' 
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' 
#' 
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
# 
# 
# pirna <- "/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.fa"
# genome <- "/home/danis31/Desktop/Temp_docs/fasta/biomartr_genome/chromosomes.fa"
# 
# 
# 
# PAC_pirna <- function(PAC, mismatches=3, pirna="piRNA", genome="mis0_genome",
#                       fasta_list=NULL, gtf_list=NULL, nuc_bias=c(1,10), cluster=TRUE, threads=1, ...){
# 
#   ## Setup
#   seqs <- rownames(PAC$Anno)
#   if(!file.exists(genome)){
#         test <- PAC$Anno[, grepl("_genome$|^genome$", colnames(PAC$Anno))]
# 
# 
#   ## Generate main pirna annotations from fasta or Anno
#   if(file.exists(pirna)){
#     cat("Input pirna was an existing file. Will treat it as a \nfasta reference and make a denovo reannotation using bowtie. \nSee ?map_reanno or ?vingette for details.\n")
#     outpath <- tempfile(pattern = "", fileext = "")
#     err <- try(map_reanno(PAC, ref_paths=list(pirna=pirna), output_path=outpath, type="external", mismatches=3,
#                    import = list(coord=FALSE, report="minimum", reduce=NULL), threads=threads, ...), silent = TRUE)
#     if(!is.null(err)){
#       err2 <- try(map_reanno(PAC, ref_paths=list(pirna=pirna), output_path=outpath, type="internal", mismatches=3,
#                    import = list(coord=FALSE, report="minimum", reduce=NULL), threads=threads, ...), silent = TRUE)
#     if(!is.null(err2)){
#       stop(paste0("\nFunction map_reanno failed. Possible reasons: \n\tNo bowtie installation\n\tBad fasta reference\n\nLast log says:\n", err2))
#       }
#     }
#     reanno <- make_reanno(outpath, PAC=PAC, mis_fasta_check=TRUE, threads=threads)
#     pirna_anno <- add_reanno(reanno, bio_search=list(pirna="pirna"), type="biotype", bio_perfect=FALSE, mismatches = 3)
#     rm(reanno)
#     unlink(outpath, recursive = TRUE)
#   } else {
#     pirna <- unlist(strsplit(pirna, "\\$"))
#     pirna <- PAC$Anno[,pirna[length(pirna)]]
#     cat("Using user defined column in PAC$Anno as main pirna classification.\n")
#   }
# 
#  ## Generate main genome annotations from fasta or Anno
#   if(file.exists(genome)){
#     cat("\nInput genome was an existing file. Will treat it as a \nfasta reference and make a denovo reannotation using bowtie. \nSee ?map_reanno or ?vingette for details.\n")
#     outpath <- tempfile(pattern = "", fileext = ".out")
#     err <- try(map_reanno(PAC, ref_paths=list(genome=genome), output_path=outpath, type="external", mismatches=3,
#                    import = list(coord=TRUE, report="full", reduce=NULL), threads=threads), silent = TRUE)
#     if(!is.null(err)){
#       outpath <- tempfile(pattern = "", fileext = "")
#       err2 <- try(map_reanno(PAC, ref_paths=list(genome=genome), output_path=outpath, type="internal", mismatches=3,
#                    import=list(coord=TRUE, report="full", reduce=NULL), threads=threads), silent = TRUE)
#     if(!is.null(err2)){
#       stop(paste0("\nFunction map_reanno failed. Possible reasons: \n\tNo bowtie installation\n\tBad fasta reference\n\tNo bowtie index (see ?Rbowtie::bowtie_build)\n\nLast log says:\n", err2))
#       }
#     }
#     reanno <- make_reanno(outpath, PAC=PAC, mis_fasta_check=TRUE, threads=threads)
#     anno_genome <- add_reanno(reanno, type="genome", mismatches = 3, genome_max="all")
# 
# 
#     rm(reanno)
#     unlink(outpath, recursive = TRUE)
#   }else{
# 
# 
# 
# 
#     genome <- unlist(strsplit(genome, "\\$"))
#     genome <- PAC$Anno[,genome[length(genome)]]
#     cat("Using user defined column in PAC$Anno as main pirna classification.")



    
