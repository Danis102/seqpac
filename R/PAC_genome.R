#' Map PAC to genome
#'
#' Map PAC object sequences to genome coordinates using primarily using
#' functions in the Rsubreads and Biostrings packages.
#' 
#' Given a PAC object and an indexed genome this function will add genome
#' alignments to the PAC$Anno object or save a genome alignments as a dataframe
#' ordered as the PAC.
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
#' @param genome
#' 
#' @param data_custom
#' @param data_custom
#' @param graphs Logical whether graphs should be reported or not.
#'   
#' @return List of ggplot2 plots and the data used for generating the plots. Use
#'   ls.str() to explore each level.
#'   
#' @examples
#' 
#' library(seqpac)
#' 
#' 
#' ########################################
#' ##### Index your genome fasta file #####
#' ########################################
#' ## Download your genome from Ensembl, UCSC, NCBI etc.
#' ## In this example we use a reference genome from Ensembl.
#' ## Ensembl always provide gtf files for easy annotation.
#' ## For aligning sequences to a genome we use the align function in the Rsubread package. 
#' ## Genome alignments with Rsubread always need an genome index.
#' ## The index is generated once using the buildindex function.
#' ## Important, setting TH_subread=<very high> assures that all repeats are mapped.
#' 
#' ## We perfer to make seperate alignments against full chromosomes and scaffold/haplotypes
#' ## Therefore we start by splitting the ref seq by chromosome and then save them as two seperate fasta files
#' 
#' fast_orig <- Biostrings::readDNAStringSet("/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz", format="fasta", use.names=TRUE)
#' fast_chr <- fast_orig[grepl("dna:chromosome", names(fast_orig)),] 
#' fast_other <- fast_orig[!grepl("dna:chromosome", names(fast_orig))]
#' Biostrings::writeXStringSet(fast_chr, filepath="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa", format="fasta")
#' Biostrings::writeXStringSet(fast_other, filepath="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/scaffolds/fast_other.fa", format="fasta")
#' 
#' ## Now we make indexes for Rsubread
#' ## We set TH_subread to a high number to include as many repeats as possible.
#' ## Remember however that this may dramatically affect performance for some species (take a long time).
#'
#'  
#' Rsubread::buildindex(basename="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr", 
#'            reference="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa", 
#'            gappedIndex=FALSE, memory= 32000, TH_subread=100000000)
#'   
#' 
#' ########################################
#' ##### Align PAC to genome          #####
#' ########################################
#' ## PAC_genome is a wrapper for the align and subjunc 
#' ## functions in the Rsubread package, designed for easy 
#' ## integration into seqpac.
#' 
#' 
#' 
#' 
#' 
#' #@export
library(seqpac)
load("/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/PAC_merge_10counts100.Rdata")
PAC <- PAC_merge_10counts100

slct <- as.character(unique(PAC$Pheno$Groups2))
PAC_sub <- PAC_filter(PAC, pheno_target=list("Groups2", slct[c(5,6,1,3)]))
v_list <- PAC_filtsep(PAC_sub, type = "rpm", threshold = 10, coverage = 100, pheno_target=list("Groups2"))
PAC_sub <- PAC_filter(PAC_sub, anno_target= unique(do.call("c", as.list(v_list))))

PAC_merge=FALSE
maxMismatches=3
maxIndels=1
nthreads=8
genome="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa"




PAC_genome <- function(x){PAC, genome, PAC_merge=FALSE, maxMismatches=3, maxIndels=0, nthreads=1)  
                            
                            tmp_nam_in  <- tempfile(pattern = "", fileext = ".fa")
                            tmp_nam_out  <- tempfile(pattern = "", fileext = ".sam")
                            
                            fst <- Biostrings::DNAStringSet(rownames(PAC$Anno))
                            names(fst) <- seq_along(fst)
                            Biostrings::writeXStringSet(fst, filepath=tmp_nam_in, format="fasta")
                            
                            indx_path <- paste0(dirname(genome),"/", gsub(".fa$|.fasta$", "", basename(genome)))  
                           
                            sam_lst <- as.list(rep(NA, times=maxMismatches+1))
                            
                            ## Find all perfect hits
                            Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
                                                      unique=FALSE, nBestLocations=16, indels=0, complexIndels=FALSE,
                                                      input_format="fasta", maxMismatches=0, nthreads=nthreads)
                            sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=100000))))))
                            if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
                            sub <- fst[names(fst) %in% sam[sam$X2 == "4",]$X1,]
                            Biostrings::writeXStringSet(sub, filepath=tmp_nam_in, format="fasta")
                            sam_lst[[1]] <- sam[!sam$X2 == "4",]
                            names(sam_lst)[1] <- paste0("mis0_ind0")
                            
                            
                            sq_thr_mis <- 0:maxMismatches
                            sq_thr_indel <- 0:maxIndels
                            
                            lapply(as.list(sq_thr_mis), function(x){sq_thr_indel+x})
                              
                            out_nam_mis  <- tempfile(pattern = "mis", fileext = ".sam")
                            out_nam_ind  <- tempfile(pattern = "ind", fileext = ".sam")
                                     
                            
                            
                                            
                            for(mis in sq_thr_mis){
                                 for(ind in sq_thr_indel){
                                      
                              
                                      Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=out_nam, output_format="SAM",
                                                      unique=FALSE, nBestLocations=16, indels=ind, complexIndels=FALSE,
                                                      input_format="fasta", maxMismatches=mis, nthreads=nthreads)
                                      }
                                      
                                      sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=10000))))))
                                      
                                      if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
                                      sub <- fst[names(fst) %in% sam[sam$X2 == "4",]$X1,]
                                      Biostrings::writeXStringSet(sub, filepath=tmp_nam_in, format="fasta")
                                      sam_lst[[paste0("mis_", i)]] <- sam[!sam$X2 == "4",]
                                 
                                      }
                            
                                        
                                        
                            data.table::setDT(sam_0)
                            
                            
                            Rsubread::align(index=genome_index, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
                                            unique=TRUE, nBestLocations=16, indels=indels, complexIndels=FALSE,
                                            input_format="fasta", maxMismatches=maxMismatches, nthreads=nthreads)
                            
                            Rsubread::align(index=genome_index, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
                                            unique=TRUE, nBestLocations=16, indels=indels, complexIndels=FALSE,
                                            input_format="fasta", maxMismatches=maxMismatches, nthreads=nthreads)
                              
                            
                            sam_0 <- data.table::setDT(sam_0)
                            test <- split(sam_0, factor(sam_0$X1))
                            test  <- lapply(split(sam, factor(sam$X1)), function(x){
                                                                          cigr <- unique(x$X6)
                                                                          stopifnot(length(unique(x$X6)) == 1)  
                                                                          data.table::data.table
                                                                          <- paste(x, collapse=

                            