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
#' fast_orig <- Biostrings::readDNAStringSet(
#'              "<your path to genome fasta>", format="fasta", use.names=TRUE)
#' fast_chr <- fast_orig[grepl("dna:chromosome", names(fast_orig)),] 
#' fast_other <- fast_orig[!grepl("dna:chromosome", names(fast_orig))]
#' Biostrings::writeXStringSet(
#'     fast_chr, filepath="<your path to new chromosome fasta>", format="fasta")
#' Biostrings::writeXStringSet(
#'     fast_other, filepath="<your path to new other fasta>", format="fasta")
#'     
#' ## Now we make indexes for Rsubread
#' ## We set TH_subread to a high number to include as many repeats as possible.
#' ## Remember however that this may dramatically affect performance for some 
#' ## species (take a long time).
#'
#'  
#' Rsubread::buildindex(
#'    basename="<example: dm6_ensembl_release_101/fasta/chr/fast_chr>", 
#'    reference="<example:dm6_ensembl_release_101/fasta/chr/fast_chr.fa>", 
#'    gappedIndex=FALSE, memory= 32000, TH_subread=10000000000)
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
# library(seqpac)
# load("/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/PAC_merge_10counts100.Rdata")
# PAC <- PAC_merge_10counts100
# 
# slct <- as.character(unique(PAC$Pheno$Groups2))
# PAC_sub <- PAC_filter(PAC, pheno_target=list("Groups2", slct[c(5,6,1,3)]))
# v_list <- PAC_filtsep(PAC_sub, type = "rpm", threshold = 10, coverage = 100, pheno_target=list("Groups2"))
# PAC_sub <- PAC_filter(PAC_sub, anno_target= unique(do.call("c", as.list(v_list))))
# 
# PAC_merge=FALSE
# maxMismatches=3
# maxIndels=3
# nthreads=8
# genome="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa"
# 
# 
# 
# 
# PAC_genome <- function(x){PAC, genome, PAC_merge=FALSE, maxMismatches=3, maxIndels=3, nthreads=1)  
#                             
#                             ## Setup ##
#                             
# 
#                             cat("Please note that this mapping function was developed for sRNA without junctions such as exon/exon boundaries.") 
#                             
#                             indx_path <- paste0(dirname(genome),"/", gsub(".fa$|.fasta$", "", basename(genome)))      
#                             tmp_nam_in  <- tempfile(pattern = "", fileext = ".fa")
#                             tmp_nam_out  <- tempfile(pattern = "", fileext = ".sam")
#                             
#                             fst <- Biostrings::DNAStringSet(rownames(PAC$Anno))
#                             if(max(nchar(paste(fst))) > 100 ){warning("Noticed that you sequences were longer than 100 nt.\nPlease use the align or Subjunc functions in the Rsubread package\nin case your intention was to map long RNA junctions (e.g. exon/exon boundaries).")}   
#                             names(fst) <- seq_along(fst)
#                             Biostrings::writeXStringSet(fst, filepath=tmp_nam_in, format="fasta")
# 
#                             sam_lst <- as.list(rep(NA, times=maxMismatches+1))
#                             
#                             ## Find all perfect hits ## 95,595
#                             # Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                             #                           unique=FALSE, nBestLocations=16, indels=0, complexIndels=FALSE,
#                             #                           input_format="fasta", maxMismatches=0, nthreads=nthreads, TH1=0, TH2=0, DP_GapExtPenalty=-1, DP_MismatchPenalty=-1)
#                             # sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=100000))))))
#                             # sam_perfect <- sam[!grepl("S|D|I|H|X|\\=|\\*", sam$X6),]  #Sort out starting/ending substitutions from sam (all should only contain M cigar)
#                             # sam_left <- sam[grepl("S|D|I|H|X|\\=|\\*", sam$X6),]
#                             # if(any(!gsub("\\d", "", unique(sam_perfect$X6)) == "M")){stop("Perfectly matched sam-file contains unidentified cigar-string entery.")}    
#                             # if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                             #   
#                             # seq_left <- fst[names(fst) %in% sam_left$X1,]
#                             # stopifnot(length(unique(sam_left$X1)) %in% length(names(seq_left)))
#                             # Biostrings::writeXStringSet(seq_left, filepath=tmp_nam_in, format="fasta")
#                             # 
#                             # sam_lst[[1]] <- sam_perfect
#                             # names(sam_lst)[1] <- paste0("div0")
#                             
#                             
#                             ## Generate mis alignment series
#                             mis_lst <- list(NULL)
#                             for(i in 0:maxMismatches){
#                                       Rbowtie::bowtie(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                                       unique=FALSE, nBestLocations=16, indels=0, complexIndels=FALSE, nsubreads=20,
#                                                       input_format="fasta", maxMismatches=i, nthreads=nthreads, minFragLength = min(fst@ranges@width), TH1=2)
#                                       
#                                   
#                                       sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=10000))))))
#                                       #sam[sam$X1 == 639,]
#                                       #sam[sam$X1 == 1611,]
#                                       sam_perfect <- sam[!grepl("S|D|I|H|X|\\=|\\*", sam$X6),]  #Sort out starting/ending substitutions from sam (all should only contain M cigar)
#                                       sam_left <- sam[grepl("S|D|I|H|X|\\=|\\*", sam$X6),]
#                                       
#                                       if(any(!gsub("\\d", "", unique(sam_perfect$X6)) == "M")){stop("Perfectly matched sam-file contains unidentified cigar-string entery.")}    
#                                       if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                                       
#                                       seq_left <- fst[names(fst) %in% sam_left$X1,]
#                                       if(!is.null(length(seq_left))){
#                                           stopifnot(length(unique(sam_left$X1)) %in% length(names(seq_left)))
#                                           Biostrings::writeXStringSet(seq_left, filepath=tmp_nam_in, format="fasta")
#                                           
#                                       }
#                                       mis_lst[[i+1]] <- sam_perfect
#                                       names(mis_lst)[i+1] <- paste0("m", i, "i0")
#                                       }
#                             
#                             
#                             
#                             
#                          
#                             ## Generate mis alignment series
#                             mis_lst <- list(NULL)
#                             for(i in 0:maxMismatches){
#                                       Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                                       unique=FALSE, nBestLocations=16, indels=0, complexIndels=FALSE, nsubreads=20,
#                                                       input_format="fasta", maxMismatches=i, nthreads=nthreads, minFragLength = min(fst@ranges@width), TH1=2)
#                                       
#                                   
#                                       sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=10000))))))
#                                       #sam[sam$X1 == 639,]
#                                       #sam[sam$X1 == 1611,]
#                                       sam_perfect <- sam[!grepl("S|D|I|H|X|\\=|\\*", sam$X6),]  #Sort out starting/ending substitutions from sam (all should only contain M cigar)
#                                       sam_left <- sam[grepl("S|D|I|H|X|\\=|\\*", sam$X6),]
#                                       
#                                       if(any(!gsub("\\d", "", unique(sam_perfect$X6)) == "M")){stop("Perfectly matched sam-file contains unidentified cigar-string entery.")}    
#                                       if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                                       
#                                       seq_left <- fst[names(fst) %in% sam_left$X1,]
#                                       if(!is.null(length(seq_left))){
#                                           stopifnot(length(unique(sam_left$X1)) %in% length(names(seq_left)))
#                                           Biostrings::writeXStringSet(seq_left, filepath=tmp_nam_in, format="fasta")
#                                           
#                                       }
#                                       mis_lst[[i+1]] <- sam_perfect
#                                       names(mis_lst)[i+1] <- paste0("m", i, "i0")
#                                       }
#                             
#                             ## Terminal substitutions alignment series
#                             Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                             unique=FALSE, nBestLocations=16, indels=0, complexIndels=FALSE,
#                                             input_format="fasta", maxMismatches=i, nthreads=nthreads, minFragLength = min(fst@ranges@width), TH1=2)
#                             sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=10000))))))
#                             if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                             cig_splt <- strsplit(sam$X6, "(?<=[[:alpha:]])", perl = TRUE)
#                             log_num <- unlist(lapply(cig_splt, function(x){if(length(x)==2|length(x)==3){sum(as.numeric(gsub("S", "", x[grepl("S", x)])))}else{0}}))
#                             
#                             termsub_lst <- list(NULL)
#                             for(i  in 1:maxMismatches){
#                                         termsub_lst[[i]] <- sam[log_num == i,]
#                                         names(termsub_lst)[i] <- paste0("term_S",i)
#                             }
#                             sam_left <- sam[!log_num %in% 1:maxMismatches,]
#                             seq_left <- fst[names(fst) %in% sam_left$X1,]
#                             if(!is.null(length(seq_left))){
#                                           stopifnot(length(unique(sam_left$X1)) %in% length(names(seq_left)))
#                                           Biostrings::writeXStringSet(seq_left, filepath=tmp_nam_in, format="fasta")
#                                       }
#                            
#                             ## Indel alignment series
#                             indel_lst <- list(NULL)
#                             for(i in 0:){
#                                       Rsubread::align(index=indx_path, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                                       unique=FALSE, nBestLocations=16, indels=1, complexIndels=FALSE,
#                                                       input_format="fasta", maxMismatches=0, nthreads=nthreads, minFragLength = 20, TH1=2)
#                               
#                                       sam <- suppressWarnings(suppressMessages(readr::read_delim(tmp_nam_out, delim="\t", col_names = FALSE, skip = sum(grepl("^@", readLines(tmp_nam_out, n=10000))))))
#                                       as.data.frame(sam[!sam$X6 == "*",])
#                                       
#                             fast_orig <- Biostrings::readDNAStringSet(genome, format="fasta", use.names=TRUE)
# 
#                                       
#                                        test <- lapply(as.list(sam$X10), function(x){grepl(x, paste(fast_orig))})
#                                        test <- 
#                                          
#                                          lapply(as.list(sam$X10), function(x){grepl(x, paste(fast_orig))})
#                                        
#                                       sam_perfect <- 
#                                         
#                                        test <- sam[!grepl("\\*", sam$X6),]  #Sort out starting/ending substitutions from sam (all should only contain M cigar)
#                                       sam_left <- sam[grepl("S|D|I|H|X|\\=|\\*", sam$X6),]
#                                       
#                                       if(any(!gsub("\\d", "", unique(sam_perfect$X6)) == "M")){stop("Perfectly matched sam-file contains unidentified cigar-string entery.")}    
#                                       if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                                       
#                                       seq_left <- fst[names(fst) %in% sam_left$X1,]
#                                       if(!is.null(length(seq_left))){
#                                           stopifnot(length(unique(sam_left$X1)) %in% length(names(seq_left)))
#                                           Biostrings::writeXStringSet(seq_left, filepath=tmp_nam_in, format="fasta")
#                                           
#                                       }
#                                       mis_lst[[i+1]] <- sam_perfect
#                                       names(mis_lst)[i+1] <- paste0("m", i, "i0")
#                                       }
#                             
#                                                 
#                             log_3 <- unlist(lapply(cig_splt, function(x){if(length(x)==2){gsub("S", "", x[grepl("S", x)])}else{0}}))  
#                                       
#                             
# 
# 
# 
#                             ## Generate mis and indel combinational hierarchy
#                             tab_vect <- expand.grid(0:maxMismatches, 0:maxIndels)
#                             sms <- rowSums(tab_vect) 
#                             
#                             tab_vect <- tab_vect[sms <= max(c(maxMismatches, maxIndels)),] # Remove instances where sums of indels and mis becomes more than max allowed 
#                             sms_upd <- rowSums(tab_vect) 
#                             hier_comb_lst <- split(tab_vect, sms_upd)
#                             names(hier_comb_lst) <- paste0("div", names(hier_comb_lst))
#                             if(any(names(hier_comb_lst) == "div0")){ hier_comb_lst <- hier_comb_lst[!names(hier_comb_lst) == "div0"]}
#                             
#                             
#                             
#                             
#                                       unique(sam_mis$X1[
#                                       div_lst[[i]] 
#                                       names(div_lst)[i] <- paste0("m",i, "i0")
#                                       
#                                       
#                             
#                             sam_lst[[1]] <- sam_perfect
#                             names(sam_lst)[1] <- paste0("div0")
#                                       
#                                   }
#                                   
#                                   
#                                   lapply(div_lst, function(t){table(!grepl("S|D|I|H|X|\\=|\\*", t$X6))})
#                                   
#                                   lapply(div_lst, function(x){x[!grepl("S|D|I|H|X|\\=|\\*", x$X6),]})
#                                   
#                                   
#                                   temp <- div_lst$m1i2[!grepl("S|D|I|H|X|\\=|\\*", div_lst$m1i2$X6),]
#                                   test2 <- lapply(div_lst, function(x){test <- x[!grepl("S|D|I|H|X|\\=|\\*", x$X6),]
#                                                               return(test[!test$X1 %in% temp$X1,])})    
#                                   
# 
#                                   lapply(div_lst, function(x){x[x$X1==63,]})
#                                   
#                                   table(div_lst$m1i2$X6)
#                                   
#                                   div_lst$m1i2[grepl("2S3M1D24M", div_lst$m1i2$X6),]
#                                   
#                                   
#                                   lapply(div_lst, function(x){x[x$X1==43075,]})
#                                   
#                                   div_lst$m3i0 %in% m2i1$m2i1
#                                   
#                                   
#                                   table(names(seq_left) == 63)
#                                   table(sam_perfect$X1 == 63)
#                                   sam_perfect
#                                   
#                                    <- strsplit(x$X6, "(?<=[[:alpha:]])", perl = TRUE)
#                                   strsplit("( x$X6)
#                                   X6
#                                   
#                                     tTTGAACGCATATCGCAGTCCATGCTG
#                                   
#                                   GGTTGAACGCATATCGCAGTCCATGCTG
#                                   
#                                   GGTTGAACGCATATCGCAGTCCATGCTG
#                                   SSD
#                                   tTGAA...
#                                   
#                                   !!! Ah, allowing deletions complicates things
#                                   !!! Check example below:
#                                   
#                                   div_lst[[1]][div_lst[[1]]$X6 == "2S2M1D24M",]
#                                   div_lst[[2]][div_lst[[2]]$X6 == "2S2M1D24M",]
#                                   div_lst[[1]][div_lst[[1]]$X1 == "54634",]
#                                   div_lst[[2]][div_lst[[2]]$X1 == "54634",]
#                                   
#                                   
#                                   
#                                   lapply(div_lst, function(t){t$X2 == 4})    
#                                   lapply(div_lst, function(t){table(t$X2)})
#                                   lapply(div_lst, function(t){table(t$X6)}) 
#                                   lapply(div_lst, function(t){length(unique(t$X10))}) 
#                                   
#                                   fst[names(fst) == 54634]
#                                   
#                                   
#                                   
#                                   
#                                   
#                                   
#                                       
#                                       if(any(!grepl("I", sam$X11))){stop(paste0("The columns in Rsubread::align sam output file did not add up.\nCheck the temp output file at:\n", tmp_nam_out))}
#                                       sub <- fst[names(fst) %in% sam[sam$X2 == "4",]$X1,]
#                                       
#                                       
#                                       
#                                       Biostrings::writeXStringSet(sub, filepath=tmp_nam_in, format="fasta")
#                                       
#                                  
#                                       }
#                             
#                                         
#                                         
#                             data.table::setDT(sam_0)
#                             
#                             
#                             Rsubread::align(index=genome_index, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                             unique=TRUE, nBestLocations=16, indels=indels, complexIndels=FALSE,
#                                             input_format="fasta", maxMismatches=maxMismatches, nthreads=nthreads)
#                             
#                             Rsubread::align(index=genome_index, readfile1=tmp_nam_in, output_file=tmp_nam_out, output_format="SAM",
#                                             unique=TRUE, nBestLocations=16, indels=indels, complexIndels=FALSE,
#                                             input_format="fasta", maxMismatches=maxMismatches, nthreads=nthreads)
#                               
#                             
#                             sam_0 <- data.table::setDT(sam_0)
#                             test <- split(sam_0, factor(sam_0$X1))
#                             test  <- lapply(split(sam, factor(sam$X1)), function(x){
#                                                                           cigr <- unique(x$X6)
#                                                                           stopifnot(length(unique(x$X6)) == 1)  
#                                                                           data.table::data.table
#                                                                           <- paste(x, collapse=

                            