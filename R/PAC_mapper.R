#' Advanced sequences mapping of a PAC object 
#'
#' \code{PAC_mapper} Mapping sequences against a reference.
#'
#' Given a PAC object and the path to a fasta reference file, this function will
#' map sequences in PAC and mapping extract information.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object.
#'
#' @param mismatches=0
#'
#' @param ref_path Character indicating the path to the fasta (.fa) reference
#'   file.
#'
#' @param mismatches Integer indicating the number of mismatches that should be
#'   allowed in the mapping.
#'
#' @param threads Integer indicating the number of parallel processes that
#'   should be used.
#'
#' @param par_type Character indicating if "PSOCK" or "FORK" should be parsed to
#'   parallel::makeCluster (default="PSOCK").
#'
#' @return Stacked list, where each object on the highest level contains:
#'                    (Object 1) Reference name and sequence. 
#'                    (Object 2) Dataframe showing the mapping results of
#'                               each quiary sequence that mapped to Object 1.
#'
#' @examples
#' 
#' ### Load data ###
#' reanno_path="/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/Processed_Pipeline31_05-03-20/R_files/"
#' load(file=paste0(reanno_path, "PAC_filt_rpm10in25.Rdata"))
#' 
#' 
#' 
#' table(PAC_filt$Anno$Biotypes_tRNA)
#' PAC_filt_tRNA <- PAC_filter(PAC_filt, anno_target=list("Biotypes_tRNA", c("Mt_tRNA", "tRNA")), subset_only=TRUE)
#' table(PAC_filt_tRNA$Anno$Biotypes_tRNA)
#' 
#' ref_path <- "/data/Data_analysis/Genomes/Pigs/Sports/Sus_scrofa/tRNA_reanno/tRNA_mature.fa"
#' ref_path <- "/data/Data_analysis/Genomes/Pigs/Sports/Sus_scrofa/tRNA/tRNA_mature.fa"

#' PAC_map <- PAC_mapper(PAC=PAC_filt_tRNA, ref_path, threads=8, mismatches=0)
#' 
#' table(rownames(PAC_filt_tRNA$Anno) %in%  unique(do.call("c", lapply(PAC_map, function(x){ rownames(x[[2]])}))))
#' 
#' save(extended_anno_lst, file=paste0(reanno_path, "/Extended_anno_tsRNA.Rdata"))
#' 
#' 
#' @export

PAC_mapper <- function(PAC, ref_path, mismatches=0, threads=1, par_type="PSOCK"){
                          require(plyr)
                          require(doParallel)
                          require(foreach)
                          require(Biostrings)
                          require(stringr)

                          # Setup parallel computing
                          cl <- parallel::makeCluster(threads, type = par_type)
                          doParallel::registerDoParallel(cl)
                          Anno_frag  <- PAC$Anno
                          query_strings <- as.list(as.character(rownames(Anno_frag)))

                          # Read files
                          cat("First reading full length rRNA from file ...", paste(Sys.time()), "\n\n")
                          full <- Biostrings::readDNAStringSet(ref_path)
                          nams_full <- names(full) # Names are lost in the next step
                          full <- Biostrings::DNAStringSet(paste("NNN", full, "NNN", sep=""))
                          names(full) <- nams_full

                          ## Aligning using parallel processing
                          cat("Now aligning", nrow(Anno_frag), "fragments over", length(full), "reference sequences using", length(cl), "threads (may take a few minutes) ...    ", paste(Sys.time()), "\n\n")
                          len <- length(full)
                          ## Run in paralell
                          fin_lst <- foreach(i=1:length(full), .packages=c("Biostrings", "stringr"), .final = function(x){names(x) <- names(full); return(x)}) %dopar% {
                                    seq_ref <- full[i]
                                    hits_vec <- do.call("c", lapply(query_strings, function(x){grepl(x, seq_ref)}))
                                    aligned_lst <- list(NA)
                                    for(t in 1:length(query_strings)){
                                                      if(t %in% which(hits_vec)){
                                                                      y <- as.data.frame(vmatchPattern(query_strings[[t]],  seq_ref, max.mismatch=mismatches, fixed=FALSE))
                                                                      y <- y[,c(1,3:5)]
                                                      }else{
                                                                      y <-data.frame(group=NA, start=NA, end=NA, width=NA)
                                                                      y$group <- nrow(y)
                                                                   }
                                                      colnames(y) <- c("n_hits", "Align_start", "Align_end", "Align_width")
                                                      aligned_lst[[t]] <- y
                                                      names(aligned_lst)[t] <- query_strings[[t]]
                                                      }
                                  	    target_lst  <- list(Ref_seq=seq_ref, Alignments=do.call("rbind", aligned_lst))
                                  	    target_lst[[2]] <- target_lst[[2]][!is.na(target_lst[[2]]$Align_start),]
                                  	    if(nrow(target_lst[[2]])<1){
                                  	               target_lst[[2]][1,]  <- rep("no_hits", each=4)
                                  	            }
                                        return(target_lst)
                                  	    }
                          cat("Finished at", paste(Sys.time()), "\n\n")
                          stopCluster(cl)
                          return(fin_lst)
                          }
