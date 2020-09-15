#' Advanced sequence mapping of a PAC object 
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
#' @param ref Character indicating the path to the fasta (.fa) reference
#'   file or a DNAStringSet with already loaded reference sequences.
#'
#' @param mismatches Integer indicating the number of mismatches that should be
#'   allowed in the mapping.
#'
#' @param N_up Character indicating a sequence that should be added to the
#'   reference at the 5' end prior to mapping. A wild card nucleotides "NNN"
#'   (any of C, T, G, A) can for example be added for mapping non-perfect
#'   reference hits. No nucleotides are added by default.
#'   
#' @param N_down Character. Same as N_up but indicating a sequence that should
#'   be added to the reference at the 3' end. Useful for tRNA analysis where the
#'   reference do not contain pre-processed tRNA. Setting N_down="NNN" or "CCA"
#'   (in many species CCA is added to mature tRNA) will allow mapping against
#'   the mature tRNA. No nucleotides are added by default.
#'   
#' @param threads Integer indicating the number of parallel processes that
#'   should be used.
#'
#' @param par_type Character indicating if "PSOCK" or "FORK" should be parsed to
#'   parallel::makeCluster (default="PSOCK").
#'   
#'   
#' @param report_string Logical whether an alignment string that shows in
#'   character where sequences align against the reference. Works well with
#'   tRNA, but makes the Alignments object difficult to work with when longer
#'   references are used (default=FALSE).
#'   
#' @return Stacked list, where each object on the highest level contains:
#'                    (Object 1) Reference name and sequence. 
#'                    (Object 2) Dataframe showing the mapping results of
#'                               each quiary sequence that mapped to Object 1.
#'
#' @examples
#' 
#' 
#' 
#' #-------------------------------------------------------------------------#
#' ### For coverage plots ###
#' #-------------------------------------------------------------------------#
#' 
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#' 
#' PAC_filt <- PAC_filter(PAC_all, size=c(16,70), threshold=10, coverage=5, type="counts", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#' PAC_filt <- PAC_rpm(PAC_filt)
#' PAC_filt <- PAC_filter(PAC_filt, size=c(16,70), threshold=10, coverage=5, type="rpm", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#' 
#' ## Remove corrupt samples and make means 
#' Ph_trg <- as.character(PAC_filt$Pheno$Sample[!PAC_filt$Pheno$Sample %in% "Inx24_200130_S12"])
#' PAC_filt <- PAC_filter(PAC_filt, pheno_target= list("Sample", Ph_trg))
#' 
#' PAC_filt$Pheno$Groups <- paste(do.call("rbind", strsplit(as.character(PAC_filt$Pheno$SampleProject), "_" ))[,1], PAC_filt$Pheno$Method, PAC_filt$Pheno$Method_tag, PAC_filt$Pheno$Tag, sep="_")
#' 
#' ## Make summaries
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "means", pheno_target=list("Groups", unique(PAC_filt$Pheno$Groups)))
#' 
#' ## Mapping
#' map_rRNA <- PAC_mapper(PAC_filt, ref="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/rRNA_reanno/drosophila_rRNA_all.fa", threads=12)
#' 
#' All_plots <- lapply(as.list(Smry_trg_all), function(x){PAC_covplot(PAC_filt, map_rRNA, summary_target = list("means_[Groups]", x), xseq=FALSE, style="line", colour="red")})
#'
#' cowplot::plot_grid(All_plots[[1]][[7]], All_plots[[2]][[7]], All_plots[[3]][[7]], 
#'                   All_plots[[4]][[7]], All_plots[[5]][[7]], All_plots[[6]][[7]],
#'                   All_plots[[7]][[7]], All_plots[[8]][[7]], All_plots[[9]][[7]], 
#'                   All_plots[[10]][[7]], All_plots[[11]][[7]], All_plots[[12]][[7]],
#'                   All_plots[[13]][[7]], All_plots[[14]][[7]], All_plots[[15]][[7]],
#'                   nrow = 5, ncol = 3)
#'
#'
#' #-------------------------------------------------------------------------#
#' #' ### For mismap ###
#' #-------------------------------------------------------------------------#
#' 
#' load(file="/home/danis31/OneDrive/Programmering/Programmering/Pipelines/Drosophila/Pipeline_3.1/seqpac/dm_test_PAC.Rdata")
#' 
#' ref_path <- "/data/Data_analysis/Genomes/Drosophila/dm6/tRNA/tRNA.fa"
#' full <- Biostrings::readDNAStringSet(ref_path)
#' ref <- full[grepl("Glu-CTC-3-1|Lys-CTT-1-1 ", names(full))]
#' 
#' 
#' #PAC_filt <- PAC_filter(PAC_all, threshold=5, coverage=4, type="counts", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#' 
#' map_refs <- PAC_mapper(PAC_filt, ref=ref, threads=8, mismatches=3)
#' 
#' 
#'
#' 
#' 
#' @export

PAC_mapper <- function(PAC, ref, mismatches=0, threads=1, par_type="PSOCK", N_up="", N_down="", report_string=FALSE){
                          require(foreach)
                          
                          # Setup
                          Anno_frag  <- Biostrings::DNAStringSet(rownames(PAC$Anno))
                          query_strings <- as.list(as.character(rownames(PAC$Anno)))

                          # Read files
                          if(class(ref)=="character"){
                                          cat("Reading full length rRNA from file ...")
                                          full <- Biostrings::readDNAStringSet(ref)
                                          }else{full <- ref}  
                          nams_full <- names(full) # Names are lost in the next step
                          full <- Biostrings::DNAStringSet(paste(N_up, full, N_down, sep=""))
                          names(full) <- nams_full

                          ## Aligning using parallel processing
                          cat("Now aligning", length(Anno_frag), "fragments over", length(full), "reference sequences using", threads, "threads (may take a few minutes) ...    ", paste(Sys.time()), "\n")
                          len <- length(full)
                          
                          cl <- parallel::makeCluster(threads, type = par_type)
                          doParallel::registerDoParallel(cl)
                          
                          ## Parallelize references 
                          # fin_lst <- foreach(i=1:length(full), .packages=c("Biostrings", "stringr"), .final = function(x){names(x) <- names(full); return(x)}) %dopar% {
                          #           seq_ref <- full[i]
                          #           hits_vec <- do.call("c", lapply(query_strings, function(x){grepl(x, seq_ref)}))
                          #           aligned_lst <- list(NA)
                          #           for(t in 1:length(query_strings)){
                          #                             if(t %in% which(hits_vec)){
                          #                                             y <- as.data.frame(vmatchPattern(query_strings[[t]],  seq_ref, max.mismatch=mismatches, fixed=FALSE))
                          #                                             y <- y[,c(1,3:5)]
                          #                             }else{
                          #                                             y <-data.frame(group=NA, start=NA, end=NA, width=NA)
                          #                                             y$group <- nrow(y)
                          #                                          }
                          #                             colnames(y) <- c("n_hits", "Align_start", "Align_end", "Align_width")
                          #                             aligned_lst[[t]] <- y
                          #                             names(aligned_lst)[t] <- query_strings[[t]]
                          #                             }
                          #         	    target_lst  <- list(Ref_seq=seq_ref, Alignments=do.call("rbind", aligned_lst))
                          #         	    target_lst[[2]] <- target_lst[[2]][!is.na(target_lst[[2]]$Align_start),]
                          #         	    if(nrow(target_lst[[2]])<1){
                          #         	               target_lst[[2]][1,]  <- rep("no_hits", each=4)
                          #         	            }
                          #               return(target_lst)
                          # }
                          
                          ## Parallelize sequences 
                          fin_lst <- list(NA)
                          for(i in 1:length(full)){ 
                                    cat(paste0("\nAligning against:\n ", names(full)[i], "\n Start ", Sys.time()))
                                    seq_ref <- full[i]
                                    aligned_lst <- foreach(t=1:length(query_strings), .packages=c("Biostrings", "stringr"), .final = function(x){names(x) <- paste0(Anno_frag); return(x)}) %dopar% {
                                                      y <- as.data.frame(Biostrings::vmatchPattern(query_strings[[t]],  seq_ref, max.mismatch=mismatches, fixed=FALSE))
                                                      
                                                      if(!nrow(y)< 1){
                                                        y <- y[,c(1,3:5)]
                                                        y$group <- nrow(y)
                                                        colnames(y) <- c("n_hits", "Align_start", "Align_end", "Align_width")
                                                        return(y)
                                                      }else{
                                                        return(NULL)}
                                      }
                                    aligned   <- do.call("rbind", aligned_lst[unlist(lapply(aligned_lst, function(x){!is.null(x)}))])                  
                                  	target_lst  <- list(Ref_seq=seq_ref, Alignments=aligned)
                                  	if(is.null(target_lst[[2]])){target_lst[[2]] <- data.frame(n_hits="no_hits", Align_start="no_hits", Align_end="no_hits", Align_width="no_hits")}
                                    fin_lst[[i]] <- target_lst
                                    names(fin_lst)[i] <- names(full)[i]
                                    cat(paste0("\n Done ", Sys.time()))
                          }

                          
                          
                          if(report_string==TRUE){
                          fin_lst <- lapply(fin_lst, function(x){
                                                          x$Ref_seq -> ref
                                                          x$Alignments -> algn
                                                          n_ref <- nchar(as.character(ref))
                                                          algn_lst <- split(algn, factor(row.names(algn), levels=row.names(algn)))
                                                          positions_lst <- foreach(j=1:length(algn_lst), .final=function(y){names(y) <- names(algn_lst);return(y)})  %dopar% {
                                                                                  ref=ref
                                                                                  n_ref=n_ref
                                                                                  algn_str <- paste(strrep("-", times=(algn_lst[[j]]$Align_start)-1), rownames(algn_lst[[j]]), strrep("-", times= n_ref-(algn_lst[[j]]$Align_end)), sep="")
                                                                                  return(algn_str)
                                                                                  }
                                                          df <- cbind(algn, data.frame(Align_string=as.character(paste(do.call("c", positions_lst))), stringsAsFactors=FALSE))
                                                          return(list(Ref_seq=ref, Alignments=df))
                                                          })
                          }
                          
                          parallel::stopCluster(cl)
                          return(fin_lst)
                          }
                          