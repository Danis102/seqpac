#' Makes annotation from R imported reannotation mapping
#'
#' \code{import_reanno} makes a reannotation list.
#'
#' Given the path to the imported R reannotation files
#' (reanno_mis0/1/2/3/4/5.Rdata) generated by import_reanno, this function will
#' summarize the reannotation files into one output that matches the order of
#' sequences in a PAC object.
#'
#' @family PAC reannotation
#'
#' @seealso \url{https://github.com/junchaoshi/sports1.0} for download and
#'   documentation about running Sports.
#'   \url{http://bowtie-bio.sourceforge.net/index.shtml} for information about
#'   Bowtie \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param reanno_path Path to a directory where reannotation .Rdata files can be
#'   found.
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names.
#' @param mis_fasta_check Logical TRUE/FALSE if checking against anno_misX.fa
#'   should be done availble in the same folder as the reannoration files.
#'
#' @param threads Integer indicating the number of paralell workers to be used.
#'
#' @return List of dataframes with additional information from reannotation
#'   using bowtie. If \emph{mis_fasta_check} is specified the function will look
#'   for a \emph{anno_misX.fa} (X = the file with the highest number) previoulsy
#'   generated by the reannotation workflow. This file is used to double check
#'   so that no sequences are missing before and after reannotation.
#'
#' @examples
#' reanno_path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Jan_IOR_200130/Data/Single/Processed_Pipeline31_26-02-20/R_files/"
#' 
#' load(paste0(reanno_path, "PAC_master.Rdata"))
#' Full_anno <- make_reanno(reanno_path, PAC=PAC_master, mis_fasta_check=TRUE, threads=10)    # Complete use
#' identical(rownames(PAC_master$Anno), rownames(Full_anno$Overview)) 
#'
#'
#' @export

make_reanno <- function(reanno_path, PAC, mis_fasta_check=FALSE, threads=1){
                  files <- list.files(reanno_path, pattern="Full_reanno_mis0|Full_reanno_mis1|Full_reanno_mis2|Full_reanno_mis3|Full_reanno_mis4|Full_reanno_mis5", full.names = TRUE)
                  seqs <- (seq(1:length(files)))-1
                  reanno_lst <- list(NA)
                  for(i in 1:length(files)){load(files[i])
                                            reanno_lst[[i]] <- reanno
                                            names(reanno_lst)[i] <- paste0("mis", seqs[i])
                                            }
              cat("\nReorganizing and matching reannotation files with PAC ...\n")
                  PAC_seq <- rownames(PAC$Anno)
                  reanno_lst_match <- lapply(reanno_lst, function(x){
                                        match_lst  <- parallel::mclapply(x,  mc.cores=threads, function(y){
                                                              y$seq <- as.character(y$seq)
                                                              y$ref_hits <- as.character(y$ref_hits)
                                                              anno_match <- y[match(PAC_seq, y$seq), ]
                                                              anno_match$seq[is.na(anno_match$seq)] <- PAC_seq[is.na(anno_match$seq)]
                                                              anno_match$mis_n <- gsub("misNA", "mis0", anno_match$mis_n)
                                                              stopifnot(identical(PAC_seq, anno_match$seq))
                                                              return(anno_match)
                                                              })
                                        return(match_lst)
                                        })
              cat("\nGenerating the overview file ...\n")
                  stopifnot(any(do.call("c", lapply(reanno_lst_match, function(t){identical(names(reanno_lst_match[[1]]), names(t))}))))
                  stopifnot(any(do.call("c", lapply(reanno_lst_match, function(t){  do.call("c", lapply(t, function(g){identical(reanno_lst_match[[1]][[1]]$seq, g$seq)}))}))))
                  bio_cat <- length(reanno_lst_match[[1]])
                  df_fin <- data.frame(matrix(NA, ncol=bio_cat, nrow=length(PAC_seq)), row.names=PAC_seq)
                  for (bio in 1:bio_cat){
                                df <- do.call("cbind", lapply(reanno_lst_match, function(x){return(x[[bio]]$mis_n)}))
                                vect <- apply(df, 1, function(x){ return(gsub("NA", "", paste(x, collapse="")))})
                                df_fin[, bio] <- vect 
                                }
                  colnames(df_fin) <- names(reanno_lst_match[[1]])                   
                  df_fin[df_fin == ""] <- "_"
                  vect_mis <- do.call("paste", as.list(df_fin)) 
                  df_fin$Any_hit <- ifelse(vect_mis == paste0(rep("_", times=bio_cat), collapse=" ") , "No_anno", "Hit") 
                  df_fin$Mis0_hit <- ifelse(grepl("mis0", vect_mis) , "Hit", "No_hit")
              if(mis_fasta_check==TRUE){
                  cat("\nChecking the anno_mis_fasta file was specified by user.\n")
                  cat("Will try to read this file from same directory as the reannotion files.\n")
                  anno_mis_fls <- list.files(reanno_path, pattern = "anno_mis\\d.fa")
                  ns <- max(as.integer(gsub("anno_mis|.fa", "", anno_mis_fls)))
                  file_nam <- paste0("anno_mis", ns, ".fa")
                  noAnno_fasta <- Biostrings::readDNAStringSet(paste0(reanno_path,"/", file_nam))
                  logi_no_anno <- df_fin$Any_hit=="No_anno"
                  logi_olap <-  rownames(df_fin)[df_fin$Any_hit=="No_anno"] %in% gsub("NO_Annotation_", "", names(noAnno_fasta))
                  cat("Of the ", length(logi_no_anno[logi_no_anno==TRUE]), "missing sequences in the reannotation files\n")
                  cat(paste0(length(logi_olap[logi_olap==TRUE]), " (", round(length(logi_olap[logi_olap==TRUE])  / length(logi_no_anno[logi_no_anno==TRUE])*100, digits=2), "%) were found in ", file_nam, ".\n"))
                        if(!length(logi_no_anno[logi_no_anno==TRUE])-length(logi_no_anno[logi_no_anno==TRUE])==0){warning(paste0("Not all missing annotations were found in ", file_nam, ". This indicates that something has went wrong in the reannotation workflow.\n"))
                        }else{cat("Good! This is how it should be...\n")}
                  }
              return(list(Overview=df_fin, Full_anno=reanno_lst_match))
              cat("Done!\n")
              }
