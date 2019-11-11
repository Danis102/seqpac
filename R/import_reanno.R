#' Imports annotation from reannotation mapping
#'
#' \code{import_reanno} imports bowtie reannotation files.
#'
#' Given the path to bowtie mapping output, \code{import_reanno} will attempt to
#' read these files into R and generate a list of unordered data.frames where
#' each row summarizes all annotations for a given sequence. \code{check_reanno}
#' will check if sequences are listed in any of the imported reannotation
#' data.frames. \code{add_reanno} will add the new annotations to a PAC-list
#' object.
#'
#' @section Important: Re-annotation must have been done using Bowtie default
#'   output settings, where output files specifically have been named '.txt'.
#'   The function will not work with other formats such as SAM/BAM formats.
#'
#' @section Important: The basenames of the bowtie output will also be used to
#'   annotated what type of annotation. Example: bowtie -v 3 -a -f
#'   '<piRBase.fasta>' '<master_anno.fasta>' piRNA.txt Will annotate sequences
#'   as piRNA since the txt-file output was named 'piRNA'
#'
#' @family PAC reannotation
#'
#' @seealso \url{https://github.com/junchaoshi/sports1.0} for download and
#'   documentation about running Sports.
#'   \url{http://bowtie-bio.sourceforge.net/index.shtml} for information about
#'   Bowtie \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names.
#' @param bowtie_path Path to a directory where bowtie output files can be
#'   found.
#' @param report Character vector "biotype" or "all"
#'
#' @param threads Integer stating the number of parallell jobs.
#'
#' @return Data.frame with additional information from reannotation using
#'   bowtie. If \emph{report="biotype"}, the function will report biotype as the
#'   file basename of a given reference used in the bowtie reannotation. If
#'   \emph{report="all"} will give the annotation as full names used in the
#'   fasta file used as reference in the bowtie reannotation. Caution:
#'   References with lots of redudancy, such as piRBase fasta files, will
#'   automatically run with \emph{report="biotype}).
#'
#' @examples
#' reanno <- import_reanno(bowtie_path, threads=8)
#' bowtie_path <- "/data/Data_analysis/Projects/Drosophila/Specific_projects/Test_temp/Processed_Sports_09-10-19/R_files"
#' bowtie_path <- "/data/Data_analysis/Projects/Drosophila/Specific_projects/Test_temp/Processed_Pipeline_3.1_11-11-19/R_files"
#' @export
import_reanno <- function(bowtie_path, base="*.txt", threads=1, report="biotype"){
                            require(data.table)
                            require(parallel)
                            require(pbmcapply)
                            require(plyr)
                            require(stringr)
                            files <- list.files(bowtie_path, pattern = base, full.names=TRUE)
                          ## Check bowtie format (8 columns; "IIIIIII" present in column 6; column 4 is an integer)
                                      row1 <- lapply(as.list(files), function(f){
                                                    x <- try(read.delim(f, nrows=1, header=FALSE), silent = TRUE)
                                                    if(inherits(x, "try-error"))
                                                    return(data.frame(V1="No_hits"))
                                                    else
                                                    return(x)
                                                    })
                                      form_logi <- lapply(row1, function(x){
                                                    if(!ncol(x) == 8){return(FALSE)
                                                    }else{form_logi[[i]] <- return(sum(c(grepl("IIIIIII", as.character(x[,6])), is.integer(x[,4]))) == 2)
                                                    }
                                                    })
                          ## Give some feedback
                            cat("\nFound ", length(files), " files of which ", sum(do.call("c", form_logi)), " had the correct bowtie format.\n")
                                  if(any(do.call("c", lapply(row1, function(x){as.character(x[1,1]) == "No_hits"})))){warning("Some bowtie outputs were empty indicating no hits at all.\n  These will be missing in the result.")}
                            cat("Now importing bowtie files to R (may take some time) ...\n")
                          ## Read files iver
                            files <- files[do.call("c", form_logi)]
                            bowtie_out <- parallel::mclapply(as.list(files), mc.cores = threads, function(x){
                                                                            data.table::setDTthreads(threads=1)
                                                                            bow_out <- data.table::fread(x, header=FALSE, select = c(3,5,8), data.table=FALSE)
                                                                            bow_out$V8 <- as.character(bow_out$V8)
                                                                            bow_out$V8[nchar(as.character(bow_out$V8)) ==0] <- NA
                                                                            if(grepl("piRBase", x)|nrow(bow_out)>10000000|report=="biotype"){
                                                                                                    if(grepl("piRBase", x)|nrow(bow_out)>10000000){
                                                                                                      warning("The bowtie mapping was done against piRBase or was extraordinary large.\nDue to extreme redudancy risk, only the biotype category will be reported, not individual reference ids.")}
                                                                                                    uni  <- unique(bow_out$V5)
                                                                                                    mis <- unique(bow_out$V8)
                                                                                                    if(any(is.na(uni))){
                                                                                                        mis <- paste0("_(NA)")
                                                                                                    }else{
                                                                                                        mis <- stringr::str_count(mis, ":")
                                                                                                        mis <- as.character(unique(mis))
                                                                                                        stopifnot(length(mis)=="1")
                                                                                                        mis <- paste0("_(mis", mis, ")")
                                                                                                        }
                                                                                                    bow_out <- data.frame(seq=uni, ref_hits= mis)
                                                                            }else{
                                                                                                    bow_out <- data.frame(seq=bow_out$V5 , ref_hits=paste0(bow_out$V3, "_(", bow_out$V8, ")"))
                                                                                                    bow_out <- plyr::ddply(bow_out, .(seq), function(x) paste(unique(x[,"ref_hits"]), collapse="; "))
                                                                                                    }
                                                                            colnames(bow_out)[2] <- "ref_hits"
                                                                            return(bow_out)}
                                                                            )
                            names(bowtie_out) <- gsub(paste0(base), "",  basename(files))
                            names(bowtie_out) <- gsub("Re-anno_", "", names(bowtie_out))
                            return(bowtie_out)
}

