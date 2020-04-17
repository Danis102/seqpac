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
#'   
#' @param bowtie_path Path to a directory where bowtie output files can be
#'   found.
#'   
#' @param report Character vector indicating what to report "minimum" or "full"
#'   (default="minimum")
#'  
#' @param reduce Character indicating a reference name (without file type
#'   extension) that should be exempted from report="full" and instead will
#'   generate a minimum report.
#'   
#' @param threads Integer stating the number of parallell jobs.
#'
#' @return List of data frames with additional information from reannotation
#'   files generated with bowtie -v -a. If \emph{report="minimum"}, the function
#'   will report each hit only as the number of mismatches for a given reference
#'   file. If \emph{report="full"} the full annotation reported in the fasta
#'   file used as reference in the bowtie reannotation will be reported. If a
#'   reference name is speficied in \emph{reduce}, this reference is excempted
#'   from \emph{report="full"} and is instead reported as
#'   \emph{report="minimum"}
#'   
#'   Caution: Large references with lots of redudancy (such as pirBase in some
#'   species) will case massive character strings if \emph{report="full"} with
#'   no restrictions. Specifing such references in
#'   \emph{reduce=<reference_names>} will circumvent this problem.
#'
#'   
#'   
#' @examples
#' bowtie_path <- "/data/Data_analysis/Projects/Drosophila/Specific_projects/Mar_Diet_4_6_2019/Processed_Pipeline31_10-03-20/R_files"
#' bowtie_path <- "/data/Data_analysis/Projects/Drosophila/Other/IOR/Jan_IOR_200130/Data/Single/Processed_Pipeline31_05-03-20/R_files"

#' reanno1 <- import_reanno(bowtie_path, report="full",  threads=1)
#' reanno2 <- import_reanno(bowtie_path, report="full",  threads=10)
#' reanno3 <- import_reanno(bowtie_path, report="full", reduce="piRNA", threads=1)
#' reanno4 <- import_reanno(bowtie_path, report="full", reduce="piRNA", threads=8)
#' reanno5 <- import_reanno(bowtie_path, report="minimum")
#' 
#' lapply(reanno5, head)
#' lapply(reanno5, names)
#' 
#' @export
  
import_reanno <- function(bowtie_path, threads=1, report="minimum", reduce="piRNA"){
                            require(data.table)
                            require(parallel)
                            require(plyr)
                            require(stringr)
                            base <- ".txt"
                            files <- list.files(bowtie_path, pattern = base, full.names=TRUE)
                            options(scipen=999)
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
                          ## Read files
                            files <- files[do.call("c", form_logi)]
                            data.table::setDTthreads(threads)
                            bowtie_out_lst <- list(NA)
                            for (i in 1:length(files)){
                                                    cat(paste0("\nImport and reorganize ", basename(files)[i], "\n"))
                                                    nam <- gsub(paste0(base, "|_piRBase"), "",  basename(files)[i])
                                                    file_len <- R.utils::countLines(files[i])[[1]] 
                                                    if(file_len>1000000){
                                                                  cat("\n|-------> Large file may take a bit longer ...\n") 
                                                                  chunk_lens <- lapply(as.list(1:ceiling(file_len/1000000)), function(x){  
                                                                                          fin <- 1000000*x
                                                                                          if(fin<=file_len){return(c(nrow=1000000, skip=(fin-1000000)))}
                                                                                          if(fin>file_len){
                                                                                                           nrw <- 1000000-(fin-file_len)
                                                                                                           skp <- file_len-nrw
                                                                                                           return(c(nrow=nrw, skip=skp))}
                                                                                          })
                                                                  read_lst <- list(NA)
                                                                  for(it in 1:length(chunk_lens)){
                                                                          bow <- data.table::fread(files[i], header=FALSE, select = c(3,5,8), data.table=TRUE, nrow = chunk_lens[[it]][1], skip= chunk_lens[[it]][2])
                                                                          bow$V8 <- as.character(bow$V8)
                                                                          read_lst[[it]] <- bow
                                                                          rm(bow)
                                                                          gc()
                                                                          Sys.sleep(0.001)
                                                                          cat(paste0("\r|-------> Reading chunks ", round(it/length(chunk_lens)*100), "%"))
                                                                          flush.console()
                                                                          }
                                                                  
                                                                  bow_out <- data.table::rbindlist(read_lst)
                                                                  rm(read_lst)
                                                   }else{
                                                                  bow_out <- data.table::fread(files[i], header=FALSE, select = c(3,5,8), data.table=TRUE)
                                                                  bow_out$V8 <- as.character(bow_out$V8)
                                                                  }
                                                    stopifnot(nrow(bow_out)==file_len)              
                                                    uni  <- unique(bow_out$V5)
                                                    mis <- unique(bow_out$V8)
                                                    n_mis <- stringr::str_count(mis, ":")
                                                    n_mis <- as.character(unique(n_mis))
                                                    if(is.na(n_mis)){n_mis <- 0}
                                                    stopifnot(length(n_mis)=="1")
                                                    n_mis <- paste0("mis", n_mis)
                 
                                                if(report=="minimum"){
                                                                bowtie_out_lst[[i]] <- data.table::data.table(seq=uni, mis_n=n_mis, mis_where="mini_report", ref_hits=nam)
                                                                names(bowtie_out_lst)[i] <- nam
                                                                }
                                                    
                                                if(report=="full"){
                                                            if(nam %in% reduce){
                                                                cat(paste0("\n|-------> ", nam, " was specified as reduced; minimum information will be extracted ..."))
                                                                bowtie_out_lst[[i]] <- data.table::data.table(seq=uni, mis_n=n_mis, mis_where="mini_report", ref_hits=nam)
                                                                names(bowtie_out_lst)[i] <- nam
                                                            }else{
                                                                bow_splt <- split(bow_out, bow_out$V5) 
                                                                cat("\n|-------> Compiling data ...")
                                                                gc(reset=TRUE)
                                                                compile_lst <- lapply(bow_splt, function(x){
                                                                                      uni_mis <- unique(x$V8)
                                                                                      uni_mis <- unique(do.call("c", stringr::str_split(uni_mis, ",")))
                                                                                      uni_mis <- uni_mis[order(as.integer(gsub( ":.*$", "", uni_mis )))]
                                                                                      if(is.na(uni_mis)){uni_mis <- "mis0"}
                                                                                      uni_mis <- paste(uni_mis, collapse="|")
                                                                                      hits <- paste(unique(x$V3), collapse=";")
                                                                                      fin <- data.table::data.table(mis_n=n_mis, mis_where=uni_mis, ref_hits=hits)
                                                                                      return(fin)
                                                                                      })
                                                                 rm(bow_splt)
                                                                 bow_fin <- data.table::rbindlist(compile_lst, idcol=TRUE)
                                                                 colnames(bow_fin)[1] <- "seq" 
                                                                 stopifnot(any(bow_fin$seq %in% uni))
                                                                 bowtie_out_lst[[i]] <- bow_fin 
                                                                 names(bowtie_out_lst)[i] <- nam
                                                            }
                                                          }
                                          cat(paste0("\n|-------> Done ", nam, "!\n"))
                                          cat(paste0("\n---------------------------------------"))
                                          }
                            return(bowtie_out_lst)
                            cat("All done!")
                                  }
                          