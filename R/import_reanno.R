#' Imports annotation from reannotation mapping
#'
#' This function imports bowtie output and summarizes the content.
#'
#' Given the path to alignment outputs from bowtie, \code{import_reanno} will
#' attempt to read these files into R and generate a list of unordered
#' data.frames where each row summarizes all annotations for a given sequence.
#' It is called by \code{map_reanno} to generate summarized data from bowtie
#' output saved as .Rdata files. These files can then be re-imported into R
#' where information is organized using \code{make_reanno}
#' and finally added to the annotation table (PAC$Anno) of the original PAC-list
#' object using \code{add_reanno}.
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
#' @seealso  \url{http://bowtie-bio.sourceforge.net/index.shtml} for information
#'   about Bowtie and for Rbowtie:
#'   \url{https://www.bioconductor.org/packages/release/bioc/html/Rbowtie.html}.
#'   \url{https://github.com/Danis102} for updates on the current package.
#'   
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names.
#'   
#' @param bowtie_path Path to a directory where bowtie output files can be
#'   found.
#'   
#' @param report Character vector indicating what to report "minimum" or "full"
#'   (default="minimum").   
#' 
#' @param coord Logical whether or not mapping coordinates should be reported
#'   when report="full".
#'   
#' @param reduce Character indicating a reference name (without file type
#'   extension) that should be exempted from report="full" and instead will
#'   generate a minimum report.
#'   
#' @param threads Integer stating the number of parallell jobs.
#'
#' @return List of data frames with additional information from reannotation
#'   files generated with bowtie. If \emph{report="minimum"}, the function
#'   will report each hit only as the number of mismatches for a given reference
#'   file. If \emph{report="full"} the full name reported in the fasta
#'   file used as reference in the bowtie reannotation will be reported. If a
#'   reference name is specified in \emph{reduce}, this reference is excempted
#'   from \emph{report="full"} and is instead reported as
#'   \emph{report="minimum"}.
#'   
#'   Caution: Large references with lots of redudancy (such as pirBase in some
#'   species) will cause massive character strings if \emph{report="full"} with
#'   no restrictions. Specifing such references in
#'   \emph{reduce=<reference_names>} will circumvent this problem.
#'
#' @examples
#' bowtie_path <- "/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/reanno"
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

import_reanno <- function(bowtie_path, threads=1, coord=FALSE, report="minimum", reduce=NULL){
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
          }else{form_logi[[x]] <- return(sum(c(grepl("IIIIIII", as.character(x[,6])), is.integer(x[,4]))) == 2)
          }
        })
  
  ## Give some feedback
  cat("\nFound ", length(files), "txt-files of which ", sum(do.call("c", form_logi)), " had the correct bowtie format.\n")
  if(any(do.call("c", lapply(row1, function(x){as.character(x[1,1]) == "No_hits"})))){warning("Some bowtie outputs were empty indicating no hits at all.\n  These will be missing in the result.")}
  cat(paste0("\n******************************************************"))
  
  ## Entering import loop
  files <- files[do.call("c", form_logi)]
  data.table::setDTthreads(threads)
  bowtie_out_lst <- list(NA)
  for (t in 1:length(files)){
    cat(paste0("\nImport and reorganize ", basename(files)[t], "\n"))
    nam <- gsub(paste0(base, "|_piRBase"), "",  basename(files)[t])
    if(coord==TRUE){bow_out <- data.table::fread(files[t], header=FALSE, select = c(3,4,1,8), data.table=TRUE)}
    if(coord==FALSE){bow_out <- data.table::fread(files[t], header=FALSE, select = c(3,1,8), data.table=TRUE)}
    bow_out$V8 <- as.character(bow_out$V8)

    uni  <- unique(bow_out$V5)
    mis <- unique(bow_out$V8)
    n_mis <- stringr::str_count(mis, ":")
    n_mis <- as.character(unique(n_mis))
    if(is.na(n_mis)){n_mis <- 0}
    stopifnot(length(n_mis)=="1")
    n_mis <- paste0("mis", n_mis)
    
    ## Generate report from imported bowtie files
    if(report=="minimum"){
      bowtie_out_lst[[t]] <- suppressPackageStartupMessages(data.table::data.table(seq=uni, mis_n=n_mis, mis_where="mini_report", ref_hits=nam))
      names(bowtie_out_lst)[t] <- nam
    }
    
    if(report=="full"){
      if(nam %in% reduce){
        cat(paste0("\n|-------> ", nam, " was specified as reduced; minimum information will be extracted ..."))
        bowtie_out_lst[[t]] <- suppressPackageStartupMessages(data.table::data.table(seq=uni, mis_n=n_mis, mis_where="mini_report", ref_hits=nam))
        names(bowtie_out_lst)[t] <- nam
      }else{
        cat("\n|-------> Compiling data for full report (may take a while)...")
        bow_splt <- split(bow_out, bow_out$V5) 
        rm(bow_out)
        gc(reset=TRUE)
        
    #### Compile everything with multi-threading 
        require("foreach", quietly = TRUE)
        chk_size <- ceiling(length(bow_splt)/100) # foreach combine every 100 instances
        chnks1 <-as.integer(seq(from=1, to=length(bow_splt), by=chk_size))
        chnks2 <-as.integer(seq(from=0, to=length(bow_splt), by=chk_size))
        chnks2 <- c(chnks2[-1], length(bow_splt))
        chnks_rng <- list(chnks1, chnks2)

        doParallel::registerDoParallel(threads) # Do not use parallel::makeClusters!!!
        bowtie_out_lst[[t]] <- foreach(s=1:length(chnks_rng[[1]]), .inorder = FALSE, .combine = "rbind", .export= c("chnks_rng", "bow_splt"), .packages=c("data.table")) %dopar% {
              compile_lst <- lapply(bow_splt[chnks_rng[[1]][s]:chnks_rng[[2]][s]], function(x){
                  uni_mis <- unique(x$V8)
                  uni_mis <- unique(do.call("c", stringr::str_split(uni_mis, ",")))
                  uni_mis <- uni_mis[order(as.integer(gsub( ":.*$", "", uni_mis )))]
                  if(any(is.na(uni_mis))){uni_mis <- "mis0"}
                  uni_mis <- paste(uni_mis, collapse="|")
                  if(coord==TRUE){hits <- paste(unique(paste(x$V3, x$V4, sep=":")), collapse="|")}
                  if(coord==FALSE){hits <- paste(unique(x$V3), collapse="|")}
                  fin <- suppressPackageStartupMessages(data.table::data.table(mis_n=n_mis, mis_where=uni_mis, ref_hits=hits))
                  return(fin)
                  })
              bow_fin <- suppressPackageStartupMessages(data.table::rbindlist(compile_lst, idcol=TRUE))
              return(bow_fin)
        }
      doParallel::stopImplicitCluster()
      names(bowtie_out_lst)[t] <- nam 
      }
    }
    cat(paste0("\n|-------> Done ", nam, "!\n"))
    cat(paste0("\n******************************************************"))
  }
  return(bowtie_out_lst)
  cat("All done!")
}

