#' Advanced sequence mapping of a PAC object 
#'
#' \code{PAC_mapper} Map sequences against a small reference.
#'
#' Given a PAC object and the path to a fasta reference file, this function will
#' map sequences in the PAC using a 'backdoor' into the reanno workflow.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object.
#'
#' @param ref Character indicating the path to the fasta (.fa) reference file or
#'   a DNAStringSet with already loaded reference sequences. If a Bowtie index
#'   is missing for the reference, PAC_mapper will attempt to temporally
#'   generate such index automatically. Thus, large large references are
#'   discouraged. Please use the original reanno workflow for large references.
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
#' @param report_string Logical whether an alignment string that shows in
#'   character where sequences align against the reference. Works well with
#'   tRNA, but makes the Alignments object difficult to work with when longer
#'   references are used (default=FALSE).
#'
#' @param multi Character indicating how to deal with multimapping. If
#'   \code{multi="keep"}, query sequences that maps multiple times to the same
#'   reference sequence will be reported >1 times in the output (indicated by
#'   .1, .2, .3 etc. in the reported sequence name). If \code{multi="remove"}
#'   (default), then all multimapping sequences will be removed, resulting in 1
#'   row for each query sequence that maps to the target reference sequence. The
#'   function will always give a warning if a query sequence maps to multiple
#'   sites within a reference sequence. However, this function discriminate
#'   multimapping only within a reference sequence. Thus, if the fasta input
#'   contains multiple reference sequences, a query sequence may be reported in
#'   multiple references sequences.
#'   
#' @param override Logical whether or not the map_reanno function should prompt
#'   you for a question if there are files in the temporary path. As default,
#'   override=FALSE will prevent deleting large files by accident, but requires
#'   an interactive R session. Setting override=TRUE may solve non-interactive
#'   problems.
#'   
#'   
#' @return Stacked list, where each object on the highest level contains:
#'                    (Object 1) Reference name and sequence. 
#'                    (Object 2) Data.frame showing the mapping results of
#'                               each query sequence that mapped to Object 1.
#'
#' @examples
#' 
#' 
#' ## Load PAC-object data ###
#'  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#' 
#' ## Make summaries and extract rRNA
#' pac <- PAC_summary(pac, norm = "cpm", type = "means", 
#'                    pheno_target=list("stage", unique(pac$Pheno$stage)))
#'                    
#'                    
#' pac_rRNA <- PAC_filter(pac, anno_target = list("Biotypes_mis0", "rRNA"))
#' pac_tRNA <- PAC_filter(pac, anno_target = list("Biotypes_mis0", "tRNA"))
#'
#'
#' ## Give path to a fasta reference (with or without bowtie index)
#' #  (Here we use an rRNA fasta included in seqpac) 
#' 
#' ref_rRNA <- system.file("extdata/rrna", "rRNA.fa", 
#'                          package = "seqpac", mustWork = TRUE)
#'                          
#' ref_tRNA_no_index <- system.file("extdata/trna_no_index", "tRNA_copy.fa", 
#'                          package = "seqpac", mustWork = TRUE)                         
#'
#'                                                                                                                                  
#' ## Map using PAC-mapper                          
#' map_rRNA <- PAC_mapper(pac_rRNA, mismatches=0, 
#'                         threads=1, ref=ref_rRNA, override=TRUE)
#'                                                                                                    
#' ## Now try a fasta with no bowtie index using PAC-mapper                                                                 
#' map_tRNA <- PAC_mapper(pac_tRNA, mismatches=0, 
#'                         threads=1, ref=ref_tRNA_no_index, override=TRUE)                        
#'                                                 
#' ## Plot rRNA according to embryonic stage using PAC_covplot                       
#' cov_rRNA<- PAC_covplot(pac_rRNA, map_rRNA, 
#'                         summary_target = list("cpmMeans_stage"), 
#'                         xseq=FALSE, style="line", 
#'                         color=c("red", "black", "blue"))
#'                         
#' cowplot::plot_grid(cov_rRNA[[1]], cov_rRNA[[2]], cov_rRNA[[3]], 
#'                     cov_rRNA[[4]], nrow=2, ncol=2)
#'
#'
#'
#' ## Plot tRNA using xseq=TRUE gives you reference sequence as X-axis:
#' # (OBS! Long reference will not )                     
#' cov_tRNA <- PAC_covplot(pac_tRNA, map_tRNA, 
#'                         summary_target = list("cpmMeans_stage"), 
#'                         xseq=TRUE, style="line", 
#'                         color=c("red", "black", "blue"))
#'
#' cov_tRNA[[1]]
#'                     
#' ## Explore the map-object                    
#' head(map_tRNA[[1]])
#' names(map_tRNA)
#' map_tRNA[[1]]
#'
#' ## Check wish reached decent number
#' # (OBS! This is a very down sampled dataset)
#' logi_hi <- unlist(lapply(map_tRNA, function(x){nrow(x$Alignments) > 10 }))
#' logi_lo <- unlist(lapply(map_tRNA, function(x){nrow(x$Alignments) > 2 }))
#' 
#' table(logi_hi)  
#' names(map_tRNA[logi_hi])
#' 
#' table(logi_lo)  
#' names(map_tRNA[logi_lo])
#'     
#' targets <- c("Ala-AGC-1-1", "Lys-CTT-1-13","Ser-AGA-2-2")
#' 
#' cov_tRNA_sub <- PAC_covplot(pac_tRNA, map_tRNA, 
#'                         summary_target = list("cpmMeans_stage"),
#'                         map_target = targets,
#'                         xseq=TRUE, style="line", 
#'                         color=c("red", "black", "blue"))                  
#'                                                       
#' cowplot::plot_grid(plotlist= cov_tRNA_sub) 
#' 
#' 
#' ###########################################################
#' ## Analyze range types with map_rangetype and PAC_trna functions
#' #
#' ## Download ss object from GtRNAdb 
#' # dest_path <- file.path("/some/path/to/destination/file/trna.tar.gz")
#' # web_path <- "http://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.tar.gz"
#' # download.file(url=web_path, destfile=dest_path)
#' # untar(dest_path, exdir= dirname(dest_path), files = "dm6-tRNAs-confidence-set.ss")
#' # ss_file <- "/some/path/to/dm6-tRNAs-confidence-set.ss"
#' #
#' # Classify fragments according to loop cleavage (small loops are omitted)       
#' # map_object_ss <- map_rangetype(<your_map_object>, type="ss", 
#' #                               ss=ss_file, min_loop_width=4)          
#' #
#' ## Remove reference tRNAs with no hits
#' # map_object_ss <-  map_object_ss[!unlist(lapply(map_object_ss, function(x){
#' #                     x[[2]][1,1] == "no_hits"
#' #                      }))]
#' # map_object_ss[[2]]
#' #
#' #
#' #
#' ###########################################################
#' ## Function classifying 5'-tRF, 5'halves, i-tRF, 3'-tRF, 3'halves
#' #
#' ## Set tolerance for classification as a terminal tRF
#' # tolerance <- 5  # 2 nucleotides from start or end of full-length tRNA)
#' #
#' ## Apply the tRNA_class function and make a tRNA type column
#' # pac_trna <- tRNA_class(pac_trna, map=map_object_ss, terminal=tolerance)
#' # head(pac_trna$Anno)
#' #
#' ## Now use PAC_trna to generate some graphs based on grand means
#' # trna_result <- PAC_trna(pac_trna, norm="cpm", filter = NULL,
#' #   join = TRUE, top = 15, log2fc = TRUE,
#' #   pheno_target = list("stage", c("Stage1", "Stage3")), 
#' #   anno_target_1 = list("type"),
#' #   anno_target_2 = list("class"))
#' #
#' # cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Grand_means,
#' #                    trna_result$plots$Log2FC_Anno_1,
#' #                    trna_result$plots$Percent_bars$Grand_means,
#' #                    nrow=1, ncol=3)
#' #
#' ## By setting join = FALSE you will get group means
#' # trna_result <- PAC_trna(pac_trna, norm="cpm", filter = NULL,
#' #   join = FALSE, top = 15, log2fc = TRUE,
#' #   pheno_target = list("stage", c("Stage1", "Stage3")), 
#' #   anno_target_1 = list("type"),
#' #   anno_target_2 = list("class"))
#' #
#' # cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Stage1,
#' #                    trna_result$plots$Expression_Anno_1$Stage3,
#' #                    trna_result$plots$Log2FC_Anno_1,
#' #                    trna_result$plots$Percent_bars$Stage1,
#' #                    trna_result$plots$Percent_bars$Stage3,
#' #                    nrow=1, ncol=5)
#'       
#'                    
#' @export

PAC_mapper <- function(PAC, ref, mismatches=0, multi="remove", 
                       threads=1, N_up="", N_down="", report_string=FALSE, 
                       override=FALSE){


## Setup
  j <- NULL
  ## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  ## Setup reference  
  if(methods::is(ref, "DNAStringSet")){
    cat("\nImporting reference from DNAStringSet ...")
    full <- ref
  }else{
    if(file.exists(ref)){
      cat("\nReading reference from fasta file ...")
      full <- Biostrings::readDNAStringSet(ref)
    }else{
      if(methods::is(ref, "character")){
        cat("\nTry to import reference from character vector ...")
        full <- Biostrings::DNAStringSet(ref)
      }else{
        stop("\nUnrecognizable reference format.",
             "\nPlease check your reference input.")
     }
   }
  }
  nams_full <- names(full) # Names are lost in the next step
  if(nchar(paste0(N_up, N_down)) >0){
    full <- Biostrings::DNAStringSet(paste(N_up, full, N_down, sep=""))
    names(full) <- nams_full
  }
  
## Setup temp folder and convert to windows format
  outpath <-  paste0(tempdir(), "/", "seqpac")
  ref_path <-  paste0(tempdir(), "/ref/reference.fa")
  if(grepl("win|WIN|Win", Sys.info()["sysname"])){
      outpath <- gsub("\\", "/", outpath, fixed=TRUE)
      }
  suppressWarnings(dir.create(outpath, recursive = TRUE))
  if(grepl("win|WIN|Win", Sys.info()["sysname"])){
      ref_path <- gsub("\\", "/", ref_path, fixed=TRUE)
      }
  suppressWarnings(dir.create(dirname(ref_path), recursive = TRUE))
  Biostrings::writeXStringSet(full, filepath=ref_path, format="fasta")
  
## Make bowtie index if not available
  # If file input check bowtie index; save results in check_file
  check_file <- FALSE
  if(is.character(ref)){
    if(file.exists(ref)){
       if(!length(list.files(dirname(ref), pattern=".ebwt"))>=2){
         check_file <- TRUE
       }
    }
  }
  if(check_file == FALSE){
    ref_path <- ref 
    cat("\nBowtie indexes found. Will try to use them...")
  }else{
    if(nchar(paste0(N_up, N_down)) >0|
     grepl("DNAString", class(full))|
     check_file == TRUE){
      cat("\nNo bowtie indexes.")
      cat("\nWill try to reindex references ...")
      Rbowtie::bowtie_build(ref_path, outdir=dirname(ref_path), 
                          prefix = "reference", force = TRUE)  
    }
  }
## Make reanno object  
  map_reanno(PAC, ref_paths=list(reference=ref_path), output_path=outpath, 
             type="internal", threads=threads, mismatches=mismatches,  
             import="genome", keep_temp=FALSE, override=override)
  map <- make_reanno(outpath, PAC=PAC, mis_fasta_check = TRUE, output="list")
  stopifnot(length(map$Full_anno$mis0) == 1)

## Reorganize reanno object to a PAC_mapper object
  align_lst <- lapply(map$Full_anno, function(x){
    x <- x[[1]][!is.na(x[[1]]$ref_hits),]
    splt_x <- strsplit(x[[4]], "(?<=;\\+\\||;-\\|)", perl=TRUE)
    names(splt_x) <- x$seq
    lst_align <- lapply(splt_x, function(y){
      y <- gsub("\\|$", "", y)
      temp <- do.call("rbind", strsplit(y, ";"))
      nams <- temp[,1]
      start_align <- as.numeric(gsub("start=", "", temp[,2]))
      strnd_align <- temp[,3]
      strnd_align <- ifelse(strnd_align=="+", "sense", "antisense")
      df <- data.frame(ref_name = nams,
                       ref_strand = "*",
                       align_start = start_align,
                       align_strand = strnd_align)
      })
    df_align <- do.call("rbind", lst_align)
  })
  for(i in 1:length(align_lst)){
    align_lst[[i]]$seqs <- gsub("\\.\\d+", "", rownames(align_lst[[i]]))
    align_lst[[i]]$mismatch <- names(align_lst)[i]
  }
  # Rbind and fix names  
  nam_mis <- paste(paste0(names(align_lst), "."), collapse="|")
  align <- do.call("rbind", align_lst)
  rownames(align) <-  gsub(nam_mis, "", rownames(align)) 
  align_splt <- split(align, align$ref_name)
  align <- lapply(align_splt, function(x){
    rownames(x) <- NULL
    dup_tab <- table(x$seqs)
    nam_multi <- names(dup_tab)[dup_tab > 1]
    # Sequences mapping multiple times are removed or kept
    if(length(nam_multi) > 0){
      if(multi=="remove"){
        warning("\nSome sequences mapped >1 to the same reference.",
                "\nSince multi='remove' these sequences will be removed:",
                immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        x <- x[!x$seqs %in% nam_multi,]
        rownames(x) <- x$seqs
      }
      if(multi=="keep"){
        warning("\nSome sequences mapped >1 to the same reference.",
                "\nSince multi='keep', these sequences will be represented",
                "\nmultiple times in the mapping (psuedoreplication):",
                immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        splt <- split(x, x$seqs)
        splt <- lapply(splt, function(y){
          rownames(y) <- NULL
          return(y)
          })
        x <- do.call("rbind", splt)
      }
    }else{
       rownames(x) <- x$seqs
    }
    ifelse(x$align_strand=="sense", "+", "-")
    df <- data.frame(Mismatch=gsub("mis", "", x$mismatch),
                     Strand= ifelse(x$align_strand=="sense", "+", "-"),
                     Align_start=x$align_start, 
                     Align_end=x$align_start+nchar(x$seqs)-1, 
                     Align_width=nchar(x$seqs))
    rownames(df) <- rownames(x)
    return(df)
  })
  
# Fix bowtie names and match with original reference
  splt_nam <- strsplit(names(full), " ")
  splt_nam <- unlist(lapply(splt_nam, function(x){x[1]}))
  nam_match <- match(splt_nam, names(align))
  align_lst <- align[nam_match]
  stopifnot(identical(names(align_lst)[!is.na(names(align_lst))],  
                      splt_nam[!is.na(names(align_lst))]))
  
# Add full length reference
  names(align_lst)[is.na(names(align_lst))] <- splt_nam[is.na(names(align_lst))]
  fin_lst <- list(NULL)
  for(i in 1:length(align_lst)){
    if(is.null(align_lst[[i]])){
      align_lst[[i]] <- data.frame(Mismatch="no_hits", Strand="no_hits", 
                                   Align_start="no_hits", 
                                   Align_end="no_hits", Align_width="no_hits")
    }
    fin_lst[[i]] <- list(Ref_seq=full[i], Alignments=align_lst[[i]])
    names(fin_lst)[i] <- names(align_lst)[i] 
  }

# Generate alignment string
  doParallel::registerDoParallel(threads) 
  `%dopar%` <- foreach::`%dopar%`
  
  if(report_string==TRUE){
    if(multi=="keep"){
      warning("\nOption multi='keep', is not compatible with report_string=TRUE",
              "\nAlignment string will not be returned.")
    }else{
      ref_lgn <- lapply(fin_lst, function(x){Biostrings::width(x$Ref_seq)})
      ref_lgn  <- max(do.call("c", ref_lgn))
      if(ref_lgn>500){
         warning("\nOption report_string=TRUE is only compatible with",
                 "\nreference < 500 nt. Alignment string will not be returned.")
      }else{
        fin_lst <- lapply(fin_lst, function(x){
          if(x$Alignments[1,1] =="no_hits"){
            x$Alignments <- cbind(x$Alignments, 
                                  data.frame(Align_string="no_hits"))
            x$Alignments <- data.frame(lapply(x$Alignments, as.character), 
                                       stringsAsFactors=FALSE)
            return(x)
          }else{
            ref <- x$Ref_seq
            algn <- x$Alignments 
            n_ref <- nchar(as.character(ref))
            algn_lst <- split(algn, factor(row.names(algn), 
                                           levels=row.names(algn)))
            positions_lst <- foreach::foreach(j=1:length(algn_lst), 
                                              .final=function(y){
                                                names(y) <- names(algn_lst)
                                                return(y)})  %dopar% {
              ref=ref
              n_ref=n_ref
              sq <- rownames(algn_lst[[j]])
              if(algn_lst[[j]]$Strand == "-"){
                 sq <- intToUtf8(rev(utf8ToInt(sq)))
              }
              algn_str <- paste(strrep("-", times=(algn_lst[[j]]$Align_start)-1), 
                                sq, 
                                strrep("-", times= n_ref-(algn_lst[[j]]$Align_end)), 
                                sep="")
            return(algn_str)
          }
          df <- cbind(algn, 
                      data.frame(
                        Align_string=as.character(
                          paste(do.call("c", positions_lst))), 
                        stringsAsFactors=FALSE))
          return(list(Ref_seq=ref, Alignments=df))
        }
      })
     }
    }
  }
  doParallel::stopImplicitCluster()
  class(fin_lst) <- c("list", "seqpac_map")
  return(fin_lst)
}
                          
