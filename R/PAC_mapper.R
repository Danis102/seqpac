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
#' @param mapper Character indicating what align engine to be used. If
#'   mapper="reanno" (default) PAC_mapper will create temporary files and
#'   generate missing bowtie indexes for a reannotion workflow (see
#'   \code{?map_reanno}. Important, for PAC_mapper to successfully tap into the
#'   reannotation workflow, the names of ref (fasta or DNAStringSet) cannot
#'   contain any ";" characters. If mapper="vmatch", PAC_mapper will instead use
#'   the \code{vmatchPattern} function in the Biostrings package. While the
#'   mapper="reanno" is by far the fastest option, mapper="vmatch" does not
#'   require bowtie indexes.
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
#' @return Stacked list, where each object on the highest level contains:
#'                    (Object 1) Reference name and sequence. 
#'                    (Object 2) Dataframe showing the mapping results of
#'                               each quiary sequence that mapped to Object 1.
#'
#' @examples
#' 
#' 
#' # # More details on the examples can be found in the vignette.
#' # 
#' # library(seqpac)
#' # load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' # 
#' # ###########################################################
#' # ### tRNA analysis in seqpac 
#' # ##----------------------------------------
#' # 
#' # # First create an annotation blanc PAC with group means
#' # pac$Anno <- pac$Anno[,1, drop=FALSE]
#' # pac_trna <- PAC_summary(pac, norm = "cpm", type = "means", pheno_target=list("stage"), merge_pac = TRUE)
#' # 
#' # # Then reannotate only tRNA using the PAC_mapper function
#' # ref <- "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.fa"
#' # map_object <- PAC_mapper(pac_trna, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=TRUE)
#' # 
#' # 
#' # ###########################################################
#' # ## Coverage plot of tRNA using PAC_covplot
#' # 
#' # # Single tRNA targeting a summary dataframe 
#' # PAC_covplot(pac_trna, map=map_object, summary_target= list("cpmMeans_stage"), map_target="tRNA-Ala-AGC-1-1_chr3R:17657145-17657217_(+)")
#' # 
#' # # Find tRNAs with many fragments
#' # n_tRFs <- unlist(lapply(map_object, function(x){nrow(x[[2]])}))
#' # selct <- (names(map_object)[n_tRFs>1])[c(1, 16, 25, 43)]
#' # cov_plt <- PAC_covplot(pac_trna, map=map_object, summary_target= list("cpmMeans_stage"), map_target=selct)
#' # cowplot::plot_grid(plotlist=cov_plt, nrow=2, ncol=2)
#' # 
#' # 
#' # ###########################################################
#' # ## Analyze range types with map_rangetype and PAC_trna functions
#' # 
#' # # Download ss object from GtRNAdb 
#' # dest_path <- file.path(
#' #              "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.tar.gz")
#' #
#' # download.file(
#' #    url="http://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.tar.gz", 
#' #    destfile=dest_path)
#' # untar(dest_path, exdir= dirname(dest_path), 
#' #    files = "dm6-tRNAs-confidence-set.ss")
#' # ss_file <- 
#' # "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/dm6-tRNAs-confidence-set.ss"
#' # 
#' # # Classify fragments according to loop cleavage (small loops are omitted)       
#' # map_object_ss <- map_rangetype(map_object, type="ss", 
#' #                                ss=ss_file, min_loop_width=4)# Gives warning         
#' # 
#' # # Remove reference tRNAs with no hits
#' # map_object_ss <-  map_object_ss[!unlist(lapply(map_object_ss, 
#'                                   function(x){x[[2]][1,1] == "no_hits"}))]
#' # map_object_ss[[2]]
#' # 
#' # 
#' # ###########################################################
#' # # Function classifying 5'-tRF, 5'halves, i-tRF, 3'-tRF, 3'halves
#' # 
#' # # Set tolerance for classification as a terminal tRF
#' # tolerance <- 5  # 2 nucleotides from start or end of full-length tRNA)
#' # 
#' # # Apply the tRNA_class function and make a tRNA type column
#' # pac_trna <- tRNA_class(pac_trna, map=map_object_ss, terminal=tolerance)
#' # pac_trna$Anno$type <- paste0(pac_trna$Anno$decoder, pac_trna$Anno$acceptor)
#' # head(pac_trna$Anno)
#' # 
#' # # Now use PAC_trna to generate some graphs based on grand means
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
#' # # By setting join = FALSE you will get group means
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
#' @export

PAC_mapper <- function(PAC, ref, mapper="reanno", mismatches=0, 
                       threads=1, N_up="", N_down="", report_string=FALSE){

############################################  
## Setup
  j <- NULL 
  ## Setup reference  
  if(class(ref)=="DNAStringSet"){
    cat("\nImporting reference from DNAStringSet ...")
    full <- ref
  }else{
    if(file.exists(ref)){
      cat("\nReading reference from fasta file ...")
      full <- Biostrings::readDNAStringSet(ref)
    }else{
      if(class(ref)=="character"){
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
  
## Setup temp folder  
  outpath <- R.utils::filePath(tempdir(), "seqpac", removeUps=FALSE)
  # Convert to windows format
  if(grepl("win|WIN|Win", Sys.info()["sysname"])){
      outpath <- gsub("\\", "/", outpath, fixed=TRUE)
      }
  suppressWarnings(dir.create(outpath, recursive = TRUE))  
  ref_path <- R.utils::filePath(outpath, "ref", "reference.fa")
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
  if(nchar(paste0(N_up, N_down)) >0|
     grepl("DNAString", class(full))|
     check_file == TRUE){
    cat("\nNo bowtie indexes were found existing or extended reference.")
    cat("\nReindexing references ...")
    Rbowtie::bowtie_build(ref_path, outdir=dirname(ref_path), 
                          prefix = "reference", force = TRUE)  
  }else{
    ref_path <- ref 
  }
  
############################################  
## PAC mapper using the reanno workflow
  if(mapper=="reanno"){
    # Make reanno object  
    map_reanno(PAC, ref_paths=list(reference=ref_path), output_path=outpath, 
               type="internal", threads=threads, mismatches=mismatches,  
               import="genome", keep_temp=FALSE)
    map <- make_reanno(outpath, PAC=PAC, mis_fasta_check = TRUE)
    stopifnot(length(map$Full_anno$mis0) == 1)
    # Reorganize reanno object to a PAC_mapper object
    tRNA_test <- 0 
    align <- lapply(map$Full_anno, function(x){
      x <- x[[1]][!is.na(x[[1]]$ref_hits),]
      splt_x <- strsplit(x[[4]], "(?<=;\\+\\||;-\\|)", perl=TRUE)
      names(splt_x) <- x$seq
      lst_align <- lapply(splt_x, function(y){
        y <- gsub("\\|$", "", y)
        temp <- do.call("rbind", strsplit(y, ";"))
        nams <- temp[,1]
        type_test1 <- grepl("tRNA", nams)
        type_test2 <- grepl("_\\(\\+\\)|_\\(-\\)", nams)
        if(sum(type_test1, type_test2) == 2 ){
          tRNA_test <- tRNA_test+1
          strnd_ref <- do.call("rbind", strsplit(nams, "\\d_\\("))[,2]
          strnd_ref <- gsub("\\)", "", strnd_ref)
          start_align <- as.numeric(gsub("start=", "", temp[,2]))
          strnd_align <- temp[,3]
          strnd_align <- ifelse(strnd_align=="+", "sense", "antisense")
          df <- data.frame(ref_name = nams,
                           ref_strand = strnd_ref,
                           align_start = start_align,
                           align_strand = strnd_align)
        }else{
          start_align <- as.numeric(gsub("start=", "", temp[,2]))
          strnd_align <- temp[,3]
          strnd_align <- ifelse(strnd_align=="+", "sense", "antisense")
          df <- data.frame(ref_name = nams,
                           ref_strand = "*",
                           align_start = start_align,
                           align_strand = strnd_align)
          }          
        })
      df_align <- do.call("rbind", lst_align)
    })
    for(i in 1:length(align)){
      align[[i]]$seqs <- gsub("\\.\\d+", "", rownames(align[[i]]))
      align[[i]]$mismatch <- names(align)[i]
    }
    align <- do.call("rbind", align)
    align_splt <- split(align, align$ref_name)
    align <- lapply(align_splt, function(x){
      dup_tab <- table(x$seqs)
      nam_multi <- names(dup_tab)[dup_tab > 1]
      ## Sequences mapping multiple times are removed
      if(length(nam_multi) > 0){
        warning("\nSequences were removed since they mapped >1",
                ",\nto the same reference sequence:",  immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        x <- x[!x$seqs %in% nam_multi,]
      }
      df <- data.frame(n_hits=1, Align_start=x$align_start, 
                       Align_end=x$align_start+nchar(x$seqs)-1, 
                       Align_width=nchar(x$seqs))  
      rownames(x) <- x$seqs
      rownames(df) <- x$seqs
      return(df)
    })
    
    # Fix bowtie names and match with original reference
    splt_nam <- strsplit(names(full), " ")
    splt_nam <- unlist(lapply(splt_nam, function(x){x[1]}))
    nam_match <- match(splt_nam, names(align))
    align_lst <- align[nam_match]
    stopifnot(identical(names(align_lst)[!is.na(names(align_lst))],  
                        splt_nam[!is.na(names(align_lst))]))
    
    ## Add full length reference
    names(align_lst)[is.na(names(align_lst))] <- splt_nam[is.na(names(align_lst))]
    fin_lst <- list(NULL)
    for(i in 1:length(align_lst)){
      if(is.null(align_lst[[i]])){
        align_lst[[i]] <- data.frame(n_hits="no_hits", Align_start="no_hits", 
                                     Align_end="no_hits", Align_width="no_hits")
      }
      fin_lst[[i]] <- list(Ref_seq=full[i], Alignments=align_lst[[i]])
      names(fin_lst)[i] <- names(align_lst)[i] 
    }
  }
      
############################################
## Old PAC_mapper:
  if(mapper=="old"){
    # Setup
    Anno_frag  <- Biostrings::DNAStringSet(rownames(PAC$Anno))
    query_strings <- as.list(as.character(rownames(PAC$Anno)))
    
    ## Aligning using parallel processing
    cat("Now aligning", length(Anno_frag), "fragments over", length(full), 
        "reference sequences using", threads, 
        "threads (may take a few minutes) ...    ", paste(Sys.time()), "\n")
    len <- length(full)
    
    ## Parallelize sequences
    fin_lst <- list(NA)
    for(i in 1:length(full)){ 
      cat(paste0("\nAligning against:\n ", names(full)[i], "\n Start ", 
                 Sys.time()))
      seq_ref <- full[i]
      
      aligned_lst <- foreach::foreach(
        t=1:length(query_strings), 
        .packages=c("Biostrings", "stringr"), 
        .final = function(x){
          names(x)<- paste0(Anno_frag); return(x)}) %dopar% {
            
        y <- as.data.frame(
          Biostrings::vmatchPattern(query_strings[[t]],  
                                    seq_ref, max.mismatch=mismatches, 
                                    fixed=FALSE))
        
        if(!nrow(y)< 1){
          y <- y[,c(1,3:5)]
          y$group <- nrow(y)
          colnames(y) <- c("n_hits", "Align_start", "Align_end", "Align_width")
          return(y)
        }else{
          return(NULL)}
      }
      aligned   <- do.call("rbind", 
                           aligned_lst[unlist(lapply(aligned_lst, function(x){
                             !is.null(x)}))])                  
      target_lst  <- list(Ref_seq=seq_ref, Alignments=aligned)
      if(is.null(target_lst[[2]])){
        target_lst[[2]] <- data.frame(n_hits="no_hits", Align_start="no_hits", 
                                      Align_end="no_hits", 
                                      Align_width="no_hits")
        }
      fin_lst[[i]] <- target_lst
      names(fin_lst)[i] <- names(full)[i]
      cat(paste0("\n Done ", Sys.time()))
    }
  }
    
############################################
## Both old and reanno
  doParallel::registerDoParallel(threads) 
  `%dopar%` <- foreach::`%dopar%`
  if(report_string==TRUE){
      fin_lst <- lapply(fin_lst, function(x){
        if(x$Alignments[1,1] =="no_hits"){
          x$Alignments <- cbind(x$Alignments, 
                                data.frame(Align_string="no_hits"))
          x$Alignments <- data.frame(lapply(x$Alignments, as.character), 
                                     stringsAsFactors=FALSE)
          
          #x$Alignments <- apply(x$Alignments, 2, as.character) 
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
            algn_str <- paste(strrep("-", times=(algn_lst[[j]]$Align_start)-1), 
                              rownames(algn_lst[[j]]), 
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
    
    doParallel::stopImplicitCluster()
    return(fin_lst)
    }
                          