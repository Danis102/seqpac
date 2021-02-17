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
#' 
#' #-------------------------------------------------------------------------#
#' ### For coverage plots ###
#' #-------------------------------------------------------------------------#
#' 
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#' 
#' 
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
#' All_plots <- lapply(as.list(Smry_trg_all), function(x){
#'                             PAC_covplot(PAC_filt, map_rRNA, 
#'                             summary_target = list("means_[Groups]", x), 
#'                             xseq=FALSE, style="line", colour="red")})
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
#' ref <- full[grepl("Glu-CTC-3-1|Lys-CTT-1-1", names(full))]
#' 
#' 
#' #PAC_filt <- PAC_filter(PAC_all, threshold=5, coverage=4, type="counts", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#' 
#' map_refs <- PAC_mapper(PAC_filt, ref=ref, threads=8, mismatches=3)
#' 
#' 
#' ref <- "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.fa"
#' 
#' 
#' @export

PAC_mapper <- function(PAC, ref, mapper="reanno", mismatches=0, threads=1, N_up="", N_down="", report_string=FALSE){

############################################  
## Setup
  doParallel::registerDoParallel(threads) 
  `%dopar%` <- foreach::`%dopar%`
  
## Setup reference
  if(file.exists(ref)){
    cat("\nReading reference from file ...")
    full <- Biostrings::readDNAStringSet(ref)
  }else{
    cat("\nUsing input dataframe as reference ...")
    full <- ref
    }  
  nams_full <- names(full) # Names are lost in the next step
  if(nchar(paste0(N_up, N_down)) >0){
    full <- Biostrings::DNAStringSet(paste(N_up, full, N_down, sep=""))
    names(full) <- nams_full
  }
  
## Setup temp folder  
  outpath <- R.utils::filePath(tempdir(), "seqpac", removeUps=FALSE)
  suppressWarnings(dir.create(outpath, recursive = TRUE))  
  ref_path <- R.utils::filePath(outpath, "ref", "reference.fa")
  suppressWarnings(dir.create(dirname(ref_path), recursive = TRUE))
  Biostrings::writeXStringSet(full, filepath=ref_path, format="fasta")
  
  # Make bowtie index if not available
  if(nchar(paste0(N_up, N_down)) >0|!file.exists(ref)|!length(list.files(dirname(ref), pattern=".ebwt"))>=2){
    cat("\nNo bowtie indexes were found existing or extended reference.")
    cat("\nReindexing references ...")
    Rbowtie::bowtie_build(ref_path, outdir=dirname(ref_path), prefix = "reference", force = TRUE)  
  }else{
    ref_path <- ref 
  }
  
############################################  
## PAC mapper using the reanno workflow
  if(mapper=="reanno"){
    # Make reanno object  
    map_reanno(PAC, ref_paths=list(reference=ref_path), output_path=outpath, type="internal", threads=threads, mismatches=mismatches,  import="genome", keep_temp=FALSE)
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
        warning("Sequences were removed since they mapped >1 to the same reference sequence:",  immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        x <- x[!x$seqs %in% nam_multi,]
      }
      df <- data.frame(n_hits=1, Align_start=x$align_start, Align_end=x$align_start+nchar(x$seqs)-1, Align_width=nchar(x$seqs))  
      rownames(x) <- x$seqs
      rownames(df) <- x$seqs
      return(df)
    })
    
    # Fix bowtie names and match with original reference
    splt_nam <- strsplit(names(full), " ")
    splt_nam <- unlist(lapply(splt_nam, function(x){x[1]}))
    nam_match <- match(splt_nam, names(align))
    align_lst <- align[nam_match]
    stopifnot(identical(names(align_lst)[!is.na(names(align_lst))],  splt_nam[!is.na(names(align_lst))]))
    
    ## Add full length reference
    names(align_lst)[is.na(names(align_lst))] <- splt_nam[is.na(names(align_lst))]
    fin_lst <- list(NULL)
    for(i in 1:length(align_lst)){
      if(is.null(align_lst[[i]])){
        align_lst[[i]] <- data.frame(n_hits="no_hits", Align_start="no_hits", Align_end="no_hits", Align_width="no_hits")
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
    cat("Now aligning", length(Anno_frag), "fragments over", length(full), "reference sequences using", threads, "threads (may take a few minutes) ...    ", paste(Sys.time()), "\n")
    len <- length(full)
    
    ## Parallelize sequences
    fin_lst <- list(NA)
    for(i in 1:length(full)){ 
      cat(paste0("\nAligning against:\n ", names(full)[i], "\n Start ", Sys.time()))
      seq_ref <- full[i]
      aligned_lst <- foreach::foreach(t=1:length(query_strings), .packages=c("Biostrings", "stringr"), .final = function(x){names(x) <- paste0(Anno_frag); return(x)}) %dopar% {
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
      if(is.null(target_lst[[2]])){
        target_lst[[2]] <- data.frame(n_hits="no_hits", Align_start="no_hits", Align_end="no_hits", Align_width="no_hits")
        }
      fin_lst[[i]] <- target_lst
      names(fin_lst)[i] <- names(full)[i]
      cat(paste0("\n Done ", Sys.time()))
    }
  }
    
############################################
## Both old and reanno
  if(report_string==TRUE){
      fin_lst <- lapply(fin_lst, function(x){
        if(x$Alignments[1,1] =="no_hits"){
          x$Alignments <- cbind(x$Alignments, data.frame(Align_string="no_hits"))
          x$Alignments <- data.frame(lapply(x$Alignments, as.character), stringsAsFactors=FALSE)
          
          #x$Alignments <- apply(x$Alignments, 2, as.character) 
          return(x)
        }else{
          x$Ref_seq -> ref
          x$Alignments -> algn
          n_ref <- nchar(as.character(ref))
          algn_lst <- split(algn, factor(row.names(algn), levels=row.names(algn)))
          positions_lst <- foreach::foreach(j=1:length(algn_lst), .final=function(y){names(y) <- names(algn_lst);return(y)})  %dopar% {
            ref=ref
            n_ref=n_ref
            algn_str <- paste(strrep("-", times=(algn_lst[[j]]$Align_start)-1), rownames(algn_lst[[j]]), strrep("-", times= n_ref-(algn_lst[[j]]$Align_end)), sep="")
            return(algn_str)
          }
          df <- cbind(algn, data.frame(Align_string=as.character(paste(do.call("c", positions_lst))), stringsAsFactors=FALSE))
          return(list(Ref_seq=ref, Alignments=df))
        }
      })
    }
    
    doParallel::stopImplicitCluster()
    return(fin_lst)
    }
                          