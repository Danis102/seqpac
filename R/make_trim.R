#' Make trimmed and quality filtered fastq files 
#'
#'\code{make_trim} uses parallel processing to generate trimmed reads.
#'
#' Given an file path to sequence files in fastq format (extension .fastq.gz or
#' .fastq) this function will first trim the sequences from adaptor sequence,
#' and then apply sequence size and phred score quality filters.
#' 
#' @family PAC generation
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current package.
#'
#' @param input Character indicating the path to a directory containing fastq
#'   formated sequence files either compressed (.fastq.gz) or text (.fastq).
#'
#' @param output Character indicating the path to a directory where the
#'   compressed trimmed fastq files should be saved (as basename trim.fastq.gz).
#'
#' @param indels Logical whether indels should be allowed in the mismatch
#'   between read and adaptor sequences. Note that indel trimming will only be
#'   done on fragments that have failed to be trimmed by other means. Default=TRUE.
#'   
#' @param adapt_3_set Character vector with three inputs named 'type', 'min' and
#'   'mismatch' that controls 3' trimming. Options are:
#'      type=
#'       "hard_trim" Trims all upto the very last nucleotide.
#'       "hard_rm"   Same as 'hard_trim' but removes untrimmed sequences.
#'       "hard_save" Same as 'hard_trim' but saves all untrimmed sequences in new fastq file.
#'       "soft_trim" No trimming on the last 3' nucleotides equal to 1/2 of 'min'. 
#'       "soft_rm"   Same as 'soft_trim' but removes untrimmed sequences.
#'       "soft_save" Same as 'soft_trim' but saves all untrimmed sequences in new fastq file.
#'      min=
#'        Integer that controls the short adaptor/poly-G/concatamer trimming and
#'        soft trimming. When min=10, a short version of the adaptor is
#'        generated containing the 10 first nucleotides. The short version of
#'        the adaptor therefore controls the minimum number of overlaps between
#'        adaptor/poly-G. Short adaptor trimming will only occur on untrimmed
#'        reads that have failed trimming with full adaptor. The min option also
#'        controls how many terminal nucleotides to save from trimming when
#'        type is set to 'soft' (1/2 of min).
#'      mismatch=
#'        Numeric controling the percen of mismatch. For instance, if min=10 and
#'        mismatch=0.1, then 1 mismatch is allowed in the minimum overlap.
#'   
#' @param adapt_5_set Same as \code{adapt_5_set} but controls the behavior of 5'
#'   trimming when a sequence is provided in \code{adapt_5}.
#'                     
#' @param polyG Character vector with three inputs named 'type', 'min' and
#'   'mismatch' that controls poly G trimming. This trimming might be necessary
#'   for fastq files from two channel illumina sequencers (e.g. NextSeq and
#'   NovaSeq) where no color signals are interpreted as 'G', but where read
#'   failure sometimes also results in long stratches of 'no-signal-Gs'.
#'
#'   Works similar to \code{adapt_3_set} but instead of an adaptor sequence a
#'   poly 'G' string is constructed from the 'min' input option. Thus if
#'   polyG=c(type="hard_rm", min=10, mismatch=0.1) then sequences containing
#'   'GGGGGGGGGGNNNNN...etc' with 10 percent mismatch will be removed from the
#'   output fastq file.
#'   
#' @param concat Integer setting the threshold for trimmiing concatamere-like
#'   adaptors. Even when adaptor synthesis and ligation are strictly controlled,
#'   concatamer-like ("di-") adaptor sequences are formed. When \code{concat} is
#'   an integer, \code{make_trim} will search all trimmed sequnces for
#'   additional adaptor sequence using a shorter version of the provided adaptor
#'   in \code{adapt_3}. The length of the short adaptor is controlled by
#'   \code{adapt_3_set$min}. If an additional, shorter, adaptor sequence is
#'   found in a read, trimming will only occur if the resulting sequence is >=
#'   \code{concat} shorter in length than the original trimmed sequence. Thus if
#'   \code{concat}=12, a read with a NNNNN-XXXXX-YYYYYYYYY composition, where N
#'   is 'real' nucleotides and X/Y are two independent adaptor sequences, the
#'   trimming will result in NNNNN-XXXXX. If instead \code{concat}=5 then
#'   trimming will result in NNNNN. Note that setting \code{concat} too low will
#'   result in trimming of real nucleotides that just happend to share sequence
#'   with the adaptor. As default \code{concat}=12, which have been carefully
#'   evaluated in relation to the risk of trimming real sequence. If
#'   \code{concat}=NULL, concatamer-like adaptors trimming will not be done.                           
#'                                            
#' @param seq_range Numeric vector with two inputs named 'min' and 'max', that
#'   controls the sequence size filter. For example, if
#'   \code{seq_range=c(min=15, max=50)} the function will extract sequences in
#'   the range between 15-50 nucleotides after trimming. As defualt,
#'   \code{seq_range=c(min=NULL, max=NULL)} and will retain all trimmed sequences.
#'   
#' @param quality Numeric vector with two inputs named 'threshold' and 'percent'
#'   that controls the fastq quality filter. For example, if
#'   \code{quality=c(threshold=20, percent=0.8)} (default) the function will
#'   extract sequences that pass a fastq phred score >=20 in 80% of the
#'   nucleotides after trimming. If \code{quality=c(threshold=NULL,
#'   percent=NULL)} then all sequences will be retained (not default!).     
#'
#' @param threads  Integer stating the number of parallell jobs. reading
#'   multiple fastq files drains memory fast, using up to 10Gb per fastq file.
#'   To avoid crashing the system due to memory shortage, make sure that each
#'   thread on the machine have at least 10 Gb of memory availabe, unless your
#'   fastq files are very small or very large files. Use
#'   \code{parallel::detectcores()} to see available threads on the machine.
#'   Never ever use all threads!
#'
#' @return Exports a compressed trimmed fastq file to output folder (file name:
#'   basename+trim.fastq.gz). If any type="save" has been set, an additional
#'   fastq file (file name includes 'REMOVED') will be exported. In addition,
#'   overview statisics will be returned as a data frame.
#'
#' @examples
#' library(seqpac)
#' 
#' input = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Drosophila_Sep_IOR_190912/Data/Double/Long/Merged_fastq/"
#' input = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Jan_IOR_200130/Data/Single/Merged_fastq/"
#' input = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Drosophila_Sep_IOR_190912/Data/Single/Merged_fastq"  # Large POOH
#' input = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Drosophila_Sep_IOR_190912/Data/Single"
#' 
#' input <- system.file("extdata", package = "seqpac", mustWork = TRUE)
#' output <- "/home/danis31/Desktop/Temp_docs/temp"
#' 
#' prog_report  <-  make_trim(input=input, output=output, threads=6, check_mem=TRUE,
#'                      adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA", 
#'                      polyG=c(type="hard_trim", min=20, mismatch=0.1),
#'                      seq_range=c(min=14, max=70),
#'                      quality=c(threshold=20, percent=0.8))
#'
#' prog_report  <-  make_trim(input=input, output=output, threads=5, 
#'                      adapt_3_set=c(type="hard_save", min=10, mismatch=0.1), adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA", 
#'                      polyG=c(type="hard_trim", min=20, mismatch=0.1),
#'                      seq_range=c(min=14, max=70),
#'                      quality=c(threshold=20, percent=0.8))                                            
#'                      
#' 
#' 
#' adapt_3_set=c(type="hard_trim", min=10, mismatch=0.1)
#' adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA"
#' adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA"
#' polyG=c(type="hard_trim", min=20, mismatch=0.1)
#' seq_range=c(min=12, max=60)
#' quality=c(threshold=20, percent=0.8) 
#' threads=7
#' check_mem=FALSE
#' indels=TRUE
#' concat=12
#' i <- 1
#' 
#' 
#'
#' 
#' @export

make_trim <- function(input, output, indels=TRUE, concat=12, check_mem=TRUE, threads=1,
                      polyG=c(type=NULL, min=NULL, mismatch=NULL),
                      adapt_3_set=c(type="trim", min=10, mismatch=0.1), adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA", 
                      adapt_5_set=c(type=NULL, min=NULL, mismatch=NULL), adapt_5=NULL,
                      seq_range=c(min=NULL, max=NULL),
                      quality=c(threshold=20, percent=0.8)){
                        
############################   
#### General setup #########
  fls <- list.files(input, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
  
  # Memory check
  if(check_mem==TRUE){
      cat("\nChecking trimming memory usage with benchmarkme...")
          mem <- as.numeric(benchmarkme::get_ram()/1000000000)
          mem <- round(mem*0.931, digits=1)
          cor <- parallel::detectCores()
          if(threads > cor*0.8){
            cat(paste0("\n--- Input will use >80% of availble cores (", threads, "/", cor, ")."))
          }
          cat("\n--- Heavy trimming processes needs approx. x10 the fastq file size available in RAM memory.")
          if(threads>length(fls)){
            par_ns <- length(fls)
          }else{
            par_ns <- threads
          }
          worst <- sum(sort(file.info(fls)$size, decreasing = TRUE)[1:par_ns])/1000000000
          worst <- round(worst*12, digits=1)
          cat(paste0("\n--- Worst scenario maximum system burden is estimated to ", worst, " GB of approx. ", mem," GB RAM available."))
          if(worst*1.1 > mem){ 
             warning("\nIn worst scenario trimming may generate an impact close to what the system can stand.\nPlease free more memory per thread if function fails.", immediate. = TRUE)
          }
          if(worst*0.5 > mem){ 
            stop("\nTrimming will generate an impact well above what the system can stand.\nIf you still wish to try, please set check_mem=FALSE.")
          }
          cat(paste0("\n--- Trimming check passed.\n"))
    }

  # Make new output file name
  nam <- gsub("\\.gz$", "", basename(fls))
  nam <- gsub("\\.fastq$|\\.fq$|fastq$", "", nam)
  nam <- gsub("\\.$|\\.tar$", "", nam)
  nam_trim <- paste0(nam, ".trim.fastq.gz")
  # Check for output folder
  out_file <- file.path(output, nam_trim)
  out_dir <- list.files(output, pattern=nam_trim, recursive = FALSE)
  if(length(out_dir)>0){
      stop(paste0("\n  Output trimmed fastq file names are identical to existing files in output:\n  ", out_dir, "\n  Please move or delete the file in the output folder."))
  }
  
  # Make dir
  if(!dir.exists(output)){
    suppressWarnings(dir.create(output, recursive = TRUE))  
  }
  

#########################################################     
##### Enter parallel loops and reading fastq      #######
  cat("\nNow entering the parallel trimming loop (R may stop respond) ...")
  cat(paste0("\n(progress may be followed in: ", output, ")"))
  doParallel::registerDoParallel(threads)  # Do not use parallel::makeClusters!!!
  suppressPackageStartupMessages(require(foreach))
  prog_report <- foreach(i=1:length(fls), .inorder = TRUE, .export= c("fls", "out_file", "nam_trim", "nam"), .final = function(x){names(x) <- basename(fls); return(x)}) %dopar% {
      # Save lists
      sav_lst <- list(NA)
      coord_lst <- list(NA)
      trim_filt <- list(NA)
    
      fstq <- paste0(ShortRead::sread(ShortRead::readFastq(fls[[i]], withIds=FALSE))) 
      num_fact <- as.factor(fstq)
      
      rm(fstq)
      gc(reset=TRUE)
      
      sav_lst$tot_reads <- length(num_fact)
      seqs <- Biostrings::DNAStringSet(levels(num_fact)) # num_fact can be used to apply new trimming coordinates later 
      num_fact <- as.numeric(num_fact)
      lgn <- nchar(paste0(seqs))
      mxlgn <- max(lgn)
      
      if(!length(unique(lgn)) == 1){
        warning("\nDiffering read lengths prior to adapter trimming.\nHave you already performed 3-prim trimming?")
      }

      
#########################################################     
##### Make NextSeq/NovaSeq poly-G trimming/filter ####### 
      if(!is.null(polyG)){
        mn <- as.numeric(polyG[names(polyG)=="min"])
        mis <- as.numeric(polyG[names(polyG)=="mismatch"])
        tp <- polyG[names(polyG)=="type"]
        ns <- mxlgn-mn
        shrt_ns <- paste0(paste(rep("G", mn), collapse=""), paste(rep("N", ns), collapse=""))
        adapt_mis_corr <- mis*(mn/mxlgn)  # Correction for N length extension
        trim_seqs <- Biostrings::trimLRPatterns(subject=seqs, ranges=TRUE, Rfixed=FALSE,  Rpattern=shrt_ns, max.Rmismatch=adapt_mis_corr)
        coord_lst$polyG <- tibble::tibble(Biostrings::end(trim_seqs))
        sav_lst$polyG <-  c(trimmed=sum(!coord_lst$polyG==lgn), polyG)
        rm(trim_seqs)
        # Type
        if(grepl("soft", tp)){
          end_vect <- mn*0.5-1
          end_vect <- lgn - end_vect
          logi_update <- coord_lst$polyG[[1]] >= end_vect
          coord_lst$polyG[[1]][logi_update] <- lgn[logi_update]
          rm(logi_update, end_vect)
        }
        if(grepl("rm|save", tp)){
        trim_filt$polyG  <- coord_lst$polyG[[1]] == lgn # Not polyG
        }
      }
      gc(reset=TRUE)
      
#####################################
##### Make adaptor trimming 3' ######
    if(!is.null(adapt_3)){
    ## setup
      mn <- as.numeric(adapt_3_set[names(adapt_3_set)=="min"])
      mis <- as.numeric(adapt_3_set[names(adapt_3_set)=="mismatch"])
      tp <- adapt_3_set[names(adapt_3_set)=="type"]
    ## Setup adapter sequences
      adapt_lgn <- nchar(adapt_3)
      ns <- mxlgn-adapt_lgn
      adapt_ns <- paste0(adapt_3, paste(rep("N", ns), collapse=""))
      adapt_shrt <- paste(c(substr(adapt_3, 1, mn), rep("N", mxlgn-mn)), collapse="")

    ## Trim perfect
      trim <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs, ranges=TRUE, Rfixed=FALSE, with.Rindels = FALSE, Rpattern=adapt_ns, max.Rmismatch=0))
      
    ## Trim mis
      trim_mis <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs, ranges=TRUE, Rfixed=FALSE, with.Rindels = FALSE, Rpattern=adapt_ns, max.Rmismatch=mis))
      logi_diff <- (trim - trim_mis) >=5 # Only update larger differences will prevent nugging
      trim[logi_diff] <- trim_mis[logi_diff]
      rm(logi_diff, trim_mis)

    ## Trimming wobbling starts 
      if(mis>=0.1){
          if(nchar(adapt_3)< 30){
                 warning("Your trimming sequence was shorter than 30 nt.\nFor best trimming results in the terminals using \nmismatches, please provide a longer adaptor sequence.")
              }
      }
      shrtst_lgn <- 15
      sb <- nchar(adapt_3)-1
      adapt_lgn <- nchar(adapt_3)
      start_gone <- list(NULL)
      while(sb >= shrtst_lgn){
              strt <- adapt_lgn-sb+1
              start_gone[[sb]] <- paste0(c(substr(adapt_3, strt, nchar(adapt_3)), rep("N", mxlgn-(adapt_lgn-strt+1))), collapse="")
              sb <- sb-1
      }
      start_gone <- unlist(start_gone[!unlist(lapply(start_gone, is.null))])
      trim_mt <- list(NULL)
      for(z in 1:length(start_gone)){
                  shrt_ns <- start_gone[z]
                  n_Ns <- stringr::str_count(shrt_ns, "N")
                  adapt_mis_corr <- (mis*(nchar(shrt_ns)-n_Ns))/nchar(shrt_ns)
                  trim_seqs <- Biostrings::trimLRPatterns(subject=seqs, ranges=TRUE, Rfixed=FALSE, with.Rindels = FALSE, Rpattern=start_gone[z], max.Rmismatch=adapt_mis_corr)
                  trim_mt[[z]] <- tibble::tibble(Biostrings::end(trim_seqs))
      }
      trim_mt <- do.call("cbind", trim_mt)
      trim_mt <- unlist(apply(trim_mt, 1, function(x){min(unique(x))<=1}))
      trim[trim_mt] <- 0
      rm(start_gone, trim_mt)

    ## Trimming indels
      if(indels==TRUE){
         logi_long <- trim >= lgn-5
         seqs_long <-  seqs[logi_long] # Apply trimming only on long (untrimmed)
         trim_seqs_indel <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs_long, ranges=TRUE, Rfixed=FALSE, with.Rindels = TRUE, Rpattern=adapt_shrt, max.Rmismatch=0))
         trim_seqs_indel_2 <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs_long, ranges=TRUE, Rfixed=FALSE, with.Rindels = TRUE, Rpattern=adapt_ns, max.Rmismatch=2/nchar(adapt_ns)))
         trim[logi_long] <- ifelse((trim_seqs_indel - trim_seqs_indel_2) > 5, trim_seqs_indel_2, trim_seqs_indel)
         rm(trim_seqs_indel, trim_seqs_indel_2, logi_long, seqs_long) 

      }

    ## Trim concatemer-like adaptors
      if(!is.na(concat)){
         n_Ns <- stringr::str_count(adapt_shrt, "N")
         adapt_mis_corr <- (mis*(nchar(adapt_shrt)-n_Ns))/nchar(adapt_shrt)
         trim_shrt <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs, ranges=TRUE, Rfixed=FALSE, with.Rindels = FALSE, Rpattern=adapt_shrt, max.Rmismatch=adapt_mis_corr))
         logi_dim <- (trim - trim_shrt) >= concat # Only update if more/equal to concat

           chck <- substr(seqs[logi_dim], trim_shrt[logi_dim]+1, trim[logi_dim]) # Check that adaptor is present inbetween
           n_Ns <- stringr::str_count(adapt_ns, "N")
           adapt_mis_corr <- (mis*(nchar(adapt_ns)-n_Ns))/nchar(adapt_ns)
           chck <- agrepl(gsub("N", "", adapt_shrt), chck, max.distance = 0.1, fixed=TRUE)

         logi_dim[logi_dim] <- chck # Update only where check was confirmed
         trim[logi_dim] <- trim_shrt[logi_dim]  # Update trim
         rm(logi_dim, trim_shrt, chck) 
      }
    
    ## Subdivide short adaptor
      logi_long <- trim >= lgn-5
      seqs_long <-  seqs[logi_long]
      adapt_shrt1 <- paste(c(substr(adapt_3, 1, 10), rep("N", mxlgn-nchar(substr(adapt_3, 1, 10)))), collapse="")
      trim_catch1 <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs_long, ranges=TRUE, Rfixed=FALSE, with.Rindels = TRUE, Rpattern=adapt_shrt, max.Rmismatch=2/nchar(adapt_shrt)))
      
      adapt_shrt2 <- paste(c(substr(adapt_3, 11, 20), rep("N", mxlgn-nchar(substr(adapt_3, 11, 20)))), collapse="")
      trim_catch2 <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs_long, ranges=TRUE, Rfixed=FALSE, with.Rindels = TRUE, Rpattern=adapt_shrt2, max.Rmismatch=2/nchar(adapt_shrt2)))
 
      adapt_shrt3 <- paste(c(substr(adapt_3, 21, 30), rep("N", mxlgn-nchar(substr(adapt_3, 21, 30)))), collapse="")
      trim_catch3 <- Biostrings::end(Biostrings::trimLRPatterns(subject=seqs_long, ranges=TRUE, Rfixed=FALSE, with.Rindels = TRUE, Rpattern=adapt_shrt3, max.Rmismatch=2/nchar(adapt_shrt3)))
      
      logi_1 <- (trim_catch2-trim_catch1) %in% (10-3):(10+1)
      logi_2 <- (trim_catch3-trim_catch1) %in% (20-3):(20+1)
      logi_3 <- (trim_catch3-trim_catch2) %in% (10-3):(10+1)

      trim_catch <- rep(75, times=length(seqs_long)) 
      trim_catch[logi_1] <- trim_catch1[logi_1]
      trim_catch[logi_2] <- trim_catch1[logi_2]
      trim_catch[logi_3] <- (trim_catch2-9)[logi_3]
      
      trim[logi_long] <-  trim_catch
      rm(logi_long, seqs_long, trim_catch, trim_catch1, trim_catch2, trim_catch3, logi_1, logi_2, logi_3)
      
    ## Correct for biostring mis-start
      sb_read  <- substr(seqs, trim-2, trim+6) 
      corrct <- stringr::str_locate(sb_read, substr(adapt_3, 1, 5))[,"start"]
      corrct2 <- stringr::str_locate(sb_read, substr(adapt_3, 1, 3))[,"start"]
      corrct <- corrct-4
      corrct2 <- corrct2-4
      corrct[is.na(corrct)] <- corrct2[is.na(corrct)] # Only update NA
      corrct[is.na(corrct)] <- 0 # Remaining NAs becomes 0
      trim <- (trim + corrct) # Add adjusted tims sites
      rm(corrct, corrct2)
      
    ## Add to coordinate list  
      coord_lst$trim_3 <- tibble::tibble(trim) 
      rm(trim)
      gc(reset=TRUE)

    ## Correct for type
      if(grepl("soft", tp)){
          end_vect <- mn*0.5
          end_vect <- lgn - end_vect
          logi_update <- coord_lst$trim_3[[1]] >= end_vect
          coord_lst$trim_3[[1]][logi_update] <- lgn[logi_update]
          rm(logi_update, end_vect)
        }
        if(grepl("rm|save", tp)){
        trim_filt$trim_3  <- coord_lst$trim_3[[1]] == lgn # Not trim 3
        }
    }
      sav_lst$trim_3 <-  c(trimmed=sum(!trim_filt$trim_3[[1]]), adapt_3_set, adapt=adapt_3)
      
      rm(lgn)
      gc(reset=TRUE)
      
################################################
##### Apply trimming coordinates on fastq ######
   ## Apply on sequences
      coord_lst <- coord_lst[!is.na(coord_lst)]
      sav_lst <- sav_lst[!is.na(sav_lst)]
      coord_lst <- apply(tibble::as_tibble(do.call("cbind", coord_lst)), 1, min)
      trim_fastq <- substring(seqs, 1, coord_lst)
      trim_fastq <- trim_fastq[num_fact]  # Creates the same order as original fastq
      
      fstq <- ShortRead::readFastq(fls[[i]]) #
      fstq@sread <- Biostrings::DNAStringSet(trim_fastq)#
      
   ## Apply on quality
      if(!is.null(quality)){
          qual <- paste(Biostrings::quality(fstq@quality)) 
          qual <- substr(qual, start=1, stop=coord_lst[num_fact])
          fstq@quality <-  ShortRead::FastqQuality(Biostrings::BStringSet(qual)) 
          }
      rm(trim_fastq, coord_lst)
      gc(reset=TRUE)
      
########################################################  
##### Filter non-trimmed depending type           ######
     trim_filt <- trim_filt[!is.na(trim_filt)]
     trim_filt <- lapply(trim_filt, function(x){x[num_fact]})
      
     if(length(trim_filt)>0){
         save_logi <- unlist(lapply(sav_lst, function(x){grepl("save", x["type"])}))
         rm_logi <- unlist(lapply(sav_lst, function(x){grepl("rm", x["type"])}))
         sav_nam <- names(save_logi)[save_logi]
         rm_nam <- names(rm_logi)[rm_logi]
         if(length(sav_nam)>0){
           save_logi <- rowSums(do.call("cbind", trim_filt[sav_nam])) > 0
           sav_nam <- paste0(sav_nam, collapse="_")
           fil_nam_un <-  gsub("\\.trim.fastq.gz$", paste0(".REMOVED_", sav_nam, ".fastq.gz"),  out_file[i])
           ShortRead::writeFastq(fstq[save_logi], fil_nam_un, mode="w", full=FALSE, compress=TRUE)
         }
         if(length(rm_nam)>0){
           rm_logi <- rowSums(do.call("cbind", trim_filt[rm_nam])) > 0
           fstq <- fstq[!rowSums(cbind(rm_logi, save_logi)) > 0]
         }else{
           fstq <- fstq[!save_logi]
         }
     }
     rm(trim_filt, save_logi, rm_logi, num_fact) 
     gc(reset=TRUE)
     
########################
##### Size filter ######
   if(!is.null(seq_range)){
        logi_min <- Biostrings::width(fstq@sread) >= seq_range["min"]
        logi_max <- Biostrings::width(fstq@sread) <= seq_range["max"]
        sav_lst$size <- c(too_short=sum(!logi_min), too_long=sum(!logi_max), seq_range)
        fstq <- fstq[logi_min+logi_max==2]
        rm(logi_min, logi_max)
        gc(reset=TRUE)
   }

##########################################
##### Phred score quality filtering ######
     # First extract trimmed qualities
     if(!is.null(quality)){
         enc <- Biostrings::encoding(fstq@quality) # Extract phred score fastq translations
         enc_srch <- as.factor(names(enc))[paste0(enc) %in% as.character(1:quality["threshold"]-1)]
         enc_srch <- paste0("[", paste(enc_srch, collapse="") , "]")
         # Apply filter if specified
         qual_logi <- stringr::str_count(paste0(Biostrings::quality(fstq@quality)), enc_srch)
         qual_logi <-  1-(qual_logi/Biostrings::width(fstq@quality)) >= quality["percent"]
         qual_logi[is.na(qual_logi)] <- FALSE
         sav_lst$quality <-  c(removed=sum(!qual_logi), quality)
         fstq <- fstq[qual_logi]
         rm(qual_logi)
         gc(reset=TRUE)
     }
     
########################################################  
##### Save fastq and create progress reportfile(s) #####
     sav_lst$out_reads <- length(fstq)
     ShortRead::writeFastq(fstq, out_file[i], mode="w", full=FALSE, compress=TRUE)
     return(sav_lst)
       }
doParallel::stopImplicitCluster()
detach(package:foreach)
gc(reset=TRUE)

nams_prog <- as.list(names(prog_report[[1]]))
names(nams_prog) <- unlist(nams_prog)
 
prog_report <-  lapply(nams_prog, function(x){
    lst_type <- list(NULL)
    for(i in 1:length(prog_report)){
     lst_type[[i]] <- cbind(data.frame(file=names(prog_report)[i]), t(data.frame(prog_report[[i]][x])))
    }
    do.call("rbind", lst_type)
})

report_fin <- prog_report$tot_reads
rownames(report_fin) <- NULL
colnames(report_fin)[2] <- "input_reads"

for(i in 1:length(prog_report)){
  rprt <- prog_report[[i]]
  rprt_nam <- names(prog_report)[i]
  
  if(rprt_nam=="polyG"){
    trim <- as.numeric(as.character(rprt$trimmed))
    perc <- paste0("(", round(trim/report_fin$input_reads, digits=4)*100, "%)")
    set <- paste0(rprt$type,"|min", rprt$min,"|mis", rprt$mismatch)
    report_fin <-  cbind(report_fin, data.frame(polyG_set=set, polyG= paste(trim, perc)))
  }
  if(rprt_nam=="trim_3"){
    trim <- as.numeric(as.character(rprt$trimmed))
    perc <- paste0("(", round(trim/report_fin$input_reads, digits=4)*100, "%)")
    set <- paste0(rprt$type,"|min", rprt$min,"|mis", rprt$mismatch)
    report_fin <-  cbind(report_fin, data.frame(trim3_set=set, trim3= paste(trim, perc)))
  }
  if(rprt_nam=="trim_5"){
    trim <- as.numeric(as.character(rprt$trimmed))
    perc <- paste0("(", round(trim/report_fin$input_reads, digits=4)*100, "%)")
    set <- paste0(rprt$type,"|min", rprt$min,"|mis", rprt$mismatch)
    report_fin <-  cbind(report_fin, data.frame(trim5_set=set, trim5= paste(trim, perc)))
  }
  if(rprt_nam=="size"){
    shrt <- as.numeric(as.character(rprt$too_short))
    lng <- as.numeric(as.character(rprt$too_long))
    tot_size <- shrt+lng
    perc <- paste0("(", round(tot_size/report_fin$input_reads, digits=4)*100, "%)")
    set <- paste0("min", rprt$min, "|max", rprt$max)
    report_fin <-  cbind(report_fin, data.frame(size_set=set, 'size_short_long'= paste0(shrt, "/", lng, " ", perc)))
  }
  if(rprt_nam=="quality"){
    rmvd <- as.numeric(as.character(rprt$removed))
    perc <- paste0("(", round(rmvd/report_fin$input_reads, digits=4)*100, "%)")
    set <- paste0("thresh", rprt$threshold, "|perc", rprt$percent)
    report_fin <-  cbind(report_fin, data.frame(quality_set=set, 'quality_removed'= paste0(rmvd, " ", perc)))
  }
  if(rprt_nam=="out_reads"){
    out <- as.numeric(as.character(rprt[,2]))
    perc <- paste0("(", round(out/report_fin$input_reads, digits=4)*100, "%)")
    report_fin <-  cbind(report_fin, data.frame(output_reads= paste(out, perc)))
  }
}

return(report_fin)         
}
