#' Make a count table
#'
#'\code{make_counts} uses parallel processing to generate a count table.
#'
#' Given a paths to fastq this function performs low-level evidence filtering,
#' generates a counts table of sequences passing the filter and plots summary
#' statistics.
#' 
#' @family PAC generation
#'
#' @seealso  \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param input A path to a directory containing input fastq-files. The script
#'   will recursively search this directory for the .fastq|.fastq.gz extension.
#'
#' @param trimming Character indicating what type of trimming tool that should
#'   be used. If \code{trimming="seqpac"}, fastq files will be sent to the
#'   \code{\link{make_trim}} function in seqpac, while if
#'   \code{trimming="cutadapt"} fastq files will be sent to the
#'   \code{\link{make_cutadapt}} function. Note that \code{trimming="seqpac"}
#'   runs independently of external software, while \code{trimming="cutadapt"}
#'   is dependent on an externally installed version of cutadapt and
#'   fastq_quality_filter. Trimmed fastq files are stored temporarily in the
#'   systems default temporary folder. Please, run \code{\link{make_trim}} and
#'   \code{\link{make_cutadapt}} seperately for perminant storage options, or
#'   set \code{save_temp=TRUE} to avoid that \code{make_counts} will delete all
#'   temporary files.  As default trimming=NULL, which indicates that input
#'   fastq files has already been trimmed.
#'     
#' @param parse Character strings defining the command that should be parsed to
#'   \code{\link{make_trim}} or \code{\link{make_cutadapt}}. This will allow
#'   you to customize your trimming according to 3' adaptor sequence and
#'   platform standards etc. Please see examples below and the manuals for
#'   \code{\link{make_trim}} and \code{\link{make_cutadapt}} for more details.
#'   For convenience, \code{parse} also have two default mode for sRNA trimming,
#'   using Illumina and New England Biotype (neb) type small RNA adaptors.
#'   \code{make_counts} will automatically print the exact setting for each
#'   default mode. Briefly, both modes involves polyG (NextSeq/NovaSeq) trimming
#'   and 3' adaptor trimming, with a 0.1 tolerance for mismatch/indels. If
#'   parse="default_illumina", then the "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 3'
#'   adaptor is trimmed and untrimmed sequences are removed. If
#'   parse="default_neb", then "AGATCGGAAGAGCACACGTCTGAACTCCA" is trimmed and
#'   untrimmed sequences are removed. Removing untrimmed sequences is
#'   recommended for sRNA sequencing.
#'   
#' @param evidence Character vector with two inputs named 'experiment' and
#'   'sample' that controls the low-level evidence filter. Users may already at
#'   this point markly reduce the level of noise in the counts table by
#'   specifying the number of independent evidence that a specific sequence must
#'   have to be included. As default,
#'   \code{evidence=c(experiment=2, sample=1)} will include all sequences that
#'   have >=1 count in at least 2 independent fastq files. Thus 'experiment'
#'   controls the number of independent fastq evidence needed across the whole
#'   experiment. Note, however, that 'sample' does not control the number of
#'   counts needed in each sample. The evidence filter will always use >=1 count
#'   in X number of fastq files. Instead 'sample' controls when a sequence
#'   should be included despite not reaching the 'experiment' threshold. Thus if
#'   \code{evidence=c(experiment=2, sample=10)}, sequences that reach 10 counts
#'   in a single sample will also be included. If evidence=NULL all unique
#'   sequences will be saved. See 'examples' below for more examples.
#'   
#' @param plot Logical whether evidence plots should be printed and saved in the
#'   returned list (default=TRUE).
#'   
#' @param threads Integer stating the number of parallell jobs. Note, that
#'   reading multiple fastq files drains memory fast, using up to 10Gb per fastq
#'   file. To avoid crashing the system due to memory shortage, make sure that
#'   each thread on the machine have at least 10 Gb of memory availabe, unless
#'   your fastq files are very small. Use \code{parallel::detectcores()} to see
#'   available threads on the machine.
#'   
#' @param save_temp Logical whether temporary files (including trimmed fastq
#'   files) should be saved or not. Note, the function will print the path to
#'   the temprorary folder in the console.
#'   
#'   
#' @return 
#' A list containing three objects:
#'   1. counts (data frame) = count table. 
#'   2. progress_report = progress report from trimming and evidence filter.
#'   3. evidence_plots = bar graphs showing the impact of evidence filter. 
#'  
#' @examples
#' 
#' \dontrun{
#' 
#'  library(seqpac)
#'   
#' ############################################################ 
#' ### Seqpac fastq trimming with the make_trim function 
#' ### using default settings for NEBNext small RNA adaptor 
#' 
#' input = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' 
#' counts  <- make_counts(input, threads=1, parse="default_neb",
#'                        type="fastq", trimming="seqpac", plot=TRUE,
#'                        evidence=c(experiment=2, sample=1))     
#'     
#'      
#' ############################################################      
#' ### Seqpac trimming using the make_trim function
#' ### using user settings parsed to make_trim
#'  
#' # Make a parse list (see ?make_trim):  
#' parse = list(adapt_3_set=c(type="hard_save", min=10, mismatch=0.1),
#'              adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA",
#'              polyG=c(type="hard_trim", min=10, mismatch=0.1),
#'              seq_range=c(min=14, max=70),
#'              quality=c(threshold=20, percent=0.8))
#'                
#' counts  <-  make_counts(input, threads=1,
#'                         type="fastq", trimming="seqpac",
#'                         parse=parse, 
#'                         evidence=c(experiment=2, sample=1))           
#'   
#'   
#' ############################################################      
#' ### Seqpac trimming using the make_cutadapt function
#' ### (Important: Needs an external installations of cutadapt 
#' ###  and fastq_quality_filter) 
#' 
#'  input <- system.file("extdata", package = "seqpac", mustWork = TRUE)
#'  
#'  Parse for make_cutadapt is a list of 2 character string expressions.
#'  The first is parsed to cutadapt and the other to fastq_quality_filter 
#'  For parallel processes '-j 1' is recommended since seqpac will   
#'  parallelize across samples and not within. Run system("cutadapt -h") and 
#'  system("fastq_quality_filter -h") for more options.
#'  
#'  cut_prse <- paste0("-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT", 
#'                     " --discard-untrimmed --nextseq-trim=20",
#'                     " -O 10 -m 14 -M 70")
#'  
#'  parse = list(cutadapt = cut_prse,
#'               fastq_quality_filter = "-q 20 -p 80")
#'  
#'  counts  <-  make_counts(input, threads=1,
#'                         type="fastq", trimming="cutadapt",
#'                         parse=parse, 
#'                         evidence=c(experiment=2, sample=1))   
#'  
#'  
#'  ## 2 evidence over two indepenent samples, saving single sample 
#'  ## sequences reaching 10 counts
#'  test <- make_counts(input=input,  type="fastq", trimming="seqpac", 
#'                      parse="default_neb",  
#'                      evidence=c(experiment=2, sample=10))
#'  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#'  test$counts[extras==1,]  28 single sample sequences reached 10 counts
#'  
#'  
#'  ## 2 evidence over two indepenent samples, saving single sample 
#'  ## sequences reaching 3 counts 
#'  test <- make_counts(input=input,  type="fastq", trimming="seqpac", 
#'                      parse="default_neb",  
#'                      evidence=c(experiment=2, sample=3))
#'  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#'  test$counts[extras==1,] 1319 single sample sequences reached 3 counts
#'  
#'  }
#'  
#'  
#' @export

make_counts <- function(input, trimming=NULL, threads=1, 
                        plot=TRUE, parse="default_illumina", 
                        evidence=c(experiment=2, sample=1), save_temp=FALSE){
  
  anno <- value <- variable <- NULL
  
  #############################
  ###### fastq as input #######
  cat("Started at ", paste0(Sys.time()), "\n")
  gc(reset=TRUE)
    
  ## Read file system
  if(sum(!dir.exists(input))== length(input)){
    count_files <- input
  }else{
    count_files <- list.files(input, pattern ="fastq.gz\\>|fastq\\>", 
                              full.names=TRUE, recursive=TRUE)
  }
  count_files <- count_files[!grepl("Undetermined_", count_files)]
  count_files_nams <- basename(count_files)
  cat("\nInput type was set to fastq.")
  cat("\nThe following fastq files were found in the path:\n")
  print(count_files)
  
  ## Setup temporary folder
  if(!is.null(trimming)){
    output <- paste0(tempdir(), "/seqpac/")
    if(dir.exists(output)){
      out_fls  <- list.files(output)
      suppressWarnings(file.remove(out_fls)) 
    }
  }else{
    trimming <- "trimmed"
  }
  
  ##########################################################
  ###### Trimming using cutadapt
  if(trimming=="cutadapt"){
    
    cat("\nSeqpac trimming using the make_cutadapt function was specified.")
    cat("\n--- Temporary trimmed fastq are stored in: ", paste0(output))
    cat("\n--- Unless moved, these will be removed when system restarts.")
    cat("\n--- Please run 'make_cutadapt' for more options.")
    if(grepl("default", parse[[1]][1])){
      if(parse[[1]][1]=="default_neb"){
        cut_prse <- paste0("-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCA",
                      " --discard-untrimmed --nextseq-trim=20", 
                      " -O 10 -m 14 -M 70 -e 0.1")
        parse <- list(cutadapt= cut_prse,
                      fastq_quality_filter="-q 20 -p 80")
      }
      if(parse[[1]][1]=="default_illumina"){ 
        
        cut_prse <- paste0("-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
                           " --discard-untrimmed --nextseq-trim=20",
                           " -O 10 -m 14 -M 70 -e 0.1")
        parse <- list(cutadapt=cut_prse,
                      fastq_quality_filter="-q 20 -p 80")
      }
      cat("\n--- Default trimming:\n   ")
      cat("\n---------------------------------------------")
      cat("\ncutadapt command:\n")
      print(parse$cutadapt)
      cat("\n---------------------------------------------")
      cat("\nfastq_quality_filter command:\n")
      print(parse$fastq_quality_filter)
      cat("\n\n")
    }
    
    ## Generate missing make_cutadapt inputs
    trim_input <- list("cutadapt", "fastq_quality_filter")
    names(trim_input) <- unlist(trim_input)
    trim_input <- lapply(trim_input, function(x){
      if(x %in% names(parse)){
        return(parse[x])
      }else{
        return(NULL)
      }
    })
    prog_report <- make_cutadapt(input=input, output=output, 
                                 threads=threads, parse=trim_input)
  }
  
  
  ##########################################################
  ###### Trimming using make_trim()
  if(trimming=="seqpac"){
    # Setup temporary file
    cat("\nSeqpac trimming using the make_trim function was specified.")
    cat("\n--- Temporary trimmed fastq are stored in: ", paste0(output))
    cat("\n--- Unless moved, these will be removed when system restarts.")
    cat("\n--- Please run 'make_trim' separately for more storage options.\n")
    if(grepl("default", parse[[1]][1])){
      if(parse[[1]][1]=="default_neb"){
        parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                      adapt_3_set=c(type="soft_rm", min=10, mismatch=0.1),
                      adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                      seq_range=c(min=14, max=70),
                      quality=c(threshold=20, percent=0.8),
                      indels=TRUE, concat=12, check_mem=FALSE)
      }
      if(parse[[1]][1]=="default_illumina"){
        parse = list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                     adapt_3_set=c(type="soft_rm", min=10, mismatch=0.1), 
                     adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
                     seq_range=c(min=14, max=70),
                     quality=c(threshold=20, percent=0.8),
                     indels=TRUE, concat=12, check_mem=FALSE)
      }
      
      cat("\nDefault trimming:   \n")
      cat("----------------------------------------------")
      cat("\nadapt_3: \n")
      cat(parse$adapt_3, "\n")
      print(parse$adapt_3_set)
      cat("----------------------------------------------")
      cat("\npoly G: \n")
      print(parse$polyG)
      cat("----------------------------------------------")
      cat("\nsize filter:\n")
      print(parse$seq_range)
      cat("----------------------------------------------")
      cat("\nquality filter:\n")
      print(parse$quality)
      cat("----------------------------------------------\n")
    }
    
    ## Generate missing make_trim inputs for both default and parse
    trim_input <- list("polyG", "adapt_3", "adapt_3_set", "adapt_5", 
                       "adapt_5_set", "seq_range", "quality", "indels", 
                       "concat", "check_mem")
    names(trim_input) <- unlist(trim_input)
    trim_input <- lapply(trim_input, function(x){
      if(x %in% names(parse)){
        return(parse[x])
      }else{
        if(x == "indels"){
          return(TRUE)
        }
        if(x == "concat"){
          return(12)
        }
        if(x == "check_mem"){
          return(TRUE)
        }
      }
    })
    prog_report <- make_trim(input=input, output=output, threads=threads,
                             adapt_3=trim_input$adapt_3[[1]],
                             adapt_3_set=trim_input$adapt_3_set[[1]],
                             polyG=trim_input$polyG[[1]],
                             seq_range=trim_input$seq_range[[1]],
                             quality=trim_input$quality[[1]],
                             concat=trim_input$concat[[1]],
                             indels=trim_input$indels[[1]],
                             check_mem=trim_input$check_mem[[1]]
    )
  }
  
  #########################################################
  ###### Read trimmed files and save sample evidence ######
  cat("\n\nIdentifying unique sequences in trimmed fastq files ...")
  if(trimming=="trimmed"){
    fls <- count_files
  }
  if(trimming=="seqpac"){
    fls <- list.files(output, pattern="trim.fastq.gz$", 
                      full.names=TRUE, recursive=FALSE)
    fls <- fls[!grepl("REMOVED", fls)]
  }
  if(trimming=="cutadapt"){
    fls <- list.files(output, pattern="trim.fastq.gz$|trim.fastq$", 
                      full.names=TRUE, recursive=FALSE)
  }
  # Do not use parallel::makeClusters!!!
  `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(threads) 
  seq_lst   <- foreach::foreach(i=1:length(fls), .final = function(x){
    names(x) <- basename(fls); return(x)}) %dopar% {
      
      if(evidence[2]>1){
        fstq <- ShortRead::readFastq(fls[[i]], withIds=FALSE)
        tab <- table(paste0(ShortRead::sread(fstq)))
        extra <- names(tab)[tab >= evidence[2]]
        return(list(uni=names(tab), logi_extra= names(tab) %in% extra))
      }else{
        # best to use ShortRead filters, before unique and table
        fstq <-  paste0(
          ShortRead::sread(
            ShortRead::readFastq(
              fls[[i]], withIds=FALSE, 
              filter=ShortRead::occurrenceFilter(duplicates="head", max=1))))  
        return(list(uni=fstq, logi_extra=NULL))
      }
      rm(fstq)
    }
  
  doParallel::stopImplicitCluster()
  gc(reset=TRUE)
  cat("\nCompiling unique sequences ...")
  seqs <- unlist(lapply(seq_lst, function(x){x$uni}), use.names=FALSE)
  
  ###################################
  ###### Apply evidence filter ######
  # Experiemnt
  if(evidence[1]==2){
    seqs_keep <- unique(as.character(seqs[duplicated(seqs)]))
  }
  if(evidence[1]>2){
    seqs_keep <- unlist(
      lapply(seq_lst, function(x){
        tab <- table(x)
        names(tab)[tab>=evidence]
      }), use.names = TRUE)
  }
  if(evidence[1]<2){
    seqs_keep <- unique(seqs)
  }
  rm(seqs)
  
  # Sample
  if(evidence[2]>1){
    lst_extras <- lapply(seq_lst, function(x){
      extras <- x$uni[x$logi_extra]
      n_extras <- sum(!extras %in% seqs_keep)
      return(list(extras=extras, n_extras=n_extras))
    })
    all_extras <- unique(unlist(lapply(lst_extras, function(x){
      x$extras}), use.names=FALSE))
    seqs_keep <- unique(c(seqs_keep, all_extras))
  }
  
  # Save n uni seqs stats
  n_uniseq <- unlist(lapply(seq_lst, function(x){
    length(x$uni)
  }))
  rm(seq_lst)
  
  ###############################
  ###### Make count table #######
  cat(paste0("\nMaking a count table with sequences appearing in at least ", 
             evidence[1]," independent samples ..."))
  if(evidence[2]>1){
    cat(paste0("\nAlso saves sequences with less evidence but >= ", 
               evidence[2], " counts in a single sample."))
    cat(paste0(
      "\n(Please set evidence$sample == 1 if this was not intended.)"))
  }
  gc(reset=TRUE)  
  doParallel::registerDoParallel(threads) # Do not use parallel::makeClusters!!!
  reads_lst <- foreach::foreach(i=1:length(fls),.final = function(x){
    names(x) <- basename(fls); return(x)}) %dopar% {
      
      stat_mat <- as.data.frame(matrix(NA, nrow=1, ncol=3))
      colnames(stat_mat) <- c("tot_reads", "reads_pass_evidence", 
                              "uniseqs_pass_evidence")
      reads <- ShortRead::readFastq(fls[i], withIds=FALSE)
      reads <- paste0(ShortRead::sread(reads))
      stat_mat$tot_reads <- length(reads)
      reads_keep <- reads[reads %in% seqs_keep]
      rm(reads)
      
      stat_mat$reads_pass_evidence <- length(reads_keep)
      stat_mat$uniseqs_pass_evidence <- length(unique(reads_keep))
      n_reads <- as.data.frame(table(reads_keep))
      rm(reads_keep)
      
      n_reads <- n_reads[match(seqs_keep, n_reads$reads),]
      n_reads$Freq[is.na(n_reads$Freq)] <- 0
      dt <- tibble::tibble(N=n_reads[,2])
      rm(n_reads)
      return(list(counts=dt, stat=stat_mat))
      gc()
    }
  doParallel::stopImplicitCluster()
  gc(reset=TRUE)
  
  cat("\nFinalizing at ", paste0(Sys.time()), "\n")
  names(reads_lst) <- gsub("_merge|\\.fastq.gz$|fastq.gz$|\\.fastq$", 
                           "", names(reads_lst))
  names(reads_lst) <- gsub("temp_|\\.trim", "", names(reads_lst)) 
  names(reads_lst) <- gsub("-", "_", names(reads_lst))
  
  ## Clean up temporary trimming files
  if(!trimming=="trimmed"){
    if(save_temp==FALSE){
      fn <- list.files(output, full.names=TRUE, recursive=FALSE)
      if(any(file.exists(fn))){
        file.remove(fn)
      }
    }
  }
  ordCount_df <- as.data.frame(do.call("cbind", lapply(reads_lst, function(x){
    x$counts})), stringsAsFactors = FALSE)
  colnames(ordCount_df) <- names(reads_lst)
  rownames(ordCount_df) <- seqs_keep
  
  if(trimming=="trimmed"){
    prog_report <- data.frame(
      trimming=as.character(rep(
        "no_trimming_was_done", times=length(reads_lst))), 
      stringsAsFactors=FALSE)
  }
  
  
  ###################################################
  ###### Add evidence stats to progress report ######
  stat_dt <- cbind(data.frame(uni_seqs=n_uniseq, stringsAsFactors =FALSE), 
                   do.call("rbind", lapply(reads_lst, function(x){x$stat})))
  stat_dt$reads_not_pass_evidence <- stat_dt$tot_reads - stat_dt$reads_pass_evidence
  stat_dt$uniseqs_not_pass_evidence <- stat_dt$uni_seqs - stat_dt$uniseqs_pass_evidence
  rownames(stat_dt) <- colnames(ordCount_df)
  stat_dt$sample <- colnames(ordCount_df) 
  if(max(nchar(stat_dt$sample))>15){
    stat_dt$sample <- paste0("fastq_", 1:length(stat_dt$sample))
  }
  df_long <- reshape2::melt(stat_dt, id.vars="sample")
  df_long$value <- as.numeric(df_long$value)
  df_long$sample <- as.factor(df_long$sample)
  main_title <- "Reads passed that evidence"
  
  p_func <- function(df, y_lab, main_title, sub_title){
    df$variable <- factor(df$variable)
    ord <- c(which(grepl("_not_", levels(df$variable))), 
             which(!grepl("_not_", levels(df$variable)))) 
    df$variable <- factor(df$variable, 
                          levels=levels(df$variable)[ord]) 
    ggplot2::ggplot(df, ggplot2::aes(x=sample, y=value, fill=variable)) +
      ggplot2::geom_col(width=1, size=0.2, color="black") + 
      ggplot2::geom_hline(yintercept=0, col="#707177", cex=0.6) +
      ggplot2::geom_text(ggplot2::aes(label=value), 
                         position = ggplot2::position_stack(vjust = 0.5), 
                         angle = 90, color="white", size=4)+
      ggplot2::scale_fill_manual(name = NULL, labels = c("not pass", "pass"), 
                                 values=rev(c("#085B84", "#B20018")))+
      ggplot2::ylab(y_lab)+
      ggplot2::labs(title=main_title, subtitle=)+
      ggplot2::theme(
        plot.caption =  ggplot2::element_text(size=10, face= "bold"),
        axis.title.y = ggplot2::element_text(size=12, face= "bold"),
        axis.title.x = ggplot2::element_text(size=12, face="bold"),
        axis.text = ggplot2::element_text(size=10),
        axis.text.x = ggplot2::element_text(angle=45, hjust=1),
        panel.background = ggplot2::element_blank())
  }
  p1 <- p_func(df_long[df_long$variable %in% c("reads_pass_evidence", 
                                               "reads_not_pass_evidence"),], 
               y_lab="number of reads", 
               main_title=paste0("Sequences present in at least ", 
                                 evidence[1], " samples"),
               sub_title=if(evidence[2]>1){paste0("(unless >", 
                                                  evidence[2], 
                                                  "reads in a single sample)")
                 }else{
                   NULL
                   })
  p2 <- p_func(df_long[df_long$variable %in% c("uniseqs_pass_evidence", 
                                               "uniseqs_not_pass_evidence"),], 
               y_lab="number of unique sequences", 
               main_title=paste0("Sequences present in at least ", 
                                 evidence[1], " samples"),
               sub_title=if(evidence[2]>1){paste0("(unless >", evidence[2], 
                                                  "reads in a single sample)")
                 }else{
                   NULL}
               )
  plt_lst <- list(reads=p1, uniseqs=p2)
  if(plot==TRUE){
    print(cowplot::plot_grid(plotlist=plt_lst, ncol=1, nrow=2))
  }else{
    plt_lst <- "Evidence plot was omitted by user"
  }
  
  prog_report <- cbind(prog_report, 
                       stat_dt[!names(stat_dt) %in% c("uni_seqs", 
                                                      "tot_reads")])
  return(list(counts=ordCount_df, progress_report=prog_report, 
              evidence_plots=plt_lst))
}
