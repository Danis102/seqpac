

context("PAC generation\n")
library(seqpac)

test_that("Testing make_counts, make_trim, make_cutadapt...", {
## Set up environment
  
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
    }
  
  input <- system.file("extdata", package = "seqpac", mustWork = TRUE)
  input <-  list.files(input, patter="fastq.gz\\>", full.names = TRUE)
  smpl <- 1
  while(length(smpl) ==1){
      smpl <- unique(round(runif(2, min = 1, max = length(input)), digits=0))
  }
  input <- input[smpl]
  
  tmp_dir <- paste0(tempdir(), "/seqpac/")
  tmp_files <- list.files(tmp_dir, full.names=TRUE, recursive=FALSE)
  if(any(file.exists(tmp_files))){file.remove(tmp_files)}
    
  parse_cut = list(cutadapt="-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 70",
              fastq_quality_filter="-q 20 -p 80")
  
  parse_seq = list(adapt_3_set=c(type="hard_save", min=10, mismatch=0.1),
               adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA",
               polyG=c(type="hard_trim", min=10, mismatch=0.1),
               seq_range=c(min=14, max=70),
               quality=c(threshold=20, percent=0.8))
  
  if(grepl("unix", .Platform$OS.type)) {
    quiet(
      counts_cut  <-  make_counts(input, threads=3,
                            trimming="cutadapt",
                            parse=parse_cut,
                            evidence=c(experiment=2, sample=1),
                            save_temp = FALSE)
      )
  }
  
  quiet(
    counts_seq  <-  make_counts(input, threads=3,
                           trimming="seqpac",
                           parse=parse_seq, 
                           evidence=c(experiment=2, sample=1),
                           save_temp = TRUE)
  )
  
    trim_files <- paste0(tempdir(), "/seqpac/")
    trim_files <- list.files(trim_files, pattern=".trim.fastq.gz$", full.names=TRUE)
    
  quiet(
    counts_trim  <-  make_counts(trim_files, threads=3,
                           trimming = NULL,
                           evidence=c(experiment=2, sample=1))
  )

## Test counting and trimming
  ## Test both windows and linux
  seq_clss  <- unlist(lapply(lapply(counts_seq, function(x){x[[1]]}), class), use.names=FALSE)
  seq_n  <- c(nrow(counts_seq$counts)>1, 
              ncol(counts_seq$counts) == length(input),
              nrow(counts_seq$progress_report) == length(input))
  trim_clss  <- unlist(lapply(lapply(counts_trim, function(x){x[[1]]}), class), use.names=FALSE)
  trim_n  <- c(nrow(counts_trim$counts)>1, 
              ncol(counts_trim$counts) == length(input),
              nrow(counts_trim$progress_report) == length(input))

  expect_identical(names(counts_seq), c("counts","progress_report", "evidence_plots"))
  expect_identical(names(counts_trim), c("counts","progress_report", "evidence_plots"))
  

  expect_true(sum(seq_n) ==3)
  expect_true(sum(trim_n) ==3)
  

  expect_equal(sum(seq_clss %in% c("numeric", "character", "gg", "ggplot")), 4)
  expect_equal(sum(trim_clss %in% c("numeric", "character", "gg", "ggplot")), 4)


# ## Test make pheno and 
# test_that("Testing make_pheno and make_PAC ...", {
  
  pheno <- as.data.frame(do.call("rbind", strsplit(basename(input), "_|\\."))[,c(1,2,3,4)]) 
  colnames(pheno) <- c("stage", "batch", "index", "sample") 
  pheno$Sample_ID <- apply(pheno, 1, function(x){paste(x, collapse="_")}) 
  
  
  invisible(capture.output(
    pheno_1 <- make_pheno(pheno=pheno, progress_report=counts_seq$progress_report, counts=counts_seq$counts)
  ))
  write.csv(pheno, file=paste0(tmp_dir, "temp.csv"), row.names=FALSE)
  
  invisible(capture.output(
  pheno_2 <- make_pheno(pheno=paste0(tmp_dir, "temp.csv"), progress_report=counts_trim$progress_report, counts=counts_trim$counts)
  ))
  
  pac_seq <- make_PAC(pheno=pheno_1, counts=counts_seq$counts, anno=NULL)
  pac_trim <- make_PAC(pheno=pheno_2, counts=counts_trim$counts, anno=NULL)
  
  expect_true(PAC_check(pac_seq))
  expect_true(PAC_check(pac_trim))
  
   ## Test cutadapt on Linux
  if(grepl("unix", .Platform$OS.type)) {
      cut_clss  <- unlist(lapply(lapply(counts_cut, function(x){x[[1]]}), class), use.names=FALSE)
      cut_n  <- c(nrow(counts_cut$counts)>1, 
                  ncol(counts_cut$counts) == length(input),
                  nrow(counts_cut$progress_report) == length(input))
      expect_identical(names(counts_cut), c("counts","progress_report", "evidence_plots"))
      expect_true(sum(cut_n) ==3)
      expect_equal(sum(cut_clss %in% c("numeric", "character", "gg", "ggplot")), 4)
      pac_cut <- make_PAC(pheno=pheno_1, counts=counts_cut$counts, anno=NULL)
      expect_true(PAC_check(pac_cut))
  }
  
  
})
  
