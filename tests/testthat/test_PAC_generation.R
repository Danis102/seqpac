

context("PAC generation\n")
library(seqpac)

test_that("Testing make_counts, make_trim, make_cutadapt...", {
## Set up environment
  
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
    }
  
  sys_path = system.file("extdata", package = "seqpac", mustWork = TRUE)
  fq <- list.files(path = sys_path, pattern = "fastq", all.files = FALSE,
                   full.names = TRUE)
  

  
  sampler <- ShortRead::FastqSampler(fq, 20000)
  set.seed(123)
  fqs <- list(fq1=ShortRead::yield(sampler),
              fq2=ShortRead::yield(sampler),
              fq3=ShortRead::yield(sampler))

  # Now generate a temp folder were we can store the fastq files
  
  input <- paste0(tempdir(), "/seqpac_temp/")
  dir.create(input, showWarnings=FALSE)
  
  # Clean up temp folder
  closeAllConnections()
  out_fls  <- list.files(input, recursive=TRUE, 
                         full.names = TRUE)
  suppressWarnings(file.remove(out_fls))
  
  # And then write the random fastq to the temp folder
  for (i in 1:length(fqs)){
    input_file <- paste0(input, names(fqs)[i], ".fastq.qz")
    ShortRead::writeFastq(fqs[[i]], input_file, mode="w", 
                          full=FALSE, compress=TRUE)
  }
  
  # Now we can run make_counts
  # Notice that make_counts will generate another temp folder, that will 
  # be emptied on finalization. By setting save_temp=TRUE you may save the 
  # content.  
  
  parse_cut = list(cutadapt="-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 70",
              fastq_quality_filter="-q 20 -p 80")
  
  parse_seq = list(adapt_3_set=c(type="hard_save", min=10, mismatch=0.1),
               adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTA",
               polyG=c(type="hard_trim", min=10, mismatch=0.1),
               seq_range=c(min=14, max=70),
               quality=c(threshold=20, percent=0.8),
               check_mem = TRUE)
  
  if(grepl("unix", .Platform$OS.type)) {
    quiet(
      counts_cut  <-  make_counts(input, threads=2,
                            trimming="cutadapt",
                            parse=parse_cut,
                            evidence=c(experiment=2, sample=1),
                            save_temp = FALSE)
      )
  }
  
  quiet(
    counts_seq  <-  make_counts(input, threads=2,
                           trimming="seqpac",
                           parse=parse_seq, 
                           evidence=c(experiment=2, sample=1),
                           save_temp = TRUE)
  )
  
    trim_files <- paste0(tempdir(), "/seqpac/")
    trim_files <- list.files(trim_files, pattern=".trim.fastq.gz$", full.names=TRUE)
    
  quiet(
    counts_trim  <-  make_counts(trim_files, threads=2,
                           trimming = NULL,
                           evidence=c(experiment=2, sample=1))
  )
  quiet(  
    counts_null  <-  make_counts(trim_files, threads=2,
                                 trimming = NULL,
                                 evidence=NULL)
  )
    
  

## Test counting and trimming
  ## Test both windows and linux
  seq_clss  <- unlist(lapply(lapply(counts_seq, function(x){x[[1]]}), class), use.names=FALSE)
  seq_n  <- c(nrow(counts_seq$counts)>1, 
              ncol(counts_seq$counts) == length(trim_files),
              nrow(counts_seq$progress_report) == length(trim_files))
  trim_clss  <- unlist(lapply(lapply(counts_trim, function(x){x[[1]]}), class), use.names=FALSE)
  trim_n  <- c(nrow(counts_trim$counts)>1, 
              ncol(counts_trim$counts) == length(trim_files),
              nrow(counts_trim$progress_report) == length(trim_files))

  expect_identical(names(counts_seq), c("counts","progress_report", "evidence_plots"))
  expect_identical(names(counts_trim), c("counts","progress_report", "evidence_plots"))
  

  expect_true(sum(seq_n) ==3)
  expect_true(sum(trim_n) ==3)
  

  expect_equal(sum(seq_clss %in% c("numeric", "character", "gg", "ggplot")), 4)
  expect_equal(sum(trim_clss %in% c("numeric", "character", "gg", "ggplot")), 4)


# ## Test make pheno and 
# test_that("Testing make_pheno and make_PAC ...", {
  
  Sample_ID <- gsub(".fastq.qz","", colnames(counts_seq$counts))
  
  pheno <- data.frame(Sample_ID=Sample_ID,
                      Treatment=c(rep("heat", times=1), 
                                  rep("control", times=2)),
                      Batch=rep(c("1", "2", "3"), times=1))
  
  quiet(                                                  
  pheno_1 <- make_pheno(pheno=pheno, progress_report=counts_seq$progress_report, 
                        counts=counts_seq$counts) 
 )  
 write.csv(pheno, file=paste0(input, "temp.csv"), row.names=FALSE)
  

 quiet(
     pheno_2 <- make_pheno(pheno=paste0(input, "temp.csv"),
                          progress_report=counts_trim$progress_report, 
                          counts=counts_trim$counts)
)
  
  # Test make_PAC both S3 and S4
  pac_seq <- make_PAC(pheno=pheno_1, 
                      counts=counts_seq$counts, 
                      output="S4")
  pac_trim <- make_PAC(pheno=pheno_2, 
                       counts=counts_trim$counts)
  
  as(pac_seq, "list") -> pac_S3
  expect_true(PAC_check(pac_seq))
  expect_true(PAC_check(pac_trim))
  expect_equal(names(pac_seq), c("Pheno","Anno","Counts"))
  expect_equal(length(pac_seq), length(pac_S3))
  expect_equal(nrow(pac_seq), nrow(pac_S3$Counts))
  expect_equal(ncol(pac_seq), 3)
  expect_equal(rownames(pac_seq), rownames(pac_S3$Counts))
  as.PAC(pac_S3) -> pac_S4
  expect_equal(rownames(pac_seq), rownames(pac_S4))

  
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

# Clean up temp
closeAllConnections()
fls_temp  <- tempdir()
fls_temp  <- list.files(fls_temp, recursive=TRUE, 
                       full.names = TRUE)
suppressWarnings(file.remove(fls_temp))
