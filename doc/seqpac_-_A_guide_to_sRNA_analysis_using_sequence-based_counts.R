## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = TRUE, eval=FALSE----------------------------------------------
#  ## Installation
#  devtools::install_github("Danis102/seqpac", upgrade="never", build_manual=TRUE, build_vignettes=TRUE)

## ---- results = "hide", eval=FALSE--------------------------------------------
#  library(seqpac)
#  
#  ## Using default settings for NEB type adaptor
#  # NEB=New England Biolabs: NEBNext® Small RNA Library Prep Set for Illumina (E7300/E7330)
#  # For illumina type adaptors, 'parse="default_illumina"' may be used
#  
#  path_to_fastq <- system.file("extdata", package = "seqpac", mustWork = TRUE)
#  count_list <- make_counts(input=path_to_fastq, type = "fastq", trimming="seqpac",
#                            parse="default_neb")
#  
#  
#  # That was quite slow, lets speed it up a bit by increasing threads:
#  count_list <- make_counts(input=path_to_fastq, type = "fastq", trimming="seqpac",
#                            parse="default_neb", threads=9)
#  
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  # Generate counts by parsing commands to cutadapt
#  # and fastq_quality_filter.
#  # (Note, the 'default_neb' parse argument is specificly designed for sRNA trimming
#  
#  count_list <- make_counts(input=path_to_fastq, type = "fastq", trimming="cutadapt", parse="default_neb")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ## Lets plot the graphs from the default output again
#  cowplot::plot_grid(plotlist=count_list$evidence_plot, nrow=2, ncol=1)
#  

## ---- results = "hide", eval=FALSE--------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Evidence over two indepenent samples, saving single sample sequences reaching 10 counts
#  test <- make_counts(input=path_to_fastq, type="fastq", trimming="seqpac",
#                      parse="default_neb", threads=9,
#                      evidence=c(experiment=2, sample=10))
#  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#  test$counts[extras==1,]  # 28 single sample sequences reach 10 counts
#  
#  ## Evidence over two indepenent samples, saving single sample sequences reaching 3 counts
#  test <- make_counts(input=path_to_fastq,  type="fastq", trimming="seqpac",
#                      parse="default_neb", threads=9,
#                      evidence=c(experiment=2, sample=3))
#  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#  test$counts[extras==1,] # A few hundred single sample sequences reach 3 counts
#  

## ---- results = "hide", eval=FALSE--------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Generate a Pheno table using the file names of test fastq
#  path_to_fastq <- system.file("extdata", package = "seqpac")
#  pheno <- as.data.frame(do.call("rbind",
#                         strsplit(list.files(path_to_fastq, pattern="*.fastq.gz"),
#                         "_|\\."))[,c(4, 1,2,3)])
#  colnames(pheno) <- c("sample", "stage", "batch", "index")
#  pheno$Sample_ID <- colnames(count_list$counts)
#  pheno <- make_pheno(pheno=pheno,
#                      progress_report=count_list$progress_report,
#                      counts=count_list$counts)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Generate PAC object
#  pac_master <- make_PAC(pheno=pheno, anno=NULL, counts=count_list$counts)
#  names(pac_master)
#  lapply(pac_master, head)
#  
#  ## If TRUE the PAC object is ok
#  PAC_check(pac_master)

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Load the PAC master and inspect the columns in Pheno and Anno
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac.Rdata", package = "seqpac", mustWork = TRUE))
#  
#  pac_master$Pheno[,1:5]
#  head(pac_master$Anno)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Extracts all sequences between 10-80 nt in length with at least 5 counts in
#  ## 20% of all samples.
#  pac_lowfilt <- PAC_filter(pac_master, size=c(10,80), threshold=5,
#                            coverage=20, norm = "counts",
#                            pheno_target=NULL, anno_target=NULL)
#  
#  ###---------------------------------------------------------------------
#  ## Extracts all sequences with 22 nt size and the samples in Batch1 and Batch2.
#  pac_subset <- PAC_filter(pac_master, subset_only = TRUE,
#                           pheno_target=list("batch", c("Batch1", "Batch2")),
#                           anno_target=list("Size", "22"))
#  
#  ###---------------------------------------------------------------------
#  ## Extracts all sequences with >=5 counts in 100% of samples a within stage
#  filtsep <- PAC_filtsep(pac_master, norm="counts", threshold=5,
#                         coverage=100, pheno_target= list("stage"))
#  
#  pac_filt <- PAC_filter(pac_master, subset_only = TRUE,
#                         anno_target= unique(do.call("c", as.list(filtsep))))
#  
#  
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Venn diagram of all sequences with >=10 counts in 100% of samples
#  ## within each group
#  olap <- reshape2::melt(filtsep,
#                         measure.vars = c("Stage1", "Stage3", "Stage5"), na.rm=TRUE)
#  plot(venneuler::venneuler(data.frame(olap[,2], olap[,1])))
#  

## ---- results=FALSE, eval=FALSE-----------------------------------------------
#  
#  ###---------------------------------------------------------------------
#  ## Example normalization in seqpac
#  pac_cpm <- PAC_norm(pac_master, norm="cpm")
#  pac_vst <- PAC_norm(pac_master, norm="vst")
#  

## ---- results=FALSE, eval=FALSE-----------------------------------------------
#  
#  ###---------------------------------------------------------------------
#  ## Filtering using cpm instead of raw counts
#  ## filter >=10 cpm in 100% of samples in >= 1 full group
#  filtsep <- PAC_filtsep(pac_cpm, norm="cpm", threshold=10, coverage=100, pheno_target= list("stage"))
#  pac_cpm_filt <- PAC_filter(pac_cpm, subset_only = TRUE,
#                         anno_target= unique(do.call("c", as.list(filtsep))))
#  
#  

## ---- results=FALSE, eval=FALSE-----------------------------------------------
#  
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt.Rdata", package = "seqpac", mustWork = TRUE))
#  head(pac_cpm_filt$Anno)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## First specify where to store all the fasta file
#  ref_path <- "/home/danis31/Desktop/Temp_docs/fasta"      # You need to change this
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## BSgenome
#  dest_path <- file.path(ref_path, "/uscs_genome/dm6_ucsc.fa")
#  BSgenome::available.genomes()
#  BiocInstaller::biocLite("BSgenome.Dmelanogaster.UCSC.dm6") # Only once
#  dm6 <- BSgenome::getBSgenome("BSgenome.Dmelanogaster.UCSC.dm6", masked=FALSE)
#  if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#  Biostrings::writeXStringSet(getSeq(dm6), filepath=dest_path, format="fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## biomartr
#  dest_path <- file.path(ref_path, "/biomartr_genome/") # Only dir
#  biomartr::is.genome.available(db = "ensembl", organism = "Drosophila melanogaster", details = TRUE)
#  file_path <- biomartr::getGenome(db="ensembl", organism = "Drosophila melanogaster", path=dest_path, gunzip=TRUE, release=101) # Bowtie don't take gzip
#  genome <- Biostrings::readDNAStringSet(filepath=list.files(dest_path, pattern=".fa$", full.names=TRUE), format="fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## biomartr ncRNA
#  dest_path <- file.path(ref_path, "/ensembl_ncRNA/") # Only dir
#  ncrna_path <- biomartr::getRNA(db="ensembl", organism = "Drosophila melanogaster", path=dest_path, release=101)
#  R.utils::gunzip(ncrna_path, destname=paste0(dest_path, "/ncrna.fa"), remove=TRUE, skip=TRUE) # Unzip to prepare for bowtie
#  ncrna <-  biomartr::read_rna(file = paste0(dest_path, "/ncrna.fa"))
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Download GtRNAdb
#  dest_path <- file.path(ref_path, "/GtRNAdb/trna.fa")
#  if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#  download.file(url="http://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.fa", destfile=dest_path)
#  trna <- Biostrings::readDNAStringSet(filepath=dest_path, format="fasta")

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Download mirbase
#  # 1. Download pre-miRNA data
#  # 2. Extract compressed fasta
#  # 3. Read as RNA
#  # 4. Extract species D. melanogaster
#  # 4. Convert to DNA
#  # 5. Overwrite RNA fasta with DNA fasta
#  # (sequenced reads are always in DNA)
#  # (we use pre-miRNA to catch all miRNA)
#  dest_path <- file.path(ref_path, "/mirbase/mirna.fa.gz")
#  if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#  download.file(url="ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz", destfile=dest_path)
#  R.utils::gunzip(dest_path, remove=TRUE, skip=TRUE)
#  mirna <- Biostrings::readRNAStringSet(filepath=gsub(".gz", "", dest_path), format="fasta")
#  mirna <- mirna[grepl("Drosophila melanogaster", names(mirna)),]
#  mirna <- Biostrings::DNAStringSet(mirna)
#  Biostrings::writeXStringSet(mirna, filepath=gsub(".gz", "", dest_path), format="fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Similar but MirGeneDB
#  # Since mirBase is rarely updated you might want to use MirGeneDB
#  dest_path <- file.path(ref_path, "/mirgenedb/mirna_pre_mirbase.fa.gz")
#  if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#  download.file(url="https://www.mirgenedb.org/static/data/dme/dme-pre.fas", destfile=dest_path)
#  mirna2 <- Biostrings::readRNAStringSet(filepath=dest_path, format="fasta")
#  mirna2 <- Biostrings::DNAStringSet(mirna2)
#  Biostrings::writeXStringSet(mirna2, filepath=dest_path, format="fasta")
#  
#  ## There are some overlaps
#  table(paste0(mirna2) %in% paste0(mirna))                          # Not much perfect overlap
#  table(grepl(paste(paste0(mirna2), collapse="|"), paste0(mirna)))  # Non-perfect overlaps better
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Extract scaffold in genome files
#  file_path <- list.files(paste0(ref_path, "/biomartr_genome/"), pattern="toplevel.fa", full.names=TRUE)
#  genome <- Biostrings::readDNAStringSet(filepath=file_path[1], format="fasta")
#  chrom <- genome[!grepl("scaffold", names(genome)),] # Remove scaffolds
#  Biostrings::writeXStringSet(chrom, filepath=paste0(dirname(file_path), "/chromsomes.fa") , format="fasta")
#  scaff <- genome[grepl("scaffold", names(genome)),] # Keep scaffolds
#  Biostrings::writeXStringSet(scaff, filepath=paste0(dirname(file_path), "/scaffolds.fa"), format="fasta")
#  
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Fixing the trna names
#  file_path <- paste0(ref_path, "/GtRNAdb/trna.fa")
#  trna <- Biostrings::readDNAStringSet(filepath=file_path, format="fasta")
#  names(trna) <- gsub("Drosophila_melanogaster_", "", names(trna))              # Remove species
#  mat <- do.call("rbind", strsplit(names(trna), " "))                           # Make a name matrix
#  names(trna) <-  paste(mat[,1], mat[,ncol(mat)-1], mat[,ncol(mat)], sep="_")   # Save the important as one single string
#  Biostrings::writeXStringSet(trna, filepath=file_path, format="fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Rearrange the ncRNA fasta names
#  file_path <- paste0(ref_path, "/ensembl_ncRNA/ncrna.fa")
#  ncrna <- Biostrings::readDNAStringSet(filepath=file_path, format="fasta")
#  mat <- do.call("rbind", strsplit(names(ncrna), " "))                        # Make a matrix of the names
#  mat <- mat[,1:7]                                                            # Pick only the 1st columns
#                                                                              # (some species contains multiple columns with same input)
#  col_bio <- grepl("gene_biotype:", mat[1,])                                  # Locate gene biotype column
#  col_coord <- grepl("chromosome:|scaffold:", mat[1,])                        # Locate coordinate column
#  identical(nrow(mat), sum(grepl("gene_biotype:", mat[,col_bio])))            # Did you catch all?
#  identical(nrow(mat), sum(grepl("chromosome:|scaffold:", mat[,col_coord])))
#  
#  new_names <- paste(mat[,ncol(mat)], mat[,1], mat[,col_bio], mat[,col_coord], sep="_") # Pick columns of your choice
#  names(ncrna) <- gsub("gene_symbol:|chromosome:BDGP6.28:|gene_biotype:", "",  new_names) # Clean up
#  Biostrings::writeXStringSet(ncrna, filepath=file_path, format="fasta")
#  
#  

## ---- eval=FALSE--------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Locate, fix names and add mito tRNA
#  trna_logi <- grepl("tRNA_mitochondrion_genome", names(ncrna))     # Locate mito tRNA
#  table(trna_logi) # Should be 22
#  mt_trna_ensembl <- ncrna[trna_logi,]
#  mt_trna <- gsub("mt:tRNA:", "MT-tRNA-", names(mt_trna_ensembl))
#  mt_trna <- gsub(":-1$", "_(-)", mt_trna)
#  mt_trna <- gsub(":1$", "_(+)", mt_trna)
#  mt_trna <- gsub("_FBtr", "-FBtr", mt_trna)
#  mt_trna <- gsub(":", "-", mt_trna)
#  mt_trna <- gsub("_tRNA_mitochondrion_genome-", "_chrMT:", mt_trna)
#  
#  names(mt_trna_ensembl) <- mt_trna
#  trna <- c(trna, mt_trna_ensembl)
#  Biostrings::writeXStringSet(trna, filepath=paste0(ref_path, "/GtRNAdb/trna.fa"), format="fasta")
#  

## ---- eval = FALSE------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Examples of generating bowtie indexes
#  
#  ref_path <- "/home/danis31/Desktop/Temp_docs/fasta"
#  Rbowtie::bowtie_build(paste0(ref_path, "/biomartr_genome/chromsomes.fa"),
#                        outdir=paste0(ref_path, "/biomartr_genome/"),
#                        prefix="chromosomes", force=TRUE)
#  Rbowtie::bowtie_build(paste0(ref_path, "/biomartr_genome/scaffolds.fa"),
#                        outdir=paste0(ref_path, "/biomartr_genome/"),
#                        prefix="scaffolds", force=TRUE)
#  Rbowtie::bowtie_build(paste0(ref_path, "/GtRNAdb/trna.fa"),
#                        outdir=paste0(ref_path, "/GtRNAdb/"),
#                        prefix="trna", force=TRUE)
#  Rbowtie::bowtie_build(paste0(ref_path, "/ensembl_ncRNA/ncrna.fa"),
#                        outdir=paste0(ref_path, "/ensembl_ncRNA/"),
#                        prefix="ncrna", force=TRUE)
#  Rbowtie::bowtie_build(paste0(ref_path, "/mirbase/mirna.fa"),
#                        outdir=paste0(ref_path, "/mirbase/"),
#                        prefix="mirna", force=TRUE)
#  

## ---- eval = FALSE------------------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Get gene gtf using biomartr
#  dest_path <- file.path(ref_path, "/gtf/") # Only dir
#  gtf_path <- biomartr::getGTF(db="ensembl", organism = "Drosophila melanogaster", path=dest_path)
#  gtf <- tibble::as_tibble(rtracklayer::readGFF(gtf_path))
#  
#  ###---------------------------------------------------------------------
#  ## Get repeatMasker gtf using biomartr (doesn't work for all species)
#  dest_path <- file.path(ref_path, "/repeatMasker/") # Only dir
#  rm_path <- biomartr::getRepeatMasker(db="refseq", organism = "Drosophila melanogaster", path=dest_path)
#  gtf <- tibble::as_tibble(rtracklayer::readGFF(gtf_path))
#  
#  ###---------------------------------------------------------------------
#  ## Get repeatMasker table and manually turn it into gtf using rtracklayer
#  # Table names can be found at:
#  # https://genome.ucsc.edu/cgi-bin/hgTables
#  dest_path <- file.path(ref_path, "/repeatMasker/repeatMasker.gtf") # Full file path
#  if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#  session <- rtracklayer::browserSession("UCSC")
#  rtracklayer::genome(session) <- "dm6"
#  rm_tab <- tibble::as_tibble(rtracklayer::getTable(rtracklayer::ucscTableQuery(session, track="RepeatMasker", table="rmsk")))
#  gr <- GenomicRanges::GRanges(seqnames=rm_tab$genoName, IRanges::IRanges(rm_tab$genoStart, rm_tab$genoEnd), strand=rm_tab$strand)
#  GenomicRanges::mcols(gr)$type <- "repeat"
#  GenomicRanges::mcols(gr)$source <- "repeatMasker_dm6"
#  GenomicRanges::mcols(gr)$repName <- rm_tab$repName
#  GenomicRanges::mcols(gr)$repClass <- rm_tab$repClass
#  GenomicRanges::mcols(gr)$repFamily <- rm_tab$repFamily
#  rtracklayer::export(gr, dest_path, format="gtf")
#  

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt.Rdata", package = "seqpac", mustWork = TRUE))
#  pac <- pac_cpm_filt
#  

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Genome mapping using external bowtie
#  
#  # Provide paths to bowtie indexed fasta references as a list
#  # OBS! You must change the paths.
#  
#  outpath_genome <- "/home/danis31/Desktop/Temp_docs/reanno_genome"
#  
#  ref_paths <- list(chromosomes="/home/danis31/Desktop/Temp_docs/fasta/biomartr_genome/chromosomes.fa",
#                    scaffolds="/home/danis31/Desktop/Temp_docs/fasta/biomartr_genome/scaffolds.fa")
#  
#  map_reanno(PAC=pac, ref_paths=ref_paths, output_path=outpath_genome,
#              type="external", mismatches=3, import="genome", threads=8, keep_temp=TRUE)
#  

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  ###---------------------------------------------------------------------
#  ## sRNA mapping using internal bowtie
#  
#  # Provide paths to bowtie indexed fasta references as a list
#   # OBS! You must change ths paths.
#  
#  outpath_biotype <- "/home/danis31/Desktop/Temp_docs/reanno_biotype"
#  
#  
#  ref_paths <- list(miRNA="/home/danis31/Desktop/Temp_docs/fasta/mirbase/mirna.fa",
#                     Ensembl="/home/danis31/Desktop/Temp_docs/fasta/ensembl_ncRNA/ncrna.fa",
#                     tRNA="/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.fa")
#  
#  map_reanno(pac, ref_paths=ref_paths, output_path=outpath_biotype,
#              type="internal", mismatches=3,  import="biotype", threads=8, keep_temp=TRUE)
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Generate a reanno object with make_reanno
#  
#  reanno_genome <- make_reanno(outpath_genome, PAC=pac, mis_fasta_check = TRUE)
#  reanno_biotype <- make_reanno(outpath_biotype, PAC=pac, mis_fasta_check = TRUE)
#  
#  # Note, setting mis_fasta_check=TRUE will double check that the number of
#  # sequences that failed to recieve an alignment in the last mismatch cycle
#  # agrees with the number sequences in the reanno object without an annotation.
#  # (these sequences are stored in mis_fasta_x.txt where x is max mismatches+1)
#  
#  # List structure
#  str(reanno_genome, max.level = 3, give.attr = FALSE)
#  str(reanno_biotype, max.level = 3, give.attr = FALSE)
#  
#  # Simple pie charts using the Overview table
#  pie(table(reanno_genome$Overview$chromosomes))
#  pie(table(reanno_genome$Overview$scaffolds))
#  pie(table(reanno_biotype$Overview$Any))
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  ###---------------------------------------------------------------------
#  ### Genomic coordinates using add_reanno
#  
#  # Output separate tibble
#  anno_genome <- add_reanno(reanno_genome, type="genome", genome_max=10, mismatches=3)
#  
#  # Output merged with provided PAC object
#  pac <- add_reanno(reanno_genome, type="genome", genome_max=10, mismatches=3, merge_pac=pac)
#  
#  # Example of original reference name annotations
#  head(reanno_genome$Full_anno$mis0$chromosomes)
#  
#  # Finished genome annotation
#  head(anno_genome)
#  head(pac$Anno)
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Classify sequences using add_reanno
#  
#  # Lets start by exploring the names in the original fasta reference:
#  
#  # Explore reference search terms
#  ref_path <- "/home/danis31/Desktop/Temp_docs/fasta"
#  
#  Ensembl_ncrna <- names(Biostrings::readDNAStringSet(paste0(ref_path,"/ensembl_ncRNA/ncrna.fa")))
#  Ensembl_ncrna <- do.call("rbind", strsplit(Ensembl_ncrna, " "))[,1]        # Extract up to 1st white space
#  Ensembl_ncrna <- do.call("rbind", strsplit(Ensembl_ncrna, "\\:"))[,1]      # Expressions up to ":"
#  table(Ensembl_ncrna)
#  
#  trna <- names(Biostrings::readDNAStringSet(paste0(ref_path,"/GtRNAdb/trna.fa")))
#  table(do.call("rbind", strsplit(trna, "\\_"))[,1])  # Divide
#  table(do.call("rbind", strsplit(trna, "\\-"))[,1])  # Further division
#  
#  mirna <- names(Biostrings::readDNAStringSet(paste0(ref_path,"/mirbase/mirna.fa")))
#  mirna <- do.call("rbind", strsplit(mirna, " "))[,1] # Extract up to 1st white space
#  table(do.call("rbind", strsplit(mirna, "\\-"))[,1]) # Divide
#  
#  
#  # Lets try two search term lists directed against each reference and written as
#  # 'regular expressions'
#  bio_search_1 <- list(
#                   Ensembl=c("lncRNA", "pre_miRNA", "rRNA", "snoRNA", "snRNA", "tRNA"),
#                   miRNA="dme-",
#                   tRNA =c("^tRNA", "MT")      # ^= Regular expression for start of string
#                  )
#  
#  bio_search_2 <- list(
#                   Ensembl=c("lncRNA", "miRNA", "pre_miRNA", "rRNA", "snoRNA",
#                             "snRNA", "tRNA", "Uhg", "7SLRNA", "asRNA", "hpRNA",
#                             "RNaseMRP","RNaseP", "sbRNA", "scaRNA", "sisRNA",
#                             "snmRNA", "snoRNA", "snRNA","Su\\(Ste\\)"), # ()= Brackets must be escaped
#                   miRNA="dme-",
#                   tRNA =c("^tRNA", "MT")
#                  )
#  
#  # Throws an error because perfect matching is required:
#  anno_temp <- add_reanno(reanno_biotype, bio_search=bio_search_1, type="biotype", bio_perfect=TRUE, mismatches = 3)
#  
#  # References with no search term hits are classified as "Other":
#  anno_temp <- add_reanno(reanno_biotype, bio_search=bio_search_1, type="biotype", bio_perfect=FALSE, mismatches = 3)
#  
#  # Better search terms gives perfect matching and merge with PAC
#  pac <- add_reanno(reanno_biotype, bio_search=bio_search_2, type="biotype", bio_perfect=TRUE, mismatches = 3, merge_pac=pac)
#  table(pac$Anno$mis0_bio)
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Hierarchical classification with simplify_reanno
#  
#  table(pac$Anno$mis0_bio)
#  
#  # Set the hierarchy (remember: order senstive; rRNA >>>>> lncRNA)
#  hierarchy <- list(rRNA="Ensembl_rRNA",
#                     Mt_tRNA="tRNA:MT",
#                     tRNA="Ensembl_tRNA|tRNA__tRNA",
#                     miRNA ="^miRNA|Ensembl_miRNA|Ensembl_pre_miRNA",
#                     snoRNA="Ensembl_snoRNA",
#                     lncRNA="Ensembl_lncRNA"
#                    )
#  
#  # No mistmach allowed
#  pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=0, bio_name="Biotypes_mis0", merge_pac=TRUE)
#  
#  # Up to 3 mismatches allowed
#  pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=3, bio_name="Biotypes_mis3", merge_pac=TRUE)
#  
#  # Example of original reference name annotations
#  head(reanno_biotype$Full_anno$mis0$Ensembl)
#  
#  # Simplified hierarchical classification
#  head(pac$Anno)
#  table(pac$Anno$Biotypes_mis0)
#  table(pac$Anno$Biotypes_mis3)
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#  
#  ###---------------------------------------------------------------------
#  ## PAC_summary in seqpac
#  
#  # Make means of counts over stages and return a data.frame
#  tab <- PAC_summary(pac, norm = "counts", type = "means",
#                     pheno_target=list("stage"), merge_pac=FALSE)
#  head(tab[[1]])
#  
#  # When merge_pac=TRUE the summarized table is added to the PAC$summary folder
#  pac_test <- PAC_summary(pac, norm = "counts", type = "means",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  
#  summary(pac)       # Structure of PAC before PAC_summary
#  summary(pac_test)  # Structure of PAC after PAC_summary
#  names(pac_test$summary)
#  head(pac_test$summary$countsMeans_stage)
#  
#  # You may want to use normalized counts
#  pac_test <- PAC_summary(pac_test, norm = "cpm", type = "means",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  
#  # Maybe only include a subset of the samples
#  pac_test <- PAC_summary(pac_test, norm = "cpm", type = "means",
#                          pheno_target=list("batch", c("Batch1", "Batch2")),
#                          merge_pac=TRUE)
#  
#  # Generate standard errors
#  pac_test <- PAC_summary(pac_test, norm = "cpm", type = "se",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  
#  # log2FC
#  pac_test <- PAC_summary(pac_test, norm = "cpm", type = "log2FC",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  
#  # log2FC generated from a grand mean over all samples
#  pac_test <- PAC_summary(pac_test, norm = "cpm", type = "log2FCgrand",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  
#  # All summarize tables have identical rownames that can be merged
#  names(pac_test$summary)
#  lapply(pac_test$summary, function(x){
#    identical(rownames(x), rownames(pac_test$summary[[1]]))
#    })
#  head(do.call("cbind", pac_test$summary))
#  

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#  
#  ###---------------------------------------------------------------------
#  ## Differential expression in seqpac
#  
#  # Simple model testing stages against using Wald test with local fit (default)
#  table(pac$Pheno$stage)
#  output_deseq <- PAC_deseq(pac, model= ~stage, threads=6)
#  
#  # More complicated, but still graphs will be generated from 'stage' since it is first in model
#  output_deseq <- PAC_deseq(pac, model= ~stage + batch)
#  
#  # Using pheno_target we can change focus
#  output_deseq <- PAC_deseq(pac, model= ~stage + batch, pheno_target=list("batch")) # No batch effect
#  
#  # With pheno_target we can also change the direction fo the comparision change focus
#  output_deseq <- PAC_deseq(pac, model= ~stage, pheno_target=list("stage", c("Stage3", "Stage1"))) # Zygotic transcription has not started
#  output_deseq <- PAC_deseq(pac, model= ~stage, pheno_target=list("stage", c("Stage5", "Stage3"))) # Start of zygotic transcritption
#  output_deseq <- PAC_deseq(pac, model= ~stage, pheno_target=list("stage", c("Stage5", "Stage1")))
#  
#  ## In the output you find PAC merged results, target plots and output_deseq
#  names(output_deseq)
#  tibble::as_tibble(output_deseq$result)
#  
#  

## ---- results = "hide", eval=FALSE--------------------------------------------
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#  
#  ###---------------------------------------------------------------------
#  ## PCA analysis in seqpac
#  
#  # As simple as possible
#  output_pca <- PAC_pca(pac)
#  names(output_pca)  # Two folders: graphs and pca results
#  
#  # Using pheno_target
#  output_pca <- PAC_pca(pac, pheno_target =list("stage"))
#  
#  # Using pheno_target with sample labels
#  output_pca <- PAC_pca(pac, pheno_target =list("stage"), label=pac$Pheno$sample)
#  
#  # Plotting sequences instead
#  output_pca <- PAC_pca(pac, type = "anno", anno_target =list("Biotypes_mis0"))
#  

## ---- include = TRUE, eval=FALSE----------------------------------------------
#  ###---------------------------------------------------------------------
#  ## Nucleotide bias in seqpac
#  
#  # Using master pac plotting 1st nt bias (default)
#  load(system.file("extdata", "drosophila_sRNA_pac.Rdata", package = "seqpac", mustWork = TRUE))
#  output_nbias <- PAC_nbias(pac_master)
#  cowplot::plot_grid(plotlist=output_nbias$Histograms)
#  
#  # Same but using filtered data
#  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#  output_nbias <- PAC_nbias(pac)
#  cowplot::plot_grid(plotlist=output_nbias$Histograms)
#  
#  # Only miRNA (Oops, heavy T-bias on 1st nt; are they piRNA?)
#  table(pac$Anno$Biotypes_mis0)
#  output_nbias <- PAC_nbias(pac, anno_target = list("Biotypes_mis0", "miRNA") )
#  cowplot::plot_grid(plotlist=output_nbias$Histograms)
#  
#  # Switch to 10:th nt bias
#  output_nbias <- PAC_nbias(pac, position=10, anno_target = list("Biotypes_mis0", "miRNA"))
#  cowplot::plot_grid(plotlist=output_nbias$Histograms)
#  
#  # Summarized over group cpm means
#  pac_test <- PAC_summary(pac, norm = "cpm", type = "means",
#                          pheno_target=list("stage"), merge_pac=TRUE)
#  output_nbias <- PAC_nbias(pac_test, summary_target = list("cpmMeans_stage") )
#  cowplot::plot_grid(plotlist=output_nbias$Histograms)
#  
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  ###---------------------------------------------------------------------
#  ## Biotype size distribution
#  
#  # Divide stacked bars by biotype with no mismach allowed
#  output_sizedist_1 <- PAC_sizedist(pac, anno_target = list("Biotypes_mis0"))
#  cowplot::plot_grid(plotlist=c(output_sizedist_1$Histograms), ncol=3, nrow=3)
#  
#  # Divide stacked bars by biotype with allowing up to 3 mismaches
#  output_sizedist_2 <- PAC_sizedist(pac, anno_target = list("Biotypes_mis3"))
#  cowplot::plot_grid(plotlist=c(output_sizedist_2$Histograms), ncol=3, nrow=3)
#  
#  # anno_target is order sensitive, thus can take care of color order issues:
#  ord_bio <- as.character(unique(pac$Anno$Biotypes_mis0))
#  ord_bio <- ord_bio[c(1,5,2,4,3,6,7)]
#  output_sizedist_1 <- PAC_sizedist(pac, anno_target = list("Biotypes_mis0", ord_bio))
#  output_sizedist_2 <- PAC_sizedist(pac, anno_target = list("Biotypes_mis3", ord_bio))
#  cowplot::plot_grid(plotlist=c(output_sizedist_1$Histograms[1:3], output_sizedist_2$Histograms[1:3]), ncol=3, nrow=2)
#  cowplot::plot_grid(plotlist=c(output_sizedist_1$Histograms[4:6], output_sizedist_2$Histograms[4:6]), ncol=3, nrow=2)
#  cowplot::plot_grid(plotlist=c(output_sizedist_1$Histograms[7:9], output_sizedist_2$Histograms[7:9]), ncol=3, nrow=2)
#  
#  ## Note: #######################################################################
#  # 1. miRNA is clearly associated with the correct size (21-23) nt, which are   #
#  # dramatically increased in Stage 5 when zygotic transcription has started.    #
#  # 2. piRNA was deliberately left out from the fasta references. Note, however, #
#  # that there is a broad peak with no annotions between 20-30 nt in Stage 1,    #
#  # which also showed a T-bias at the first nt. These are likely piRNA.          #
#  ################################################################################
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  ###---------------------------------------------------------------------
#  ## Stacked bars in seqpac
#  
#  # Choose an anno_target and plot samples (summary="samples")
#  PAC_stackbar(pac, anno_target=list("Biotypes_mis0"))
#  
#  # 'no_anno' and 'other' will always end on top not matter the order
#  ord_bio <- as.character(sort(unique(pac$Anno$Biotypes_mis3)))
#  p1 <- PAC_stackbar(pac, anno_target=list("Biotypes_mis0", ord_bio))
#  p2 <- PAC_stackbar(pac, anno_target=list("Biotypes_mis0", rev(ord_bio)))
#  cowplot::plot_grid(plotlist=list(p1, p2))
#  # (Hint: if you want them to appear not on top, rename them)
#  
#  # Reorder samples by pheno_targets
#  PAC_stackbar(pac, pheno_target=list("batch"), summary="samples", anno_target=list("Biotypes_mis0"))
#  
#  # Summarized over pheno_target
#  # (as default PAC_stackbar orders by pheno_target but plots all samples, unless summary="pheno")
#  PAC_stackbar(pac, anno_target=list("Biotypes_mis0"), summary="pheno", pheno_target=list("stage"))
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  ###---------------------------------------------------------------------
#  ## Pie chart in seqpac
#  
#  # Choose an anno_target and plot samples (summary="samples"; default)
#  output_pie <- PAC_pie(pac, anno_target=list("Biotypes_mis0"))
#  cowplot::plot_grid(plotlist=output_pie)
#  
#  output_pie <- PAC_pie(pac, anno_target=list("Biotypes_mis3"))
#  cowplot::plot_grid(plotlist=output_pie, ncol=3, scale = 0.8)
#  
#  
#  # Make ordered pie charts of grand mean percent of all samples
#  ord_bio <- as.character(sort(unique(pac$Anno$Biotypes_mis3)), unique(pac$Anno$Biotypes_mis0))
#  output_pie_1 <- PAC_pie(pac, anno_target=list("Biotypes_mis0", ord_bio), summary="all")
#  output_pie_2 <- PAC_pie(pac, anno_target=list("Biotypes_mis3", ord_bio), summary="all")
#  cowplot::plot_grid(plotlist=c(output_pie_1, output_pie_2), nrow=2)
#  
#  # Rotate
#  PAC_pie(pac, anno_target=list("Biotypes_mis0"), summary="all", angle=180)
#  PAC_pie(pac, anno_target=list("Biotypes_mis0"), summary="all", angle=40)
#  
#  # Compare biotype mapping with or without mismaches and group by PAC$Pheno
#  # (Notice: all stages gains a lot of 'hidden' lncRNA by allowing mismatches)
#  ord_bio <- as.character(sort(unique(pac$Anno$Biotypes_mis0)))  # Make sure that both get same biotype order
#  output_mis0 <- PAC_pie(pac, pheno_target=list("stage"), summary="pheno", anno_target=list("Biotypes_mis0", ord_bio))
#  output_mis3 <- PAC_pie(pac, pheno_target=list("stage"), summary="pheno", anno_target=list("Biotypes_mis3", ord_bio))
#  cowplot::plot_grid(plotlist=c(output_mis0, output_mis3), labels = names(c(output_mis0, output_mis3)), nrow=2)
#  
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  library(seqpac)
#  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#  
#  ###---------------------------------------------------------------------
#  ## tRNA analysis in seqpac
#  
#  # First create an annotation blanc PAC with group means
#  pac$Anno <- pac$Anno[,1, drop=FALSE]
#  pac_trna <- PAC_summary(pac, norm = "cpm", type = "means", pheno_target=list("stage"), merge_pac = TRUE)
#  
#  # Then reannotate only tRNA using the PAC_mapper function
#  ref <- "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.fa"  # Your path to tRNA fasta
#  map_object <- PAC_mapper(pac_trna, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=TRUE)
#  
#  ### Hint:
#  # By adding N_up ad N_down you can make sure that modified fragments (like 3'
#  # -CAA in mature tRNA are included).
#  
#  ###---------------------------------------------------------------------
#  ## Coverage plot of tRNA using PAC_covplot
#  
#  # Single tRNA targeting a summary dataframe
#  PAC_covplot(pac_trna, map=map_object, summary_target= list("cpmMeans_stage"), map_target="tRNA-Ala-AGC-1-1_chr3R:17657145-17657217_(+)")
#  
#  # Find tRNAs with many fragments
#  n_tRFs <- unlist(lapply(map_object, function(x){nrow(x[[2]])}))
#  table(n_tRFs)
#  names(map_object)[n_tRFs>2]
#  selct <- (names(map_object)[n_tRFs>1])[c(1, 16, 25, 43)]
#  cov_plt <- PAC_covplot(pac_trna, map=map_object, summary_target= list("cpmMeans_stage"), map_target=selct)
#  cowplot::plot_grid(plotlist=cov_plt, nrow=2, ncol=2)
#  
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  
#  ###---------------------------------------------------------------------
#  ## Analyze range types with map_rangetype and PAC_trna functions
#  
#  # Download ss object from GtRNAdb
#  dest_path <- file.path("/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/trna.tar.gz")
#  download.file(url="http://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-tRNAs.tar.gz", destfile=dest_path)
#  untar(dest_path, exdir= dirname(dest_path), files = "dm6-tRNAs-confidence-set.ss")
#  ss_file <- "/home/danis31/Desktop/Temp_docs/fasta/GtRNAdb/dm6-tRNAs-confidence-set.ss"
#  
#  # Classify fragments according to loop cleavage (small loops are omitted)
#  map_object_ss <- map_rangetype(map_object, type="ss", ss=ss_file, min_loop_width=4)
#  
#  # Remove reference tRNAs with no hits
#  map_object_ss <-  map_object_ss[!unlist(lapply(map_object_ss, function(x){x[[2]][1,1] == "no_hits"}))]
#  map_object_ss[[2]]
#  
#  
#  ###---------------------------------------------------------------------
#  ## Function classifying 5'-tRF, 5'halves, i-tRF, 3'-tRF, 3'halves
#  
#  # Does all tRNAs have 3 loops?
#  table(unlist(lapply(map_object_ss, function(x){unique(x$Alignments$n_ssloop)})))
#  
#  # Set tolerance for classification as a terminal tRF
#  tolerance <- 5  # 2 nucleotides from start or end of full-length tRNA)
#  
#  ### Important:
#  # We set N_up and N_down to "NNN" in the PAC_mapper step. To make sure
#  # that we have a tolerance that include the original tRNA sequence
#  # we set terminal= 2+3 (5).
#  
#  ## tRNA classifying function
#  # Apply the tRNA_class function and make a tRNA type column
#  pac_trna <- tRNA_class(pac_trna, map=map_object_ss, terminal=tolerance)
#  pac_trna$Anno$type <- paste0(pac_trna$Anno$decoder, pac_trna$Anno$acceptor)
#  head(pac_trna$Anno)
#  

## ---- message=FALSE, include = TRUE, eval=FALSE-------------------------------
#  
#  ###---------------------------------------------------------------------
#  ## Plotting tRNA types
#  
#  # Now use PAC_trna to generate some graphs based on grand means
#  trna_result <- PAC_trna(pac_trna, norm="cpm", filter = NULL,
#    join = TRUE, top = 15, log2fc = TRUE,
#    pheno_target = list("stage", c("Stage1", "Stage3")),
#    anno_target_1 = list("type"),
#    anno_target_2 = list("class"))
#  
#  names(trna_result)
#  names(trna_result$plots)
#  names(trna_result$plots$Expression_Anno_1)
#  
#  cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Grand_means,
#                     trna_result$plots$Log2FC_Anno_1,
#                     trna_result$plots$Percent_bars$Grand_means,
#                     nrow=1, ncol=3)
#  
#  # By setting join = FALSE you will get group means
#  trna_result <- PAC_trna(pac_trna, norm="cpm", filter = NULL,
#    join = FALSE, top = 15, log2fc = TRUE,
#    pheno_target = list("stage", c("Stage1", "Stage3")),
#    anno_target_1 = list("type"),
#    anno_target_2 = list("class"))
#  
#  names(trna_result$plots$Expression_Anno_1)
#  
#  cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Stage1,
#                     trna_result$plots$Expression_Anno_1$Stage3,
#                     trna_result$plots$Log2FC_Anno_1,
#                     trna_result$plots$Percent_bars$Stage1,
#                     trna_result$plots$Percent_bars$Stage3,
#                     nrow=1, ncol=5)
#  

