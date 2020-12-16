# # Benchmarking seqpac
# 
# library(rbenchmark)
# library(seqpac)
# 
# 
# #### Download SRA ######
# ## Avital:
# sra_acc <- as.character(read.csv("/data/Data_analysis/Projects/Drosophila/SRA_Download/Avital_etal_Embryo_stages/SRR_Acc_List.txt", header=FALSE)[,1])
# sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname )
# SRAdb::getSRAfile(in_acc=sra_acc, sra_con, destDir = getwd(), fileType = 'fastq', srcType = 'ftp', makeDirectory = FALSE, method = 'curl', ascpCMD = NULL )
# 
# 
# 
# SRR5935244, SRR5935245, SRR5935247, SRR5935249, SRR5935251, SRR5935253, SRR5935255, SRR5935257
# SRR5935259	SRR5935261	SRR5935263	SRR5935265	SRR5935266	SRR5935268
# 
# 
# 
# ##### Kang et al with Illumina adaptor #######
# input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq"
# output = "/home/danis31/Desktop/Temp_docs/temp/test/bench_temp"
# 
# make_cutadapt_alone <- function(input, output, parse=NULL, threads=1){
# 
#               ## Setup
#               cat("\nRunning make_cutadapt ...")
#               cat("\n--- cutadapt and fastq_quality_filter must be correctly installed.")
#               fls <- list.files(input, pattern ="fastq.gz\\>|fastq\\>", full.names=TRUE, recursive=TRUE)
# 
#               # Make dir
#               if(!dir.exists(output)){
#                 suppressWarnings(dir.create(output))
#               }
# 
#               # Make output names and check output folder
#               nam <- gsub("\\.gz$", "", basename(fls))
#               nam <- gsub("\\.fastq$|\\.fq$|fastq$", "", nam)
#               nam <- gsub("\\.$|\\.tar$", "", nam)
#               nam_trim <- paste0(nam, ".trim.fastq.gz")
#               # Check for output folder
#               out_file <- file.path(output, nam_trim)
#               out_dir <- list.files(output, pattern=nam_trim, recursive = FALSE)
#               if(length(out_dir)>0){
#                   stop(paste0("\n  Output trimmed fastq file names are identical to existing files in output:\n  ", out_dir, "\n  Please move or delete the file in the output folder:\n  ", output))
#               }
# 
#               ## cutadapt and fastq_quality_filter
#               prog_report <- list(NULL)
#               for(i in 1:length(fls)){
#                       log_lst <- list(NULL, NULL)
#                       names(log_lst) <- c("cutadapt", "fastq_quality_filter")
#                       spl_nam <- nam_trim[i]
#                       temp_out <- gsub("trim.fastq.gz$", "temp.fastq", out_file[i])
#                       if(!is.null(parse[[1]])){
#                          log_lst[[1]] <- system(paste0("cutadapt ", parse[[1]], " -o ", temp_out, " ", fls[i]), intern = TRUE)
#                       }
#                       if(!is.null(parse[[2]])){
#                          log_lst[[2]] <- system(paste0("fastq_quality_filter ", parse[[2]], " -v -i ", temp_out, " -o ", out_file[i], " -z"), intern = TRUE)
#                       }
# 
#                       if(!is.null(parse[[1]])){
#                            file.remove(temp_out)
#                           }
#                       prog_report[[1]] <- log_lst
#                 }
#               doParallel::stopImplicitCluster()
#               cat("\n--- Finished generating trimmed temporary files.")
#               return(prog_report)
#               gc(reset=TRUE)
#               }
# ##########################################################################
# ##### Benchmarking Kang et al with Illumina adaptor #######
# bn <- rbenchmark::benchmark("seqpac" = {
#            prog_report  <-  make_trim(input=input, output=output, threads=7, check_mem=FALSE, indels=TRUE, concat=12,
#                      adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA",
#                      polyG=c(type="hard_trim", min=20, mismatch=0.1),
#                      seq_range=c(min=14, max=45),
#                      quality=c(threshold=20, percent=0.8))
# 
#            fn <- list.files(output, full.names=TRUE, recursive=FALSE)
#            file.remove(fn)
#           },
#           "cutadapt" = {
#            prog_report  <-  make_cutadapt(input, output, threads=7, parse= list(
#                      cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45",
#                      fastq_quality_filter="-q 20 -p 80"))
#            fn <- list.files(output, full.names=TRUE, recursive=FALSE)
#            file.remove(fn)
#           },
#           "cutadapt_alone" = {
#            prog_report  <-  make_cutadapt_alone(input, output, threads=7, parse= list(
#                      cutadapt="-j 7 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45",
#                      fastq_quality_filter="-q 20 -p 80"))
#            fn <- list.files(output, full.names=TRUE, recursive=FALSE)
#            file.remove(fn)
#           },
#           replications = 9,
#           columns = c("test", "replications", "elapsed",
#                     "relative", "user.self", "sys.self", "user.child", "sys.child"))
# 
# 
# ##########################################################################
# ##### Kang et al with Illumina adaptor #######
# input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq/"
# 
# output = "/home/danis31/Desktop/Temp_docs/temp/test/cutadapt"
# prog_report  <-  make_cutadapt(input, output, threads=7, parse= list(
#                      cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45 -e 0.1",
#                      fastq_quality_filter="-q 20 -p 80"))
# 
# output = "/home/danis31/Desktop/Temp_docs/temp/test/seqpac"
# prog_report  <-   make_trim(input=input, output=output, threads=7, indels=TRUE, concat=12, check_mem=FALSE,
#                      adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA",
#                      polyG=c(type="hard_trim", min=20, mismatch=0.1),
#                      seq_range=c(min=14, max=45),
#                      quality=c(threshold=20, percent=0.8))
# 
# ##########################################################################
# ######### Seqpac drosophila NEBNext adator dataset ##############
# input <- system.file("extdata", package = "seqpac", mustWork = TRUE)
# 
# output = "/home/danis31/Desktop/Temp_docs/temp/test2/cutadapt/"
# prog_report  <-  make_cutadapt(input, output, threads=7, parse= list(
#                      cutadapt="-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCA --discard-untrimmed --nextseq-trim=20 -O 10 -m 12 -M 45 -e 0.1",
#                      fastq_quality_filter="-q 20 -p 80"))
# 
# output = "/home/danis31/Desktop/Temp_docs/temp/test2/seqpac"
# prog_report  <-  make_trim(input=input, output=output, threads=7, indels=TRUE, concat=12, check_mem=FALSE,
#                      adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
#                      polyG=c(type="hard_trim", min=20, mismatch=0.1),
#                      seq_range=c(min=14, max=45),
#                      quality=c(threshold=20, percent=0.8))
# 
# 
# 
# ##################################################################
# 
# input_seq <- "/home/danis31/Desktop/Temp_docs/temp/test2/seqpac"
# input_cut <- "/home/danis31/Desktop/Temp_docs/temp/test2/cutadapt"
# #input_orig <- "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq"
# input_orig <- system.file("extdata", package = "seqpac", mustWork = TRUE)
# 
# fls_seq <- list.files(input_seq, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
# fls_cut <- list.files(input_cut, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
# fls_orig  <- list.files(input_orig, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
# fastq_lst <- list(NULL)
# fastq_lst[[1]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_seq[[1]], withIds=FALSE)))
# fastq_lst[[2]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_cut[[1]], withIds=FALSE)))
# fastq_lst[[3]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_orig[[1]], withIds=FALSE)))
# 
# table(unique(fastq_lst[[1]]) %in% fastq_lst[[2]])
# table(unique(fastq_lst[[2]]) %in% fastq_lst[[1]])
# 
# seq <- fastq_lst[[1]][ !fastq_lst[[1]] %in% unique(fastq_lst[[2]])]
# cut <- fastq_lst[[2]][ !fastq_lst[[2]] %in% unique(fastq_lst[[1]])]
# 
# head(cut)
# head(seq)
# 
# tar_seq <- "AACACTAAGCGGTGGATCACTCGGC"
# fastq_lst[[1]][fastq_lst[[1]] == tar_seq]
# fastq_lst[[2]][fastq_lst[[2]] == tar_seq]
# 
# fastq_lst[[1]][grepl(paste0("^", tar_seq), fastq_lst[[1]])]
# fastq_lst[[2]][grepl(paste0("^", tar_seq), fastq_lst[[2]])]
# 
# fastq_lst[[3]][grepl(paste0("^", tar_seq), fastq_lst[[3]])]
# 
# 
# 
# tar_seq <- paste0("TATTTGAACGAGAAACCTGTAACCAACTCTCAACTGATGAGATCGAAGAGCACACGTCTGAACTCCAGTCACTAG")
# which(grepl(tar_seq, fastq_lst[[3]])) #7043
# 
# 
# AACACTAAGCGGTGGATCACTCGGC
# AACACTAAGCGGTGGATCACTCGGCAGAACGGAAGAGCACACGACTGAACACCAGACACTAGCTTAACTCGTATG
#                          AGATCGGAAGAGCACACGTCTGAACTCCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# 
# 
# 
# 
# 
# 
# ###################################################################################################################################
# ###################################################################################################################################
# ###################################################################################################################################
# ##########################################################################
# ##### Kang et al count table #######
# input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq/"
# 
# ##########
# ## Cutadapt
# parse = list(cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45 -e 0.1",
#              fastq_quality_filter="-q 20 -p 80")
# counts_cut  <- make_counts(input, threads=7, parse=parse,
#                         type="fastq", trimming="cutadapt", plot=TRUE,
#                         evidence=c(experiment=2, sample=1))
# anno_cut  <- make_anno(counts_cut, type="counts")
# pheno_cut <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace", type="manual", counts=counts_cut[[1]], progress_report=counts_cut[[2]])
# pac_cut <- make_PAC(pheno_input=pheno_cut, anno_input=anno_cut, counts_input=counts_cut[[1]])
# PAC_check(pac_cut)
# save(pac_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")
# 
# #########
# ## Seqpac
# parse = list(adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
#              polyG=c(type="hard_trim", min=20, mismatch=0.1),
#              seq_range=c(min=14, max=45),
#              quality=c(threshold=20, percent=0.8),
#              indels=TRUE, concat=12, check_mem=FALSE)
# counts_seq  <- make_counts(input, threads=7, parse=parse,
#                         type="fastq", trimming="seqpac", plot=TRUE,
#                         evidence=c(experiment=2, sample=1))
# anno_seq  <- make_anno(counts_seq, type="counts")
# pheno_seq <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace", type="manual", counts=counts_seq[[1]], progress_report=counts_seq[[2]])
# 
# pac_seq <- make_PAC(pheno_input=pheno_seq, anno_input=anno_seq, counts_input=counts_seq[[1]])
# PAC_check(pac_seq)
# #save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
# 
# 
# 
# 
# ###################################################################################################################################
# ##########################################################################
# ##### Kang et al reanno genome  #######
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")
# 
# output_path = "/home/danis31/Desktop/Temp_docs/temp"
# 
# ref_paths_gen <- list(dm6="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa",
#                       hg38="/data/Data_analysis/Genomes/Humans/Ensembl/CRCh38.101/Homo_sapiens.GRCh38.dna.toplevel.fa")
# 
# ### seqpac
# map_reanno(pac_seq, type = "internal", output_path, ref_paths=ref_paths_gen, mismatches = 3, threads = 7, import="genome")
# reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_seq,  mis_fasta_check = TRUE, threads = 7)
# #save(reanno_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_seqpac.Rdata")
# #load(file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_seqpac.Rdata")
# pac_seq <- add_reanno(reanno=reanno_seq, mismatches = 3,  merge_pac=pac_seq, type = "genome", genome_max = 10)
# save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
# # Sequences in total: 491968
# #	mis0	dm6:  	Ref_hits= 57901     (11.77%)
# #	mis0	hg38: 	Ref_hits= 229082    (46.56%)
# #	mis1	dm6:  	Ref_hits= 22014     (4.47%)
# #	mis1	hg38: 	Ref_hits= 88590     (18.01%)
# #	mis2	dm6:  	Ref_hits= 13616     (2.77%)
# #	mis2	hg38: 	Ref_hits= 54194     (11.02%)
# #	mis3	dm6:  	Ref_hits= 15112     (3.07%)
# #	mis3	hg38: 	Ref_hits= 36839     (7.49%)
# 
# ### cutadapt
# map_reanno(pac_cut, type = "internal", output_path, ref_paths=ref_paths_gen, mismatches = 3, threads = 7, import="genome")
# reanno_cut <- make_reanno(reanno_path=output_path, PAC=pac_cut,  mis_fasta_check = TRUE, threads = 7)
# #save(reanno_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_cutadapt.Rdata")
# #load(file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_cutadapt.Rdata")
# pac_cut <- add_reanno(reanno=reanno_cut, mismatches = 3,  merge_pac=pac_cut, type = "genome", genome_max = 10)
# #save(pac_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadapt.Rdata")
# # Sequences in total: 479491
# #	mis0	dm6:  	Ref_hits= 57167     (11.92%)
# #	mis0	hg38: 	Ref_hits= 223113    (46.53%)
# #	mis1	dm6:  	Ref_hits= 21765     (4.54%)
# #	mis1	hg38: 	Ref_hits= 86814     (18.11%)
# #	mis2	dm6:  	Ref_hits= 13410     (2.8%)
# #	mis2	hg38: 	Ref_hits= 53352     (11.13%)
# #	mis3	dm6:  	Ref_hits= 14955     (3.12%)
# #	mis3	hg38: 	Ref_hits= 36183     (7.55%)
# 
# ###################################################################################################################################
# ##########################################################################
# ##### Kang et al reanno genome  #######
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")
# 
# output_path = "/home/danis31/Desktop/Temp_docs/temp"
# 
# ref_paths_bio <- list(dm6_miRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/miRNA/miRBase_21-dme.fa",
#                   dm6_Ensembl="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/Ensembl/Drosophila_melanogaster.BDGP6.ncrna.fa",
#                   dm6_tRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/tRNA_reanno/tRNA_mature.fa",
#                   dm6_piRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.fa",
# 
#                   hg38_miRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/miRBase_21/miRBase_21-hsa.fa",
#                   hg38_Ensembl="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/Ensembl/Homo_sapiens.GRCh38.ncrna.fa",
#                   hg38_tRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/GtRNAdb/hg19-tRNAs_mature.fa",
#                   hg38_piRNA="/data/Data_analysis/Genomes/Humans/pirBase/hsa.fa")
# 
# 
# bio_search <- list(
#                  dm6_Ensembl=c("lincRNA", "miRNA", "rRNA", "pre_miRNA", "snoRNA", "snRNA", "_tRNA", "mt_tRNA"),
#                  dm6_miRNA="dme-",
#                  dm6_piRNA = "piR-dme|piRNA",
#                  dm6_tRNA =c("^tRNA", "mt-tRNA"),
# 
#                  hg38_miRNA="hsa",
#                  hg38_Ensembl=c("lncRNA", "miRNA", "rRNA", "snoRNA", "snRNA", "tRNA"),
#                  hg38_tRNA=c("tRNA"),
#                  hg38_piRNA="piR-hsa|piRNA")
# 
# hierarchy <- list(rRNA="rRNA",
#                    tRNA="tRNA",
#                    miRNA ="miRNA",
#                    snoRNA="snoRNA",
#                    snRNA="snRNA",
#                    lncRNA="lncRNA|lincRNA",
#                    piRNA="piRNA"
#                   )
# 
# 
# ### seqpac
# map_reanno(pac_seq, type = "internal", output_path, ref_paths=ref_paths_bio, mismatches = 3, threads = 7, import="biotype")
# reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_seq,  mis_fasta_check = TRUE, threads = 7)
# #save(reanno_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_biotype_seqpac.Rdata")
# 
# pac_seq <- add_reanno(reanno_seq, bio_search=bio_search, type="biotype", bio_perfect=FALSE, mismatches = 0, merge_pac=pac_seq)
# pac_seq <- simplify_reanno(pac_seq, hierarchy=hierarchy, mismatches = 0, bio_name = "Biotype", merge_pac = TRUE)
# 
# 
# save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac_anno.Rdata")
# 
# 
# ###################################################################################################################################
# ###################################################################################################################################
# ###################################################################################################################################
# ##########################################################################
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac_anno.Rdata")
# load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadapt.Rdata")
# 
# ## Wenn-diagram 
# input <- list(seqpac=rownames(pac_seq$Anno), cutadapt=rownames(pac_cut$Anno))
# PAC_venn <- function(input, type="venn"){
#     n_input <- length(input)
#     nams <- names(input)
#     names(input) <- LETTERS[1:n_input]
#     all_seqs <- unique(unlist(input))
#   # For sav logi table:
#     comb_logi <- do.call("cbind", lapply(1:n_input, function(x){
#       logi <- all_seqs %in% input[[x]]
#      }))
#     colnames(comb_logi) <- nams
#     comb_logi <- cbind(data.frame(seqs=all_seqs, comb_logi))
#   # For venn diagram:  
#     comb_ven <- do.call("cbind", lapply(1:n_input, function(x){
#       logi <- all_seqs %in% input[[x]]
#       ifelse(logi, names(input)[x], NA)
#     }))
#     comb_ven <- table(apply(comb_ven, 1, function(x){paste(x , collapse="&")}))
#     names(cnts) <- gsub("NA&|&NA", "", names(cnts))
#     intg <- as.integer(cnts)
#     names(intg) <- names(cnts)
#   
#   if(type=="euler"){
#     ven <- venneuler::venneuler(intg) 
#     ven$labels <- nams
#     return(list(venn=ven, olap=comb_logi)) 
#   }
#   if(type=="venn"){
#     ven <- ggvenn::ggvenn(input)
#   }
#   return(list(venn=ven, olap=comb_logi)) 
# }
# 
# ## Nbias plots seqpac
# nbias_seq <- PAC_nbias(pac_seq)
# cowplot::plot_grid(plotlist=nbias_seq$Histograms)
# nbias_cut <- PAC_nbias(pac_cut)
# cowplot::plot_grid(plotlist=nbias_cut$Histograms)
# 
# ## Size dist seqpac
# ord <- c("no_anno", "other", "miRNA", "tRNA", "snoRNA", "piRNA" , "snRNA", "lncRNA", "rRNA")
# sizebio_seq <- PAC_sizedist(pac_seq, anno_target=list("Biotype", ord))
# cowplot::plot_grid(plotlist=sizebio_seq$Histograms)
# 
# 
# ## Species genome mapping pie
# 
# pac_seq$Anno$mis0_genome_olap <- paste(ifelse(pac_seq$Anno$dm6 == "mis0", "fly", "no"), ifelse(pac_seq$Anno$hg38 == "mis0", "human", "no"))
# plts_mis0 <- PAC_pie(pac_seq, pheno_target=list("Sample_ID"), anno_target=list("mis0_genome_olap"))
# cowplot::plot_grid(plotlist=plts_mis0, nrow=3, ncol=3, scale=0.8)
# 
# 
# pac_seq$Anno$mis3_genome_olap <- paste(ifelse(pac_seq$Anno$dm6 %in% c("mis0","mis1", "mis2", "mis3"), "fly", "no"), ifelse(pac_seq$Anno$hg38 %in% c("mis0","mis1", "mis2", "mis3"), "human", "no"))
# plts_mis3 <- PAC_pie(pac_seq, pheno_target=list("Sample_ID"), anno_target=list("mis3_genome_olap", c("no no", "fly no", "no human", "fly human")))
# cowplot::plot_grid(plotlist=plts_mis3, nrow=3, ncol=3, scale=0.8)
# 
# 
# ###################################################################################################################################
# ###################################################################################################################################
# ###################################################################################################################################
# ##########################################################################
# ## PAC_tRNA
# 




