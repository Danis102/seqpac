# Benchmarking seqpac

library(rbenchmark) library(seqpac)


#### Download SRA ###### ## Avital: sra_acc <-
as.character(read.csv("/data/Data_analysis/Projects/Drosophila/SRA_Download/Avital_etal_Embryo_stages/SRR_Acc_List.txt",
                      header=FALSE)[,1]) sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname )
SRAdb::getSRAfile(in_acc=sra_acc, sra_con, destDir = getwd(), fileType =
                    'fastq', srcType = 'ftp', makeDirectory = FALSE, method = 'curl', ascpCMD =
                    NULL )



SRR5935244, SRR5935245, SRR5935247, SRR5935249, SRR5935251, SRR5935253,
SRR5935255, SRR5935257
SRR5935259	SRR5935261	SRR5935263	SRR5935265	SRR5935266	SRR5935268



##### Kang et al with Illumina adaptor #######

input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq"
output = "/home/danis31/Desktop/Temp_docs/temp/test/bench_temp"

make_cutadapt_alone <- function(input, output, parse=NULL, threads=1){
  ## Setup
  cat("\nRunning make_cutadapt ...")
  cat("\n--- cutadapt and fastq_quality_filter must be correctly installed.")
  fls <- list.files(input, pattern ="fastq.gz\\>|fastq\\>", full.names=TRUE, recursive=TRUE)

# Make dir if(!dir.exists(output)){ suppressWarnings(dir.create(output)) }

# Make output names and check output folder
  nam <- gsub("\\.gz$", "", basename(fls))
  nam <- gsub("\\.fastq$|\\.fq$|fastq$", "", nam)
  nam <- gsub("\\.$|\\.tar$", "", nam)
  nam_trim <- paste0(nam, ".trim.fastq.gz")
# Check for output folder
  out_file <- file.path(output, nam_trim)
  out_dir <- list.files(output, pattern=nam_trim, recursive = FALSE)
  if(length(out_dir)>0){
    stop(paste0("\n  Output trimmed fastq file names are identical to existing files in output:\n  ", out_dir, "\n  Please move or delete the file in the output folder:\n  ", output))
    }
## cutadapt and fastq_quality_filter
  prog_report <- list(NULL) for(i in 1:length(fls)){
    log_lst <- list(NULL, NULL)
    names(log_lst) <- c("cutadapt", "fastq_quality_filter")
    spl_nam <- nam_trim[i]
    temp_out <- gsub("trim.fastq.gz$", "temp.fastq", out_file[i])
    if(!is.null(parse[[1]])){
      log_lst[[1]] <- system(paste0("cutadapt ", parse[[1]], " -o ", temp_out, " ", fls[i]), intern = TRUE)
      }
    if(!is.null(parse[[2]])){
      log_lst[[2]] <- system(paste0("fastq_quality_filter ", parse[[2]], " -v -i ", temp_out, " -o", out_file[i], " -z"), intern = TRUE)
      }
    if(!is.null(parse[[1]])){
      file.remove(temp_out) } prog_report[[1]] <- log_lst
    }
    doParallel::stopImplicitCluster()
    cat("\n--- Finished generating trimmed temporary files.")
    return(prog_report) gc(reset=TRUE) }



##########################################################################
##### Benchmarking Kang et al with Illumina adaptor ####### 
bn <-
rbenchmark::benchmark(
  "seqpac" = {
    prog_report  <-  make_trim(input=input, output=output, threads=7, check_mem=FALSE, indels=TRUE, concat=12,
                             adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1),
                             adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA",
                             polyG=c(type="hard_trim",min=20, mismatch=0.1),
                             seq_range=c(min=14, max=45), quality=c(threshold=20, percent=0.8))

  fn <- list.files(output, full.names=TRUE, recursive=FALSE)
  file.remove(fn) },

  "cutadapt" = {
    prog_report  <-  make_cutadapt(input, output, threads=7,
                                   parse=list(cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45",
                                   fastq_quality_filter="-q 20 -p 80"))
    fn <- list.files(output, full.names=TRUE, recursive=FALSE)
    file.remove(fn) },

  "cutadapt_alone" = {
    prog_report  <- make_cutadapt_alone(input, output, threads=7,
                                        parse= list( cutadapt="-j 7 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45",
                                        fastq_quality_filter="-q 20 -p 80"))
    fn <- list.files(output, full.names=TRUE, recursive=FALSE)
    file.remove(fn) },

  replications = 9, columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self", "user.child", "sys.child"))


##########################################################################
##### Plot benchmark results

res_bench <- read.delim("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/long_benchmark.csv", header=TRUE, sep=",")

library(ggplot2)

ggplot(res_bench, aes(x=treatment, y=s, fill=treatment)) +
   							stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
							stat_summary(geom = "bar", fun = mean, position = "dodge") +
							geom_jitter(col = "black", size=4, position=position_jitter(w=0.2, h=0.2)) +
                            scale_fill_manual(values=c("#094A6B", "grey", "blue"))+
							scale_color_manual(values=c("black", "black", "black"))+
							theme_bw()+
							theme(legend.position="none", axis.line.x =element_blank())


##########################################################################
##### Kang et al with Illumina adaptor ####### 
input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq/"

output = "/home/danis31/Desktop/Temp_docs/temp/test/cutadapt"
prog_report  <- make_cutadapt(input, output, threads=7,
                              parse= list(cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45 -e 0.1",
                                          fastq_quality_filter="-q 20 -p 80"))

output = "/home/danis31/Desktop/Temp_docs/temp/test/seqpac"
prog_report  <- make_trim(input=input, output=output, threads=7, indels=TRUE, concat=12, check_mem=FALSE,
                          adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1),
                          adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTA",
                          polyG=c(type="hard_trim", min=20, mismatch=0.1),
                          seq_range=c(min=14, max=45), quality=c(threshold=20, percent=0.8))


##########################################################################
######### Seqpac drosophila NEBNext adator dataset ############## 
input <- system.file("extdata", package = "seqpac", mustWork = TRUE)

output = "/home/danis31/Desktop/Temp_docs/temp/test2/cutadapt/"
prog_report <-  make_cutadapt(input, output, threads=7,
                              parse= list( cutadapt="-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCA --discard-untrimmed --nextseq-trim=20 -O 10 -m 12 -M 45 -e 0.1",
                                           fastq_quality_filter="-q 20 -p 80"))

output = "/home/danis31/Desktop/Temp_docs/temp/test2/seqpac"
prog_report  <- make_trim(input=input, output=output, threads=7, indels=TRUE, concat=12,
                          check_mem=FALSE, adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1),
                          adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                          polyG=c(type="hard_trim", min=20, mismatch=0.1),
                          seq_range=c(min=14, max=45), quality=c(threshold=20, percent=0.8))



##################################################################

input_seq <- "/home/danis31/Desktop/Temp_docs/temp/test2/seqpac"
input_cut <- "/home/danis31/Desktop/Temp_docs/temp/test2/cutadapt"
#input_orig <- "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq"
input_orig <- system.file("extdata", package = "seqpac", mustWork = TRUE)

fls_seq <- list.files(input_seq, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
fls_cut <- list.files(input_cut, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)
fls_orig  <- list.files(input_orig, pattern ="fastq.gz\\>|\\.fastq\\>", full.names=TRUE, recursive=FALSE, include.dirs = FALSE)

fastq_lst <- list(NULL)
fastq_lst[[1]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_seq[[1]], withIds=FALSE)))
fastq_lst[[2]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_cut[[1]], withIds=FALSE)))
fastq_lst[[3]] <- paste0(ShortRead::sread(ShortRead::readFastq(fls_orig[[1]], withIds=FALSE)))

table(unique(fastq_lst[[1]]) %in% fastq_lst[[2]])
table(unique(fastq_lst[[2]]) %in% fastq_lst[[1]])

seq <- fastq_lst[[1]][ !fastq_lst[[1]] %in% unique(fastq_lst[[2]])]
cut <- fastq_lst[[2]][ !fastq_lst[[2]] %in% unique(fastq_lst[[1]])]

head(cut)
head(seq)

tar_seq <- "AACACTAAGCGGTGGATCACTCGGC"
fastq_lst[[1]][fastq_lst[[1]] == tar_seq]
fastq_lst[[2]][fastq_lst[[2]] == tar_seq]

fastq_lst[[1]][grepl(paste0("^", tar_seq), fastq_lst[[1]])]
fastq_lst[[2]][grepl(paste0("^", tar_seq), fastq_lst[[2]])]

fastq_lst[[3]][grepl(paste0("^", tar_seq), fastq_lst[[3]])]



tar_seq <- paste0("TATTTGAACGAGAAACCTGTAACCAACTCTCAACTGATGAGATCGAAGAGCACACGTCTGAACTCCAGTCACTAG")
which(grepl(tar_seq, fastq_lst[[3]])) #7043

AACACTAAGCGGTGGATCACTCGGC
AACACTAAGCGGTGGATCACTCGGCAGAACGGAAGAGCACACGACTGAACACCAGACACTAGCTTAACTCGTATG
AGATCGGAAGAGCACACGTCTGAACTCCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
##########################################################################
##### Kang et al count table ####### 

input = "/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/fastq/"

##########
## Cutadapt

parse = list(cutadapt="-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --discard-untrimmed --nextseq-trim=20 -O 10 -m 14 -M 45 -e 0.1",
             fastq_quality_filter="-q 20 -p 80")

counts_cut  <- make_counts(input, threads=7, parse=parse, type="fastq", trimming="cutadapt", plot=TRUE, evidence=c(experiment=2, sample=1))
anno_cut  <- make_anno(counts_cut, type="counts")
pheno_cut <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace",
                        type="manual", counts=counts_cut[[1]], progress_report=counts_cut[[2]])

pac_cut <- make_PAC(pheno_input=pheno_cut, anno_input=anno_cut, counts_input=counts_cut[[1]]) PAC_check(pac_cut)
# save(pac_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")

#########
## Seqpac
parse = list(adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1),
             adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
             polyG=c(type="hard_trim", min=20, mismatch=0.1),
             seq_range=c(min=14, max=45), quality=c(threshold=20, percent=0.8),
             indels=TRUE, concat=12, check_mem=FALSE)

counts_seq  <- make_counts(input, threads=7, parse=parse, type="fastq",
                           trimming="seqpac", plot=TRUE, evidence=c(experiment=2, sample=1))
anno_seq <- make_anno(counts_seq, type="counts")
pheno_seq <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace",
                        type="manual", counts=counts_seq[[1]], progress_report=counts_seq[[2]])

pac_seq <- make_PAC(pheno_input=pheno_seq, anno_input=anno_seq, counts_input=counts_seq[[1]]) PAC_check(pac_seq)
# save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")


###################################################################################################################################
##########################################################################
##### Kang et al reanno genome  #######
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
load("/data/Data_analysis/Projeimportcts/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")

output_path = "/home/danis31/Desktop/Temp_docs/temp"

ref_paths_gen <- list(dm6="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/fasta/chr/fast_chr.fa",
                      hg38="/data/Data_analysis/Genomes/Humans/Ensembl/CRCh38.101/Homo_sapiens.GRCh38.dna.toplevel.fa")

### seqpac
map_reanno(pac_seq, type = "internal", output_path, ref_paths = ref_paths_gen,
           mismatches = 3, threads = 7, import="genome")
reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_seq, mis_fasta_check = TRUE, threads = 7)
#save(reanno_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_seqpac.Rdata")
#load(file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_seqpac.Rdata")

pac_seq <- add_reanno(reanno=reanno_seq, mismatches = 3,  merge_pac=pac_seq, type = "genome", genome_max = 10)
save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")


# Sequences in total: 491968 #
# mis0	dm6:  	Ref_hits= 57901     (11.77%)
# mis0	hg38: 	Ref_hits= 229082    (46.56%)
# mis1	dm6:  	Ref_hits= 22014     (4.47%)
# mis1	hg38: 	Ref_hits= 88590     (18.01%)
# mis2	dm6:  	Ref_hits= 13616     (2.77%)
# mis2	hg38: 	Ref_hits= 54194     (11.02%)
# mis3	dm6:    Ref_hits= 15112     (3.07%)
# mis3	hg38: 	Ref_hits= 36839     (7.49%)

### cutadapt
map_reanno(pac_cut, type = "internal", output_path, ref_paths=ref_paths_gen,
           mismatches = 3, threads = 7, import="genome")
reanno_cut <- make_reanno(reanno_path=output_path, PAC=pac_cut, mis_fasta_check = TRUE, threads = 7)

#save(reanno_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_cutadapt.Rdata")
#load(file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_genome_cutadapt.Rdata")

pac_cut <- add_reanno(reanno=reanno_cut, mismatches = 3,  merge_pac=pac_cut, type = "genome", genome_max = 10)
#save(pac_cut, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadapt.Rdata")

# Sequences in total: 479491
# mis0	dm6:  	Ref_hits= 57167     (11.92%)
# mis0	hg38: 	Ref_hits= 223113    (46.53%)
# mis1	dm6:  	Ref_hits= 21765     (4.54%)
# mis1	hg38: 	Ref_hits= 86814     (18.11%)
# mis2	dm6:  	Ref_hits= 13410     (2.8%)
# mis2	hg38: 	Ref_hits= 53352     (11.13%)
# mis3	dm6:    Ref_hits= 14955     (3.12%)
# mis3	hg38: 	Ref_hits= 36183     (7.55%)


###################################################################################################################################
##########################################################################
##### Kang et al reanno genome  #######
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac.Rdata")
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadpt.Rdata")

output_path = "/home/danis31/Desktop/Temp_docs/temp"

ref_paths_bio <- list(dm6_miRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/miRNA/miRBase_21-dme.fa",
                      dm6_Ensembl="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/Ensembl/Drosophila_melanogaster.BDGP6.ncrna.fa",
                      dm6_tRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/tRNA_reanno/tRNA_mature.fa",
                      dm6_piRNA="/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.fa",

                      hg38_miRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/miRBase_21/miRBase_21-hsa.fa",
                      hg38_Ensembl="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/Ensembl/Homo_sapiens.GRCh38.ncrna.fa",
                      hg38_tRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/GtRNAdb/hg19-tRNAs_mature.fa",
                      hg38_piRNA="/data/Data_analysis/Genomes/Humans/pirBase/hsa.fa")


bio_search <- list( dm6_Ensembl=c("lincRNA", "miRNA", "rRNA", "pre_miRNA", "snoRNA", "snRNA", "_tRNA", "mt_tRNA"),
                    dm6_miRNA="dme-",
                    dm6_piRNA = "piR-dme|piRNA",
                    dm6_tRNA =c("^tRNA", "mt-tRNA"),

                    hg38_miRNA="hsa",
                    hg38_Ensembl=c("lncRNA", "miRNA", "rRNA", "snoRNA", "snRNA", "tRNA"),
                    hg38_tRNA=c("tRNA"),
                    hg38_piRNA="piR-hsa|piRNA")

hierarchy <- list(rRNA="rRNA", tRNA="tRNA", miRNA ="miRNA", snoRNA="snoRNA", snRNA="snRNA", lncRNA="lncRNA|lincRNA", piRNA="piRNA" )

hierarchy <- list(piRNA="piRNA", rRNA="rRNA", tRNA="tRNA", miRNA ="miRNA", snoRNA="snoRNA", snRNA="snRNA", lncRNA="lncRNA|lincRNA")


### seqpac
#map_reanno(pac_seq, type = "internal", output_path, ref_paths=ref_paths_bio, mismatches = 3, threads = 7, import="biotype")
#reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_seq, mis_fasta_check = TRUE, threads = 7)
#save(reanno_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_biotype_seqpac.Rdata")
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/reanno_biotype_seqpac.Rdata")

pac_seq <- add_reanno(reanno_seq, bio_search=bio_search, type="biotype", bio_perfect=FALSE, mismatches = 0, merge_pac=pac_seq)
pac_seq <- simplify_reanno(pac_seq, hierarchy=hierarchy, mismatches = 0, bio_name = "Biotype_hir2", merge_pac = TRUE)


#save(pac_seq, file="/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac_anno.Rdata")


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
##########################################################################
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac_anno.Rdata")
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_cutadapt.Rdata")

## Wenn-diagram
input <- list(seqpac=rownames(pac_seq$Anno), cutadapt=rownames(pac_cut$Anno))
PAC_venn <- function(input, type="venn"){
  n_input <- length(input)
  nams <- names(input)
  names(input) <-LETTERS[1:n_input]
  all_seqs <- unique(unlist(input))
  # For sav logi table:
  comb_logi <- do.call("cbind", lapply(1:n_input, function(x){
    logi <- all_seqs %in% input[[x]] }))
  colnames(comb_logi) <- nams 
  comb_logi <- cbind(data.frame(seqs=all_seqs, comb_logi))
  # For venn diagram:
  comb_ven <- do.call("cbind", lapply(1:n_input, function(x){
    logi <- all_seqs %in% input[[x]]
    ifelse(logi, names(input)[x], NA)
    }))
  comb_ven <- table(apply(comb_ven, 1, function(x){
    paste(x , collapse="&")
    }))
  names(cnts) <- gsub("NA&|&NA", "", names(cnts))
  intg <- as.integer(cnts)
  names(intg) <- names(cnts)
  if(type=="euler"){
    ven <- venneuler::venneuler(intg)
    ven$labels <- nams
    return(list(venn=ven, olap=comb_logi))
    }
  if(type=="venn"){
    ven <- ggvenn::ggvenn(input)
  }
  return(list(venn=ven, olap=comb_logi))
  }

## Nbias plots seqpac
nbias_seq <- PAC_nbias(pac_seq)
cowplot::plot_grid(plotlist=nbias_seq$Histograms)
bias_cut <- PAC_nbias(pac_cut) cowplot::plot_grid(plotlist=nbias_cut$Histograms)

## Size dist seqpac 
ord <- c("no_anno", "other", "miRNA", "tRNA", "snoRNA", "piRNA" , "snRNA", "lncRNA", "rRNA")
sizebio_seq <- PAC_sizedist(pac_seq, anno_target=list("Biotype", ord))
sizebio_seq <- PAC_sizedist(pac_seq, anno_target=list("Biotype_hir2", ord))
cowplot::plot_grid(plotlist=sizebio_seq$Histograms)


## Species genome mapping pie
pac_seq$Anno$mis0_genome_olap <- paste(ifelse(pac_seq$Anno$dm6 == "mis0", "fly", "no"), ifelse(pac_seq$Anno$hg38 == "mis0", "human", "no"))
plts_mis0 <- PAC_pie(pac_seq, pheno_target=list("Sample_ID"), anno_target=list("mis0_genome_olap"))
cowplot::plot_grid(plotlist=plts_mis0, nrow=3, ncol=3, scale=0.8)


pac_seq$Anno$mis3_genome_olap <- paste(ifelse(pac_seq$Anno$dm6 %in% c("mis0","mis1", "mis2", "mis3"), "fly", "no"),
                                       ifelse(pac_seq$Anno$hg38 %in% c("mis0","mis1", "mis2", "mis3"), "human", "no"))

plts_mis3 <- PAC_pie(pac_seq, pheno_target=list("Sample_ID"),
                     anno_target=list("mis3_genome_olap", c("no no", "fly no", "no human", "fly human")))

cowplot::plot_grid(plotlist=plts_mis3, nrow=3, ncol=3, scale=0.8)


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
########################################################################## ##
# PAC_tRNA
library(seqpac)

#######################
## Fix tRNA reference
file_path <- paste0("/data/Data_analysis/Genomes/Humans/GtRNAdb/hg38/hg38-tRNAs/hg38-tRNAs.fa")
trna <- Biostrings::readDNAStringSet(filepath=file_path, format="fasta")

names(trna) <- gsub("Homo_sapiens_", "", names(trna))              # Remove species
mat <- do.call("rbind", strsplit(names(trna), " "))                           # Make a name matrix
names(trna) <-  paste(mat[,1], mat[,ncol(mat)-1], mat[,ncol(mat)], sep="_")   # Save the important as one single string
file_path <- paste0("/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/Ensembl/Homo_sapiens.GRCh38.ncrna.fa")

ncrna <- Biostrings::readDNAStringSet(filepath=file_path, format="fasta")
mat <- do.call("rbind", strsplit(names(ncrna), " "))                        # Make a matrix of the names
mat <- mat[,1:7]

trna_logi <- grepl("Mt_tRNA_MT", names(ncrna)) # Locate mito tRNA
table(trna_logi) # Should be 22
mt_trna_ensembl <- ncrna[trna_logi,]
mt_trna <- gsub("Mt_tRNA:", "MT-tRNA-", names(mt_trna_ensembl))
mt_trna <- gsub(":-1$", "_(-)", mt_trna)
mt_trna <- gsub(":1$", "_(+)", mt_trna)
mt_trna <- gsub("_ENST", "-ENST", mt_trna)
mt_trna <- gsub(":", "-", mt_trna)
mt_trna <- gsub("_Mt_tRNA_MT-", "_chrMT:", mt_trna)

names(mt_trna_ensembl) <- mt_trna
trna <- c(trna, mt_trna_ensembl)
Biostrings::writeXStringSet(trna, filepath="/data/Data_analysis/Genomes/Humans/GtRNAdb/hg38/hg38-tRNAs/hg38-tRNAs_seqpac_fixed.fa", format="fasta")

## Make bowtie index
ref_path <- "/data/Data_analysis/Genomes/Humans/GtRNAdb/hg38/hg38-tRNAs/hg38-tRNAs_seqpac_fixed.fa"
Rbowtie::bowtie_build(ref_path, outdir=dirname(ref_path), prefix="hg38-tRNAs_seqpac_fixed", force=TRUE)


##########################################################
#### PAC generation Cancer cell lines dataset
path_to_fastq <- "/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/fastq/ena_files"

# Need to run "hard_trim" (not default_NEB) since only 50 nt reads
parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                        adapt_3_set=c(type="hard_trim", min=10, mismatch=0.1),
                        adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                        seq_range=c(min=14, max=70),
                        quality=c(threshold=20, percent=0.8),
                        indels=TRUE, concat=12, check_mem=FALSE)

count_list <- make_counts(input=path_to_fastq, type = "fastq", trimming="seqpac", 
                          parse=parse, threads=12, save_temp = TRUE)
anno <- make_anno(count_list, type = "counts", stat = TRUE)
pheno <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/pheno.csv",
                    counts=count_list$counts,
                    progress_report=count_list$progress_report)
pac  <- make_PAC(pheno_input=pheno, anno_input=anno, counts_input=count_list$counts)
table(pac$Pheno$cell_type)

# Check nbias
nbias_out <- PAC_nbias(pac)
nbias_out$Histograms[[1]]
cowplot::plot_grid(plotlist=nbias_out$Histograms[1:12], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[13:24], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[25:36], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[37:42], nrow=3, ncol=4)

#save(pac, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_master.Rdata")
#load(file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_master.Rdata")

pheno_ord <- c("C33A", "C41", "Caski", "HeLa", "HT-3", "SiHa", "SW756", "SCC-090", "SCC-1", "SCC-154", "SCC-4", "SCC-47", "UPCI-017", "UPCI-068")

pac_filt  <- PAC_filter(pac, size = c(16,50), threshold=10, coverage=50, pheno_target=list("cell_type", pheno_ord), stat=TRUE)
pac_filt  <- PAC_norm(pac_filt, norm = "cpm")
pac_filt  <- PAC_norm(pac_filt, norm = "vst")

pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("cell_type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), label=pac_filt$Pheno$cell_type, graphs=TRUE)

#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")


## Different filters
# pac_filt_hi <- PAC_filtsep(pac_filt, norm="counts", coverage=100, threshold=10, pheno_target=list("cell_type"))
# pac_filt_hi <- unique(do.call("c", as.list(pac_filt_hi)))
# pac_filt_hi <- PAC_filter(pac_filt, subset_only = TRUE, anno_target= pac_filt_hi) 
# pca_out <- PAC_pca(pac_filt_hi, norm="vst", pheno_target=list("type"), graphs=TRUE)
# 
# temp <- PAC_filtsep(pac, norm="counts", coverage=100, threshold=10, pheno_target=list("cell_type"))
# temp <- unique(do.call("c", as.list(temp)))
# temp <- PAC_filter(pac, subset_only = TRUE, anno_target= temp) 
# temp  <- PAC_norm(temp, norm = "vst")
# pca_out <- PAC_pca(temp, norm="vst", pheno_target=list("type"), graphs=TRUE)

## Check highly expressed rRNA
pac_filt$Counts[rownames(pac_filt$Counts) == "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGC",] # rRNA
pac_filt$Anno[rownames(pac_filt$Anno) == "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGC",]
temp <- reanno_seq$Full_anno
lapply(temp, function(x){x$hg38_Ensembl[x$hg38_Ensembl$seq == "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGC",]})


###########################################################
## Generate Reanno Cancer cell lines dataset
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")
output_path = "/home/danis31/Desktop/Temp_docs/temp"

### Genome
ref_paths_gen <- list(hg38="/data/Data_analysis/Genomes/Humans/Ensembl/CRCh38.101/Homo_sapiens.GRCh38.dna.toplevel.fa",
                      M.hyorhinis="/data/Data_analysis/Genomes/Other/Mycoplasma hyorhinis ATCC/GCF_000383515.1_ASM38351v1_genomic.fa")

map_reanno(pac_filt, type = "internal", output_path, ref_paths = ref_paths_gen,
           mismatches = 3, threads = 12, import="genome")
reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_filt, mis_fasta_check = TRUE)
pac_filt <- add_reanno(reanno=reanno_seq, mismatches = 3,  merge_pac=pac_filt, type = "genome", genome_max = 10)

pac_filt$Anno[ rownames(pac_filt$Anno) == "TCCCTGTGGTCTAGTGGTTAGGATTCGGCGCA",]

# Rbowtie::bowtie_build("/data/Data_analysis/Genomes/Other/Mycoplasma orale/GCF_900660435.1_50465_D02-3_genomic.fa", 
#                       outdir="/data/Data_analysis/Genomes/Other/Mycoplasma orale/",
#                       prefix="GCF_900660435.1_50465_D02-3_genomic", force=TRUE)


### Biotype
ref_paths_bio <- list(hg38_miRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/miRBase_21/miRBase_21-hsa.fa",
                      hg38_Ensembl="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/Ensembl/Homo_sapiens.GRCh38.ncrna.fa",
                      hg38_tRNA="/data/Data_analysis/Genomes/Humans/Sports_2019-06-19/GtRNAdb/hg19-tRNAs_mature.fa",
                      hg38_piRNA="/data/Data_analysis/Genomes/Humans/pirBase/hsa.fa")
bio_search <- list( hg38_miRNA="hsa",
                    hg38_Ensembl=c("lncRNA", "miRNA", "rRNA", "snoRNA", "snRNA", "tRNA"),
                    hg38_tRNA=c("tRNA"),
                    hg38_piRNA="piR-hsa|piRNA")
hierarchy <- list(rRNA="rRNA", tRNA="tRNA", miRNA ="miRNA", snoRNA="snoRNA", snRNA="snRNA", lncRNA="lncRNA|lincRNA", piRNA="piRNA" )
map_reanno(pac_filt, type = "internal", output_path, ref_paths=ref_paths_bio, mismatches = 3, threads = 12, import="biotype")
reanno_seq <- make_reanno(reanno_path=output_path, PAC=pac_filt, mis_fasta_check = TRUE)

pac_filt <- add_reanno(reanno_seq, bio_search=bio_search, type="biotype", bio_perfect=FALSE, mismatches = 3, merge_pac=pac_filt)
pac_filt <- simplify_reanno(pac_filt, hierarchy=hierarchy, mismatches = 0, bio_name = "Biotype_mis0", merge_pac = TRUE)
pac_filt <- simplify_reanno(pac_filt, hierarchy=hierarchy, mismatches = 3, bio_name = "Biotype_mis3", merge_pac = TRUE)

#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")


###########################################################
## gtf type annotation of rRNA
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")

## Use refseq gtf (more RNA annotations)
# Found at: refSeq https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/

# gtf <- tibble::as_tibble(rtracklayer::readGFF("/data/Data_analysis/Genomes/Humans/NCBI/GRCh38_latest_genomic.gtf"))
# gtf_genes <- gtf[gtf$type == "gene",]
# gtf_rrna <- gtf_genes[grepl("^RNA", gtf_genes$description),]
# gtf_rrna <- gtf_rrna[grepl("ribosomal", gtf_rrna$description),]
# table(gtf$gene_biotype)
# table(gtf_genes$gene_biotype)
# table(gtf_rrna$gene_biotype)
# gtf_rrna$description
# rtracklayer::export(gtf_rrna, "/data/Data_analysis/Genomes/Humans/NCBI/rRNA_extract.gtf", format="gtf")

## Load converter file and extract available chrom
# conv <- read.delim("/data/Data_analysis/Genomes/Humans/Convert_refSeq_UCSC_Enembl_hg38.csv", header=TRUE, sep="\t")
# conv <- conv[!nchar(as.character(conv$Ensembl))==0,]    # Removes not avaialble chromosomes in Ensembl
# identical(conv$UCSC, conv$UCSC2)                        # Check UCSC chrom columns connecting refseq with Ensembl
# table(as.character(gtf_rrna$seqid) %in% as.character(conv$refSeq)) # 53 Not avaialble
# availb <- as.character(gtf_rrna$seqid) %in% as.character(conv$refSeq)
# gtf_rrna <- gtf_rrna[availb,]
# gtf_rrna$seqid  <- as.character(conv$Ensembl[match(as.character(gtf_rrna$seqid), as.character(conv$refSeq))]) # Update seqid in 
# rtracklayer::export(gtf_rrna, "/data/Data_analysis/Genomes/Humans/NCBI/rRNA_extract_chrom_conv.gtf", format="gtf")

# Run PAC_gtf using refseq gtf
gtf_rrna <- tibble::as_tibble(rtracklayer::readGFF("/data/Data_analysis/Genomes/Humans/NCBI/rRNA_extract_chrom_conv.gtf"))
gtf_anno_mis0 <- PAC_gtf(PAC=pac_filt, genome="/data/Data_analysis/Genomes/Humans/Ensembl/CRCh38.101/Homo_sapiens.GRCh38.dna.toplevel.fa", 
          mismatches=0, gtf_other=list(rrna=gtf_rrna), target_other=list(rrna=c("gene", "gene_id", "gene_biotype", "description")), threads=12, return="simplify" )

gtf_anno_mis3 <- PAC_gtf(PAC=pac_filt, genome="/data/Data_analysis/Genomes/Humans/Ensembl/CRCh38.101/Homo_sapiens.GRCh38.dna.toplevel.fa", 
          mismatches=3, gtf_other=list(rrna=gtf_rrna), target_other=list(rrna=c("gene", "gene_id", "gene_biotype", "description")), threads=12, return="simplify" )

head(pac_filt$Anno)
table(pac_filt$Anno$mis0_bio)
identical(rownames(pac_filt$Anno), as.character(gtf_anno_mis0$seqs))
identical(rownames(pac_filt$Anno), as.character(gtf_anno_mis3$seqs))

# Update new biotype mis variable and then simplify with hierarchry
pac_filt_ext <- pac_filt

new_var <- paste(pac_filt_ext$Anno$mis0_bio, ifelse(!is.na(gtf_anno_mis0$gene_biotype), "hg38_refSeq_rRNA", ""), sep=";")
pac_filt_ext$Anno$mis0_bio <- gsub(";$", "", new_var)
     
new_var <- paste(pac_filt_ext$Anno$mis3_bio, ifelse(!is.na(gtf_anno_mis3$gene_biotype), "hg38_refSeq_rRNA", ""), sep=";")
pac_filt_ext$Anno$mis3_bio <- gsub(";$", "", new_var)

pac_filt_ext <- simplify_reanno(pac_filt_ext, hierarchy=hierarchy, mismatches = 0, bio_name = "Biotype_mis0_refseq", merge_pac = TRUE)
pac_filt_ext <- simplify_reanno(pac_filt_ext, hierarchy=hierarchy, mismatches = 3, bio_name = "Biotype_mis3_refseq", merge_pac = TRUE)

#save(pac_filt_ext, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10_refseq_rRNA.Rdata")


###########################################################
## Make sizedist biotype plots Cancer cell lines dataset
library(seqpac)
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10_refseq_rRNA.Rdata")
pac_filt <- pac_filt_ext

# Fix subgroups
pac_filt$Pheno$subtype <- ifelse(pac_filt$Pheno$cell_type %in% "SCC-154", "SCC-154",
                                 ifelse(pac_filt$Pheno$cell_type %in% "SCC-4", "SCC-4", "all_other"))
pac_filt$Pheno$subtype2 <- paste(pac_filt$Pheno$subtype, pac_filt$Pheno$type, sep="_")

pac_filt <- PAC_summary(pac_filt, norm = "cpm", type = "means", pheno_target=list("cell_type"), merge_pac = TRUE)
pac_filt <- PAC_summary(pac_filt, norm = "cpm", type = "means", pheno_target=list("type"), merge_pac = TRUE)
pac_filt <- PAC_summary(pac_filt, norm = "cpm", type = "means", pheno_target=list("subtype"), merge_pac = TRUE)

## plot size dist 
out_sizedist_mis0 <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis0"), summary_target=list("cpmMeans_cell_type"))
cowplot::plot_grid(plotlist=out_sizedist_mis0$Histograms)

out_sizedist_mis3 <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis3"), summary_target=list("cpmMeans_cell_type"))
cowplot::plot_grid(plotlist=out_sizedist_mis3$Histograms)

out_sizedist_mis0_ext <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis0_refseq"), summary_target=list("cpmMeans_cell_type"))
cowplot::plot_grid(plotlist=out_sizedist_mis0_ext$Histograms)

out_sizedist_mis3_ext <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis3_refseq"), summary_target=list("cpmMeans_cell_type"))
cowplot::plot_grid(plotlist=out_sizedist_mis3_ext$Histograms)

# Fix subtype2 mixing means of all with individual samples of SCC4 and SCC154
pac_filt <- PAC_summary(pac_filt, norm = "cpm", type = "means", pheno_target=list("subtype2", c("all_other_cell", "all_other_EXO", "all_other_MV")), merge_pac = TRUE)
names(pac_filt$summary)
head(pac_filt$summary$cpmMeans_subtype2)

temp <-  pac_filt$norm$cpm[,rownames(pac_filt$Pheno)[!grepl("all", pac_filt$Pheno$subtype2)]]
colnames(temp) <- pac_filt$Pheno$subtype2[!grepl("all", pac_filt$Pheno$subtype2)]
head(temp)
pac_filt$summary$cpmMeans_subtype2 <- cbind(pac_filt$summary$cpmMeans_subtype2, temp)

## plot size dist for subtype 2
ord <-  c("no_anno", "other",  "snoRNA",  "snRNA", "rRNA", "tRNA", "piRNA",  "miRNA")
out_sizedist_mis0_ext <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis0_refseq", ord), summary_target=list("cpmMeans_subtype2"))
cowplot::plot_grid(plotlist=out_sizedist_mis0_ext$Histograms)



## Plot percent mycoplasma
pac_filt$Anno$genome_hit <- ifelse(grepl("mis",  pac_filt$Anno$hg38_genome), "human",
                              ifelse(grepl("mis",  pac_filt$Anno$M.hyorhinis_genome), "mycoplasma",
                                     ifelse(grepl("mis",  pac_filt$Anno$Myco.orale_genome), "mycoplasma", "no_hit")))
                                     


ord <-  c("no_hit", "mycoplasma", "human")
PAC_stackbar(pac_filt, anno_target= list("genome_hit", ord),  pheno_target=list("type"), color=c("#BCBCBD", "#9D0014","#094A6B"))

#######################
#### Remove no-annotated and filter on range (Cancer cell lines dataset)
# group means
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10_refseq_rRNA.Rdata")
pac_filt <- pac_filt_ext

pac_filt$Anno$genome_hit <- ifelse(grepl("mis",  pac_filt$Anno$hg38_genome), "human",
                              ifelse(grepl("mis",  pac_filt$Anno$M.hyorhinis_genome), "mycoplasma",
                                     ifelse(grepl("mis",  pac_filt$Anno$Myco.orale_genome), "mycoplasma", "no_hit")))

pure_filt <- 
  as.data.frame(
    cbind(grepl("mis",  pac_filt$Anno$hg38_genome),
          grepl("mis",  pac_filt$Anno$M.hyorhinis_genome),
          grepl("mis",  pac_filt$Anno$Myco.orale_genome)))

pure_filt <- apply(pure_filt, 1,  function(x){paste(x, collapse="_")}) == "TRUE_FALSE_FALSE"
pac_filt$Anno$pure<- ifelse(pure_filt==TRUE, "pure_human", "not")
pac_filt  <- PAC_filter(pac_filt, size = c(16,45), anno_target=list("pure", "pure_human"))

# Test
ord <-  c("no_hit", "mycoplasma", "human")
PAC_stackbar(pac_filt, anno_target= list("genome_hit", ord),  pheno_target=list("type"), color=c("#BCBCBD", "#9D0014","#094A6B"))

# Remove old vst norm
pac_filt$norm <- pac_filt$norm[!names(pac_filt$norm) %in% "vst"]
pac_filt <- PAC_norm(pac_filt, norm = "vst")

pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="cpm", pheno_target=list("type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), label=pac_filt$Pheno$cell_type, graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="cpm", pheno_target=list("type"), label=pac_filt$Pheno$cell_type, graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="counts", pheno_target=list("type"), label=pac_filt$Pheno$cell_type, graphs=TRUE)
#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_pure.Rdata")


###########################################################
## Coverage plot of rRNA Cancer cell lines dataset
library(seqpac)
#load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_pure.Rdata") # Range filtered will not work.
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10_refseq_rRNA.Rdata")
pac_filt <- pac_filt_ext
ref <- "/data/Data_analysis/Genomes/Humans/rRNA/rRNA_Peaks1-4_45S.fa"
map_object <- PAC_mapper(pac_filt, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)


# Remove contaminated samples and make summaries
sampl <- as.character(unique(pac_filt$Pheno$cell_name))
sampl <- sampl[!sampl %in%  c("SCC-154", "SCC-4")]
pac_filt_sub <- PAC_filter(pac_filt, pheno_target=list("cell_name", sampl))
pac_filt_sub <- PAC_summary(pac_filt_sub, norm="cpm", pheno_target=list("type"))

# Plot coverage plots
out_covp_sub <- PAC_covplot(pac_filt_sub, map=map_object, summary_target=list("cpmMeans_type"), xseq = FALSE)
out_covp_sub[1]
cowplot::plot_grid(plotlist=out_covp[2:5], nrow=4, ncol=1)
out_covp_sub[2]
out_covp_sub[6]

# Only HeLa
pac_sub <- PAC_filter(pac_filt, pheno_target=list("cell_type", "HeLa"))
pac_sub$summary$cpm <- pac_sub$norm$cpm
names(pac_sub$summary$cpm) <- pac_sub$Pheno$name
out_covp <- PAC_covplot(pac_sub, map=map_object, summary_target=list("cpm"), xseq = FALSE)
out_covp[1]
cowplot::plot_grid(plotlist=out_covp[2:5], nrow=4, ncol=1)
out_covp[2]
out_covp[3]
out_covp[4]
out_covp[5]


######################################
## Patient cervical cancer Snoek et al 2018 dataset
library(seqpac)
input <- "/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/fastq"
output <- "/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/trimmed_seqpac"

# They used Illumina TrueSeq; SRA run selector says average size 50 nt (use hard_trim)
parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                        adapt_3_set=c(type="hard_trim", min=10, mismatch=0.1),
                        adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
                        seq_range=c(min=14, max=70),
                        quality=c(threshold=20, percent=0.8),
                        indels=TRUE, concat=12, check_mem=TRUE)

count_list <- make_counts(input=input, type = "fastq", trimming="seqpac", 
                          parse=parse, threads=1, save_temp = TRUE)

save(count_list, file="/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/R/counts.Rdata")
anno <- make_anno(count_list, type = "counts", stat = TRUE)
pheno <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/pheno.csv",
                    counts=count_list$counts,
                    progress_report=count_list$progress_report)
pac  <- make_PAC(pheno_input=pheno, anno_input=anno, counts_input=count_list$counts)
#save(pac, file="/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/R/pac_master.Rdata")

table(pac$Pheno$histology)
table(pac$Pheno$source_name)
table(pac$Anno$Size)

pac_filt  <- PAC_filter(pac, size = c(14,70), threshold=3, coverage=50, pheno_target=list("histology", c("CIN3", "normal")), stat=TRUE)
pac_filt  <- PAC_norm(pac_filt, norm = "cpm")
pac_filt  <- PAC_norm(pac_filt, norm = "vst")
pac_filt  <- PAC_summary(pac_filt, norm = "cpm", pheno_target=list("histology"))

nbias_out <- PAC_nbias(pac_filt, summary_target= list("cpmMeans_histology", c("CIN3", "normal")))
cowplot::plot_grid(plotlist=nbias_out$Histograms)

#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/R/pac_filt.Rdata")

library(seqpac)
load(file="/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/R/pac_filt.Rdata")

## Get fragemtns  from HeLa experiment
ref <- "/data/Data_analysis/Genomes/Humans/rRNA/rRNA_Peaks1-4_45S.fa"
map_object <- PAC_mapper(pac_filt, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)


pac_filt$summary$cpm <- pac_filt$norm$cpm

out_covp <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_histology", c("CIN3", "normal")), xseq = FALSE)
out_covp[[1]]
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), graphs=TRUE)

nchar("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA") # Main fragment

table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", rownames(pac_filt$Anno)))
table(grepl("TCGTGATCGATGTGGTGACG", rownames(pac_filt$Anno)))
table(grepl("GTGATCGATG", rownames(pac_filt$Anno))) # The fragment is completely missing


library(seqpac)
load(file="/data/Data_analysis/Projects/Human/SRA_download/Snoek_2018_Cervical_cancer_patients/R/pac_filt.Rdata")

## Get fragemtns  from HeLa experiment
ref <- "/data/Data_analysis/Genomes/Humans/rRNA/rRNA_Peaks1-4_45S.fa"
map_object <- PAC_mapper(pac_filt, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)
map_object[[2]]


out_covp <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_tissue_type", c("cervical_tumor", "normal_cervical")), xseq = FALSE)
out_covp[[1]] # Much variability 
out_covp[[2]]
out_covp[[3]]
out_covp[[4]]
out_covp[[5]]
out_covp[[6]]




######################################
## Patient cervical Xu et al dataset
library(seqpac)

#### Initial check
## Control
fls <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_XXX_cervical_patient_both_RNA_sRNA/fastq/ena_files/sRNA_seq/read_1/SRR11095741/SRR11095741_1.fastq.gz" 
fstq <- paste0(ShortRead::sread(ShortRead::readFastq(fls, withIds=FALSE))) 
table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", fstq)) # TRUE 78
table(grepl("TCGTGATCGATGTGGTGACG", fstq)) # TRUE 1423

## Tumor
# sRNA:
fls_2 <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_XXX_cervical_patient_both_RNA_sRNA/fastq/ena_files/sRNA_seq/read_1/SRR11095745/SRR11095745_1.fastq.gz"
#fls_2 <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_XXX_cervical_patient_both_RNA_sRNA/fastq/ena_files/sRNA_seq/read_2/SRR11095745/SRR11095745_2.fastq.gz" # the paired read found 0
fstq_2 <- paste0(ShortRead::sread(ShortRead::readFastq(fls_2, withIds=FALSE))) 
table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", fstq_2)) # TRUE 107
table(grepl("TCGTGATCGATGTGGTGACG", fstq_2))  # TRUE 5916
# long RNA
fls_2_long <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_XXX_cervical_patient_both_RNA_sRNA/fastq/ena_files/RNA_seq/SRR11095733/SRR11095733_1.fastq.gz"
fstq_2_long <- paste0(ShortRead::sread(ShortRead::readFastq(fls_2_long, withIds=FALSE))) 
table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", fstq_2_long)) # TRUE 0
table(grepl("TCGTGATCGATGTGGTGACG", fstq_2_long))  # TRUE 0

### NEB next adaptor:
# TGAGATGAAGCACTGTAGCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTGAAAAAAACCCCCCCGTGCCCCCGCCCGCCCCCCCCCCCACCGCCCCCTCTCCCCCCTCCCCTCCCC
#                      AGATCGGAAGAGCACACGTCTGAACTCCA

input <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/fastq/ena_files/sRNA_seq/read_1"
output <- "/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/trimmed_seqpac_read_1_sRNA"

# Only  
parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                        adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1),
                        adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                        seq_range=c(min=14, max=70),
                        quality=c(threshold=20, percent=0.8),
                        indels=TRUE, concat=12, check_mem=TRUE)

count_list <- make_counts(input=input, type = "fastq", trimming="seqpac", 
                          parse=parse, threads=3, save_temp = TRUE)

#save(count_list, file="/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/R/counts_list.Rdata")
anno <- make_anno(count_list, type = "counts", stat = TRUE)
pheno <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/pheno_sRNA_read_1.csv",
                    counts=count_list$counts,
                    progress_report=count_list$progress_report)
pac  <- make_PAC(pheno_input=pheno, anno_input=anno, counts_input=count_list$counts)
#save(pac, file="/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/R/pac_master.Rdata")

table(pac$Pheno$tissue_type)
table(pac$Anno$Size)

pac_filt  <- PAC_filter(pac, size = c(16,70), threshold=3, coverage=50, pheno_target=list("tissue_type", c("cervical_tumor", "normal_cervical")), stat=TRUE)
pac_filt  <- PAC_norm(pac_filt, norm = "cpm")
pac_filt  <- PAC_norm(pac_filt, norm = "vst")
pac_filt  <- PAC_summary(pac_filt, norm = "cpm", pheno_target=list("tissue_type"))
#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/R/pac_filt.Rdata")

nbias_out <- PAC_nbias(pac_filt, summary_target= list("cpmMeans_tissue_type", c("cervical_tumor", "normal_cervical")))
cowplot::plot_grid(plotlist=nbias_out$Histograms)

library(seqpac)
load(file="/data/Data_analysis/Projects/Human/SRA_download/Xu_et_at_2020_cervical_patient_both_RNA_sRNA/R/pac_filt.Rdata")

## Peak 1 from HeLa experiment
ref <- "/data/Data_analysis/Genomes/Humans/rRNA/rRNA_Peaks1-4_45S.fa"
map_object <- PAC_mapper(pac_filt, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)
map_object[[2]]


out_covp <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_tissue_type", c("cervical_tumor", "normal_cervical")), xseq = FALSE)
out_covp[[1]] # Much variability 
out_covp[[2]]
out_covp[[3]]
out_covp[[4]]
out_covp[[5]]
out_covp[[6]]


##  Plots error bars 
# Create new variable
logi_peak1<- rownames(pac_filt$Anno) %in% rownames(map_object[[2]]$Alignments)
table(logi_peak1)

pac_filt$Anno$rRNA_peak_1 <- ifelse(logi_peak1, "peak1", "not_peak")
pac_peaks <- PAC_filter(pac_filt, anno_target=list("rRNA_peak_1", c("peak1")))

dt_wide <- aggregate(pac_peaks$norm$cpm, list(pac_peaks$Anno$rRNA_peak_1), sum)
nam_peaks <-  dt_wide[,1]
dt_wide <- as.data.frame(t(dt_wide[,-1]))
colnames(dt_wide) <- nam_peaks 
identical(rownames(dt_wide), rownames(pac_peaks$Pheno))
dt_wide <- cbind(dt_wide, pac_peaks$Pheno$tissue_type)
colnames(dt_wide)[2] <- "type"
library(ggplot2)
ggplot(dt_wide, aes(x=type, y=peak1, fill=type)) +
   							stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
							stat_summary(geom = "bar", fun = mean, position = "dodge") +
							geom_jitter(col = "black", size=4, position=position_jitter(w=0.2, h=0.2)) +
							scale_fill_manual(values=c("#094A6B", "grey", "blue"))+
							scale_color_manual(values=c("black", "black", "black"))+
							theme_bw()+
							theme(legend.position="none", axis.line.x =element_blank())






######################################
## Patient leukemia Veneziano et al 2019 dataset
library(seqpac)

#### Initial check
table(grepl("TCGTGATCGATGTGGTGACG", fstq_2))  # TRUE 5916
# long RNA
fls <- "/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/fastq/ena_files/SRR10282933/SRR10282933.fastq.gz"
fstq <- paste0(ShortRead::sread(ShortRead::readFastq(fls, withIds=FALSE))) 
table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", fstq)) # TRUE 0
table(grepl("TCGTGATCGATGTGGTGACG", fstq))  # TRUE 4490


fls <- "/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/trimmed_fastq/SRR10282918.trim.fastq.gz"
fstq <- paste0(ShortRead::sread(ShortRead::readFastq(fls, withIds=FALSE))) 
table(grepl("CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGA", fstq)) # TRUE 0
table(grepl("TCGTGATCGATGTGGTGACG", fstq)) 

# NEBNext adaptor
# CTTCGAGGAACTGTAGGCACCATCAATCCCCCCTAAGTAAGATCGGAAGA
#                                         AGATCGGAAGAGCACACGTCTGAACTCCA

input <- "/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/fastq"
output <- "/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/trimmed_fastq"

# Only  
parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                        adapt_3_set=c(type="hard_trim", min=10, mismatch=0.1),
                        adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                        seq_range=c(min=14, max=70),
                        quality=c(threshold=20, percent=0.8),
                        indels=TRUE, concat=12, check_mem=TRUE)

count_list <- make_counts(input=input, type = "fastq", trimming="seqpac", 
                          parse=parse, threads=4, save_temp = TRUE)

count_list <- make_counts(input=output, type = "fastq", trimming=NULL, threads=3)

#save(count_list, file="/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/R/counts_list.Rdata")
anno <- make_anno(count_list, type = "counts", stat = TRUE)
pheno <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/pheno.csv",
                    counts=count_list$counts,
                    progress_report=count_list$progress_report)
pac  <- make_PAC(pheno_input=pheno, anno_input=anno, counts_input=count_list$counts)
#save(pac, file="/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/R/pac_master.Rdata")

table(pac$Pheno$Diagnosis)
table(pac$Anno$Size)

pac_filt  <- PAC_filter(pac, size = c(16,70), threshold=3, coverage=50, pheno_target=list("Diagnosis", c("Chronic lymphocytic leukemia", "Healthy")), stat=TRUE)
pac_filt  <- PAC_norm(pac_filt, norm = "cpm")
pac_filt  <- PAC_norm(pac_filt, norm = "vst")
pac_filt  <- PAC_summary(pac_filt, norm = "cpm", pheno_target=list("Diagnosis"))

nbias_out <- PAC_nbias(pac_filt, summary_target= list("cpmMeans_Diagnosis", c("Chronic lymphocytic leukemia", "Healthy")))
cowplot::plot_grid(plotlist=nbias_out$Histograms)

#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/R/pac_filt.Rdata")

library(seqpac)
load(file="/data/Data_analysis/Projects/Human/SRA_download/Veneziano_et_al_2019_PNAS_leukemia/R/pac_filt.Rdata")

## Get fragemtns  from HeLa experiment
ref <- "/data/Data_analysis/Genomes/Humans/rRNA/rRNA_Peaks1-4_45S.fa"
map_object <- PAC_mapper(pac_filt, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)


out_covp <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_Diagnosis", c("Chronic lymphocytic leukemia", "Healthy")), xseq = FALSE)
out_covp[[1]] # Much variability 
out_covp[[2]]
out_covp[[3]]
out_covp[[4]]
out_covp[[5]]

load(file="/data/Data_analysis/Projects/Human/nov_HeLa_human_201111/R/pac_filt_pure.Rdata")
map_object_pure <- PAC_mapper(pac_pure, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=FALSE)

##  Plots error bars 


# Create new variable
logi_peak1<- rownames(pac_filt$Anno) %in% rownames(map_object_pure[[2]]$Alignments)
table(logi_peak1)

pac_filt$Anno$rRNA_peak_1 <- ifelse(logi_peak1, "peak1", "not_peak")
pac_peaks <- PAC_filter(pac_filt, anno_target=list("rRNA_peak_1", c("peak1")))

dt_wide <- aggregate(pac_peaks$norm$cpm, list(pac_peaks$Anno$rRNA_peak_1), sum)
nam_peaks <-  dt_wide[,1]
dt_wide <- as.data.frame(t(dt_wide[,-1]))
colnames(dt_wide) <- nam_peaks 
identical(rownames(dt_wide), rownames(pac_peaks$Pheno))
dt_wide <- cbind(dt_wide, pac_peaks$Pheno$tissue_type)
colnames(dt_wide)[2] <- "type" 
ggplot(dt_wide, aes(x=type, y=peak1, fill=type)) +
   							stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
							stat_summary(geom = "bar", fun = mean, position = "dodge") +
							geom_jitter(col = "black", size=4, position=position_jitter(w=0.2, h=0.2)) +
							scale_fill_manual(values=c("#094A6B", "grey", "blue"))+
							scale_color_manual(values=c("black", "black", "black"))+
							theme_bw()+
							theme(legend.position="none", axis.line.x =element_blank())




















path_to_fastq <- "/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/fastq/ena_files"

# Need to run "hard_trim" (not default_NEB) since only 50 nt reads
parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                        adapt_3_set=c(type="hard_trim", min=10, mismatch=0.1),
                        adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                        seq_range=c(min=14, max=70),
                        quality=c(threshold=20, percent=0.8),
                        indels=TRUE, concat=12, check_mem=FALSE)

count_list <- make_counts(input=path_to_fastq, type = "fastq", trimming="seqpac", 
                          parse=parse, threads=12, save_temp = TRUE)
anno <- make_anno(count_list, type = "counts", stat = TRUE)
pheno <- make_pheno("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/pheno.csv",
                    counts=count_list$counts,
                    progress_report=count_list$progress_report)
pac  <- make_PAC(pheno_input=pheno, anno_input=anno, counts_input=count_list$counts)
table(pac$Pheno$cell_type)

# Check nbias
nbias_out <- PAC_nbias(pac)
nbias_out$Histograms[[1]]
cowplot::plot_grid(plotlist=nbias_out$Histograms[1:12], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[13:24], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[25:36], nrow=3, ncol=4)
cowplot::plot_grid(plotlist=nbias_out$Histograms[37:42], nrow=3, ncol=4)

#save(pac, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_master.Rdata")
#load(file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_master.Rdata")

pheno_ord <- c("C33A", "C41", "Caski", "HeLa", "HT-3", "SiHa", "SW756", "SCC-090", "SCC-1", "SCC-154", "SCC-4", "SCC-47", "UPCI-017", "UPCI-068")

pac_filt  <- PAC_filter(pac, size = c(16,50), threshold=10, coverage=50, pheno_target=list("cell_type", pheno_ord), stat=TRUE)
pac_filt  <- PAC_norm(pac_filt, norm = "cpm")
pac_filt  <- PAC_norm(pac_filt, norm = "vst")

pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("cell_type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt, norm="vst", pheno_target=list("type"), label=pac_filt$Pheno$cell_type, graphs=TRUE)

#save(pac_filt, file="/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")











###########################################################
## Coverage plot of tRNA using PAC_covplot  Cancer cell lines dataset

# Single tRNA targeting a summary dataframe
PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_source_name"), map_target="tRNA-Glu-CTC-1-1_chr1:146035692-146035763_(+)")

# Find tRNAs with many fragments
n_tRFs <- unlist(lapply(map_object, function(x){nrow(x[[2]])}))
table(n_tRFs)
names(map_object)[n_tRFs>30]
selct <- (names(map_object)[n_tRFs>30])[c(1, 8, 9, 11, 16, 22, 23, 24, 29, 30)]
cov_plt <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_source_name"), map_target=selct)
cov_plt <- PAC_covplot(pac_filt, map=map_object, summary_target= list("cpmMeans_Cell_type"), map_target=selct)
cowplot::plot_grid(plotlist=cov_plt, nrow=2, ncol=5)

trna_seqs <- unique(unlist(lapply(map_object, function(x){rownames(x$Alignments)})))
pac_filt_trna <- PAC_filter(pac_filt, anno_target=trna_seqs)

pca_out <- PAC_pca(pac_filt_trna, norm="vst", pheno_target=list("Cell_type"), graphs=TRUE)
pca_out <- PAC_pca(pac_filt_trna, norm="vst", pheno_target=list("source_name"), graphs=TRUE)



## Divided by subtype
out_sizedist_mis0_sub <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis0"), summary_target=list("cpmMeans_group")) 
cowplot::plot_grid(plotlist=out_sizedist_mis0_sub$Histograms, nrow=3, ncol=3)

## Without no_anno
test <- PAC_filter(pac_filt, anno_target=list("Any_bio", "Hit"))
test_mis0 <- PAC_sizedist(test, norm = NULL, anno_target= list("Biotype_mis0"), summary_target=list("cpmMeans_group")) 
cowplot::plot_grid(plotlist=test_mis0$Histograms)
test_mis3 <- PAC_sizedist(test, norm = NULL, anno_target= list("Biotype_mis3"), summary_target=list("cpmMeans_group")) 
cowplot::plot_grid(plotlist=test_mis3$Histograms)


head(test$Anno)

Any_bio

out_sizedist_mis3 <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis3"), summary_target=list("cpmMeans_group"))

# Remove no anno
no_anno_rm <-  as.character(unique(pac_filt$Anno$Biotype_mis0))
no_anno_rm <- no_anno_rm[!no_anno_rm=="no_anno"]
out_sizedist_mis0 <- PAC_sizedist(pac_filt, norm = NULL, anno_target= list("Biotype_mis0", no_anno_rm), summary_target=list("cpmMeans_source_name")) 
cowplot::plot_grid(plotlist=out_sizedist_mis0$Histograms)

test <- cbind(pac_filt$Anno, pac_filt$summary$cpmMeans_source_name)
test <- test[order(test$'SCC-154', decreasing = TRUE),]
head(test, n=20)
test2 <- test[order(test$'HeLa', decreasing = TRUE),]
head(test2)

# Whare are all rRNA?
test[grepl("rRNA", test$mis0_bio),]


GCCGATTTAGCTCAGCGGTAGAGCAGCTG                     
GCCGATTTAGCTCAGCGGTAGAGCAGCTGG                       
GCCG AT TTAGCTCAGC GG T AGAGCA G CT A
GCCG TC TTAGCTCAGCTGG C AGAGCA A CT G

AGATCGGAAGAGCACACGTCTGAACTCCA #NEB adapt
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC #Illumina adapt
                     
GCCGTCTTAGCTCAGCTGGCAGAGCAACTG # Top expressed
CAGTTGCTCTGCCAGCTGAGCTAAGACGGC # Revcomp top expressed 

GCCGTCTTAGCTCAGCTGGCAGAGCAACTG
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT # Top tRNA HeLa 
GCATTGGTGGTTCAGTGGTAGAATTCTCGCCT

GCCGATTTAGCTCAGCGGTAGAGCAGCTG
TGCTTGGACTACATATGGTTGAGGGTTGTA

GCCGATTTAGCTCAGCGGTAGAGCAGCTG
GCCGTCTTAGCTCAGCTGGCAGAGCAACTG
GCCGATTTAGCTCAGCGGTAGAGCAGCTGG
GCCGATTTAGCTCAGCGGTAGAGCAGCTA
GGCTCTGTAGCTCAGTAGGTAGAGCAACG
GCGCTTGTAGCTCAATTGGACAGAGTGTTTG




###########################################################
## Check rRNA fragment peak
load("/data/Data_analysis/Projects/Human/SRA_download/Kang_2018_KTH_miRTrace/R/pac_master_seqpac_anno.Rdata")
nbias_seq <- PAC_nbias(pac_seq)
cowplot::plot_grid(plotlist=nbias_seq$Histograms)

test <- cbind(pac_seq$Anno, pac_seq$Counts)
test <- test[order(test$'SRR7687078', decreasing = TRUE),]

table(rownames(pac_seq$Anno) ==  "GCCGATTTAGCTCAGCGGTAGAGCAGCTG")

head(test)

###########################################################
## In depth tRNA analysis of cancer cell line dataset
library(seqpac)
load("/data/Data_analysis/Projects/Human/SRA_download/Tong_et_al_2020_Cancer_cell_lines_exosomes/R/pac_filt_20_in_10.Rdata")
pac_trna <- pac_filt
pac_trna$Anno <- pac_trna$Anno[,1, drop=FALSE]

ref <- "/data/Data_analysis/Genomes/Humans/GtRNAdb/hg38/hg38-tRNAs/hg38-tRNAs_seqpac_fixed.fa"
map_object_mis0 <- PAC_mapper(pac_trna, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=TRUE)
map_object_mis3 <- PAC_mapper(pac_trna, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=3, threads=8, report_string=TRUE)

# Remove unmapped tRNA ref and check 3 loops criteria?
ss_file <- "/data/Data_analysis/Genomes/Humans/GtRNAdb/hg38/hg38-tRNAs/hg38-tRNAs-detailed.ss"
map_object_ss <- map_rangetype(map_object_mis0, type="ss", ss=ss_file, min_loop_width=5)
map_object_ss <-  map_object_ss[!unlist(lapply(map_object_ss, function(x){x[[2]][1,1] == "no_hits"}))]
nloop <- unlist(lapply(map_object_ss, function(x){unique(x$Alignments$n_ssloop)}))
table(nloop) # Most have 3

test <- map_rangetype(map_object_ss[!nloop==3], type="ss", ss=ss_file, min_loop_width=4) # change loop width
nloop_test <- unlist(lapply(test, function(x){unique(x$Alignments$n_ssloop)}))
table(nloop_test)  # Strange: some jumps directly from 2 to 4 loops; better remove them  
map_object_ss <- map_object_ss[nloop==3]

# Run tRNA_class function 
pac_trna <- tRNA_class(pac_trna, map=map_object_ss, terminal=5) # 2 nucleotides from start or end of full-length tRNA)
pac_trna$Anno$type <- paste0(pac_trna$Anno$decoder, pac_trna$Anno$acceptor)

# Run PAC_trna
trna_result <- PAC_trna(pac_trna, norm="cpm", filter = 100,
  join = TRUE, top = 15, log2fc = TRUE,
  pheno_target = list("cell_name", c("HeLa", "SCC4")), 
  anno_target_1 = list("type"),
  anno_target_2 = list("class"))


cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Grand_means,
                   trna_result$plots$Log2FC_Anno_1,
                   trna_result$plots$Percent_bars$Grand_means,
                   nrow=1, ncol=3)

# By setting on subtype
pac_trna$Pheno$subtype <- ifelse(pac_trna$Pheno$cell_name %in% c("SCC154", "SCC4"), "subtype_1", "subtype_2")

trna_result <- PAC_trna(pac_trna, norm="cpm", filter = 100,
  join = TRUE, top = 15, log2fc = TRUE,
  pheno_target = list("subtype"), 
  anno_target_1 = list("type"),
  anno_target_2 = list("class"))
cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Grand_means,
                   trna_result$plots$Log2FC_Anno_1,
                   trna_result$plots$Percent_bars$Grand_means,
                   nrow=1, ncol=3)


# tRNA pca
out_pca <- PAC_pca(pac_trna, norm="vst", pheno_target=list("vesicle_type")) # Cells are very different
out_pca <- PAC_pca(pac_trna, norm="vst", pheno_target=list("cell_name")) # But then SCC154 and SCC4 stands out
out_pca <- PAC_pca(pac_trna, norm="vst", pheno_target=list("subtype"))

temp <- PAC_filter(pac_trna, anno_target=list("class", c("3'-tRF", "3'-half")))
temp <- PAC_filter(pac_trna, pheno_target=list("vesicle_type", c("EXO", "MV")))
temp <- PAC_filter(pac_trna, pheno_target=list("vesicle_type", c("cell")))
out_pca <- PAC_pca(temp, norm="vst", pheno_target=list("subtype"))



out_pca <- PAC_pca(temp, norm="vst", pheno_target=list("cell_name"))



# pca using all
out_pca <- PAC_pca(pac_filt, norm="vst", pheno_target=list("vesicle_type")) # Cells are very different

out_pca <- PAC_pca(pac_filt, norm="vst", pheno_target=list("vesicle_type"), labels=pac_filt$Pheno$cell_name)
out_pca <- PAC_pca(pac_filt, norm="vst", pheno_target=list("cell_name")) # But then SCC154 and SCC4 stands out

temp <- PAC_filter(pac_filt, pheno_target=list("vesicle_type", c("EXO", "MV")))
test <- cbind(temp$Anno, temp$summary$cpmMeans_source_name)
row_nam <- rownames(test[order(test$'SCC-154', decreasing = TRUE),])[1]  # strange RNA
temp$Pheno$temp <- as.numeric(temp$norm$vst[row_nam,])
pca_res <- FactoMineR::PCA(t(as.matrix(temp$norm$vst)), graph=FALSE)
dat <- data.frame(pc1=pca_res$ind$coord[, 1],
                  pc2=pca_res$ind$coord[, 2],
                  pc3=pca_res$ind$coord[, 3])
dat$col  <- temp$Pheno$temp

ggplot(data = dat, aes(x = pc2, y = pc3, color = col)) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_vline(xintercept = 0, lty = 2) +
      geom_point(alpha = 0.8, size = 2)

      

