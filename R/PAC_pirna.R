#' piRNA analysis of PAC object
#'
#' Analyzing and plotting piRNAs between groups.
#'
#' Piwi interacting RNA are often classified according to what type of genomic
#' element they silence. Sileincing targets may differ from cell-type to
#' cell-type, which makes piRNA analysis challenging. Given a PAC object, this
#' function automates the mapping and annotation of piRNA.
#' 
#' @family PAC analysis
#'   
#' @seealso \url{https://github.com/Danis102} for updates.
#'
#' @param PAC PAC-list object. The Anno object needs to contain a piRNA column
#'   generated during the \emph{PAC reannotation} workflow
#'   (\code{map_reanno},\code{import_reanno},\code{make_reanno},
#'   \code{add_reanno}.
#' 
#' @param mismatches 
#' 
#' @param pirna
#' 
#' @param genome
#' 
#' @param fasta_list
#' 
#' @param gtf_list
#' 
#' @param nuc_bias Integer vector indicating which nucleotide to perform
#'   nucluotide bias analysis on. As default, nuc_bias = c(1,10). Nucleotide
#'   bias analysis is then performed on the 1st and 10th nucleotide.
#'   
#' @param graphs Logical whether graphs should be reported or not.
#'  
#' @param threads Integer indiciating the number of parallel processes.   
#'     
#' @return List of ggplot2 plots and the data used for generating the plots. Use
#'   ls.str() to explore each level.
#'   
#' @examples
#' 
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' 
#' #######################################################################
#' ### Convert of UCSC to Ensembl (Tanks to Devon Ryan and co-workers) ### 
#' 
#' # Go to https://github.com/dpryan79/ChromosomeMappings
#' # Locate your genome of choice
#' # Download from raw.githubusercontent:
#' dm6_conv <- readr::read_tsv("https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/BDGP6_UCSC2ensembl.txt", col_names=FALSE)
#' names(dm6_conv) <- c("ucsc", "ensembl")
#' 
#' # Load local UCSC formated file and make a new converted vector 
#' gtf_repeat <- "/home/danis31/Desktop/Temp_docs/fasta/repeatMasker/repeatMasker.gtf"
#' repeat_gtf <- rtracklayer::readGFF(gtf_repeat)
#' convec <- dm6_conv$ensembl[match(as.character(repeat_gtf$seqid),   as.character(dm6_conv$ucsc))]
#' 
#' # Test that the convercsion went well
#' test <- cbind(convec, as.character(repeat_gtf$seqid))
#' test2 <- unique(paste( test[,1], test[,2], sep="|"))
#' table(table(do.call("rbind", strsplit(test2, "\\|"))[,1]) == 1) # All TRUE
#'
#' # Exchange the chromosomal names, convert to genomic ranges and save your fil as gtf 
#' repeat_gtf$seqid <- convec
#' gr <- GenomicRanges::GRanges(seqnames=repeat_gtf$seqid, IRanges::IRanges(repeat_gtf$start, repeat_gtf$end), strand=repeat_gtf$strand)
#' GenomicRanges::mcols(gr) <- data.frame(type="repeat", source="repeatMasker_dm6_ucsc",
#'                                        repName = repeat_gtf$repName, repClass = repeat_gtf$repClass,
#'                                        repFamily = repeat_gtf$repFamily)
#' rtracklayer::export(gr, file.path(dirname(gtf_repeat), "repeatMasker_ensembl.gtf"), format="gtf")
#' 
#' 
#' ################################################## 
#' ### Prepare pirna_meta from pirbase data file  ###
#' # - Download pirbase data from pirbase. Make sure you get a table where
#' #   Pubmed ID for each unique experiment is specified in a seperate column. 
#' # - You may also use your own file, with a custom 2nd column. 
#' # - Resulting data.table should have 2 columns named: "name", "n_evidence". 
#' # - The name column must either match fasta reference names or sequence names in PAC.
#' # - The 2nd column will be embedded in PAC_pirna output.   
#' # - Noter: Reading and preparing pirbase data may take several min. 
#'   
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Humans/pirBase/hg38/piR_hsa.txt.gz"), delim="\t", col_names = TRUE)
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.txt.gz"), delim="\t", col_names = TRUE)
#' 
#' datset <- strsplit(pirbase_dat$pubmed, " ")
#' n_evidence  <- unlist(lapply(datset, length))
#' pirna_meta <- tibble::tibble(name=pirbase_dat$name, n_evidence=n_evidence)
#' table(pirna_meta$n_evidence)
#' 
#' ##############################################################################
#' ### Get repeatMasker table and manually turn it into gtf using rtracklayer ###
#' # Table names can be found at:
#' # https://genome.ucsc.edu/cgi-bin/hgTables
#' 
#' dest_path <- file.path(ref_path, "/repeatMasker/repeatMasker.gtf") # Full file path
#' if(!file.exists(dirname(dest_path))){dir.create(dirname(dest_path))}
#' session <- rtracklayer::browserSession("UCSC")
#' rtracklayer::genome(session) <- "dm6"
#' rm_tab <- tibble::as_tibble(rtracklayer::getTable(rtracklayer::ucscTableQuery(session, track="RepeatMasker", table="rmsk")))
#' gr <- GenomicRanges::GRanges(seqnames=rm_tab$genoName, IRanges::IRanges(rm_tab$genoStart, rm_tab$genoEnd), strand=rm_tab$strand)
#' GenomicRanges::mcols(gr)$type <- "repeat"
#' GenomicRanges::mcols(gr)$source <- "repeatMasker_dm6"
#' GenomicRanges::mcols(gr)$repName <- rm_tab$repName
#' GenomicRanges::mcols(gr)$repClass <- rm_tab$repClass
#' GenomicRanges::mcols(gr)$repFamily <- rm_tab$repFamily  
#' rtracklayer::export(gr, dest_path, format="gtf")
#' 
#' 
#' #####################
#' ### Run PAC_pirna ###
#' # 
#' library(seqpac)
#' load("/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/PAC_merge_10counts100.Rdata")
#' load(system.file("extdata", "drosophila_sRNA_pac_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' pirna <- "/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.fa"
#' genome <- "/home/danis31/Desktop/Temp_docs/fasta/biomartr_genome/chromosomes.fa"
#' pirna <- "pirna"
#' genome <- "mis0_chromosomes_genome"
#'
#' pirna_meta=NULL
#'
#' gtf_repeat <- "/home/danis31/Desktop/Temp_docs/fasta/repeatMasker/repeatMasker_ensembl.gtf"
#' gtf_protein <- "/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/gtf/Drosophila_melanogaster.gtf"
#' 
#' PAC=pac
#' threads=10
#' mismatches=3
#' norm="rpm"
#' pheno_target= list("type")
#' cluster_gap=1000
#' pirna_meta=pirna_meta
#' 
#' #export

PAC_pirna <- function(PAC, norm="counts", genome=NULL, mismatches=3, stranded=FALSE, cluster=NULL,
                      pheno_target=NULL, pirna=NULL, pirna_meta=NULL, gtf_repeat=NULL, gtf_protein=NULL,
                      threads=1){

##### Setup general ####################################
  seqs <- rownames(PAC$Anno)
  
  if(!is.null(pirna_meta)|!is.null(cluster)){
    cat("\nRunning this function using a 'pirna_meta' input against \nfasta reference IDs or using a 'cluster' input may reduce \nperformance since a full piRNA mapping report must be generated.\nIf this becomes a problem reduce the number of sequences in your\nPAC object or run 'pirna_meta=NULL' and 'cluster=NULL' (default).\n") 
    rprt <- "full"
  }else{
    rprt <- "minimum" 
  }

##### Setup pirna annotation and run reanno if necessary ####################################
  # If user know the columns:
  if(!file.exists(pirna) & !is.null(pirna)){
   anno_pirna <- PAC$Anno[, pirna, drop=FALSE]
  # If user do not know the columns:
  }else{
    # If user don't provide fasta reference:
    if(!file.exists(pirna) & is.null(pirna)){
          anno_pirna <- PAC$Anno
    }
    # If user provide fasta reference:
    if(file.exists(pirna)){
      cat("\nInput pirna was an existing file. Will treat it as a \nfasta reference and make a denovo reannotation using bowtie. \nSee ?map_reanno or ?vingette for details.\n")
      outpath <- tempfile(pattern = "", fileext = "")
      err <- try(map_reanno(PAC, ref_paths=list(pirna=pirna), output_path=outpath, type="external", mismatches=mismatches,
                     import = list(coord=FALSE, report=rprt, reduce=NULL), threads=threads), silent = TRUE)
      if(!is.null(err)){
        err2 <- try(map_reanno(PAC, ref_paths=list(pirna=pirna), output_path=outpath, type="internal", mismatches=mismatches,
                     import = list(coord=FALSE, report=rprt, reduce=NULL), threads=threads), silent = TRUE)
      if(!is.null(err2)){
        stop(paste0("\nFunction map_reanno failed. Possible reasons: \n\tNo bowtie installation\n\tBad fasta reference\n\nLast log says:\n", err2))
        }
      }
      reanno <- make_reanno(outpath, PAC=PAC, mis_fasta_check=TRUE, threads=threads)
      anno_pirna <- add_reanno(reanno, bio_search=list(pirna="pirna"), type="biotype", bio_perfect=FALSE, mismatches = mismatches)
      
      unlink(outpath, recursive = TRUE)
    }
    # Extract pirna column if user don't provide columns:
    nam <- colnames(anno_pirna)
    logi_col1 <- grepl("^piRNA|^pirna", nam)
    prefix <- suppressWarnings(do.call("rbind", strsplit(nam[logi_col1], "_"))[,1])
    if(any(duplicated(prefix))){
      stop("\nFound more than 1 possible pirna annotation.\nAutomatic column identification only works when the results \nof one fasta pirna reference is reported in PAC$Anno.\nPlease, specify in 'pirna=' the exact columns reporting the alignments \nagainst 1 reference genome (column prefix: mis0_, mis1_, mis2_ etc.) \nor proved a path to a bowtie indexed fasta pirna file.")
    }
    pirna_col <- nam[logi_col1]
    anno_pirna <- anno_pirna[,logi_col1, drop=FALSE]
    anno_pirna[anno_pirna == "_"] <- "noAlign"
    }
  anno_pirna$xyzx <- as.character(ifelse(anno_pirna[[1]] %in% c("noAlign", "_"), "not_pirna", "pirna"))
  pirna_hit_col <- paste0("pirna_hits_mis", mismatches)
  colnames(anno_pirna)[colnames(anno_pirna) =="xyzx"] <- pirna_hit_col

  ##### Handle pirna_meta  ####################################            
    if(!is.null(pirna_meta)){
      pirna_meta <- tibble::as_tibble(pirna_meta, .name_repair="minimal")
      chck_str <- unique(unlist(strsplit(paste(pirna_meta[[1]][1:5], collapse=""), "")))
      
      # pirna_meta against PAC
      if(sum(chck_str %in% c("N", "T", "C", "A", "G")) == length(chck_str)){ 
           if(any(!as.character(pirna_meta[[1]]) %in% rownames(PAC$Anno))){
           warning("The 1st ,'ID', column in 'pirna_meta' does not match rownames in 'PAC$Anno'.\nPlease, double check all sequence names in 'pirna_meta'.")
        }
        pirna_meta <- pirna_meta[match(pirna_meta[[1]], rownames(PAC$Anno)),]
      }else{
      # pirna_meta against pirna reference
        pirna_id  <- dplyr::bind_cols(lapply(reanno$Full_anno, function(x){
             x$pirna[names(x$pirna) == "ref_hits"]
          }), .name_repair = "minimal")
        
        pirna_id <- tibble::tibble(pirna_id=apply(pirna_id, 1, function(x){
          paste(x, collapse="|")}))
        pirna_id[[1]] <- gsub(":sense|:antisense", "", pirna_id[[1]])
        
        pirna_id_splt <- strsplit(pirna_id[[1]], "\\|")
        uni_id <- unique(unlist(pirna_id_splt))
        uni_id <- uni_id[!uni_id == "NA"]
         
        if(any(!uni_id %in% pirna_meta[[1]])){
           stop("The 1st ,'ID', column in pirna_meta does not match piRNA reference names.\nPlease, make sure that names in 'pirna_meta' are compatilble with fasta names.")
        }
        
        
        
        
        !!!! HÄR ÄR JAG 
        Fixa så att pirna_meta för reference får length = PAC och ordered as PAC
        
        
        
        
        
        
        pirna_meta <- pirna_meta[match(pirna_meta[[1]], rownames(PAC$Anno)),]
        }
      }
       
            
            
  
    rm(reanno)            

  ##### Run PAC_gtf and generate anno ####################################
  full <- PAC_gtf(PAC, genome=genome, return="all", mismatches=mismatches, gtf_repeat=gtf_repeat, gtf_protein=gtf_protein, stranded=TRUE, threads=threads)
  
  if(is.null(pirna_meta){
    PAC$Anno <- cbind(data.frame(Length=PAC$Anno$Length), anno_pirna, full$simplify)
  }else{
    PAC$Anno <- cbind(data.frame(Length=PAC$Anno$Length, pirna_meta=pirna_meta), anno_pirna, full$simplify)
  }
  rownames(PAC$Anno) <- seqs
  
  PAC$Anno$repeats <- paste(PAC$Anno[,pirna_hit_col], ifelse(!is.na(full$simplify$repName), "repeats", "not_repeats"), sep="|")
  PAC$Anno$protein <- paste(PAC$Anno[,pirna_hit_col], ifelse(grepl("^protein_coding|\\|protein_coding", full$simplify$gene_biotype), "protein", "not_protein"), sep="|")

  ##### Create group means  ####################################
  if(!norm=="counts"){
    if(!any(names(PAC) == "norm")){
       PAC <- PAC_norm(PAC, type=norm, PAC_merge=TRUE)
       warning("\nThere were no normalized table named '", norm, "' in PAC$norm.\nWill try to generate normalized data with the existing PAC.\nNote, that generating normalized values from filtered \ndata may be incorrect.")
    }else{
       if(!names(PAC$norm) == norm){
         warning("\nThere were no normalized table named '", norm, "' in PAC$norm.\nWill try to generate normalized data with the existing PAC.\nNote, that generating normalized values from filtered \ndata may be incorrect.")
         PAC <- PAC_norm(PAC, type=norm, PAC_merge=TRUE)
       }
     }
   }
  PAC <- PAC_summary(PAC, norm=norm, type="means", pheno_target=pheno_target, PAC_merge = TRUE)

  ##### Create clusters  ####################################
  full_pirna <- full$full[anno_pirna[[2]] == "pirna"]
  lapply(full_pirna, function(x){
                  gr <- GenomicRanges::GRanges(seqnames=x$seqid,  ranges=IRanges(start=x$start, x$end), strand=x$strand)
                  <- GenomicRanges::resize(gr, width=cluster_gap, fix="start", use.names=TRUE)
  GenomicRanges::findOverlaps(gK.reg.ext , gr_repeats, minoverlap=1L) ->olap_repeats

  tt.all -> res_tab

colnames(res_tab)
res_tab[c(4,1:2,8)] -> man
colnames(man) <- c("SNP","CHR","BP","P")
as.numeric(man$BP) -> man$BP          # Need to make numeric
man$CHR <- gsub("chr", "", man$CHR)   # Need to get rid of "chr"
as.numeric(man$CHR) -> man$CHR        # Need to make numeric
head(man)
top<- as.character(man[res_tab$adj.P.Val_MvalC < 0.1,1])
temp<- man[1:(round(nrow(man)*0.01)),]


# Top 5% (rank=13632; p<0.02787317)
qqman::manhattan(man, ylim = c(0, 8), cex=1.5, ylab="", xlab="", pwd=2, col = c("blue4", "gray3"), highlight=top, suggestiveline = F, genomewideline = -log10(res_tab$P.Value_MvalC[round(nrow(res_tab)*0.01)]))
res_tab$P.Value_MvalC[round(nrow(res_tab)*0.01)]



  ##### Overview graphs  ####################################
invisible(capture.output(
    pie_lst_1$mis <- PAC_pie(PAC, pheno_target=pheno_target, anno_target=list(pirna_col))
    PAC_pirna <- PAC_filter(PAC, subset_only=TRUE, anno_target=list(pirna_hit_col, "pirna"))
    PAC_not_pirna <- PAC_filter(PAC, subset_only=TRUE, anno_target=list(pirna_hit_col, "not_pirna"))
    PAC_lst <- list(pirna=PAC_pirna, not_pirna=PAC_not_pirna)
    pie_lst_2 <- lapply(PAC_lst, function(x){
              pie_lst <- list(NULL)
              pie_lst$rep <- PAC_pie(x, pheno_target=pheno_target, anno_target=list("repeats"))
              pie_lst$repClass <- PAC_pie(x, pheno_target=pheno_target, anno_target=list("repClass"))
              pie_lst$prot <- PAC_pie(x, pheno_target=pheno_target, anno_target=list("protein"))
              pie_lst$protBio <- PAC_pie(x, pheno_target=pheno_target, anno_target=list("gene_biotype"))
              }))
    list(all=pie_lst_1, 
  ))


  ##### Nuc bias graphs  ####################################
  anno_trgt_lst <- list(NULL, "repeats", "protein")
  names(anno_trgt_lst) <- c(all="all", repeats="repeats", "protein")

  # 1st nuc bias
  nbias_pos1 <- lapply(anno_trgt_lst, function(x){
        if(is.null(x)){
          mis <- paste0("mis", 0:mismatches)
          res_lst <- list(NULL)
          res_lst[[1]] <- PAC_nbias(PAC, position=1, summary_target=list(paste0(norm,"Means_", pheno_target[[1]])), anno_target=NULL)
          names(res_lst)[1] <- "all"
          for(i in 1:length(mis)){
                PAC$Anno$temp <- as.character(PAC$Anno[,pirna_col] %in% c(mis[1:i]))
                res_lst[[i+1]] <- PAC_nbias(PAC, position=1, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list("temp", "TRUE"))
                names(res_lst)[i+1]<- mis[i]
          }
          res_lst$not_pirna <- PAC_nbias(PAC, position=1, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list("temp", "FALSE"))
        }else{
          if(x=="repeats"){
            ann_targ <- c("pirna|repeats", "pirna|not_repeats", "not_pirna|repeats", "not_pirna|not_repeats")
          }
          if(x=="protein"){
            ann_targ <- c("pirna|protein", "pirna|not_protein", "not_pirna|protein", "not_pirna|not_protein")
          }
          res_lst <- list(NULL)
          for(i in 1:length(ann_targ)){
            res_lst[[i]] <- PAC_nbias(PAC, position=1, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list(x, ann_targ[i]))
            names(res_lst)[i] <- ann_targ[i]
          }
        }
        res_lst  <- lapply(res_lst, function(x){x[["Histograms"]]})
        return(res_lst)
  })

  # 10th nuc bias
nbias_pos10 <- lapply(anno_trgt_lst, function(x){
        if(is.null(x)){
          mis <- paste0("mis", 0:mismatches)
          res_lst <- list(NULL)
          res_lst[[1]] <- PAC_nbias(PAC, position=10, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=NULL)
          names(res_lst)[1] <- "all"
          for(i in 1:length(mis)){
                PAC$Anno$temp <- as.character(PAC$Anno[,pirna_col] %in% c(mis[1:i]))
                res_lst[[i+1]] <- PAC_nbias(PAC, position=10, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list("temp", "TRUE"))
                names(res_lst)[i+1]<- mis[i]
          }
          res_lst$not_pirna <- PAC_nbias(PAC, position=10, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list("temp", "FALSE"))
        }else{
          if(x=="repeats"){
            ann_targ <- c("pirna|repeats", "pirna|not_repeats", "not_pirna|repeats", "not_pirna|not_repeats")
          }
          if(x=="protein"){
            ann_targ <- c("pirna|protein", "pirna|not_protein", "not_pirna|protein", "not_pirna|not_protein")
          }
          res_lst <- list(NULL)
          for(i in 1:length(ann_targ)){
            res_lst[[i]] <- PAC_nbias(PAC, position=10, summary_target=list(paste0(norm,"Means_",pheno_target[[1]])), anno_target=list(x, ann_targ[i]))
            names(res_lst)[i] <- ann_targ[i]
          }
        }
        res_lst  <- lapply(res_lst, function(x){x[["Histograms"]]})
        return(res_lst)
  })


  
  
  
  