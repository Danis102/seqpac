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
#' @param PAC PAC-list object.
#'   
#' @param pirna Character, either indicating a file path to a bowtie indexed
#'   piRNA fasta reference (see vignette for more information) or the name of a
#'   column in PAC$Anno holding piRNA mapping information (single column
#'   containing 'mis0', 'mis1', 'mis2', 'mis3', '_' mismatch classifications).
#'   In case a fasta file exist the function will use the reannotation workflow
#'   to map PAC sequences against the fasta reference.
#'   
#' @param pirna_meta Data frame or tibble with exactly 2 columns with character
#'   strings. The first column should hold either the names of features in the
#'   pirna fasta reference. In such case the second column should contain
#'   additional annotations for each pirna fasta feature. In the examples below
#'   are an example on how a pirbase data file can be used for additional depths
#'   in the analysis. If the first column contains a nuclotide sequence (C, G, T, A, N
#'   characters) the additional annotations will be matched against the sequence
#'   names in PAC instead.
#'   
#' @param genome Character, either indicating a file path to a bowtie indexed
#'   reference genome in fasta format (see vignette for more information) or the
#'   name of a column in PAC$Anno holding genome mapping coordinates (multiple
#'   column with the prefix 'mis0_', 'mis1_', 'mis2_', 'mis3_' etc). In case a
#'   fasta file exist the function will use the reannotation workflow to map PAC
#'   sequences against the fasta reference.
#' 
#' @param mismatches Integer indicating the number of mismatches to be allowed
#'   in the analysis.
#' 
#' @param chrom_size Character vector. In case a fasta reference is not provided
#'   in \code{genome} then the user can provide the lengths of the chromosomes
#'   for optimal manhattan plotting.
#'   
#' @param stranded Logical whether or not matches between gtf coordinates and
#'   genome coordinates should be strand sensitive. If stranded=FALSE (default),
#'   hits on the opposite strand will be allowed and reported reported as '+/-',
#'   beside +/+ and  -/- hits. If stranded=TRUE, only +/+ and  -/- hits will be
#'   reported.
#'   
#' @param gtf_repeat Character indicating the file path to a repeatMasker
#'   formated gtf file (containing 'repFamily', 'repClass' and 'repName' column
#'   names when imported using rtracklayer::readGFF). Can also directly be
#'   provided as a tibble data frame.
#'
#' @param gtf_protein Character indicating the file path to an Ensembl formated
#'   gtf file (containing 'type', 'gene_name', 'gene_biotype', 'exon_number',
#'   'strand' column names when imported using rtracklayer::readGFF). Can also
#'   directly be provided as a tibble data frame.
#'   
#' @param norm Character indicating what type of normalization that should be
#'   used for plotting. Note, in case there are no match between \code{norm} and
#'   a table in \code{PAC$norm}, the function will automatically attempt to generate
#'   normalized values using the function \code{PAC_norm}.
#'   
#' @param pheno_target List with: 1st object being character vector of target
#'   column(s) in Pheno, 2nd object being a character vector of the target
#'   group(s) in the target column (1st object). Used for subsetting data and
#'   create group divided plots.
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
#' # Humans:
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Humans/pirBase/hg38/piR_hsa.txt.gz"), delim="\t", col_names = TRUE)
#' # Flies:
#' pirbase_dat <- readr::read_delim(gzfile("/data/Data_analysis/Genomes/Drosophila/dm6/sports/Drosophila_melanogaster/piRNA_piRBase/piR_dme.txt.gz"), delim="\t", col_names = TRUE)
#' pirna_meta <- tibble::tibble(name=pirbase_dat$name, pubmed_id=pirbase_dat$pubmed)
#' rm(pirbase_dat)
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
#' cluster=list(max_gap=1000, min_n=2) 
#' pirna_meta=pirna_meta
#' chrom_size=NULL
#' 
#' 
#' pirna_analysis <- PAC_pirna(PAC, pirna=pirna, pirna_meta=pirna_meta, genome=genome, chrom_size=NULL, mismatches=3, stranded=FALSE,
#'                       gtf_repeat=gtf_repeat, gtf_protein=gtf_protein, norm=norm, pheno_target=pheno_target, threads=1)
#'  
#'      
#' @export

PAC_pirna <- function(PAC, pirna=NULL, pirna_meta=NULL,  genome=NULL, chrom_size=NULL, mismatches=3, stranded=FALSE,
                       gtf_repeat=NULL, gtf_protein=NULL, norm="counts", pheno_target=NULL, threads=1){

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
           cat("\nID column in pirna_meta holds nucleotide sequences and \nwill therefore be matched with PAC$Anno sequence(row) names.")
           if(any(!as.character(pirna_meta[[1]]) %in% rownames(PAC$Anno))){
           warning("The 1st ,'ID', column in 'pirna_meta' does not match rownames in 'PAC$Anno'.\nPlease, double check all sequence names in 'pirna_meta'.")
        }
        pirna_meta <- pirna_meta[match(pirna_meta[[1]], rownames(PAC$Anno)),]
      }else{
      # pirna_meta against pirna reference
        cat("\nProcessing pirna_meta against fasta reference ... \n")
        pirna_id  <- dplyr::bind_cols(lapply(reanno$Full_anno, function(x){
             x$pirna[names(x$pirna) == "ref_hits"]
          }), .name_repair = "minimal")

        pirna_id <- tibble::tibble(pirna_id=apply(pirna_id, 1, function(x){
          paste(x, collapse="|")
          }))
        pirna_id[[1]] <- gsub(":sense|:antisense", "", pirna_id[[1]])

        pirna_id_splt <- strsplit(pirna_id[[1]], "\\|")
        uni_id <- unique(unlist(pirna_id_splt))
        uni_id <- uni_id[!uni_id == "NA"]

        meta_logi <- pirna_meta[[1]] %in% uni_id
        if(!sum(meta_logi) == length(uni_id)){
           stop("The 1st ,'ID', column in pirna_meta does not match piRNA reference names.\nPlease, make sure that names in 'pirna_meta' are compatilble with fasta names.")
        }

        pirna_meta <- pirna_meta[meta_logi,]
        splt_meta <- strsplit(pirna_meta[[2]], " |\\||;")
        hits_meta <- lapply(pirna_id_splt, function(x){
            hits <- pirna_meta[pirna_meta[[1]] %in% x,]
            hits <- unique(unlist(strsplit(hits[[2]], " |\\||;")))
            return(tibble::tibble(x=paste0(hits, collapse="|"),  y=length(hits)))
        })
        hits_meta <- do.call("rbind", hits_meta)
        col_1 <- paste0("pirna_",names(pirna_meta)[2])
        col_2 <- paste("pirna_n", names(pirna_meta)[2], sep="_")
        names(hits_meta) <- c(col_1, col_2)
        }
      }
  rm(reanno, splt_meta, meta_logi)

  ##### Run PAC_gtf and generate anno ####################################
  full <- PAC_gtf(PAC, genome=genome, return="all", mismatches=mismatches, gtf_repeat=gtf_repeat, gtf_protein=gtf_protein, stranded=TRUE, threads=threads)

  if(is.null(pirna_meta)){
    PAC$Anno <- cbind(data.frame(Length=PAC$Anno$Length), anno_pirna, full$simplify)
  }else{
    PAC$Anno <- cbind(data.frame(Length=PAC$Anno$Length), anno_pirna, hits_meta, full$simplify)
  }
  rownames(PAC$Anno) <- seqs

  PAC$Anno$repeats <- paste(PAC$Anno[,pirna_hit_col], ifelse(!is.na(full$simplify$repName), "repeats", "not_repeats"), sep="|")
  PAC$Anno$protein <- paste(PAC$Anno[,pirna_hit_col], ifelse(grepl("^protein_coding|\\|protein_coding", full$simplify$gene_biotype), "protein", "not_protein"), sep="|")



  ##### Create group means  ####################################
  if(!norm=="counts"){
    if(!any(names(PAC) == "norm")){
       PAC <- PAC_norm(PAC, type=norm, PAC_merge=TRUE)
       cat("\n")
       warning("\nThere were no normalized table named '", norm, "' in PAC$norm.\nWill try to generate normalized data with the existing PAC.\nNote, generating normalized values from filtered \ndata may be incorrect.")
    }else{
       if(!names(PAC$norm) == norm){
         cat("\n")
         warning("\nThere were no normalized table named '", norm, "' in PAC$norm.\nWill try to generate normalized data with the existing PAC.\nNote, generating normalized values from filtered \ndata may be incorrect.")
         PAC <- PAC_norm(PAC, type=norm, PAC_merge=TRUE)
       }
     }
   }
  PAC <- PAC_summary(PAC, norm=norm, type="means", pheno_target=pheno_target, PAC_merge = TRUE)
  sum_nam_pheno  <- paste0(norm,"Means_", pheno_target[[1]])

  PAC <- suppressWarnings(PAC_summary(PAC, norm=norm, type="means", pheno_target=NULL, PAC_merge = TRUE))
  sum_nam_all  <- paste0(norm,"Means_All")

  PAC <- suppressWarnings(PAC_summary(PAC, norm=norm, type="log2FC", pheno_target=pheno_target, PAC_merge = TRUE))
  sum_nam_log2FC  <- paste0("Log2FC_", pheno_target[[1]])

  stopifnot(PAC_check(PAC))


  ##### Create clusters using GenomicRanges::findOverlaps ####################################
  full_pirna <- full$full[anno_pirna[[2]] == "pirna"]

  sum_lst <- list(X=sum_nam_all, Y=sum_nam_log2FC)
  names(sum_lst) <- c(sum_nam_all, sum_nam_log2FC)

  tibb_lst <- lapply(sum_lst, function(x){
    sm <- PAC$summary[[x]]
    sm <- sm[anno_pirna[[2]] == "pirna",]
    tibb <- list(NULL)
    for(i in 1:length(full_pirna)){
       tibb[[i]] <- full_pirna[[i]][c("seqid","start","end","strand")]
       tibb[[i]]$id <- names(full_pirna)[i]
       tibb[[i]]$xyz <- sm[i]
    }
    tibb <- do.call("rbind", tibb)
    names(tibb)[names(tibb)=="xyz"] <- x
    return(tibb[!is.na(tibb[[1]]),])
  })
  stopifnot(identical(tibb_lst[[1]]$id, tibb_lst[[2]]$id))

  tibb  <- dplyr::bind_cols(tibb_lst[[1]], tibb_lst[[2]][length(tibb_lst[[2]])])
  gr <- GenomicRanges::GRanges(seqnames=tibb$seqid, ranges=IRanges::IRanges(start=tibb$start, tibb$end), strand="*")
  gr <- GenomicRanges::resize(gr, width=(IRanges::width(gr)+(2*cluster$max_gap)), fix="center", use.names=FALSE)
  gr_cluster <- GenomicRanges::reduce(gr)
  olap <- GenomicRanges::findOverlaps(gr_cluster, gr, minoverlap=1L)
  olap <- tibble::tibble(from=olap@from, to=olap@to)
  splt_cluster <- split(olap, as.factor(olap$from))
  splt_logi <- lapply(splt_cluster, function(x){length(x)<cluster$min_n})
  splt_cluster <- splt_cluster[!unlist(splt_logi)]
  lst_cluster <- list(NULL)
  for(i in 1:length(splt_cluster)){
               clster <-  tibb[splt_cluster[[i]]$to,]
               srnd <- sort(paste(unique(clster$strand)), decreasing = TRUE)
               lst_cluster[[i]] <- clster
               names(lst_cluster)[i] <- paste0(unique(clster$seqid),":", min(clster$start),"-", max(clster$end), ":", paste(srnd, collapse=""))
    }

  ##### Plot manhattan of clusters and mean rpm/log2FC ####################################
  if(is.null(chrom_size) && file.exists(genome)){
    genm <- Biostrings::readDNAStringSet(genome)
    chrom_size <- Biostrings::width(genm)
    names(chrom_size) <- do.call("rbind", strsplit(names(genm), " "))[,1]
  }
  # Check max value at start
  chr <- as.character(NULL)
  for(i in 1:nrow(tibb)){
      if(tibb$strand[i]=="-"){
        coord <- tibb$end[i]
      }else{
        coord <- tibb$start[i]
      }
      chr[i] <- paste(tibb$seqid[i], coord, sep=":")
  }

  mam_pl_lst<- lapply(sum_lst, function(x){
    coord_val <- tibble::tibble(coord=chr, value=tibb[[x]])
    coord_val <- split(coord_val, as.factor(coord_val$coord))
    coord_val <- unlist(lapply(coord_val, function(x){unique(max(x$value))}))
    coord_nam <- do.call("rbind", strsplit(names(coord_val), ":"))
    chrm <- coord_nam[,1]
    chrm_uni <- sort(unique(coord_nam[,1]))

    for(i in 1:length(chrm_uni)){
      chrm[chrm==chrm_uni[i]] <- paste0(i,"XWDC")
    }
    chrm <- gsub("XWDC", "", chrm)
    chrm_uni[grepl("mitochond", chrm_uni)] <- "MT"
    man  <- data.frame(chrom=as.numeric(chrm), position=as.numeric(coord_nam[,2]), value=as.numeric(coord_val))
    man_dat <- man %>%
      dplyr::group_by(chrom) %>%
      dplyr::summarise(chr_len=max(position)) %>%
      dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
      dplyr::select(-chr_len) %>%
      dplyr::left_join(man, ., by=c("chrom"="chrom")) %>%
      dplyr::arrange(chrom, position) %>%
      dplyr::mutate( position_cum=position+tot)

    axs <- man_dat %>%
      group_by(chrom) %>%
      summarise(center=( max(position_cum) + min(position_cum) ) / 2 )

     if(grepl("Log2|log2", x)){
       ylb <- paste0("Log2FC (", norm, ")")
     }else{
       ylb <- paste0("mean ", norm)
       man_dat$value <- log2(man_dat$value)
     }

    man_plots <- ggplot(man_dat, aes(x=position_cum, y=value)) +
                geom_point( aes(color=as.factor(chrom)), alpha=0.8, size=1.3) +
                scale_color_manual(values = rep(c("#094A6B", "#9D0014"), 50 )) +
                scale_x_continuous( label = chrm_uni[axs$chrom], breaks= axs$center ) +
                scale_y_continuous(expand = c(0, 1) ) +
                ggtitle("Manhattan piRNA", subtitle ="(Note: due to multimapping pseudoreplication is likely)") +
                ylab(ylb) +
                xlab("chromosomes") +
                theme_bw() +
                theme(
                  legend.position="none",
                  panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()
            )
    return(man_plots)
  })

  ##### Overview graphs  ####################################

    pie_lst_1 <- list(NULL)
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
              })
   pie_res <- list(all=pie_lst_1, divided=pie_lst_2)



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

  ##### Put humpty-dumty together  ####################################

return(list(data=list(PAC=PAC,  pie=pie_res, manhattan=mam_pl_lst, nbais=list(positon_1=nbias_pos1, Postion_10=nbias_pos10))))
}
  
  