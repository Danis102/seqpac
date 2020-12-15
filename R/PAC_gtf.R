#' Annotate against a GTF file
#'
#' This function will annotate a PAC object using input from a GTF/GFF file.
#'
#' Given a PAC object and a gtf formated annotation file(s), this function will
#' attempt to annotate sequences mapped to a reference genome against genomic
#' coordinates in the gtf file(s). In case no genomic mapping coordinates are
#' available in the PAC object, the function provides a backdoor into the
#' \emph{PAC reannotation} workflow, where genome mapping is performed using
#' Bowtie.
#' 
#' @family PAC analysis
#'   
#' @seealso \url{https://github.com/Danis102} for updates.
#'
#' @param PAC PAC-list object. The Anno object needs to contain genome
#'   coordinates in the following columns: "chr","start","end", and "strand".
#'   These coordinates are eaisly obtained by using \code{PAC_genome}
#'   
#' @param genome Character indicating the path to a reference genome file in
#'   fasta format. Note that a bowtie index, created by for example
#'   Rbowtie::bowtie-build,  must be contained in the same folder and having the
#'   same basename as the fasta.
#'   
#' @param mismatches Integer indicating the number of allowed mismatches for the
#'   genome mapping.
#'   
#' @param return Character indicating what information to return. 
#' 
#'    If return="simplify" (default), a table containing with one unique
#'    sequence per row will be returned. Multiple hits between genomic
#'    coordinates and gtf coordinates will be merged and only unique annotations
#'    will be reported.
#'    
#'    If return="full", a list will be returned containing 1 table reporting all
#'    annotations for each genomic coordinate of each unique sequence.
#'    
#'    If return="all", both a simplified table and a full annotation list will
#'    be returned as a list.
#'    
#'    If return="merge", a simplified table will be merged with the Anno table
#'    of the provided PAC object, and an updated PAC object containing the new
#'    annotations will be returned.
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
#' @param gtf_other Named list of characters, indicating file path(s) to other
#'   gtf files with differing formats. Can also directly be provided as a tibble
#'   dataframe in a named list.
#'   
#' @param target_other Named list of character vectors indicating target columns
#'   in the listed gtf files in \emph{gtf_other}. Important, the listed objects
#'   must have the same length and names as in \emph{gtf_other}. The vector
#'   indicates the column names as if imported by rtracklayer::readGFF.
#'   
#' @param threads Integer indicating the number of parallel processes.
#'   
#' @return List, tibble dataframe or updated PAC object. See \emph{return} for
#'   more information.
#'   
#'   
#' @examples
#'
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' ##############################################################
#' ## Simplified repeatmasker annotation with genomic mapping
#'   
#' genome <- "/home/danis31/Desktop/Temp_docs/fasta/biomartr_genome/chromosomes.fa"
#' gtf_repeat <- "/home/danis31/Desktop/Temp_docs/fasta/repeatMasker/repeatMasker_ensembl.gtf"
#' 
#' # Only returns tibble:
#' repeat_simple <- PAC_gtf(pac, genome=genome, return="simplify", gtf_repeat=gtf_repeat, threads=10)
#' 
#' # Merge with PAC$Anno dataframe: 
#' pac_merge <- PAC_gtf(pac, genome=genome, return="merge", gtf_repeat=gtf_repeat, threads=10)        
#' 
#' 
#' ##############################################################
#' ## Full output previously mapped columns up to 3 mismatches
#' 
#' # Generates an error because genome mapping was done
#' # with add_reanno(genome_max=10):
#' genome_col <- colnames(pac$Anno)[grepl("chromosomes_genome", colnames(pac$Anno))]

#' 
#' repeat_full <- PAC_gtf(pac, genome=genome_col, return="full", gtf_repeat=gtf_repeat, threads=10) 
#' 
#' # Works because PAC_gtf automatically maps the genome with add_reanno(genome_max="all")
#' genome_col <- colnames(pac_merge$Anno)[grepl("^genome|mis\\d_genome", colnames(pac_merge$Anno))]
#' repeat_full <- PAC_gtf(pac_merge, genome=genome_col, return="full", gtf_repeat=gtf_repeat, threads=10)
#' 
#' # return="full" returns all annotation for all each coordinate
#' repeat_full[800:820]
#' head(repeat_full["TGCGGAAGGATCATTA_mis0"])
#'      
#' 
#' ##############################################################
#' ## With additional gtfs including custom other
#' gtf_repeat <- "/home/danis31/Desktop/Temp_docs/fasta/repeatMasker/repeatMasker_ensembl.gtf"
#' gtf_protein <- "/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/gtf/Drosophila_melanogaster.gtf"
#' 
#' # Note, target_other points to columns in each gtf listed in gtf_other 
#' # and the list objects must therefore have the same names:
#' gtf_other=list(rep=gtf_repeat, prot=gtf_protein)
#' target_other=list(rep="repFamily", prot=c("type", "gene_id")) 
#' 
#' genome_col <- colnames(pac_merge$Anno)[grepl("^genome|mis\\d_genome", colnames(pac_merge$Anno))]
#' many_simply <- PAC_gtf(pac_merge, genome=genome_col, return="simplify", mismatches=3,
#'                        gtf_repeat=gtf_repeat, gtf_protein=gtf_protein,  
#'                        gtf_other=gtf_other, target_other=target_other, threads=10)
#' 
#' # With perfect alignments (0 mismatches)
#' many_simply <- PAC_gtf(pac_merge, genome=genome_col, return="simplify", mismatches=0,
#'                        gtf_repeat=gtf_repeat, gtf_protein=gtf_protein, threads=10)
#' 
#' ##############################################################
#' ## Convert of UCSC to Ensembl (Tanks to Devon Ryan and co-workers)
#' 
#' # Go to https://github.com/dpryan79/ChromosomeMappings
#' # Locate your genome of choice
#' # Download from raw.githubusercontent:
#' dm6_conv <- readr::read_tsv("https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/BDGP6_UCSC2ensembl.txt", col_names=FALSE)
#' names(dm6_conv) <- c("ucsc", "ensembl")
#' 
#' # Fix mito genome name (should not contain "dmel_")
#' dm6_conv[grepl("mito", dm6_conv$ensembl),]
#' dm6_conv$ensembl <- gsub("dmel_mitochondrion_genome", "mitochondrion_genome", dm6_conv$ensembl) 
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
#' @export

PAC_gtf<- function(PAC, genome=NULL, mismatches=3, return="simplify", stranded=FALSE,
                   gtf_repeat=NULL, gtf_protein=NULL, gtf_other=NULL, target_other=NULL, 
                   threads=1){

##### Setup general ####################################
  seqs <- rownames(PAC$Anno)

##### Setup genome and run reanno if necessary ####################################
  # If user do not know columns
  if(is.null(genome)){
      anno_genome <- tibble::as_tibble(PAC$Anno, .name_repair="minimal")
    }else{
      # If user know the columns:
      if(any(!file.exists(genome))){
          anno_genome <- tibble::as_tibble(PAC$Anno[, genome, drop=FALSE], .name_repair="minimal")
      }else{
      # If user provide fasta reference:
      if(file.exists(genome)){
          cat("\nInput genome was an existing file. Will treat it as a \nfasta reference and make a denovo reannotation using bowtie. \nSee ?map_reanno or ?vingette for details.\n")
          outpath <- tempfile(pattern = "", fileext = ".out")
          err <- try(map_reanno(PAC, ref_paths=list(genome=genome), output_path=outpath, type="external", mismatches=mismatches,
                         import ="genome", threads=threads), silent = TRUE)
          if(!is.null(err)){
            outpath <- tempfile(pattern = "", fileext = "")
            err2 <- try(map_reanno(PAC, ref_paths=list(genome=genome), output_path=outpath, type="internal", mismatches=mismatches,
                         import ="genome", threads=threads), silent = TRUE)
          if(!is.null(err2)){
            stop(paste0("\nFunction map_reanno failed. Possible reasons: \n\tNo bowtie installation\n\tBad fasta reference\n\tNo bowtie index (see ?Rbowtie::bowtie_build)\n\nLast log says:\n", err2))
            }
          }
          reanno <- make_reanno(outpath, PAC=PAC, mis_fasta_check=TRUE, threads=threads)
          anno_genome <- add_reanno(reanno, type="genome", mismatches = 3, genome_max="all")
          rm(reanno)
          unlink(outpath, recursive = TRUE)
      }
     }
    }
  # Extract genome columns
  srch <- paste0("^mis", 0:mismatches)
  nam <- colnames(anno_genome)
  logi_col1 <- grepl("_genome$", nam)
  logi_col2 <- grepl(paste(srch, collapse="|"), nam)
  nam_col <- gsub("mis0_", "", nam[grepl("^mis0_", nam)])
  mis_col <- grepl(paste0("^", nam_col), nam)
  logi_coord <- rowSums(cbind(logi_col1, logi_col2))==2
  prefix <- suppressWarnings(do.call("rbind", strsplit(nam[logi_coord], "_"))[,1])
  if(any(duplicated(prefix))){
      stop("\nFound more than 1 possible genome annotation.\nAutomatic column identification only works when the results \nof one fasta reference genome is reported in PAC$Anno.\nPlease, specify in 'genome=' the exact columns reporting the alignments \nagainst 1 reference genome (column prefix: mis0_, mis1_, mis2_ etc.) \nor proved a path to a bowtie indexed fasta genome file.")
    }
  coord_genome <- anno_genome[,logi_coord, drop=FALSE]
  mis_genome <- anno_genome[,mis_col, drop=FALSE]
  # Check that all columns are included
  logi_check  <- apply(coord_genome, 2, function(x){any(grepl("Warning>\\d", x))})
  if(any(logi_check)){
        stop("\nCoordinates in genome appear to have truncation warnings ('Warning>'). \nThis indicates that genome was mapped with reporting a limited number of \nalignments per sequence. Please, rerun the 'add_reanno' function using \nmax_genome='all', or provide the path to a fasta genome reference in 'genome'.")
      }

##### Organize mapped coordinates ####################################
  cat("\n\nReorganizing coordinates ...")
  coord_lst <- list(NULL)
  lng_all <- PAC$Anno[,colnames(PAC$Anno) %in% c("Length", "Size")]
  mis_incl <- unique(mis_genome[[1]])
  mis_incl <- mis_incl[mis_incl %in% paste0("mis", 0:mismatches)]
  for(i in 1:nrow(mis_genome)){
    mis <- mis_genome[[1]][i]
    lng <- lng_all[i]
    if(!mis %in% mis_incl){
      tb <- tibble::tibble(seqid=NA, start=NA, end=NA, strand=NA)
    }else{  
      targt <- coord_genome[i, grepl(mis, names(coord_genome))]
      splt <- unlist(strsplit(targt[[1]], "\\|"))
      tb <- tibble::as_tibble(do.call("rbind", strsplit(splt, ":")), .name_repair="minimal")
      tb <- tibble::tibble(seqid=tb[[1]], start=as.numeric(tb[[2]]), end= as.numeric(tb[[2]])+lng, strand=tb[[3]])
    }
    coord_lst[[i]] <- tb
    names(coord_lst)[i] <- mis
  }
  
##### Import gtfs ####################################
  gtf_lst <- list(repeats=gtf_repeat, protein=gtf_protein, other=gtf_other)
  logi_list <- unlist(lapply(gtf_lst, class))=="list"
  if(any(logi_list)){
     gtf_lst <- c(gtf_lst[!logi_list], unlist(gtf_lst[logi_list]))
  }

  # Import gtf
  cat("\nImport gtf files ...")
  gtf_lst <- gtf_lst[!unlist(lapply(gtf_lst, is.null))]
  gtf_import <- lapply(gtf_lst, function(x){
        if(file.exists(x)){
          gtf <- tibble::as_tibble(rtracklayer::readGFF(x), .name_repair="minimal")
        }
        if(!file.exists(x)){
          gtf <- x
        }
      return(gtf)
      })

##### Check format in repeats and protein ####################################  
    for (i in 1:length(gtf_lst)){
      nam <- names(gtf_lst)[i]
      if(nam=="repeats"){
         col_rep <- c("repName", "repClass", "repFamily", "strand")
         logi_RM<- sum(names(gtf_import[[nam]]) %in% col_rep)== 4
         if(!logi_RM){
           stop("\nInput 'gtf_repeats' does not contain repeatMasker standard columns\n('repName', 'repClass', 'repFamily'). Please, check column names \nusing for example rtracklayer::readGFF('<path_to_gtf>') or move the \ngtf path to 'gtf_other' and set target column(s) with 'target_other'.") 
         }
      }
      if(nam=="protein"){ 
         col_prot <- c("type", "gene_name", "gene_biotype", "exon_number", "strand")
         logi_prot <- sum(names(gtf_import[[nam]]) %in% col_prot)== 5
         if(!logi_prot){
           stop("\nInput 'gtf_protein' does not contain Ensembl standard columns\n('type', 'gene_name', 'exon_number', 'gene_biotype'). Please, check \ncolumn names using for example rtracklayer::readGFF('<path_to_gtf>') \nor move the gtf path to 'gtf_other' and set target column(s) \nmanually with 'target_other'.") 
         }
      }
     logi_all <- sum(names(gtf_import[[nam]]) %in% c("seqid", "start", "end", "strand"))== 4
     if(!logi_all){
           stop("\nInput gtf '", nam, "' does not contain essential columns\n('seqid', 'start', 'end', 'strand'). Please, check column \nnames using for example rtracklayer::readGFF('<path_to_gtf>').") 
     }
    }
  
##### Check all chromosomes names in gtf ####################################
  cat("\nChromosome names compatibility check ...")
  chrm_nam <- unique(unlist(lapply(coord_lst, function(x){x$seqid})))
  chrm_nam <- chrm_nam[!is.na(chrm_nam)]
  tot_chrom <- length(chrm_nam)
  
  for(i in 1:length(gtf_import)){
    chrom_in_gtf <- sum(chrm_nam %in% unique(gtf_import[[i]]$seqid))
    if(chrom_in_gtf==0){
      stop("\nYour chromosome names in 'gtf_", names(gtf_import)[i], "' input \ndid not match your genome alignments. Did you forget \nto convert between database formats (e.g. UCSC vs Ensembl vs \nNCBI)? Please, check ?PAC_gtf and vingette(sepac) to instructions \non how to harmonize your gtf files.")
      }
    if(chrom_in_gtf/tot_chrom<0.5 ){
      cat("\n")
      warning("Low overlap between chromosome names of 'gtf_", names(gtf_import)[i], "'\ninput and genome alignments was detected. Please, double \ncheck that conversion between database formats (e.g. UCSC vs \nEnsembl vs NCBI) was completed. Check ?PAC_gtf and vingette(sepac) \nfor instructions on how to harmonize your gtf files to \nthe reference genome.")
    }
      if(!chrom_in_gtf/tot_chrom==1 ){
      cat("\n")
      cat(paste0("  -- Note, not all chrosomsome names were represented in 'gtf_", names(gtf_import)[i],"'.\n     Could be wise to double check that chrosomsome names are compatible."))
      }
  }                       

##### Annotating ####################################
  cat("\nAnnotating against the gtf file(s)  ...")
  # Convert gtf to gr
  gtf_gr<- lapply(gtf_import, function(x){
    gr <- GenomicRanges::GRanges(seqnames=x$seqid, IRanges::IRanges(start=x$start, end=x$end), strand=x$strand)
    return(gr)
  })
  # Convert coordinates to gr
  coord_gr <- lapply(coord_lst, function(x){
    if(is.na(x[[1,1]])){
      gr <- NA
    }else{
      if(stranded=="TRUE"){
          gr <- GenomicRanges::GRanges(seqnames=x$seqid, IRanges::IRanges(start=x$start, end=x$end), strand=x$strand)
      }
      if(stranded=="FALSE"){
          gr <- GenomicRanges::GRanges(seqnames=x$seqid, IRanges::IRanges(start=x$start, end=x$end), strand=NA)
      }
    }
    return(gr)
  })
  
  # Extraction loop
  anno_lst <- list(NULL)
  for(i in 1:length(gtf_gr)){
    gtf_nam <- names(gtf_gr)[i]
    cat(paste0("\n   |--> Extract and compile ", gtf_nam, " ..."))  
    
    # Set up columns to extract
    if(gtf_nam=="repeats"){
       col_target <- col_rep
    }
    if(gtf_nam=="protein"){
       col_target <- col_prot
    }
    if(!gtf_nam %in% c("repeats", "protein")){
       col_target <- target_other[[gsub("other.", "", gtf_nam)]]
    }
    
    # Run overlap and extract anno 
    coord_anno <- lapply(coord_gr, function(x){
      if(paste(x[1])=="NA"){
        uni_gtf <- tibble::as_tibble(matrix(NA, nrow=1, ncol=length(col_target)), .name_repair="minimal")
        names(uni_gtf) <- col_target
      }else{
        olap <- GenomicRanges::findOverlaps(x, gtf_gr[[i]])
        gtf <- gtf_import[[i]][olap@to, col_target]
        uni_gtf <- data.frame(matrix(NA, nrow=length(x), ncol=length(col_target)))
        names(uni_gtf) <- col_target
        if(!nrow(gtf)==0){
          splt <- split(gtf, as.factor(olap@from))
          uni_anno <- lapply(splt, function(z){
              uni <- apply(z, 2, function(y){
                     y <- y[!is.na(y)]
                     paste(sort(unique(y)), collapse="|")
                    })
              return(uni)
              })
          uni_anno <- do.call("rbind", uni_anno)
          uni_gtf[rownames(uni_anno),] <- uni_anno
        }
      }
      return(uni_gtf)
    })
  anno_lst[[i]] <- coord_anno
  names(anno_lst)[i] <- gtf_nam
  }
  
##### Returning objects ####################################  
  if(return %in% c("full","all")){
    cat("\n")
    full_lst <- list(NULL)
    loop <- length(coord_lst)
    for(i in 1:loop){
      cat("\rGenerating full annotation list: ", i, "/", loop)
      flush.console() 
      tab <- coord_lst[[i]]
      for(z in 1:length(anno_lst)){
        tab <- cbind(tab, anno_lst[[z]][[i]])
      }
      full_lst[[i]] <- tibble::as_tibble(tab, .name_repair="minimal")
      names(full_lst)[i] <- paste(seqs[i], names(coord_lst)[i], sep="_")
    }
  }
    
  if(return %in% c("simplify", "all", "merge")){
    cat("\n")
    simp_lst <- list(NULL)
    loop <- nrow(coord_genome)
    for(i in 1:loop){
      cat("\rGenerating simplified annotation table: ", i, "/", loop)
      flush.console() 
      tab <- coord_genome[i,]
      for(z in 1:length(anno_lst)){
         simp_tab <- apply(anno_lst[[z]][[i]], 2, function(x){
            simp <- sort(unique(unlist(strsplit(as.character(x), "\\|"), use.names = FALSE)))
            simp <- paste(simp[!is.na(simp)], collapse="|")
            if(nchar(simp) == 0){
              simp <- NA
            }
            return(simp)
          })
        tab <- cbind(tab, t(as.data.frame(simp_tab)))
      } 
      simp_lst[[i]] <- tab
      names(simp_lst)[i] <- seqs[i]
    }
    simp_df <- tibble::as_tibble(cbind(data.frame(seqs=as.character(seqs)), mis_genome, do.call("rbind", simp_lst)), .name_repair="minimal")
   }

  if(return=="simplify"){
    return(simp_df)
    }
  if(return=="full"){
    return(full_lst)
  }
  if(return=="all"){
    return(list(simplify=simp_df, full=full_lst))
  }
  if(return=="merge"){
     stopifnot(identical(rownames(PAC$Anno), as.character(simp_df$seqs)))
     PAC$Anno <- cbind(PAC$Anno, simp_df[,-1])
     stopifnot(PAC_check(PAC))
     return(PAC)
  }
}
