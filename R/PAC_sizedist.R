#' Generates a size distribution analysis from sequences and counts in a PAC object
#'
#' \code{PAC_sizedist} analyses size distribution with bar charts, allowing for grouping of data with annotation information and summaries.
#'
#' Given a PAC object, the function will attempt to extract the ratios of
#' specific biotypes at a given size of sequences in the Anno data.frame
#' in relation to the sequence counts in Counts. Additionally, count data from Summaries made with
#' \code{\link[seqpac]{PAC_summary}} can also be used.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a count table with raw counts.
#' @param range Integer vector giving the range  in sequence lengths (default=c(min, max)). 
#' 
#' @param anno_target List with: 
#'                          1st object being character vector of target column(s) in Anno, 
#'                          2nd object being a character vector of the target biotypes(s) in the target column (1st object).
#'                          (default=NULL)
#' 
#' @param pheno_target List with: 
#'                          1st object being character vector of target column(s) in Pheno, 
#'                          2nd object being a character vector of the target group(s) in the target column (1st object).
#'                          (default=NULL)
#'                          
#' @param summary_target List with: 
#'                          1st object being character target object in PAC$summary, 
#'                          2nd object being a character vector of the target columns(s) in the target object (1st object).
#'                          (default=NULL)
#'                          
#'
#' @return A list of objects: 
#'               1st object (Histograms::Samples): Individual histograms showing the nucleotide ratios per sample over the specified range.   
#'               2nd object (Data::Samples): Data used to generate the plots.
#'               3rd object (Stacked_bars::Groups): Stacked bars showing the mean ratios of each nucleotide per group over the specified range.
#'               4th object (Error_bars::Groups): Error bar plots with mean ratio of each nucleotide per group over the specified range. 
#'
#' @examples
#' 
#'
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' PAC_filt <- PAC_norm(pac, norm="cpm")
#' PAC_filt <- PAC_summary(PAC=PAC_filt, norm = "cpm", type = "means", pheno_target=list("stage"))
#' 
#' 
#'  ord <- c("no_anno", "other", "miRNA", "tRNA", "rRNA", "snoRNA",  "lncRNA")
#' 
#' sizedist_plots <- PAC_sizedist(PAC_filt, anno_target=list("Biotypes_mis0", ord), summary_target=list("cpmMeans_stage"))
#' cowplot::plot_grid(plotlist=sizedist_plots[[1]], nrow = 3, ncol = 1)
#' 
#' sizedist_plots <- PAC_sizedist(PAC_filt, norm="counts", anno_target=list("Biotypes_mis0", ord), pheno_target=list("batch", "Batch1"))
#' cowplot::plot_grid(plotlist=sizedist_plots[[1]], nrow = 3, ncol = 1)
#' 
#' 
#' @export

PAC_sizedist <- function(PAC, norm="counts", range=NULL, anno_target, pheno_target=NULL, summary_target=NULL, colors=NULL){
  
  ## Organize input
  # anno <- PAC$Anno
  # if(!is.null(norm)){
  #   if(norm=="counts|Counts"){
  #     data <- PAC$Counts
  #     labl <- "rawCounts"
  #   }else{
  #     if(is.null(summary_target)){ 
  #       data <- PAC$norm[[norm]]
  #       labl <- norm
  #     }else{
  #       stop("\nYou have specified both 'norm' and 'summary_target'.\nPlease specify one and set the other to 'NULL'.") 
  #     }
  #   }
  # }else{
  #   data <- PAC$summary[[summary_target[[1]]]]
  #   labl <- summary_target[[1]]
  # }
  # 
  # if(length(summary_target)==1){
  #   summary_target[[2]]  <- names(PAC$summary[[summary_target[[1]]]])
  # }
  # if(!is.null(summary_target)){
  #   data <- data[, colnames(data) %in% summary_target[[2]], drop=FALSE]
  # }   
  # 
  # ## Add range filter
  # if(is.null(range)){
  #   range <- c(min(anno$Size), max(anno$Size))
  # }
  # filt <- anno$Size >= range[1] & anno$Size <= range[2] 
  # anno <- anno[filt,]
  # data <- data[filt,]
  # 
  # ## Reomve unwanted biotypes
  # if(length(anno_target)==1){ 
  #   anno_target[[2]] <- as.character(unique(anno[,anno_target[[1]]]))
  # }
  # filt2 <- anno[,anno_target[[1]]] %in% anno_target[[2]]
  # anno <- anno[filt2,]
  # data <- data[filt2,]
  # 
  # ## Remove unwanted samples
  # if(!is.null(pheno_target)){ 
  #   if(length(pheno_target)==1){ 
  #     pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))}
  #     filt3 <- PAC$Pheno[,pheno_target[[1]]] %in%  pheno_target[[2]]
  #     data <- data[,filt3,drop=FALSE]
  #     ph <- PAC$Pheno[filt3,,drop=FALSE]
  #     match_pfilt <-  order(match(ph[,pheno_target[[1]]], pheno_target[[2]]))
  #     data <- data[,match_pfilt,drop=FALSE]
  #     ph <- ph[match_pfilt,,drop=FALSE]
  # }else{ 
  #   if(!is.null(summary_target)){
  #     ph <- data.frame(colnames(data))
  #     rownames(ph) <- ph[,1] 
  #   }else{ 
  #     ph <-PAC$Pheno 
  #     }
  #   }
  
  # Prepare filtered PAC
  if(!is.null(pheno_target)){ 
    if(length(pheno_target)==1){ 
      pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))
     }
    }
  
  if(!is.null(anno_target)){ 
    if(length(anno_target)==1){ 
      anno_target[[2]] <- as.character(unique(PAC$Anno[,anno_target[[1]]]))
     }
    }
  if(is.null(range)){
    range <- c(min(PAC$Anno$Size), max(PAC$Anno$Size))
    }
  
  PAC <- seqpac::PAC_filter(PAC, size=range, pheno_target=pheno_target, anno_target=anno_target, subset_only=TRUE)
  
  stopifnot(PAC_check(PAC))
  ph <- PAC$Pheno
  anno <- PAC$Anno
  
  # Extract data
  if(is.null(summary_target)){
    if(is.null(norm)){
      data <- PAC$Counts; labl <- "Counts"
    }else{
      if(norm=="counts"){
        data <- PAC$Counts; labl <- "Counts"
      }else{
        data <- PAC$norm[[norm]]; labl <- norm
      }
    }
  }else{
    data <- PAC$summary[[summary_target[[1]]]]; labl <- summary_target
    }
  

  #### Summarize over size and biotype
  bio_fact <- factor(anno[, anno_target[[1]]], levels=anno_target[[2]])
  seq_range <- seq(range[1], range[2])
  size_fact <- factor(anno$Size, levels=seq_range)
  size_lst <- lapply(as.list(data), function(x){
    bio_agg <- aggregate(x, list(bio_fact, size_fact), sum)
    colnames(bio_agg) <- c("biotype", "size", "data")
    bio_agg_lst <- lapply(split(bio_agg, bio_agg$biotype), function(y){
      if(any(!seq_range %in% y$size )){
        y <- rbind(y, data.frame(biotype=as.character(unique(y$biotype)), size=seq_range[!seq_range %in% as.character(y$size)], data=0))
      }
      fin <- y[order(y$size),]  
      stopifnot(identical(as.character(fin$size), as.character(seq_range)))
      return(fin)
    })
    fin_df <- do.call("rbind", bio_agg_lst)
    rownames(fin_df) <- NULL
    return(fin_df)
  })
  
  #### Generate individial histograms 
  options(scipen=999)
  #### Set up colors colors ###
  if(is.null(colors)){
    colfunc <- grDevices::colorRampPalette(c("#094A6B", "#EBEBA6", "#9D0014"))
    rgb <- colfunc(sum(!anno_target[[2]] %in% c("other", "no_anno")))
    rgb_vec <- NULL
    cnt <- 0
    for(i in 1:length(anno_target[[2]])){ 
      if(anno_target[[2]][i] == "other"){
        rgb_vec[i] <- "#808080";  cnt <- cnt+1}
      if(anno_target[[2]][i] == "no_anno"){
        rgb_vec[i] <- "#C0C0C0";  cnt <- cnt+1}
      if(!anno_target[[2]][i] %in% c("other","no_anno")){
        rgb_vec[i] <- rgb[i-cnt]}
    }
  }else{
    rgb_vec <- colors
    }
  
  #### Plot individual plots ###
  histo_lst <- list(NA)
  if(!is.null(summary_target)){
    samp <- colnames(data)
  }else{
    if(is.null(pheno_target)){
      samp <- rownames(ph)
    }else{
      samp <- paste0(ph[,pheno_target[[1]]],"-", rownames(ph)) 
    }
  }
  for(i in 1:length(size_lst)){
    histo_lst[[i]] <- ggplot2::ggplot(size_lst[[i]], ggplot2::aes(x=size, y=data, fill=biotype))+
      ggplot2::geom_bar(width = 0.9, cex=0.2, colour="black", stat="identity")+
      ggplot2::geom_hline(yintercept=0, col="azure4")+
      ggplot2::xlab("Size (nt)")+
      ggplot2::ylab(paste0(labl))+
      ggplot2::labs(subtitle = samp[i])+
      ggplot2::scale_fill_manual(values=rgb_vec)+
      #coord_cartesian(ylim=c(0,10000))+
      #scale_y_continuous(breaks = seq(0, 10000, 2500))+
      #ggthemes::geom_rangeframe(aes(x=range))+   
      ggthemes::theme_tufte()+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
    names(histo_lst)[i] <- names(size_lst)[i]
  }
  return(list(Histograms=histo_lst, Data=size_lst))
}

