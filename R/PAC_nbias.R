#' Generates a nucleotide bias analysis from sequences and counts in a PAC object
#'
#' \code{PAC_nbias} analyses nucleotide bias.
#'
#' Given a PAC object the function will attempt to extract the ratios of
#' specific nucleotides at a given position in sequences in the Anno data.frame
#' in relation to the sequence counts in Counts.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a count table with raw counts.
#' @param position Integer giving the nucleotide postion from 3' to 5' position (default=1).
#' @param range Integer vector giving the range  in sequence lengths (default=c(min, max)).#' 
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
#'                          1st object being character vector of target object in PAC$summary, 
#'                          2nd object being a character vector of the target column(s) in the target summary object (1st object).
#'                          (default=NULL)
#'
#' @param colors Character vector with color codes to be parsed to ggplot2.
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
#' library(seqpac)
#' load("/home/danis31/OneDrive/Programmering/Programmering/Pipelines/Drosophila/Pipeline_3.1/seqpac/dm_test_PAC.Rdata")
#' 
#' PAC_filt <- PAC_rpm(PAC_filt)
#' PAC_filt <- PAC_summary(PAC=PAC_filt, norm = "rpm", type = "means", pheno_target=list("Method"))
#'
#' hierarchy <- list( Mt_rRNA= "12S|16S|Mt_rRNA",
#'                  rRNA="5S|5.8S|18S|28S|S45|Ensembl_rRNA|rRNA_Other",
#'                  Mt_tRNA= "tRNA_mt-tRNA",
#'                  tRNA="Ensembl_tRNA|tRNA_nuc-tRNA",
#'                  miRNA="^miRNA|Ensembl_miRNA|Ensembl_pre_miRNA",
#'                  piRNA="piRNA")
#'
#' PAC_filt <- simplify_reanno(PAC_filt, hierarchy=hierarchy, mismatches=0, bio_name="Biotypes_mis0", PAC_merge = TRUE)
#'
#' ### Plot using raw counts and order samples using pheno_target:
#' nbias_result <- PAC_nbias(PAC=PAC_filt, position=1, norm="raw", pheno_target=list("Method", c("TGIRT", "Proto")))
#' nbias_result <- PAC_nbias(PAC=PAC_filt, position=1, norm="raw", pheno_target=list("Method", "Proto"))    
#' cowplot::plot_grid(plotlist=nbias_result[[1]], nrow = 3, ncol = 2)
#' 



#' ### Plot using rpm counts and order samples using pheno_target:
#' nbias_result <- PAC_nbias(PAC=PAC_filt, position=1, norm="rpm") 
#' 
#' ### Plot using summarized data in summary_target: 
#' nbias_result <- PAC_nbias(PAC=PAC_filt, position=1, summary_target=list("rpmMeans_Method")) 
#' 
#' 
#' cowplot::plot_grid(plotlist=nbias_result[[1]], nrow = 3, ncol = 2)
#' 
#' 
#' @export


PAC_nbias <- function(PAC, position=1, norm=NULL, range=NULL, anno_target=NULL, pheno_target=NULL, summary_target=NULL, colors=NULL){
  # Prepare filtered PAC
  if(!is.null(pheno_target)){ 
    if(length(pheno_target)==1){ pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))
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
  
  
  # Extract data
  if(is.null(summary_target)){
    if(is.null(norm)){
      dat <- PAC$Counts; labl <- "Counts"
    }
    else{
      if(norm=="counts"){
        dat <- PAC$Counts; labl <- "Counts"
      }
      else{
        dat <- PAC$norm[[norm]]; labl <- norm
      }
    }
  }
  else{
    dat <- PAC$summary[[summary_target[[1]]]]; labl <- summary_target
    }
  
  # Extract nt anno
  anno <-PAC$Anno
  anno$nuc_bias <- substr(rownames(anno), start=position, stop=position)
  cat(paste0("\nCounting nucleotides"))
  combin <- c(paste0(range[1]:range[2], "_A"), paste0(range[1]:range[2], "_T"), paste0(range[1]:range[2], "_C"), paste0(range[1]:range[2], "_G"), paste0(range[1]:range[2], "_N"))
  nuc_lst <- lapply(as.list(dat), function(x){
    nuc_agg <- aggregate(x, list(factor(paste(anno$Size, anno$nuc_bias, sep="_"))), sum)
    colnames(nuc_agg) <- c("position_nuc", "counts")
    zeros <- combin[!combin %in% as.character(nuc_agg$position_nuc)]
    if(length(zeros) >0){
      nuc_agg <- rbind(nuc_agg, data.frame(position_nuc=zeros, counts=0))}
    nuc_agg_ord <- nuc_agg[match(combin, nuc_agg$position_nuc),]
    rownames(nuc_agg_ord) <- NULL
    stopifnot(identical(as.character(nuc_agg_ord$position_nuc), combin))
    splt<- do.call("rbind", strsplit(as.character(nuc_agg_ord$position_nuc), split="_"))
    return(data.frame(length=splt[,1], nucleotide=splt[,2], counts=nuc_agg_ord[,2]))
  })
  
  #### Set options and load requirements
  options(scipen=999)
  #### Set up colors colors ###
  if(is.null(colors)){
    colfunc_sports <- grDevices::colorRampPalette(c("#094A6B", "#FFFFCC", "#9D0014"))
    colors <- colfunc_sports(5)
    colors <- c("#A0A0A0",colors[c(1,2,3,5)])
  }
  
  #### Plot ###                 
  histo_lst <- list(NA)
  if(is.null(summary_target)){samp <- rownames(PAC$Pheno)}
  else{samp <- names(PAC$summary[[summary_target[[1]]]])}
  for(i in 1:length(nuc_lst)){
    nuc_lst[[i]]$nucleotide <- factor(nuc_lst[[i]]$nucleotide, levels=c("N","C","G","A","T"))
    uni_chr_len <- as.integer(unique(as.character(nuc_lst[[i]]$length)))
    nuc_lst[[i]]$length <- factor(nuc_lst[[i]]$length, levels=uni_chr_len[order(uni_chr_len)] )
    histo_lst[[i]] <- ggplot2::ggplot(nuc_lst[[i]], ggplot2::aes(x=length, y=counts, fill=nucleotide))+
      ggplot2::geom_bar(width = 0.9, cex=0.2, colour="black", stat="identity")+
      ggplot2::geom_hline(yintercept=0, col="azure4")+
      ggplot2::xlab("Size (nt)")+
      ggplot2::ylab(paste(labl))+
      ggplot2::labs(subtitle = samp[i])+
      ggplot2::scale_fill_manual(values=colors)+
      #coord_cartesian(ylim=c(0,10000))+
      #scale_y_continuous(breaks = seq(0, 10000, 2500))+
      #ggthemes::geom_rangeframe(aes(x=range))+   
      ggthemes::theme_tufte()+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
    names(histo_lst)[i] <- names(nuc_lst)[i]
  }
  return(list(Histograms=histo_lst, Data=nuc_lst))
}
