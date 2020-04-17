#' Generates a size distribution analysis from sequences and counts in a PAC object
#'
#' \code{PAC_sizedist} analyses nucleotide bias.
#'
#' Given a PAC object the function will attempt to extract the ratios of
#' specific biotypes at a given position in sequences of the Anno data.frame
#' in relation to the sequence counts in Counts.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a count table with raw counts.
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
#' @return A list of objects: 
#'               1st object (Histograms::Samples): Individual histograms showing the nucleotide ratios per sample over the specified range.   
#'               2nd object (Data::Samples): Data used to generate the plots.
#'               3rd object (Stacked_bars::Groups): Stacked bars showing the mean ratios of each nucleotide per group over the specified range.
#'               4th object (Error_bars::Groups): Error bar plots with mean ratio of each nucleotide per group over the specified range. 
#'
#' @examples
#' 
#' library(seqpac)
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#'
#' ### First apply shallow counts filter, and then generate rpm on which rpm filter is applied.  
#' PAC_filt <- PAC_filter(PAC_all, size=c(16,45), threshold=10, coverage=5, type="counts", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#' PAC_filt <- PAC_rpm(PAC_filt)
#' PAC_filt <- PAC_filter(PAC_filt, size=c(16,45), threshold=10, coverage=5, type="rpm", stat=FALSE, pheno_target=NULL, anno_target=NULL)
#'
#' hierachy <- list( Mt_rRNA= "12S|16S|Mt_rRNA",
#'                  rRNA="5S|5.8S|18S|28S|S45|Ensembl_rRNA|rRNA_Other",
#'                  Mt_tRNA= "tRNA_mt-tRNA",
#'                  tRNA="Ensembl_tRNA|tRNA_nuc-tRNA",
#'                  miRNA="^miRNA|Ensembl_miRNA|Ensembl_pre_miRNA",
#'                  piRNA="piRNA")
#'                  
#' hierachy <- list(miRNA="^miRNA|Ensembl_miRNA|Ensembl_pre_miRNA",
#'                  Mt_rRNA= "12S|16S|Mt_rRNA",
#'                  rRNA="5S|5.8S|18S|28S|S45|Ensembl_rRNA|rRNA_Other",
#'                  Mt_tRNA= "tRNA_mt-tRNA",
#'                  tRNA="Ensembl_tRNA|tRNA_nuc-tRNA",
#'                  piRNA="piRNA")
#'                  
#'
#' PAC_filt <- simplify_reanno(PAC_filt, hierachy=hierachy, mismatches=0, bio_name="Biotypes_1", PAC_merge = TRUE)
#'
#' PAC_filt$Pheno$Groups <- paste(PAC_filt$Pheno$Batch, PAC_filt$Pheno$Method, PAC_filt$Pheno$Method_tag, PAC_filt$Pheno$Tag, sep="_")
#' 
#' 
#' ## Plot only histograms
#' sizedist_result <- PAC_sizedist(PAC_filt, anno_target=list("Biotypes_1", c("no_anno", "other", "miRNA", "tRNA", "Mt_tRNA", "rRNA", "Mt_rRNA",  "piRNA")), 
#' pheno_target=list("Groups", unique(PAC_filt$Pheno$Groups)))
#' 
#' sizedist_result <- PAC_sizedist(PAC=PAC_filt, anno_target=list("Biotypes_1"))
#'       
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_POOH__","Sep_Proto_Long_POOH_tag","Sep_Proto_Short_POOH_tag") 
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_POOH__","Jan_Proto_Long_POOH_tag","Jan_Proto_Short_POOH_tag")
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_IOR__","Jan_Proto_IOR__","Sep_TGIRT_POOH__")
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_Long_IOR_notag","Jan_Proto_Long_IOR_notag","Sep_Proto_Short_IOR_notag") 
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_Short_IOR_notag","Jan_TGIRT_IOR__","Sep_TGIRT_IOR__") 
#' filt <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_Short_IOR_notag","Sep_TGIRT_POOH__","Jan_R2_R2R_TGIRT_R2_R2R__") 
#' 
#' cowplot::plot_grid(plotlist=sizedist_result[[1]][filt], nrow = 3, ncol = 3)
#' 
#' cowplot::plot_grid(plotlist=sizedist_result[[1]][37:45], nrow = 3, ncol = 3)
#' 
#' 
#' 
#' @export

PAC_sizedist <- function(PAC, range=NULL, anno_target, pheno_target=NULL, colvect=NULL){
                    anno <- PAC$Anno
										counts <- PAC$Counts 		  
                    
										## Add range filter
										if(is.null(range)){range <- c(min(anno$Length), max(anno$Length))}
										filt <- anno$Length >= range[1] & anno$Length <= range[2] 
							      anno <- anno[filt,]
										counts <- counts[filt,]
										
										## Reomve unwanted biotypes
										if(length(anno_target)==1){ anno_target[[2]] <- as.character(unique(anno[,anno_target[[1]]]))}
										filt2 <- anno[,anno_target[[1]]] %in% anno_target[[2]]
							      anno <- anno[filt2,]
										counts <- counts[filt2,]
										
										#### Summarize over size and biotype
										
										bio_fact <- factor(anno[, anno_target[[1]]], levels=anno_target[[2]])
										seq_range <- seq(range[1], range[2])
										size_fact <- factor(anno$Length, levels=seq_range)
                    size_lst <- lapply(as.list(counts), function(x){
                                              bio_agg <- aggregate(x, list(bio_fact, size_fact), sum)
                                              colnames(bio_agg) <- c("biotype", "size", "counts")
                                              bio_agg_lst <- lapply(split(bio_agg, bio_agg$biotype), function(y){
                                                                  if(any(!seq_range %in% y$size )){
                                                                  y <- rbind(y, data.frame(biotype=as.character(unique(y$biotype)), size=seq_range[!seq_range %in% as.character(y$size)], counts=0))
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
                    require(RColorBrewer)
                    require(scales)
                    require(ggthemes)
                    require(extrafont)
                    require(ggplot2)  
                    #### Set up colors colors ###
                    if(is.null(colvect)){
                              colfunc <- grDevices::colorRampPalette(c("#094A6B", "#FFFFCC", "#9D0014"))
                              rgb <- colfunc(sum(!anno_target[[2]] %in% c("other", "no_anno")))
                              rgb_vec <- NULL
                              cnt <- 0
                              for(i in 1:length(anno_target[[2]])){ 
                                                  if(anno_target[[2]][i] == "other"){rgb_vec[i] <- "#808080";  cnt <- cnt+1}
                                                  if(anno_target[[2]][i] == "no_anno"){rgb_vec[i] <- "#C0C0C0";  cnt <- cnt+1}
                                                  if(!anno_target[[2]][i] %in% c("other","no_anno")){rgb_vec[i] <- rgb[i-cnt]}
                              }
                    }else{rgb_vec <- colvect}
                    #### Plot individual plots ###
										histo_lst <- list(NA)
										for(i in 1:length(size_lst)){
										                          if(is.null(pheno_target)){ph <- names(size_lst)[i]} 
										                          else {ph <- paste(PAC$Pheno[,pheno_target[[1]]][i], names(size_lst)[i])}
 										                          histo_lst[[i]] <- ggplot(size_lst[[i]], aes(x=size, y=counts, fill=biotype))+
                                                              	    geom_bar(width = 0.9, cex=0.2, colour="black", stat="identity")+
                                                              	    geom_hline(yintercept=0, col="azure4")+
                                                                  	xlab("Length (bp)")+
 										                                                ylab("Counts")+
                                                                    labs(subtitle = ph)+
                                                                  	scale_fill_manual(values=rgb_vec)+
                                                                  	#coord_cartesian(ylim=c(0,10000))+
                                                                  	#scale_y_continuous(breaks = seq(0, 10000, 2500))+
 										                                                #ggthemes::geom_rangeframe(aes(x=range))+   
                                                                    ggthemes::theme_tufte()+
                                                                  	theme_classic()+
                                                                  	theme(axis.text.x = element_text(angle = 0))
 			                                        names(histo_lst)[i] <- names(size_lst)[i]
										                          }
 			                return(list(Histograms=histo_lst, Data=size_lst))
                      }
#  			              if(is.null(pheno_target)){ph <- "NA"} else {ph <- PAC$Pheno[,pheno_target[[1]]]}
# 										if(!is.null(pheno_target[[2]])){ph <- ph[ph %in% pheno_target[[2]]]}
