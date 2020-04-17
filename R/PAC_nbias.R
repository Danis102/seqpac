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
#' PAC_filt <- simplify_reanno(PAC_filt, hierachy=hierachy, mismatches=0, bio_name="Biotypes_1", PAC_merge = TRUE)
#'
#' PAC_filt$Pheno$Groups <- paste(do.call("rbind", strsplit(as.character(PAC_filt$Pheno$SampleProject), "_" ))[,1], PAC_filt$Pheno$Method, PAC_filt$Pheno$Method_tag, PAC_filt$Pheno$Tag, sep="_")
#' 
#' ## Plot only histograms
#' nbias_result <- PAC_nbias(PAC=PAC_filt, position=1, pheno_target=list("Groups", NULL))   
#' filt1 <- PAC_filt$Pheno$Groups %in% c("Sep_Proto_Reg__","Jan_Proto_Long_Tag","Jan_Proto_Short_Tag") 
#' cowplot::plot_grid(plotlist=nbias_result[[1]][filt1], nrow = 3, ncol = 3)
#' 
#' ## Generate graphs for a specific biotype 
#' nbias_result2 <- PAC_nbias(PAC=PAC_filt, position=1, anno_target=list("Biotypes_1", "piRNA"), pheno_target=list("Groups", unique(PAC_filt$Pheno$Groups)))   

#' cowplot::plot_grid(plotlist=nbias_result2[[1]][1:9], nrow = 3, ncol = 3)
#' cowplot::plot_grid(plotlist=nbias_result2[[1]][10:18], nrow = 3, ncol = 3)
#' cowplot::plot_grid(plotlist=nbias_result2[[1]][19:27], nrow = 3, ncol = 3)
#' cowplot::plot_grid(plotlist=nbias_result2[[1]][28:36], nrow = 3, ncol = 3)
#' cowplot::plot_grid(plotlist=nbias_result2[[1]][37:45], nrow = 3, ncol = 3)
#' 
#' ## Get group based plots
#' nbias_result2 <- PAC_nbias(PAC=PAC_filt, position=1, pheno_target=list("Groups", unique(PAC_filt$Pheno$Groups)))   
#' 
#' 
#' 
#' 
#' names(nbias_result[[3]])
#' 
#' @export

PAC_nbias <- function(PAC, position=1, range=NULL, anno_target=NULL, pheno_target=NULL){
										cat(paste0("Started: ", Sys.time(),"\n"))
										PAC -> PAC_work
										if(is.null(range)){range <- c(min(PAC_work$Anno$Length), max(PAC_work$Anno$Length))}
										cat(paste0("Filtering out position ", position , " in fragments ranging between ", range[1], "-", range[2], "bp.\n")) 
										filt <- PAC_work$Anno$Length >= range[1] & PAC_work$Anno$Length <= range[2] 
							      if(!is.null(anno_target)){
										              bio_anno <- as.character(PAC_work$Anno[, anno_target[[1]]])
										              filt_bio <- bio_anno %in% anno_target[[2]]
										              filt <- filt_bio + filt == 2
										              }
                    PAC_work$Anno <- PAC_work$Anno[filt,]
										PAC_work$Counts <- PAC_work$Counts[filt,]
										
										#### Counting nucs
										PAC_work$Anno$nuc_bias <- substr(rownames(PAC_work$Anno), start=position, stop=position)
                    cat(paste0("Counting nucleotides"))
                    combin <- c(paste0(range[1]:range[2], "_A"), paste0(range[1]:range[2], "_T"), paste0(range[1]:range[2], "_C"), paste0(range[1]:range[2], "_G"), paste0(range[1]:range[2], "_N"))
										nuc_lst <- lapply(as.list(PAC_work$Counts), function(x){
										                          nuc_agg <- aggregate(x, list(factor(paste(PAC_work$Anno$Length, PAC_work$Anno$nuc_bias, sep="_"))), sum)
										                          colnames(nuc_agg) <- c("position_nuc", "counts")
										                          nuc_agg <- rbind(nuc_agg, data.frame(position_nuc=combin[!combin %in% as.character(nuc_agg$position_nuc)], counts=0))
										                          nuc_agg_ord <- nuc_agg[match(combin, nuc_agg$position_nuc),]
										                          rownames(nuc_agg_ord) <- NULL
										                          stopifnot(identical(as.character(nuc_agg_ord$position_nuc), combin))
										                          splt<- do.call("rbind", strsplit(as.character(nuc_agg_ord$position_nuc), split="_"))
										                          return(data.frame(length=splt[,1], nucleotide=splt[,2], counts=nuc_agg_ord[,2]))
										                          })
										
										#### Generate individial histograms 
                            options(scipen=999)
                            require(RColorBrewer)
                            require(scales)
                            require(ggthemes)
                            require(cowplot)
                            require(extrafont)
                            ###############################################################################################################
                            #### Set up colors colors ###
                                      colfunc_sports <- grDevices::colorRampPalette(c("#094A6B", "#FFFFCC", "#9D0014"))
                                      rgb_vec_sports <- colfunc_sports(5)
										require(ggplot2)   
										histo_lst <- list(NA)
										for(i in 1:length(nuc_lst)){
										                          nuc_lst[[i]]$nucleotide <- factor(nuc_lst[[i]]$nucleotide, levels=c("A","T","C","G","N"))
										                          if(is.null(pheno_target)){ph <- "NA"} else {ph <- PAC_work$Pheno[,pheno_target[[1]]]}
										                          if(!is.null(pheno_target[[2]])){ph <- ph[ph %in% pheno_target[[2]]]}                        
 										                          histo_lst[[i]] <- ggplot(nuc_lst[[i]], aes(x=length, y=counts, fill=nucleotide))+
                                                              	    geom_bar(width = 0.9, cex=0.2, colour="black", stat="identity")+
                                                              	    geom_hline(yintercept=0, col="azure4")+
                                                                  	xlab("Length (bp)")+
 										                                                ylab("Nucleotide counts")+
                                                                    labs(subtitle = ph[i])+
                                                                  	scale_fill_manual(values=rgb_vec_sports)+
                                                                  	#coord_cartesian(ylim=c(0,10000))+
                                                                  	#scale_y_continuous(breaks = seq(0, 10000, 2500))+
 										                                                #ggthemes::geom_rangeframe(aes(x=range))+   
                                                                    ggthemes::theme_tufte()+
                                                                  	theme_classic()+
                                                                  	theme(axis.text.x = element_text(angle = 0))
 			                                        names(histo_lst)[i] <- names(nuc_lst)[i]
										                        }
                      return(list(Histograms=histo_lst, Data=nuc_lst))
                  }
