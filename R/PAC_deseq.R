#' A wrapper to DESeq2 that apply differential expression analysis on a PAC object  
#'
#' \code{PAC_deseq} DESeq2 analysis on PAC_object.
#'
#' Given a PAC object this function will apply a differential expression analysis using DESeq2.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC object containing a Pheno data.frame with samples as row
#'   names, an Anno data.frame with sequences as row names and a Counts table with
#'   raw counts. The Counts table must have the sample names as column names and
#'   sequences as row names.
#'
#' @param model Character vector describing the statistical model based on column names in Pheno.
#'   
#' @param main_factor Character vector (default = first term in model).
#' 
#' @param norm Logical whether to return deseq normalized values or not (default=TRUE).
#' 
#' @param histogram Logical whether to return a p-value destributions as a histogram (default=TRUE).
#' 
#' @param pheno_target (optional) List with: 
#'          1st object being a character vector of target column in Pheno 
#'          2nd object being a character vector of the target group(s) in the target Pheno column (1st object).
#'          (default=NULL)
#'          
#' @param anno_target (optional) List with: 
#'          1st object being a character vector of target column in Anno 
#'          2nd object being a character vector of the target type/biotypes(s) in the target Anno column (1st object).
#'          (default=NULL)
#'
#' @param threads Integer number of threads to run in parallel.          
#'   
#' @return A list of objects: 
#'               1st object - result table 
#'               2nd object - p-value histogram (optional)    
#' @examples
#' 
#' load(file="/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/Processed_Pipeline3_02-01-20/R_files/PAC_master_reanno.Rdata")
#' PAC <- PAC_rpm(PAC_master_reanno)
#' PAC$Pheno$Hemi_cast <- ifelse(grepl("HC", PAC$Pheno$Index), "Hemi", "Cont") 
#' PAC$Pheno$Type <- ifelse(grepl("sperm_cells", PAC$Pheno$Index), "Sperm", "Plasma")
#' pheno_target <- list("Type", "Sperm")
#' PAC  <- PAC_filter(PAC, size = c(16, 45), threshold = 1, coverage = 100, type = "rpm", stat = TRUE, pheno_target = pheno_target, anno_target = NULL)
#' pheno_target <- list("Hemi_cast", c("Hemi", "Cont"))  
#' 
#' test <- PAC_deseq(PAC, model=~Hemi_cast, main_factor="Hemi_cast", norm=TRUE, histogram=TRUE, threads=4, pheno_target=NULL, anno_target=NULL)
#' 
#' @export

PAC_deseq <- function(PAC, model, main_factor=NULL, norm=TRUE, histogram=TRUE, threads=1, pheno_target=NULL, anno_target=NULL){
                            ### Subset samples by Pheno                                                
                              							if(!is.null(pheno_target)){
                                              sub_pheno <- as.character(PAC$Pheno[pheno_target[[1]]][,1]) %in% pheno_target[[2]]
                                              if(any(names(PAC)=="rpm")){ PAC$rpm  <- PAC$rpm[,sub_pheno] }
                                              PAC$Counts   <- PAC$Counts [,sub_pheno]
                                              PAC$Pheno  <- PAC$Pheno[sub_pheno,]
                                              tab_pheno <- as.data.frame(table(sub_pheno))
                                              cat(paste0("\nPheno filter was specified, will retain: ", tab_pheno[tab_pheno[,1]==TRUE, 2], " of ", length(sub_pheno), " samples\n"))                                            ### Subset data by groups
                              							}              
                            ### Subset data by Anno
                                           if(!is.null(anno_target)){
                                              sub_anno <- as.character(PAC$Anno[anno_target[[1]]][,1]) %in% anno_target[[2]]
                                              if(any(names(PAC)=="rpm")){ PAC$rpm  <- PAC$rpm[sub_anno,] }
                                              PAC$Counts   <- PAC$Counts [sub_anno,] 
                                              PAC$Anno  <- PAC$Anno[sub_anno,]
                                              tab_anno <- as.data.frame(table(sub_anno))
                                              cat(paste0("\nAnno filter was specified, will retain: ", tab_anno[tab_anno[,1]==TRUE, 2], " of ", length(sub_anno), " seqs\n"))                                                ### Subset data by groups
                                           }
  
                            ### Prepare data
                                            anno <- PAC$Anno
                                            pheno <- PAC$Pheno
                                            if(!is.null(pheno_target)){
                                                            groups_ord <- pheno_target[[2]]
                                                            trg <- as.character(pheno[, colnames(pheno) == main_factor])
                                                            pheno[,colnames(pheno) == main_factor] <- factor(trg, levels=rev(pheno_target[[2]]))
                                                            }
                                            
                                            dds <- DESeq2::DESeqDataSetFromMatrix(countData = PAC$Counts , colData = droplevels(PAC$Pheno), design=model)
                                            dds <- DESeq2::estimateSizeFactors(dds)
                                            
                            ### DEseq function
                                            BiocParallel::register(BiocParallel::MulticoreParam(workers=threads))
                                					  deseq <- function(dds, anno, model, norm){
                                    												dds_fit <- DESeq2::DESeq(dds, fitType='local', parallel = TRUE)
                                    												res_nam <- DESeq2::resultsNames(dds_fit)
                                    												res_DESeq2 <- DESeq2::results(dds_fit, name=paste0(res_nam[2]))
                                    												cat(summary(res_DESeq2))
                                    												res_DESeq2_df <- as.data.frame(res_DESeq2[order(as.numeric(res_DESeq2$pvalue)),])
                                    												res_DESeq2_df_not_ord <- as.data.frame(res_DESeq2)
                                    												anno_filt <- anno[match(rownames(res_DESeq2_df_not_ord), rownames(anno)),]
                                    												if(identical(rownames(dds), rownames(res_DESeq2_df_not_ord))==FALSE){stop("Error! Not identical ids in result and dds.\nHave you been messing with the code?\n")}
                                    												if(identical(rownames(dds), rownames(anno_filt))==FALSE){stop("Error! Not identical ids in anno and dds.\n")}
                                    												Res_counts <- cbind(data.frame(feature_ID=rownames(dds)), res_DESeq2_df_not_ord, DESeq2::counts(dds, normalized=norm), anno_filt)
                                    												colnames(Res_counts)[colnames(Res_counts)=="log2FoldChange"] <-  paste("logFC", res_nam[2], sep="_")
                                    												return(Res_counts)
                        												}
      									#--------------------------------------------------------
      									# Apply the deseq function and add histograms:	
      									res <- deseq(dds, anno, model, norm)
      									p <- ggplot2::ggplot(data=res, aes(res$pvalue)) + 
        														ggplot2::geom_histogram(breaks=seq(0.0, 1.0, by=0.025), col="black", fill="green", alpha=1) +
      															ggplot2::labs(title="p-value destributions", x="p-value", y = "Number of features")
      											   res_lst<- list(result=res, histogram=p)
      								  return(res_lst)
      }

