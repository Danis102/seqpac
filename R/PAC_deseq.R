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
#' @param deseq_norm Logical whether to return deseq normalized values or not (default=TRUE).
#' 
#' @param histogram Logical whether to return a p-value destributions as a histogram (default=TRUE).
#' 
#' @param pheno_target (optional) List with: 
#'          1st object being a character indicating the main target column in
#'          Pheno.
#'          2nd object being a character vector of the target group(s) in the
#'          target Pheno column (1st object).
#'          
#'          Note: In PAC_deseq pheno_target controls the main comparative factor
#'          category. If for instance a target column named "groups" in
#'          PAC$Pheno contains "control" and "treatment" categories, setting
#'          pheno_target=list("groups", c("treatment", "controls") ensures that
#'          "treatment" is presented 1st in the factor levels, making for
#'          example log2FC appear as "treatment vs control". As default,
#'          pheno_target=NULL will result in the factor levels being
#'          automatically sorted in ascending order, which in the example above
#'          would result in control vs treatment log2FC.
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
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' PAC <- PAC_norm(pac, type="rpm")
#' pheno_target <- list("replicate", c("test3", "test1"))
#' pheno_target <- list("type", c("Sperm", "Ovaries"))
#' pheno_target <- list("type", c("Ovaries", "Sperm"))
#' 
#' model=~replicate
#' model=~replicate+type  
#'     
#'     
#' pheno_target <- list("replicate")
#' 
#' test <- PAC_deseq(PAC, model=model, deseq_norm=TRUE, histogram=TRUE, threads=4, pheno_target=pheno_target, anno_target=NULL)
#' 
#' @export

PAC_deseq <- function(PAC, model, deseq_norm=FALSE, histogram=TRUE, threads=1, pheno_target=NULL, anno_target=NULL){
  
  ### Subset samples and sequences
  PAC <-  PAC_filter(PAC, subset_only=TRUE, pheno_target=pheno_target, anno_target=anno_target)
  cat("\n")
  ### Prepare data
  anno <- PAC$Anno
  pheno <- PAC$Pheno
  if(!is.null(pheno_target)){
    trg <- pheno[,colnames(pheno) == pheno_target[[1]]]
    pheno[,colnames(pheno) == pheno_target[[1]]] <- factor(trg, levels=rev(pheno_target[[2]]))
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = PAC$Counts , colData = droplevels(PAC$Pheno), design=model)
  dds <- DESeq2::estimateSizeFactors(dds)
  
  ### DEseq function
  BiocParallel::register(BiocParallel::MulticoreParam(workers=threads))
  deseq <- function(dds, anno, model, norm){
    dds_fit <- DESeq2::DESeq(dds, fitType='local', parallel = TRUE)
    res_nam <- DESeq2::resultsNames(dds_fit)
    res_DESeq2 <- DESeq2::results(dds_fit, name=paste0(res_nam[2]))
    comp <- strsplit(res_DESeq2@elementMetadata@listData$description[2], ": replicate ")[[1]][2]
    cat(comp)
    cat(summary(res_DESeq2))
    res_DESeq2_df <- as.data.frame(res_DESeq2[order(as.numeric(res_DESeq2$pvalue)),])
    res_DESeq2_df_not_ord <- as.data.frame(res_DESeq2)
    anno_filt <- anno[match(rownames(res_DESeq2_df_not_ord), rownames(anno)),]
    if(identical(rownames(dds), rownames(res_DESeq2_df_not_ord))==FALSE){stop("Not identical ids in result and dds.\nHave you been messing with the code?\n")}
    if(identical(rownames(dds), rownames(anno_filt))==FALSE){stop("Not identical ids in anno and dds.\n")}
    Res_counts <- cbind(data.frame(feature_ID=rownames(dds)), res_DESeq2_df_not_ord, DESeq2::counts(dds, normalized=norm), anno_filt)
    colnames(Res_counts)[colnames(Res_counts)=="log2FoldChange"] <-  paste("logFC", res_nam[2], sep="_")
    return(Res_counts)
  }
  #--------------------------------------------------------
  # Apply the deseq function and add histograms:	
  require(DESeq2, quietly = TRUE) 
  res <- deseq(dds, anno, model, norm=deseq_norm)
  detach(package:DESeq2)
  p <- ggplot2::ggplot(data=res, aes(res$pvalue)) + 
    ggplot2::geom_histogram(breaks=seq(0.0, 1.0, by=0.025), col="black", fill="green", alpha=1) +
    ggplot2::labs(title="p-value destributions", x="p-value", y = "Number of features")
  res_lst<- list(result=res, histogram=p)
  return(res_lst)
}

