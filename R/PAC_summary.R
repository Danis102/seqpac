#' Summarizes PAC objects 
#'
#' \code{PAC_summary} summarizes data stored in a PAC object.
#'
#' Given a PAC object this function summarize data in PAC$Counts or in the norm
#' 'folder' according to a grouping columns in PAC$Pheno.
#' 
#' @family PAC analysis
#' 
#' @seealso \url{https://github.com/Danis102} for updates on the current
#' package.
#'
#' @param PAC PAC object containing a Pheno data.frame with samples as row
#'   names and a Counts table with raw counts. Optionally, the PAC object may
#'   contain normalized counts tables, saved in the norm list ('folder'). Such
#'   normalized table can be generated using the \code{\link{PAC_norm}}
#'   function.
#'   
#' @param norm Character indicating what type of data to be used. If 'counts' the raw
#'   counts in Counts will be used (default). Given any other value, the
#'   function will search for the value as a name on a data.frame stored in the
#'   normalized list-folder.
#'
#' @param type Character indicating what type of summary to be applied to the
#'   data. The function currently supports:
#'   type="means"         # Group means
#'   type="sd"            # Group standard deviation
#'   type="se"            # Group standard error of the mean
#'   type="log2FC"        # Group log2 fold changes against other groups
#'   type="log2FCgrand"   # Group log2 fold changes against a grand mean 
#'   type="percentgrand"  # Group log2 fold changes against a grand mean 

#'
#' @param pheno_target List with: 1st object being a character vector
#'   of target column in Pheno 2nd object being a character vector of the target
#'   group(s) in the target Pheno column (1st object).
#'   
#' @param rev Logical whether pairwise comparisions (e.g. log2FC) should be
#'   reversed (default=FALSE).
#'   
#' @param merge_pac Logical whether simplified annotation vector should
#'   automatically be added to the Anno object of the input PAC list object
#'   (default=TRUE). If \code{merge_pac=FALSE} a dataframe is returned. 
#'   
#' @return A PAC object with a pheno_summary folder containing the summarized
#'   data in a dataframe. The dataframe will be named according to the
#'   pheno_target, type and norm input.
#'
#' @examples
#' 
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' 
#' PAC_check(PAC_filt) # TRUE
#' 
#' # Easy to generate simple group summaries 
#' pac <- PAC_summary(pac, norm = "cpm", type = "means", pheno_target=list("stage"))       
#' pac <- PAC_summary(pac, norm = "cpm", type = "se", pheno_target=list("stage"))
#' pac <- PAC_summary(pac, norm = "cpm", type = "log2FC", pheno_target=list("stage"))
#' 
#' names(pac$summary)               # Names of individual summaries
#' head(pac$summary$cpmMeans_stage) # View individual individual summaries
#' tibble::as_tibble(pac$summary)  # View merge summaries
#' df <- as.data.frame(tibble::as_tibble(pac$summary))  # Merge multiple summaries
#' head(df)
#' 
#' 
#' # If a pheno_target is left out, a mean of all samples will be returned:
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE))
#' pac <- PAC_summary(pac, norm = "cpm", type = "mean")  
#' pac <- PAC_summary(pac, norm = "cpm", type = "percentgrand")
#' names(pac$summary)
#' tibble::as_tibble(pac$summary)   
#' 
#' @export
#' 
PAC_summary <- function(PAC, norm="counts", type="means", pheno_target=NULL, rev=FALSE, merge_pac=TRUE){
  
## Prepare and subset ################
  if(!is.null(pheno_target)){ 
    if(length(pheno_target)==1){ pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))
    }
  }
  PAC -> sav
  PAC <- PAC_filter(PAC, subset_only=TRUE, pheno_target=pheno_target, anno_target=NULL)

  ### Extract data ###
  if(norm=="counts"){
    data <- PAC$Counts
  }else{  
    if(is.null(PAC$norm[[norm]])){
      stop(paste0("There is no object called '", norm, "' in the norm list.\n  (Hint: Did you forget to normalize the data using for example PAC_rpm,\n  or would you rather run the function on raw counts using norm='counts'?)"))}  
    data <- PAC$norm[[norm]]
  }
  
  ### Subset dataset ###
  pheno <- PAC$Pheno
  pheno$All <- "All"
  if(is.null(pheno_target)){
    warning("No grouping factor was specified with pheno_target.\nCalculations are based on all samples.")
    pheno_target <- list("All","All")
  }else{
    if(length(pheno_target)==1){pheno_target[[2]] <- unique(pheno[, pheno_target[[1]]])
    }else{
      if(is.null(pheno_target[[2]])){pheno_target[[2]] <- unique(pheno[, pheno_target[[1]]])
      }}}
  indx <- pheno[, pheno_target[[1]]] %in% pheno_target[[2]]
  pheno <- pheno[indx,]
  data <- data[,indx]
  stopifnot(identical(rownames(pheno), colnames(data)))
  start_lst <- as.list(pheno_target[[2]])
  names(start_lst) <- pheno_target[[2]]
  sub_data_lst <- lapply(start_lst, function(x){ y <- data[, pheno[, pheno_target[[1]]] == x]; return(y)})
  
  ### Create pairwise combinations of pheno_target###
  
  if(length(sub_data_lst)>1){
    combn_lst <- as.list(data.frame(combn(1:length(sub_data_lst), m=2)))
    if(rev==TRUE){combn_lst <- lapply(combn_lst, function(x){c(x[2], x[1])})}
  }
  
  ### Apply log2FC to all pairwise combinations ###        
  if(type %in% c("log2FC", "Log2FC", "Log2Fc")){
    group_means <- lapply(sub_data_lst, function(x){ mns <- rowMeans(x)
    mns[mns == 0] <- 0.0000001
    return(mns)})
    log2FC_lst <- lapply(combn_lst, function(x){
      combn_data <- do.call("cbind", group_means[x])
      log2FC <- as.data.frame(log2(combn_data[,1]/combn_data[,2]))
      colnames(log2FC) <- paste(names(group_means)[x], collapse="_vs_")
      return(log2FC)
    })
    fin <- do.call("cbind", log2FC_lst)
    rownames(fin) <- rownames(data)
    df_nam <- paste0("Log2FC_", pheno_target[[1]])
  }
  if(type %in% c("log2FCgrand", "Log2FCgrand", "log2FCGrand")){
    grnd_mns <- rowMeans(data)
    grnd_mns[grnd_mns == 0] <- 0.0000001
    group_means <- lapply(sub_data_lst, function(x){ mns <- rowMeans(x)
    mns[mns == 0] <- 0.0000001
    return(mns)})
    log2FC_lst <- lapply(group_means, function(x){
      combn_data <- data.frame(group_means=x, grand_means=grnd_mns)
      log2FC <- as.data.frame(log2(combn_data[,1]/combn_data[,2]))
      #colnames(log2FC) <- paste(names(group_means)[x], collapse="_vs_grandMeans")
      return(log2FC)
    })
    fin <- do.call("cbind", log2FC_lst)
    rownames(fin) <- rownames(data)
    colnames(fin) <- paste0(names(log2FC_lst), "_vs_grandMeans")
    df_nam <- paste0("Log2FCgrand_", pheno_target[[1]])
  }
  
  ### Apply perc_change to all pairwise combinations ###       
  if(type %in% c("percent", "perc")){
    group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
    perc_diff_lst <- lapply(combn_lst, function(x){
      combn_data <- do.call("cbind", group_means[x])
      perc_diff <- as.data.frame((combn_data[,1]/combn_data[,2])*100)
      colnames(perc_diff) <- paste(names(group_means)[x], collapse="_vs_")
      return(perc_diff)
    })
    fin <- do.call("cbind", perc_diff_lst)
    rownames(fin) <- rownames(data)
    df_nam <- paste0("perc_", pheno_target[[1]])
  }
  if(type %in% c("percentgrand", "percgrand", "percentGrand", "percGrand")){
    grnd_mns <- rowMeans(data)
    group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
    perc_diff_lst <- lapply(group_means, function(x){
      combn_data <- data.frame(group_means=x, grand_means=grnd_mns)
      perc_diff <- as.data.frame((combn_data[,1]/combn_data[,2])*100)
      #colnames(perc_diff) <- paste(names(group_means)[x], collapse="_vs_")
      return(perc_diff)
    })
    fin <- do.call("cbind", perc_diff_lst)
    rownames(fin) <- rownames(data)
    colnames(fin) <- paste0(names(perc_diff_lst), "_vs_grandMeans")
    df_nam <- paste0("percGrand_", pheno_target[[1]])
  }
  
  
### means, sd and se ###  
  if(type %in% c("means", "Means", "mean", "Mean")){
    group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
    fin <- do.call("cbind", group_means)
    colnames(fin) <- names(group_means)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"Means_", pheno_target[[1]])
  }
  if(type %in% c("sd", "SD")){
    group_sd <- lapply(sub_data_lst, function(x){as.data.frame(apply(x, 1, stats::sd))})
    fin <- do.call("cbind", group_sd)
    colnames(fin) <- names(group_sd)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"SD_", pheno_target[[1]])
  }
  if(type %in% c("se", "SE")){
    group_se <- lapply(sub_data_lst, function(x){
      as.data.frame(apply(x, 1, function(y){stats::sd(y)/sqrt(length(y))}))
      }) 
    fin <- do.call("cbind", group_se)
    colnames(fin) <- names(group_se)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"SE_", pheno_target[[1]])
  }     
  
  ### Add more summary functions###
  # if(type=="perc_change"){}
  # if(type=="mean_diff"){}
  # etc
  ### Fix names and return object t###
  if(merge_pac==TRUE){
    PAC <- sav
    PAC$summary$new <- list(NULL)
    PAC$summary[[which(names(PAC$summary)=="new")]] <- as.data.frame(fin)
    names(PAC$summary)[which(names(PAC$summary)=="new")] <- df_nam
    stopifnot(PAC_check(PAC))
  return(PAC)
  }else{
    fin_lst <- list(fin)
    names(fin_lst) <- df_nam
    return(fin_lst)
  }
}



  