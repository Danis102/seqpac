#' Summarizes PAC objects 
#'
#' \code{PAC_summary} Summarizes data stored in a PAC object.
#'
#' Given a PAC object this function summarize data in Counts or in the norm
#' list-folder according to a grouping variable in Pheno.
#' 
#' @family PAC analysis
#' 
#' @seealso \url{https://github.com/Danis102} for updates on the current
#' package.
#'
#' @param PAC PAC object containing a Pheno dataframe with samples as row names,
#'   a Counts table with raw counts and a normalized list-folder containing
#'   tables with normalized Counts (e.g. rpm).
#'
#' @param norm Character indicating what type of data to be used. If 'counts' the raw
#'   counts in Counts will be used (default). Given any other value, the
#'   function will search for the value as a name on a dataframe stored in the
#'   normalized list-folder (created for example by PAC_norm).
#'
#' @param type Character indicating what type of summary to be applied to the
#'   data.
#'
#' @param pheno_target List with: 1st object being a character vector
#'   of target column in Pheno 2nd object being a character vector of the target
#'   group(s) in the target Pheno column (1st object).
#'   
#' @param rev Logical whether pairwise comparisions should be reversed
#'   (default=FALSE).
#'   
#' @param PAC_merge Logical whether simplified annotation vector should
#'  automatically be added to Anno object if a PAC list object were given as
#'  input (default=FALSE)
#'
#' @return A PAC object with a pheno_summary folder containing the summarized
#'   data in a dataframe. The dataframe will be named according to the
#'   pheno_target, type and norm input.
#'
#' @examples
#' load(file="/home/danis31/OneDrive/Programmering/Programmering/Pipelines/Drosophila/Pipeline_3.1/seqpac/dm_test_PAC.Rdata")
#' 
#' 
#' PAC_check(PAC_filt) # TRUE
#'
#' PAC_filt <- PAC_rpm(PAC_filt)
#'
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "means", pheno_target=list("Method", NULL))
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "log2FC", pheno_target=list("Method", NULL))
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "log2FCgrand", pheno_target=list("Method"))
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "log2FCgrand")
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "percentgrand")
#' 
#' PAC_test <- PAC_summary(PAC_filt, norm = "rpm", type = "log2FC", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' PAC_test <- PAC_summary(PAC_test, norm = "rpm", type = "means", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' PAC_test <- PAC_summary(PAC_test, norm = "rpm", type = "se", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' 
#' PAC_filt <- PAC_summary(PAC_filt, norm = "rpm", type = "means", pheno_target=list("Method"))
#' 
#' 
#' 
#' 
#' @export
#' 
PAC_summary <- function(PAC, norm="counts", type="means", pheno_target=NULL, rev=FALSE, PAC_merge=FALSE){

                              ### Extract data ###
                                    if(norm=="counts"){
                                            data <- PAC$Counts
                                    }else{  
                                          if(is.null(PAC$norm[[norm]])){stop(paste0("There is no object called '", norm, "' in the norm list.\n  (Hint: Did you forget to normalize the data using for example PAC_rpm,\n  or would you rather run the function on raw counts using norm='counts'?)"))}  
                                          data <- PAC$norm[[norm]]
                                    }
                              
                              ### Subset dataset ###
                                    pheno <- PAC$Pheno
                                    pheno$All <- "All"
                                    if(is.null(pheno_target)){warning("No grouping factor was specified with pheno_target.\nCalculations are based on all samples.")
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

                                              PAC$summary$new <- list(NULL)
                              ### Apply log2FC to all pairwise combinations ###        
                                    if(type=="log2FC"){
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
                                    if(type=="log2FCgrand"){
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
                                  if(type=="percent"){
                                              group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
                                              perc_diff_lst <- lapply(combn_lst, function(x){
                                                                                    combn_data <- do.call("cbind", group_means[x])
                                                                                    perc_diff <- as.data.frame((combn_data[,1]/combn_data[,2])*100)
                                                                                    colnames(perc_diff) <- paste(names(group_means)[x], collapse="_vs_")
                                                                                    return(perc_diff)
                                                                                    })
                                              fin <- do.call("cbind", perc_diff_lst)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("Percdiff_", pheno_target[[1]])
                                  }
                                 if(type=="percentgrand"){
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
                                              df_nam <- paste0("Percdiffgrand_", pheno_target[[1]])
                                  }
                                              
                                              
                                      
                                 if(type=="means"){
                                              group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
                                              fin <- do.call("cbind", group_means)
                                              colnames(fin) <- names(group_means)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0(norm,"Means_", pheno_target[[1]])
                                 }
                                 if(type=="sd"){
                                              group_sd <- lapply(sub_data_lst, function(x){as.data.frame(apply(x, 1, stats::sd))})
                                              fin <- do.call("cbind", group_sd)
                                              colnames(fin) <- names(group_sd)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("Sd_", pheno_target[[1]])
                                 }
                                 if(type=="se"){
                                              group_se <- lapply(sub_data_lst, function(x){as.data.frame(apply(x, 1, function(y){stats::sd(y)/sqrt(length(y))}))}) 
                                              fin <- do.call("cbind", group_se)
                                              colnames(fin) <- names(group_se)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("Se_", pheno_target[[1]])
                                 }     
  
                              ### Add more summary functions###
                                    # if(type=="perc_change"){}
                                    # if(type=="mean_diff"){}
                                    # etc
                              ### Fix names and return object t###
                              if(PAC_merge==TRUE){PAC$summary[[which(names(PAC$summary)=="new")]] <- fin
                                                  names(PAC$summary)[which(names(PAC$summary)=="new")] <- df_nam
                                                  return(PAC)
                              }else{
                                                  
                                                  fin_lst <- list(fin)
                                                  names(fin_lst) <- df_nam
                                                  return(fin_lst)
                              }
                        }



  