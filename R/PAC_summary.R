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
#' @param norm Character indicating what type of data to be used. If raw the raw
#'   counts in Counts will be used (default="raw"). Given any other value, the
#'   function will search for the value as a name on a dataframe stored in the
#'   normalized list-folder (created for example by PAC_rpm).
#'
#' @param type Character indicating what type of summary to be applied to the
#'   data.
#'
#' @param pheno_target List with: 1st object being a character vector
#'   of target column in Pheno 2nd object being a character vector of the target
#'   group(s) in the target Pheno column (1st object).
#'   
#' @return A PAC object with a pheno_summary folder containing the summarized
#'   data in a dataframe. The dataframe will be named according to the
#'   pheno_target, type and norm input.
#'
#' @examples
#' reanno_path="/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/Processed_Pipeline31_05-03-20/R_files/"
#' load(file=paste0(reanno_path, "PAC_filt_rpm10in25.Rdata"))
#' 
#' PAC_test <- PAC_summary(PAC_filt, norm = "rpm", type = "log2FC", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' PAC_test <- PAC_summary(PAC_test, norm = "rpm", type = "means", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' PAC_test <- PAC_summary(PAC_test, norm = "rpm", type = "se", pheno_target=list("Index", c("sperm_cells_HC", "sperm_cells_CT")))
#' 
#' 
#' @export
#' 
PAC_summary <- function(PAC, norm="raw", type="log2FC", pheno_target){

                              ### Extract data ###
                                    if(norm=="raw"){
                                            data <- PAC$Counts
                                    }else{
                                            data <- PAC$norm[[norm]]
                                    }

                              ### Subset dataset ###
                                    pheno <- PAC$Pheno
                                    indx <- pheno[, pheno_target[[1]]] %in% pheno_target[[2]]
                                    pheno <- pheno[indx,]
                                    data <- data[,indx]
                                    stopifnot(identical(rownames(pheno), colnames(data)))
                                    start_lst <- as.list(pheno_target[[2]])
                                    names(start_lst) <- pheno_target[[2]]
                                    sub_data_lst <- lapply(start_lst, function(x){ y <- data[, pheno[, pheno_target[[1]]] == x]; return(y)})

                              ### Create pairwise combinations of pheno_target###
                                    combn_lst <- as.list(data.frame(combn(1:length(sub_data_lst), m=2)))
                                    PAC$summary$new <- list(NULL)
                              ### Apply log2FC to all pairwise combinations ###        
                                    if(type=="log2FC"){
                                              group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
                                              log2FC_lst <- lapply(combn_lst, function(x){
                                                                                    combn_data <- do.call("cbind", group_means[x])
                                                                                    log2FC <- as.data.frame(log2(combn_data[,1]/combn_data[,2]))
                                                                                    colnames(log2FC) <- paste(names(group_means)[x], collapse="_vs_")
                                                                                    return(log2FC)
                                                                                    })
                                              fin <- do.call("cbind", log2FC_lst)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("log2FC_[", pheno_target[[1]], "]")
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
                                              df_nam <- paste0("perc_diff_[", pheno_target[[1]], "]")
                                  }
                                                                  ### Apply perc_change to all pairwise combinations ###       
                                 if(type=="means"){
                                              group_means <- lapply(sub_data_lst, function(x){ as.data.frame(rowMeans(x))})
                                              fin <- do.call("cbind", group_means)
                                              colnames(fin) <- names(group_means)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("means_[", pheno_target[[1]], "]")
                                 }
                                 if(type=="sd"){
                                              group_sd <- lapply(sub_data_lst, function(x){as.data.frame(apply(x, 1, stats::sd))})
                                              fin <- do.call("cbind", group_sd)
                                              colnames(fin) <- names(group_sd)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("sd_[", pheno_target[[1]], "]")
                                 }
                                 if(type=="se"){
                                              group_se <- lapply(sub_data_lst, function(x){as.data.frame(apply(x, 1, function(y){stats::sd(y)/sqrt(length(y))}))}) 
                                              fin <- do.call("cbind", group_se)
                                              colnames(fin) <- names(group_se)
                                              rownames(fin) <- rownames(data)
                                              df_nam <- paste0("se_[", pheno_target[[1]], "]")
                                 }     
  
                              ### Add more summary functions###
                                    # if(type=="perc_change"){}
                                    # if(type=="mean_diff"){}
                                    # etc
                              ### Fix names and return object t###
                              PAC$summary[[which(names(PAC$summary)=="new")]] <- fin
                              names(PAC$summary)[which(names(PAC$summary)=="new")] <- df_nam
                              return(PAC)
                        }



  