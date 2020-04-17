#' Filter a PAC object on sequence size and covarage
#'
#' \code{PAC_saturation} Perfoms an sequence diversity/saturation analysis on a
#' PAC objects.
#'
#' Given a PAC object the function will perform a sequence saturation analysis.
#' This is done by downsampling the original dataset by permutation at different
#' percentages of the original dataset. The flatter the curve at the original
#' sequence depth (100%) the more saturated is the diversity of sequences at the
#' original dataset.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#' @param PAC PAC-list object containing a Counts data.frame with sequences as
#'   row names and samples as column names.
#' @param resample Integer setting the number of permutations at each percentage
#'   step (default=10).
#' @param steps Integer the number of percentage steps between 0-100% of the
#'   original dataset (default=10).
#' @param thresh Integer vector containing two mean count thresholds that will
#'   be targeted. Default is set to c(1,500), where each new occurance (>=1) and
#'   each new occurance reaching 500 counts (>=500) will be analyzed.
#' @param cumulative (optional) Logical. FALSE make the analysis on total unique
#'   sequences when resampling at different percentages (default). TRUE make the
#'   analysis on unique sequences gained at each step compared to a given
#'   downsampled dataset given in 'start_perc' (default = 1%). While the integer
#'   in the 'resample' input sets how many (n) permutations to be made at each
#'   step, permutation with FALSE will permute n times at each percentage step,
#'   while TRUE will permute 1 time per step and then repeat this n times. The
#'   two approaches gives similar results.
#' @param start_perc (optional) Logical specifying if an coverage graph should
#'   be generated or not (default=FALSE).
#'
#' @param threads Number of cores to be used for performing the permutations.
#'
#' @return A list with two graph objects: The 1:st graph (A) shows
#'   saturation/diversity result at the 1:st threshold. The 2:nd graph (B) shows
#'   saturation/diversity result at the 2:nd threshold.
#' @examples
#'
#' library(seqpac)
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#' 
#' plot_lst  <- PAC_saturation(PAC, resample=10, steps=10, thresh=c(1,500), cumulative=TRUE, start_perc=1, threads=8)
#' cowplot::plot_grid(plot_lst$A,plot_lst$B)
#'
#' @export
#'
PAC_saturation <- function(PAC, resample=10, steps=10, thresh=c(1,500), cumulative=FALSE, start_perc=1, threads=8){
              ###################### Setting up data ######################
                              cat("\nOrganizing data...\n") 
                              df <- data.frame(seq=rownames(PAC$Counts), mean_counts=as.integer(rowMeans(PAC$Counts)))
                              rep_vect <- as.character(rep(df$seq, time=df$mean_counts))
                              size <- length(rep_vect)
                              x_vect <-  round(seq(0, 100, length.out=steps+1), digits=0)
                              x_vect <- x_vect[-1]
                              
              ###################### No start point ######################
                              if(cumulative==FALSE){
                                          ## Permuting over x_vect using thresh
                                              cat(paste0("\nPermutations will be generated using ", threads, " of ", parallel::detectCores(), " available core workers.")) 
                                              cat(paste0("\nMaking ", resample, " permutations at ", steps," equally destributed sequence depths."))
                                              cat(paste0("\nData will be saved for two coverage thresholds at >= ", thresh[1], " counts and >= ", thresh[2], " counts\n")) 
                                              pb = txtProgressBar(min = 0, max = length(x_vect), initial = 0,  style = 3, width=40)
                                              permuted_data <- list(NULL)  
                                              for (i in 1:length(x_vect)){
                                                  resampl_lst <- as.list(rep(x_vect[i], time=resample))
                                                        resampl_sub_lst <- parallel::mclapply(resampl_lst, mc.cores=threads, function(resampl){
                                                                                  resampl_data <- sample(rep_vect, size=round(size*(resampl*0.01), digit=0), replace =FALSE)
                                                                                  df_res <- as.data.frame(table(resampl_data))
                                                                                  res_thresh_A<- table(df_res$Freq>=thresh[1])
                                                                                  res_thresh_B<- table(df_res$Freq>=thresh[2])
                                                                                  df <- data.frame(A=0, B=0)
                                                                                  if(length(res_thresh_A)==2|length(res_thresh_A)==1 && !names(res_thresh_A)==FALSE){df$A  <- res_thresh_A["TRUE"]}
                                                                                  if(length(res_thresh_B)==2|length(res_thresh_B)==1 && !names(res_thresh_B)==FALSE){df$B  <- res_thresh_B["TRUE"]}
                                                                                  return(df)
                                                                        })
                                                          permuted_data[[i]] <- do.call("rbind", resampl_sub_lst)
                                                          names(permuted_data)[i] <- x_vect[i]
                                                          setTxtProgressBar(pb,i)
                                              }
                                                data_A <- reshape2::melt(do.call("cbind", lapply(permuted_data, function(x){x$A})))
                                                data_A <- rbind(data_A, data.frame(Var1=NA, Var2=0, value=0))
                                                data_A$variable <- "A"
                                                data_B <- reshape2::melt(do.call("cbind", lapply(permuted_data, function(x){x$B}))) 
                                                data_B <- rbind(data_B, data.frame(Var1=NA, Var2=0, value=0))
                                                data_B$variable <- "B"
                                                data <- rbind(data_A, data_B)
                                                names(data)[2] <- "perc"
                                          }
                                            
                             
              ###################### Generate cumulative new sequences from a start point depth ######################
                              if(cumulative==TRUE){
                                          cat(paste0("\nPermutations will be generated using ", threads, " of ", parallel::detectCores(), " available core workers.")) 
                                          cat(paste0("\nMaking ", resample, " permutations at ", steps," equally destributed sequence depths."))
                                          cat(paste0("\nData will be saved for two coverage thresholds at >= ", thresh[1], " counts and >= ", thresh[2], " counts\n")) 
                                          permuted_data <- list(NULL)
                                          pb = txtProgressBar(min = 0, max = length(x_vect), initial = 0,  style = 3, width=40)
                                          for (i in 1:resample){
                                                    resampl_lst <- as.list(x_vect)
                                                    names(resampl_lst) <- x_vect
                                                    resampl_start <- sample(rep_vect, size=round(size*(start_perc*0.01), digit=0), replace =FALSE)
                                                    start <- as.data.frame(table(resampl_start))
                                                    start_A <- as.character(start[start$Freq >= thresh[1], 1])
                                                    start_B <- as.character(start[start$Freq >= thresh[2], 1])
                                                    resampl_sub_lst <- parallel::mclapply(resampl_lst, mc.cores=threads, function(resampl){
                                                                                              resampl_data <- sample(rep_vect, size=round(size*(resampl*0.01), digit=0), replace =TRUE)
                                                                                              df_res <- as.data.frame(table(resampl_data))
                                                                                              seq_A <- as.character(df_res[df_res$Freq >= thresh[1], 1])
                                                                                              seq_B <- as.character(df_res[df_res$Freq >= thresh[2], 1])
                                                                                              seq_A_res <- table(!seq_A %in% start_A)
                                                                                              seq_B_res <- table(!seq_B %in% start_B)
                                                                                              df <- data.frame(A=0, B=0)
                                                                                              if(length(seq_A_res)==2|length(seq_A_res)==1 && !names(seq_A_res)==FALSE){df$A  <- seq_A_res["TRUE"]}
                                                                                              if(length(seq_B_res)==2|length(seq_B_res)==1 && !names(seq_B_res)==FALSE){df$B  <- seq_B_res["TRUE"]}
                                                                                              return(df)
                                                                                      })
                                                    permuted_data[[i]] <- do.call("rbind", resampl_sub_lst)
                                                    setTxtProgressBar(pb,i)
                                          }
                                          data <-  suppressMessages(reshape2::melt(cbind(data.frame(perc=rownames(do.call("cbind", permuted_data))), do.call("cbind", permuted_data))))
                                          data$perc <- as.numeric(as.character(data$perc))
                              }
                              
              ###################### Plotting results ######################
                                                options(scipen=10)
                                                cat("\nNow plotting the results...\n")
                                                A <- ggplot2::ggplot(data[data$variable=="A",], ggplot2::aes(x = perc, y = value))+
                                                                            ggplot2::geom_point()+
                                                                            ggplot2::xlim(c(0, 200))+                          
                                                                            ggplot2::stat_smooth(method = "gam", formula = y ~ s(x), size = 1, fullrange=TRUE)+
                                                                            ggplot2::ylim(c(0, max(data[data$variable=="A",]$value)*1.5))+
                                                                            ggplot2::xlab(paste0("Percent of original dataset (mean counts = ", size, ")"))+
                                                                            ggplot2::ylab(paste0("Number of unique sequences reaching threshold"))+
                                                                            ggplot2::ggtitle(paste0("Threshold >=", thresh[1], " occurences"))+
                                                                            ggplot2::geom_vline(xintercept=100, color="red")+
                                                                            ggplot2::theme_bw() 
                                                                            
                                               B <- ggplot2::ggplot(data[data$variable=="B",], ggplot2::aes(x = perc, y = value))+
                                                                            ggplot2::geom_point()+
                                                                            ggplot2::xlim(c(0, 200))+                          
                                                                            ggplot2::stat_smooth(method = "gam", formula = y ~ s(x), size = 1, fullrange=TRUE)+
                                                                            ggplot2::ylim(c(0, max(data[data$variable=="B",]$value)*1.5))+
                                                                            ggplot2::xlab(paste0("Percent of original dataset (mean counts = ", size, ")"))+
                                                                            ggplot2::ylab(paste0("Number of unique sequences reaching threshold"))+
                                                                            ggplot2::ggtitle(paste0("Threshold >=", thresh[2], " occurences"))+
                                                                            ggplot2::geom_vline(xintercept=100, color="red")+
                                                                            ggplot2::theme_bw()
                                               
                                               plot_lst <- list(A=A, B=B)
                                               return(plot_lst)
                                               print(cowplot::plot_grid(plot_lst$A,plot_lst$B))
                                               options(scipen=0)
                            }