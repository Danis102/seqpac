#' Filter a PAC object on sequence size and coverage
#'
#' \code{PAC_saturation} Perfoms an sequence diversity/saturation analysis on a
#' PAC objects.
#'
#' Given a PAC object the function will perform a sequence saturation analysis.
#' This is done by downsampling the original dataset by permutation at different
#' percentages of the original dataset. The closer the curve at the original
#' sequence depth (100%) is to a plateau phase the more saturated is the
#' diversity of sequences for the original dataset. Approaching the plateau,
#' usually means that the sequencing depth of the library have sampled the full
#' population of sequences available in the sample. Here we use an none-linear
#' least square (\code{nls}) model with a self-starter for asymptotic
#' regresssion (\code{SSasympt}) to describe the rate in which the library
#' approaches the plateau.
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
#' @param thresh Integer vector containing mean count thresholds that will
#'   be targeted. Default is set to c(1,10), where each new occurance (>=1) and
#'   each new occurance reaching 10 counts (>=10) will be analyzed.
#' @param threads Number of cores to be used for performing the permutations.
#'
#' @return A list with ggplot2 graph objects: The 1:st graph shows
#'   saturation/diversity result at the 1:st threshold. The 2:nd graph shows
#'   saturation/diversity result at the 2:nd threshold, etc.
#' @examples
#'
#' library(seqpac)
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Joint_analysis/R_analysis/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#' 
#' plot_lst  <- PAC_saturation(PAC, resample=10, steps=10, thresh=c(1,500), threads=8)
#' cowplot::plot_grid(plot_lst$A,plot_lst$B)
#'
#' @export
#'
PAC_saturation <- function(PAC, resample=10, steps=10, thresh=c(1,10), threads=4, par_type="PSOCK"){
              ###################### Setting up data ######################
                              require("foreach")
                              cat("\nOrganizing data...\n") 
                              df <- data.frame(seq=rownames(PAC$Counts), mean_counts=as.integer(rowMeans(PAC$Counts)))
                              rep_vect <- as.character(rep(df$seq, time=df$mean_counts))
                              size <- length(rep_vect)
                              x_vect <-  round(seq(0, 100, length.out=steps+1), digits=0)
                              x_vect <- x_vect[-1]
                              
                              cl <- parallel::makeCluster(threads, type = par_type)
                              doParallel::registerDoParallel(cl)
                              
              ###################### Generate data ######################
                              cat(paste0("\nPermutations will be generated using ", threads, " of ", parallel::detectCores(), " available core workers.")) 
                              cat(paste0("\nMaking ", resample, " permutations at ", steps," equally distributed sequence depths."))
                              cat(paste0("\nData will be saved for two coverage thresholds at >= ", thresh[1], " counts and >= ", thresh[2], " counts\n")) 
                              resampl_sub_lst <- foreach::foreach(i=1:length(x_vect), .inorder=FALSE) %dopar% {
                                                                            resampl_lst <- vector(mode = "list", length = resample)  
                                                                            resampl_lst <- lapply(resampl_lst, function(x){
                                                                                            res <- sample(rep_vect, size=round(size*(x_vect[i]*0.01), digit=0), replace =TRUE)
                                                                                            df_res <- as.data.frame(table(res))
                                                                                            lst <- lapply(as.list(thresh), function(z){
                                                                                                      seq <- as.character(df_res[df_res$Freq >= z, 1])
                                                                                                      seq_uni  <- as.character(unique(length(seq)))
                                                                                                      return(seq_uni)})
                                                                                            return(do.call("c", lst))})
                                                                            
                                                                            res <- as.data.frame(do.call("rbind", resampl_lst), stringsAsFactors=FALSE)
                                                                            res_lst <- as.list(res)
                                                                            df_lst <- list(NULL)
                                                                            for(t in 1:length(thresh)){df_lst[[t]] <- data.frame(value=res_lst[[t]], thresh=thresh[t], stringsAsFactors=FALSE)}         
                                                                            df <- do.call("rbind", df_lst)
                                                                            df$perc <- x_vect[i]
                                                                            return(df)
                              }
                              closeAllConnections()
                              dat <- as.data.frame(do.call("rbind", resampl_sub_lst))
                              ## Add intercept and fix classes  
                              intercpt <-  dat[!duplicated(paste0(dat$perc, dat$thresh)),]
                              intercpt[,!colnames(intercpt)=="thresh"] <- 0
                              dat <- rbind(dat, intercpt)
                              dat$thresh <- as.factor(dat$thresh)
                              dat$value <-  as.numeric(dat$value)
                              dat$perc <-  as.integer(dat$perc)

              ###################### Running nls and plotting results ######################
                              options(scipen=10)
                              cat("\nNow fitting the Asymptotic curve and plotting the results...\n")
                              plots <- lapply(as.list(thresh), function(thr){
                                                dat_sub <- dat[dat$thresh==thr,]
                                                
                                                ## Calculate Nonlinear Least Squares model ta obtain plateau rate using asymptotic regression
                                                
                                                p <- ggplot2::ggplot(dat_sub, ggplot2::aes(x = perc, y = value))+
                                                          ggplot2::geom_point()+
                                                          ggplot2::xlim(c(0, 200))+ 
                                                          ggplot2::ylim(c(0, max(dat[dat$thresh==thr,]$value)*1.5))+
                                                          ggplot2::xlab(paste0("Percent of original dataset (mean counts = ", size, ")"))+
                                                          ggplot2::ylab(paste0("Number of unique sequences reaching threshold"))+
                                                          ggplot2::ggtitle(paste0("Threshold >=", thr, " occurences"))+
                                                          ggplot2::geom_vline(xintercept=100, color="red")+
                                                          ggplot2::theme_bw()
                                                
                                                fm1 <- tryCatch(nls(value ~ SSasymp(perc, Asym, R0, lrc), data = dat_sub), error=function(x){return(NULL)})
                                                if(!is.null(fm1)){
                                                          y_predict <- predict(fm1)
                                                          good_fit <- round(cor(dat_sub$value, y_predict), digits=4)
                                                          plateau_rate <- round(exp(summary(fm1)$coefficients["lrc","Estimate"]), digits=4)
                                                          p <- p + ggplot2::stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, fullrange = TRUE)
                                                          return(p + annotate(geom = "text", x = 105, 
                                                                              y = max(dat_sub$value)*0.5, 
                                                                              label = paste0("Plateau rate=", plateau_rate, "\nGood.of.fit=", good_fit), 
                                                                              hjust = "left"))
                                                }else{          
                                                          p <- p + ggplot2::stat_smooth(method = "gam", formula = y ~ s(x), size = 1, fullrange=TRUE)
                                                          return(p + annotate(geom = "text", x = 105, 
                                                                              y = max(dat_sub$value)*0.5, 
                                                                              label = paste0("Asymptotic nls failed\nstat_smooth='gam', y ~ s(x)\nmodel was used instead"),
                                                                              hjust = "left"))
                                                          }
                                                                                    
                                                          ## Calculate linear model on slope based on size
                                                          # ggplot2::stat_smooth(method = "gam", formula = y ~ s(x), size = 1, fullrange=TRUE)+
                                                          # modl <- mgcv::gam(data=dat[dat$thresh==thr,], formula = value ~ s(perc,  bs = "cs"), method = "REML")
                                                          # y_predct <- unique(modl$fitted.values)
                                                          # y_predct <- y_predct[order(y_predct, decreasing=TRUE)]
                                                          # x_val <- unique(modl$model$perc)
                                                          # x_val <- x_val[order(x_val, decreasing=TRUE)]
                                                          # slp <- round((y_predct[2]-y_predct[1])/((size*(x_val[2]*0.01))-size), digits=2)
                                                          # return(p + annotate(geom = "text", x = 120, y = y_predct[1]*0.8, label = paste0("Size normalized\npredicted slope=", slp), hjust = "left"))
                                                }) 

                             return(plots)
                             parallel::stopCluster(cl)
                             options(scipen=0)
                            }

