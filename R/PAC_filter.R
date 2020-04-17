#' Filter a PAC object on sequence size and covarage  
#'
#' \code{PAC_filter} Filter PAC objects.
#'
#' Given a PAC object the function will extract sequences within a given size
#' interval and percent coverage across independent samples.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a Counts table with raw counts or reads per million (rpm).
#'   
#' @param size Integer vector giving the size interval, as c(min,max), that
#'   should be saved (default=c(16,45)).
#'   
#' @param treshold Integer giving the threshold in rpm or counts that needs to
#'   be reached for a sequence to be included (default=10).
#'   
#' @param coverage Integer giving the percent of independent samples that need
#'   to reach the threshold for a sequence to be included (default=100).
#'   
#' @param type Character specifying if filtering should be done "rpm" or
#'   "counts" (default="rpm").
#' 
#' @param stat (optional) Logical specifying if an coverage graph should be
#'   generated or not (default=FALSE).
#'
#' @param subset_only Logical whether only subsetting using pheno_target and/or
#'   anno_target should be done. If subset=FALSE (default) both subsetting and
#'   other filtering will be done.
#'   
#' @param pheno_target (optional) List with: 
#'          1st object being a character vector of target column in Pheno, 
#'          2nd object being a character vector of the target group(s) in the target Pheno column (1st object).
#'          (default=NULL)
#'          
#' @param anno_target (optional) List with: 
#'          1st object being a character vector of target column in Anno, 
#'          2nd object being a character vector of the target type/biotypes(s) in the target Anno column (1st object).
#'          (default=NULL) 
#'
#' @return A list of objects: 
#'               PAC object with filtered data.   
#'               (optional) A covarage plot 
#' @examples
#' 
#' path="/data/Data_analysis/Projects/Drosophila/Other/IOR/Jan_IOR_200130/R_analysis_full/"
#' load(file=paste0(path, "PAC_all.Rdata"))
#' PAC_all <- PAC_rpm(PAC_all)
#' pheno_target = list("Method", c("IOR1_proto", "IOR1_tgirt"))
#'
#' 
#' test <- PAC_filter(PAC_all, size=c(16,45), threshold=10, coverage=50, type="rpm", stat=TRUE, pheno_target=pheno_target, anno_target=NULL)
#' test <- PAC_filter(PAC_all, size=c(16,45), threshold=10, coverage=50, type="counts", stat=TRUE, pheno_target=pheno_target, anno_target=NULL)
#' 
#' 
#' @export

PAC_filter <- function(PAC, size=c(16,45), threshold=10, coverage=100, type="counts", subset_only=FALSE, stat=FALSE, pheno_target=NULL, anno_target=NULL){
                                            library(ggplot2, quietly=TRUE)
                                            options(scipen=999)
                                            strt <- nrow(PAC$Counts)
                                            ### Subset samples by Pheno                                                
                              							if(!is.null(pheno_target)){
                                              sub_pheno <- as.character(PAC$Pheno[, pheno_target[[1]]]) %in% pheno_target[[2]]
                                              if(any(names(PAC)=="norm")){PAC$norm <- lapply(as.list(PAC$norm), function(x){x[,sub_pheno]})}
                                              PAC$Counts  <- PAC$Counts[,sub_pheno]
                                              PAC$Pheno  <- PAC$Pheno[sub_pheno,]
                                              tab_pheno <- as.data.frame(table(sub_pheno))
                                              cat(paste0("\nPheno filter was specified, will retain: ", tab_pheno[tab_pheno[,1]==TRUE, 2], " of ", length(sub_pheno), " samples\n"))                                            ### Subset data by groups
                              							}  
                                                                                        
                                            ### Subset data by Anno
                                            if(!is.null(anno_target)){
                                              sub_anno <- as.character(PAC$Anno[, anno_target[[1]]]) %in% anno_target[[2]]
                                              if(any(names(PAC)=="norm")){PAC$norm <- lapply(as.list(PAC$norm), function(x){x[sub_anno,]})}
                                              PAC$Counts  <- PAC$Counts[sub_anno,]
                                              PAC$Anno  <- PAC$Anno[sub_anno,]
                                              tab_anno <- as.data.frame(table(sub_anno))
                                              cat(paste0("\nAnno filter was specified, will retain: ", tab_anno[tab_anno[,1]==TRUE, 2], " of ", length(sub_anno), " seqs\n"))                                                ### Subset data by groups
                                            }  
                                            
                                            if(subset_only==TRUE){return(PAC)
                                            }else{
                                            
                                              ### Subset data by Size
                                            if(!is.null(size)){
                                              sub_size <- PAC$Anno$Length >= size[1] & PAC$Anno$Length <= size[2] 
                                              if(any(names(PAC)=="norm")){PAC$norm <- lapply(as.list(PAC$norm), function(x){x[sub_size,]})}
                                              PAC$Counts  <- PAC$Counts[sub_size,] 
                                              PAC$Anno  <- PAC$Anno[sub_size,]
                                              tab_anno <- as.data.frame(table(sub_size))
                                              cat(paste0("\nSize filter was specified, will retain: ", tab_anno[tab_anno[,1]==TRUE, 2], " of ", length(sub_size), " seqs\n"))                                                ### Subset data by groups
                                            }
                                              
                                              ### Extract essential information 
                                            if(type=="rpm"){ df <- PAC$norm$rpm; cat("\nRPM filter was specified\n")}
                                            if(type=="counts"){ df <- PAC$Counts; cat("\nCount filter was specified\n")}    
                                            
                                            ### Check col and row names
                                            if(!identical(rownames(PAC$Anno), rownames(df))){stop("\nError: Not matching rownames in input files! (Anno vs data)")}
                                            if(!identical(rownames(PAC$Pheno), colnames(df))){stop("\nError: Not matching rownames in input files!")}
                               							
                                            ### Create rpm filter
                                            idx_filt <- data.frame(rowSums(df >= threshold)) >= round(ncol(df)*(coverage*0.01))
                                            idx_tab <- as.data.frame(table(idx_filt))

                                            ### Generate stat
                                            if(stat==TRUE){
                                                      filt_plot <- data.frame(matrix(NA, ncol=2, nrow=9))
                                        							colnames(filt_plot) <- c("filter","n_features")
                                        							filt_plot$x_graph <- as.numeric(c(1,5,10,20,30,40,50,75,100))
                                        							filt_plot$filter <- c(paste0(">1_in_", coverage, "%"), paste0(">5_in_", coverage, "%"), paste0(">10_in_", coverage, "%"), paste0(">20_in_", coverage, "%"), 
                                        																							paste0(">30_in_", coverage, "%"), paste0(">40_in_", coverage, "%"), paste0(">50_in_", coverage, "%"),
                                        																							paste0(">75_in_", coverage, "%"), paste0(">100_in_", coverage, "%"))
                                        												for (i in 1:length(filt_plot$x_graph)){ 
                                        													  tab <- as.data.frame(table(data.frame(rowSums(df >= filt_plot$x_graph[i])) >= round(ncol(df)*(coverage*0.01))))
                                        														filt_plot[i,2] <- tab$Freq[tab$Var1=="TRUE"]
                                        												}
                              							### Plot graph
                              												p <- ggplot(filt_plot, aes(x=x_graph, y=n_features, fill=x_graph))+
                              																			geom_line()+
                              																			geom_point(cex=2, fill="blue")+
                              																			geom_hline(yintercept=0)+
                              																			scale_x_discrete(limit=filt_plot$x_graph, labels= as.character(filt_plot$filter))+
                              																			ggtitle("User filter:")+					
                              																			xlab(NULL)+
                              																			theme_classic()+
                              																			theme(axis.text.x = element_text(angle = 90, hjust = 0))+
                              																			geom_vline(xintercept=threshold, col="red")+
                              												              geom_label(x=50, y=max(filt_plot[,2])*0.95, label=paste0("n input sequences: ", strt, "\nn analyzed sequences: ", nrow(df), "\nn sequences after filter: ",idx_tab[idx_tab[,1]==TRUE, 2]), show.legend = FALSE)+
                              																			theme(plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit="cm"), plot.title = element_text(color="red", size=10))
                              												Sys.sleep(0.01)
                              												print(p)
                              							### Promt for user input					
                              												cat("\n!!           !!\nUser input needed:\n") 
                              												answer <- readline(prompt="Continute with this filter? [Y/n]")
                              												#if(answer=="N" | answer=="n"){stop("Script was terminated by user.")}
                              												if(answer=="Y" | answer=="y"){}else{stop("Script was terminated by user.")}
                              												}

                                            ### Apply rpm filter
                                              if(type=="rpm"){ 
                                      							cat(paste0("\nThe chosen filters will retain: ", idx_tab[idx_tab[,1]==TRUE, 2], " of ", strt, " seqs\n"))                                               
                                                    if(any(names(PAC)=="norm")){PAC$norm <- lapply(as.list(PAC$norm), function(x){x[idx_filt,]})}
                                                    PAC$Counts  <- PAC$Counts[idx_filt,] 
                                                    PAC$Anno  <- PAC$Anno[idx_filt,]
                                              }
                                            }
                                            ## Double check
                                                   if(!identical(rownames(PAC$Anno), rownames(PAC$Counts))){stop("Error: Not matching rownames in input files! (Anno vs Counts)")}
                                                   if(!identical(rownames(PAC$Pheno), colnames(PAC$Counts))){stop("Error: Not matching rownames/colnames in input files! (Pheno vs Counts)")}
                                              if(type=="rpm"){  
                                                   if(!any(do.call("c", lapply(as.list(PAC$norm), function(x){identical(rownames(PAC$Anno), rownames(x))})))){stop("Error: Not matching rownames in input files! (Anno vs RPM)")}
                                                   if(!any(do.call("c", lapply(as.list(PAC$norm), function(x){identical(colnames(PAC$Counts), colnames(x))})))){stop("Error: Not matching rownames/colnames in input files! (Pheno vs RPM)")}
                                              }
                                    if(!subset_only==TRUE){return(PAC)}
                                      }      

