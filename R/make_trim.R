# AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT
# AGATCGGAAGAGCACGCGTCTGAACTCCA
# 19652834   13382796    12785870
# 17169481/21695341       
# make_trim
# path_to_fastq -> path
# 
# 
# fls = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Drosophila_Sep_IOR_190912/Data/Double/Long/Merged_fastq/Inx14-190912_S17_merge.fastq.gz"
# fls = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Jan_IOR_200130/Data/Single/Merged_fastq/Inx18-200130_S15_merge.fastq.gz" # medium TGIRT
# fls = "/data/Data_analysis/Projects/Drosophila/Other/IOR/Drosophila_Sep_IOR_190912/Data/Single/Merged_fastq/Inx3-190912_S13_merge.fastq.gz"  # Large POOH
# 
# input =  path_to_fastq <- system.file("extdata", package = "seqpac", mustWork = TRUE)
# path=
# adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT"
# adapt_min=5
# adapt_mis=20
# mis_inc=10
# max_inc=3
# qc_thesh=20
# qc_prop=80
# 
#  if(type=="fastq"){
#                                     cat("Started at ", paste0(Sys.time()), "\n")
#                                     ## Parallel setup
#                                     gc(reset=TRUE)
#                                     cl <- parallel::makeCluster(threads, type = par_type)
#                                     doParallel::registerDoParallel(cl)
#                             
#                                     ## Read file system
#                                     path <- input
#                                     count_files <- list.files(path, pattern ="fastq.gz\\>|fastq\\>", full.names=TRUE, recursive=TRUE)
#                                     count_files <- count_files[!grepl("Undetermined_", count_files)]
#                                     count_files_nams <- basename(count_files)
#                                     cat("\nInput type was set to fastq.\n")
#                                     cat("The following fastq files were found in the path:\n")
#                                     print(count_files)
#   
#                                     
# AGATCGTTAGATCGGAAGCGCACAAGTCTGAACTCCAGACACTTAGGCATCGCGTATGC                                  
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT                                    
#                                     
# make_trim <- function(path, adapt="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT", adapt_min=5, adapt_mis=20, mis_inc=10, max_inc=3, qc_thesh=20, qc_prop=80)  
#                                       
#                     fls <- list.files(path, pattern ="fastq.gz\\>|fastq\\>", full.names=TRUE, recursive=TRUE)
#  
#                     cl <- parallel::makeCluster(threads, type = par_type)
#                     doParallel::registerDoParallel(cl)
#                    
#                     
#                      
# 
#                     
#                     srch_1st <- substr(adapt, 1, adapt_min)
#                     srch_2nd <- substr(adapt, adapt_min+1, adapt_min+adapt_min)
#                     
#                     srch_last <- substr(adapt, 1, adapt_min*2)
#                     srch_mis <- substr(adapt, 1, adapt_mis)
#                     
#                     seq(mis_inc, mis_inc*3, mis_inc)
#                     seq_inc <-  seq(mis_inc, nchar(adapt), mis_inc)  
#                                    
#                     
#                     
#                     
# 
#                     foreach::foreach(i=1:length(fls), .packages=c("ShortRead"), .final = function(x){names(x) <- basename(fls); return(x)}) %dopar% {
#                                             
#                                             ## Read fastq and setup files
#                                             fstq <- ShortRead::readFastq(fls[[i]])
#                                             cnts <- table(paste0(ShortRead::sread(fstq)))
#                                             vect <- names(cnts)
#                                             names(vect) <- vect
#                                             lgn <- unique(nchar(vect))
#                                             vect[] <- "empty"
#                                             check_again_lst <- list(NULL) 
#                                             if(!length(lgn) == 1){stop("\nDiffering read lengths prior to adapter trimming.\nHave you already performed 3-prim trimming?")}
# 
#                                             ## Find obvious hits with first followed by another
#                                             adp_coord_1st <- stringr::str_locate(names(cnts), srch_1st)[,1]
#                                             adp_left <- substr(names(cnts), adp_coord_1st, lgn)
#                                             dub_hit <- stringr::str_locate(adp_left, srch_2nd)[,1]
#                                             filt <- dub_hit == (adapt_min+1)
#                                             filt[is.na(filt)] <- FALSE
#                                             vect[filt] <- adp_coord_1st[filt]
#                                             
#                                             
#                                             ## Find new coord if no double hit
#                                             filt_no <- !dub_hit == (adapt_min+1)
# 
#                                             adp_left_new <- substr(adp_left, 2, nchar(adp_left))
#                                             adp_coord_1st_new <- stringr::str_locate(adp_left_new, srch_1st)[,1]
#   
#                                             adp_left_new_new <- substr(adp_left_new, adp_coord_1st_new, nchar(adp_left_new))
#                                             dub_hit_new <- stringr::str_locate(adp_left_new_new, srch_2nd)[,1]
#                                             filt_new <- dub_hit_new == (adapt_min+1)
#                                             filt_new[is.na(filt_new)] <- FALSE
# 
#                                             if(!unique(vect[paste0(filt,  filt_new) == "FALSETRUE"]) == "empty"){stop("Unresolved double hit(s). This should simply not happen!")}
#                                             vect[paste0(filt,  filt_new) == "FALSETRUE"] <- (lgn - nchar(adp_left_new_new)+1)[paste0(filt,  filt_new) == "FALSETRUE"]
#                                             check_again_lst[[1]] <- paste0(filt,  filt_new) == "TRUETRUE"  # Save for later
#                                             
#                                             rm(adp_coord_1st, adp_left , dub_hit, filt, filt_no, adp_left_new, adp_coord_1st_new, adp_left_new_new, dub_hit_new, filt_new)
#                                             
#                                             ## 3prim sequential end-trimming
#                                             srch_loop <- srch_last
#                                             nt <- nchar(srch_last)
#                                             start <- lgn+1-nt
# 
#                                             filt_last_lst <- as.list(rep.int("na", times=adapt_min+1))
#                                             names(filt_last_lst) <- paste0("n", (lgn+1-nt):(lgn+1-adapt_min))
#                                             
#                                             while(start <= lgn+1-adapt_min ){
#                                                       cat(paste0("\n", start))
#                                                       last_nt <- substr(names(cnts), start, lgn)
#                                                       filt_last_lst[[paste0("n",start)]] <- grepl(srch_loop, last_nt)
#                                                       srch_loop <- substr(srch_loop, 1, nchar(srch_loop)-1)
#                                                       start <- start+1
#                                             }
#                                             for(t in 1:length(filt_last_lst)){
#                                                       prev_log <- paste0(!vect=="empty", filt_last_lst[[t]])
#                                                       check_again_lst[[length(check_again_lst)+1]] <- prev_log == "TRUETRUE"
#                                                       vect[prev_log == "FALSETRUE"] <- as.integer(gsub("n", "", names(filt_last_lst)[t]))}
#            
#                                             
#                                             ## Mismatch search
#                                             mis_lst <- list(NULL)
#                                             for(t in 1:length(seq_inc)){
#                                                               inc <- substr(adapt, 1, seq_inc[t])
#                                                               hit_mis <- agrepl(inc, names(cnts), max.distance = t-1)
#                                                               mis_lst[[t]] <- paste0(!vect=="empty", hit_mis)
#                                                               }
#                                             lapply(mis_lst, function(x){return(table(x=="FALSETRUE"))})
#                                             
#                                             test <- do.call("cbind", check_again_lst)
#                                             test <- rowSums(test) >0
#                                             
#                                             table(paste0(mis_lst[[1]], check_again_lst[[1]]) == "FALSETRUETRUE")
#                                             
#                                             names(vect)[ 
#                                             vect[mis_lst[[1]]=="FALSETRUE"]
#                                             vect[check_again_lst[[1]]=="TRUE"]
#                                             
#                                             <- seq(round(lgn*mismatch_rate, digits=0), lgn,  by=round(lgn*mismatch_rate, digits=0))
#                                             
#                                             
#                                             
#                                             
#                                             
#                                             
#                                             
#                                             
#                                             lapply(filt_last_lst, table)
#                                             
#                                             
#                                             
#                                             adp_left[filt_last_lst[[6]]]
#                                             
#                                             
#                                             
#                                             
#                                             
#                                             adp_coord_2nd <- stringr::str_locate(adp_left, srch_2nd)[,1]
#                                             
#                                             
#                                             adp_left_re <- substr(adp_left, 2, nchar(adp_left))
#                                             
#                                             adp_coord_new <- stringr::str_locate(adp_left_re, srch_1st)
#                                             adp_left_new <- substr(adp_left_re, adp_coord_new, nchar(adp_left_re))
#                                             
#                                             
#                                             dub_hit_new <- stringr::str_locate(adp_left_new, srch_2nd)
#                                             dub_hit_new[is.na(dub_hit_new[,1]),1] <- lgn
#                                             sub_filt <- dub_hit_new[,1] == (adapt_min+1)
#                                             vect[sub_filt] <- adp_coord_new[sub_filt,1]
#                                             
#                                             adp_coord_new[dub_hit[,1] == (adapt_min+1)]
#                                            
#                                             vect <-
#                                             
#                                             
#                                             ## Find long this with mismatch
#                                             
#                                             
#                                             
#                                             names(cnts)
#                                             
#                                             
#                                             test <- stringr::str_locate(adapt_mis, srch_1_5)
#                                             
#                                             
#                                             table(dub_hit[,1])
#                                             
#                                             filt_best_30  <- agrepl(srch_1_30, names(cnts), max.distance = 2)
#                                             
#                                             table((adp_coord[,1]==75) - (is.na(dub_hit[,1])) == 0)
#                                             table((adp_coord[,1]==75) + (is.na(dub_hit[,1])) == 1)
#                                             table(!adp_coord[,1]==75)
#                                             table(!is.na(dub_hit[,1]))
#                                             
#                                             srch_6_10
#                                             str_split
#                                             
#                                             
#                                             
#                                             
#                                             Här är jag!!!!
#                                               Gör så här:
#                                               1. Sök upp alla srch_1_5 med nedanstående function (str_locate)
#                                               2. Gör en split på nucleotiden framför (borde då få adaptor bara)
#                                               3. Leta använd efter adaptor med agrepl och 2 mismatch (borde få hits på alla som har adaptor)
#                                               4. Alla som har bekräftad adaptor kan sedan klippas där str_locate hittade.
# 
#                                             test5 <- stringr::str_locate(names(cnts), srch_1_5)
#                                             
#                                             
#                                             table(is.na(test30[,1]))
#                                             table(is.na(test20[,1]))
#                                             table(is.na(test10[,1]))
#                                             table(is.na(test5[,1]))
#                                             
#                                             filt_best_30  <- agrepl(srch_1_30, names(cnts), max.distance = 2)
#                                             filt_best_20  <- agrepl(srch_1_20, names(cnts), max.distance = 1)
#                                             
#                                             rest <- names(cnts)[filt_best_20+filt_best_30==0]
#                                             lgn <- unique(nchar(rest))
#                                             
# 
#                                             ## Narrowing search up to the last 5 nt at 3prim end
#                                             nt <- start
#                                             filt_last_lst <- as.list(rep.int("na", times=75-start+1))
#                                             names(filt_last_lst) <- start:lgn
#                                             while(nt < lgn ){
#                                                       cat(paste0("\n", nt))
#                                                       last_nt <- substr(rest, lgn-(nt-1), lgn)
#                                                       srch_last <- substr(adpt, 1, nt)
#                                                       if(nt<mis_inc){filt_last_lst[[paste(nt)]] <- grepl(srch_last, last_nt)}
#                                                       if(nt>=mis_inc){filt_last_lst[[paste(nt)]] <-  agrepl(srch_last, last_nt, max.distance = 1)}
#                                                       if(nt>=mis_inc*2){filt_last_lst[[paste(nt)]] <-  agrepl(srch_last, last_nt, max.distance = 2)}   
#                                                       nt <- nt+1
#                                                       
#                                             }
#                                             
#                                             sqs <- names(cnts)
#                                             nt <- start
#                                             filt_last_lst <- as.list(rep.int("na", times=75-start+1))
#                                             names(filt_last_lst) <- start:lgn
#                                             while(nt <= lgn ){
#                                                       cat(paste0("\n", nt))
#                                                       last_nt <- substr(sqs, lgn-(nt-1), lgn)
#                                                       srch_last <- substr(adpt, 1, nt)
#                                                       filt_last_lst[[paste(nt)]] <- grepl(srch_last, last_nt)
#                                                       #if(nt>=mis_inc){filt_last_lst[[paste(nt)]] <-  agrepl(srch_last, last_nt, max.distance = 1)}
#                                                       #if(nt>=mis_inc*2){filt_last_lst[[paste(nt)]] <-  agrepl(srch_last, last_nt, max.distance = 2)}   
#                                                       nt <- nt+1
#                                                       
#                                             }
#                                             
#                                             narr <- do.call("cbind", filt_last_lst)
#                                             table(rowSums(narr) >0)
#                                             
#                                             rowSums(narr)
#                                             rest_filt <- (filt_last20 + filt_last14 + filt_narr) > 0
#                                             
#                                             left <- rest[!rest_filt]
#                                             
#                                             
#                                             table(grepl("AGAT", left))
#                                             table(grepl("ACAT", left))
#                                             
#                                             AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT
#                                             
#                                             test  <- left[!agrepl(srch_1_15, left, max.distance = 1)]
#                                             
#                                             test2 <- cnts[names(cnts) %in% test]
#                                             
#                                             
#                                             ## Map with 1 mismatch
#                                             filt_best_1mis  <- agrepl(srch_1_20, mismtch, max.distance = 1)
#                                             test <- agrepl(srch_1_20, names(cnts), max.distance = 1)
#                                             test <- mismtch[filt_best_1mis]
#                                                           
#                                                               
#                                                               
#                             ####### Phred score quality filtering   
#                             ## Extract phred score fastq translations
#                             enc <- encoding(quality(fstq))
#                             enc_srch <- as.factor(names(enc))[paste0(enc) %in% as.character(1:qc_thesh-1)]
#                             enc_srch <- paste0("[", paste(enc_srch, collapse="") , "]")
#                             
#                             ## Count phred scores
#                             qual <- paste(quality(Biostrings::quality(fstq)))
#                             
#                             test <- stringr::str_count(paste(quality(Biostrings::quality(fstq))), enc_srch)
#                             
#                             test <- stringr::str_count(qual, "[#]")
#                             table(test)
#                             ShortRead::sread(fstq)[test==6]
#                             Biostrings::quality(fstq)[test==6]
#                                                            
#                                                               
#                                                               
#                                                               
#                                                               
                                                              
                                                           