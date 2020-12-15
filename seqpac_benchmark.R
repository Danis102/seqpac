######## From FastQC output #############
fls1 <- as.list(list.files(path="/home/danis31/Desktop/Temp_docs/temp/test/seqpac/", pattern="fastqc.zip", full.names =TRUE))
names(fls1) <- paste0("seqpac_", gsub(".trim_fastqc.zip", "", basename(unlist(fls1))))
fqc1 <- lapply(fls1, function(x){
         unzfls <- unzip(x, list = TRUE)
         fl <- unzfls$Name[grepl("fastqc_data.txt", unzfls$Name)]
         fqc <- readLines(unzip(x, files=fl))
         
         adapt <- fqc[(which(grepl(">>Adapter Content", fqc))+1):(which(grepl(">>Kmer Content", fqc))-2)]
         temp <- strsplit(as.character(adapt),"\t", fixed=TRUE)
         adapt <- data.frame(do.call('rbind',temp[-1]))
         colnames(adapt) <- temp[[1]]
         
         over <- fqc[(which(grepl(">>Overrepresented sequences", fqc))+1):(which(grepl(">>Adapter Content", fqc))-2)]
         temp <- strsplit(as.character(over),"\t", fixed=TRUE)
         over <- data.frame(do.call('rbind',temp[-1]))
         colnames(over) <- temp[[1]]
        
         return(list(adapt=adapt, over=over))
         })

fls2 <- as.list(list.files(path="/home/danis31/Desktop/Temp_docs/temp/test/cutadapt/", pattern="fastqc.zip", full.names =TRUE))
names(fls2) <- paste0("cutadapt_", gsub(".trim_fastqc.zip", "", basename(unlist(fls2))))
fqc2 <- lapply(fls2, function(x){
         unzfls <- unzip(x, list = TRUE)
         fl <- unzfls$Name[grepl("fastqc_data.txt", unzfls$Name)]
         fqc <- readLines(unzip(x, files=fl))
         
         adapt <- fqc[(which(grepl(">>Adapter Content", fqc))+1):(which(grepl(">>Kmer Content", fqc))-2)]
         temp <- strsplit(as.character(adapt),"\t", fixed=TRUE)
         adapt <- data.frame(do.call('rbind',temp[-1]))
         colnames(adapt) <- temp[[1]]
         
         over <- fqc[(which(grepl(">>Overrepresented sequences", fqc))+1):(which(grepl(">>Adapter Content", fqc))-2)]
         temp <- strsplit(as.character(over),"\t", fixed=TRUE)
         over <- data.frame(do.call('rbind',temp[-1]))
         colnames(over) <- temp[[1]]
        
         return(list(adapt=adapt, over=over))
         })

# Overlapping seqs
uni_seqs <- unique(unlist(lapply(list(seqpac=fqc1, cutadapt=fqc2), function(x){
    as.character(unlist(lapply(x, function(y){y$over$'#Sequence'})))})))

ord_fqcs <- lapply(list(seqpac=fqc1, cutadapt=fqc2), function(x){
      lapply(x, function(y){y$over$Count[match(uni_seqs, y$over$'#Sequence')]})})

# Correlation
lst <- list()
for(i in 1:length(ord_fqcs[[1]])){
     nam_seq<- paste(names(ord_fqcs[[1]][i]), sum(!is.na(ord_fqcs[[1]][[i]])), sep="_")
     nam_cut<- paste(names(ord_fqcs[[2]][i]), sum(!is.na(ord_fqcs[[1]][[i]])), sep="_")
     
     dt <- data.frame(log10(as.numeric(as.character(ord_fqcs[[1]][[i]]))), log10(as.numeric(as.character(ord_fqcs[[2]][[i]]))))
     colnames(dt) <- c(nam_seq, nam_cut)
     corl <- Hmisc::rcorr(as.matrix(dt))
     lst[[i]] <- ggplot2::ggplot(dt, ggplot2::aes_string(x=names(dt[1]) ,y=names(dt[2])))+
                           ggplot2::geom_point() +
                           ggplot2::geom_smooth(method = "lm", col="red") +
                           ggplot2::labs(subtitle=paste("r2=", round(corl$r[1,2]^2, digits =7), ", p=", corl$P[1,2], ", n=", corl$n[1,2], collapse=""))
     names(lst)[i] <- names(ord_fqcs[[1]][i])
     }

cowplot::plot_grid(plotlist=lst, nrow=3, ncol=3)

######## Directly from fastq #############

fls1 <- as.list(list.files(path="/home/danis31/Desktop/Temp_docs/temp/test/seqpac/", pattern="fastq.gz", full.names =TRUE))
fls2 <- as.list(list.files(path="/home/danis31/Desktop/Temp_docs/temp/test/cutadapt/", pattern="fastq.gz", full.names =TRUE))
names(fls1) <- paste0("seqpac_", gsub(".trim.fastq.gz", "", basename(unlist(fls1))))
names(fls2) <- paste0("cutadapt_", gsub(".trim.fastq.gz", "", basename(unlist(fls2))))
 

doParallel::registerDoParallel(3)  # Do not use parallel::makeClusters!!!
require(foreach)
tab_plt <- foreach(i=1:length(fls1), .inorder = TRUE, .export= c("fls1", "fls2"), .final = function(x){names(x) <- gsub("seqpac_", "", names(fls1)); return(x)}) %dopar% {
      df <- data.frame(matrix(NA, nrow=1, ncol=11))
      colnames(df) <- c("ID", "seq_tot_reads", "cut_tot_reads", "seq_uni_seqs", "cut_uni_seqs", "seq_precSeqsNotIn_cut", "cut_precSeqsNotIn_seq", "seq_nReadsNotIn_cut", "cut_nReadsNotIn_seq", "R2_olap_counts", "p_val_olap_counts")                           
      df$ID <- gsub("seqpac_", "", names(fls1)[i])
      fstq_seq <- paste0(ShortRead::sread(ShortRead::readFastq(fls1[[i]], withIds=FALSE)))
      fstq_seq <- dplyr::count(tibble::tibble(seqpac=fstq_seq), seqpac)
      fstq_cut <- paste0(ShortRead::sread(ShortRead::readFastq(fls2[[i]], withIds=FALSE)))
      fstq_cut <- dplyr::count(tibble::tibble(cutadapt=fstq_cut), cutadapt)
     
      df$seq_tot_reads <- sum(fstq_seq$n)
      df$cut_tot_reads <- sum(fstq_cut$n)
       
      df$seq_uni_seqs <- length(fstq_seq$seqpac)
      df$cut_uni_seqs <- length(fstq_cut$cutadapt)
      
      df$seq_precSeqsNotIn_cut <- sum(!fstq_seq$seqpac %in% fstq_cut$cutadapt) / df$seq_uni_seqs *100 # 0.024%
      df$cut_precSeqsNotIn_seq <- sum(!fstq_cut$cutadapt %in% fstq_seq$seqpac) / df$cut_uni_seqs *100 # 0.043%

      df$seq_nReadsNotIn_cut <- sum(fstq_seq$n[!fstq_seq$seqpac %in% fstq_cut$cutadapt]) / sum(fstq_seq$n) *100 # 0.024%
      df$cut_nReadsNotIn_seq <- sum(fstq_cut$n[!fstq_cut$cutadapt %in% fstq_seq$seqpac]) / sum(fstq_cut$n) *100 # 0.043%      
      
      ord_seqs <- fstq_seq$seqpac[fstq_seq$seqpac %in% fstq_cut$cutadapt]
      ord_cut <- fstq_cut[match(ord_seqs, fstq_cut$cutadapt),]
      ord_seq <- fstq_seq[match(ord_seqs, fstq_seq$seqpac),]
      rm(fstq_seq, fstq_cut)
      ord_cut <- ord_cut[!is.na(ord_cut$cutadapt),]
      ord_seq <- ord_seq[!is.na(ord_seq$seqpac),]
      dt <- data.frame(log10(ord_seq$n), log10(ord_cut$n))
      rm(ord_cut, ord_seq)
      colnames(dt) <- c(paste0(df$ID, "_seqpac"),  paste0(df$ID, "_cutadapt"))
      corl <- Hmisc::rcorr(as.matrix(dt))
      df$R2_olap_counts <- corl$r[1,2]^2
      df$p_val_olap_counts <- corl$P[1,2]
      plot <- ggplot2::ggplot(dt, ggplot2::aes_string(x=names(dt[1]) ,y=names(dt[2])))+
                           #ggplot2::geom_smooth(method = "lm", formula= 'y ~ x', size=0.5, col="red") +
                           ggplot2::geom_point(size=2, col="darkblue") +
                           ggplot2::labs(subtitle=paste("r2=", round(corl$r[1,2]^2, digits =7), ", p=", corl$P[1,2], ", n=", corl$n[1,2], collapse="")) +
                           ggplot2::theme_classic()
      rm(dt)
      return(list(df=df, plot=plot))
}
      
df <- do.call("rbind", lapply(tab_plt, function(x){x$df}))     
write.table(df, file="/home/danis31/OneDrive/Skriverier/Daniel - seqpac/Results_benchmark/BenchmarkTable.xls", sep="\t", row.names=FALSE)

plts <- lapply(tab_plt, function(x){x$plot})
plts[[1]]
cowplot::plot_grid(plotlist=plts, nrow=3, ncol=3)





