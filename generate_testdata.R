
input1 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index11_S3_merge.fastq.gz" #sperm
nam1 <- "sperm_wt_test1"
input2<-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index14_S2_merge.fastq.gz" #sperm
nam2 <- "sperm_wt_test2"
input3 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index17_S1_merge.fastq.gz" #sperm
nam3 <- "sperm_wt_test3"

input4 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index13_S16_merge.fastq.gz" #sperm Aub
nam4 <- "sperm_Aub_test1"
input5 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index16_S13_merge.fastq.gz" #sperm Aub
nam5 <- "sperm_Aub_test2"
input6 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index19_S24_merge.fastq.gz" #sperm Aub
nam6 <- "sperm_Aub_test3"

input7 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index20_S17_merge.fastq.gz" #ovary
nam7 <- "ovary_wt_test1"
input8 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index20_S20_merge.fastq.gz" #ovary
nam8 <- "ovary_wt_test2"
input9 <-  "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/Index10_S1_Ovaries_test3.fastq.gz" #ovary
nam9 <- "ovary_wt_test3"
nams <- c(nam1, nam2, nam3, nam4, nam5, nam6, nam7, nam8, nam9)


input1 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx3-190327_S3_merge.fastq.gz"
input2 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx4-190327_S4_merge.fastq.gz"
input3 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx5-190327_S5_merge.fastq.gz"


input4 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx13-190327_S13_merge.fastq.gz"
input5 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx14-190327_S14_merge.fastq.gz"
input6 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx15-190327_S15_merge.fastq.gz"

input7 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx20-190327_S20_merge.fastq.gz"
input8 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx22-190327_S22_merge.fastq.gz"
input9 <- "/data/Data_analysis/Projects/Drosophila/Other/Lovisa/Stage1to5_wm4h_1/Merged_fastq/Inx23-190327_S23_merge.fastq.gz"

input <- list(input1, input2, input3, input4, input5, input6, input7, input8, input9)
nams <- c("Stage1_Batch1", "Stage1_Batch2", "Stage1_Batch3", 
          "Stage3_Batch1", "Stage3_Batch2", "Stage3_Batch3", 
          "Stage5_Batch1", "Stage5_Batch2", "Stage5_Batch3")
output_path <- "/data/Data_analysis/Projects/Drosophila/Specific_projects/temp/"


lst <- list(NA)
for(i in 1:length(input)){
    fstq <- ShortRead::readFastq(input[[i]], withIds=FALSE)
    RandomNum <- round(runif(100000, 1, length(fstq)), digits=0)
    fastq_red <- fstq[RandomNum]
    output <- gsub("_merge.fastq.gz", "", basename(input[[i]])) 
    output <- gsub("-190327", "", output)
    output <- paste0(output_path, nams[i], "_", output, ".fastq.gz")
    ShortRead::writeFastq(fastq_red, output, mode="w", full=FALSE, compress=TRUE)
}





