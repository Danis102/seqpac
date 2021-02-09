#' Make Phenotype File
#'
#' \code{make_pheno} generates a Pheno data.frame object containing additional
#' metadata for unique samples.
#'
#' Given the path to a directory with a single csv file, the function will read
#' the file. If provided with a countTable (see \code{make_countTable}), it will
#' attempt to organize the the row names according to the column names in the
#' countTable.
#'
#' The function was specifically written for Illumina SampleSheetUsed.csv files
#' downloaded from the BaseSpace server using the basemount software. It will,
#' however, be compatible with any comma seperated text file with any of its
#' first lines containing a header with the first column named 'Sample_ID'.
#'
#' The function will also attempt to read a 'Samples_stat.txt' (tab seperate
#' file) file also downloaded from Illumina BaseSpace. This file will add
#' original number of reads (past filter) generated in the sequencing. If
#' provided with a Sports progress_report object (see \code{progress_report} the
#' script will attempt to add this data to the Pheno data.frame.
#'
#' @family PAC generation
#'
#' @seealso \url{https://github.com/junchaoshi/sports1.0} for download and
#'   documentation about running Sports. \url{https://github.com/Danis102} for
#'   updates on the current package.
#'
#' @param pheno_input Character vector with the path to a .csv file or a
#'   dataframe with rownames as sample names.
#'
#' @param type Character indicating what type to file to expect. If
#'   \emph{type="basespace"} the function will attempt to read nd join
#'   information from 'SampleSheetUsed.csv' and 'Samples_stat.txt' in
#'   pheno_input that have been downloaded from Illumina Basespace. If
#'   type="manual", the function will attempt to read a file named 'pheno.csv',
#'   in which the first column has been named 'Sample_ID' containing the exact
#'   sample names matching basenames of the fastq sample files.
#'
#' @param counts Data.frame object with the same column names as in
#'   Sample_ID column of the .csv file.
#'
#' @param progress_report Data.frame object with progress report.
#'
#'
#' @return data.frame
#'
#' @examples
#'
#' ### Sports input
#'   path_sport <- "/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/Processed_Pipeline31_05-03-20"
#'
#' sports_lst <- import_sports("/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/Processed_Pipeline31_05-03-20", threads=8)
#'
#' anno <- make_anno(sports_lst = sports_lst, threads = 10, filt = 2, stat = TRUE)
#' counts <- make_counts(sports_lst = sports_lst, anno = anno, threads = 10)
#' report <- progress_report(path_sport)
#'
#' pheno <- make_pheno(<your_path_to_sports_output_directory>)
#'
#' pheno_input <- "/data/Data_analysis/Projects/Pigs/Specific_projects/SRA_download/SRP135969_Sperm_Exosomes_Hemicastration/summary"
#' pheno <- make_pheno(pheno_path, countTable=counts, progress_report=report, type="manual")
#'
#' ### Simple manual
#'
#' @export
make_pheno<- function(pheno_input, type="manual", counts=NULL, progress_report=NULL){
  
  ### Read using basespace download files
  if(type=="basespace"){
    if(grepl("SampleSheetUsed.csv|SampleSheet.csv", pheno_input)){
      lines <- readLines(pheno_input, n=20)
      header <- which(grepl("\\<Sample_ID", lines))
      if(!length(header) == 1){
        stop("Error! Cannot find comma seperated header with first column named 'Sample_ID'")
      }
      pheno <- read.delim(pheno_input, skip= header-1,  header=TRUE, sep=",")
      cat("Illumina type SampleSheet.csv file was found.\n")
      ## Universal search for colnames:
      col_srch <- c("^Sample_ID$", "^Sample_Name$", "^SampleProject$|^Sample_Project$", "^Index$|^index$")
      col_nr <- unlist(lapply(col_srch, function(x){which(grepl(x, colnames(pheno)))})) 
      pheno <- pheno[, col_nr]
      ## Try to add Illumina stat
      # path_stat <- list.files(dirnames(pheno_input, pattern=".txt", full.names = TRUE)
      # if(length(path_stat) == 1){
      #               stat <- read.delim(paste0(pheno_input, "/Samples_stat.txt"), header=TRUE, sep="\t")
      #               stopifnot(any(as.character(pheno$Sample_ID) %in% as.character(stat$SampleId)))
      #               pheno <- cbind(pheno, stat[match(as.character(pheno$Sample_ID), as.character(stat$SampleId)), !colnames(stat) %in% c("SampleId", "Name", "Index")])
      # }
    }else{
      type <- "manual"
      warning("Did not find SampleSheet.csv. Will try to read pheno-file using type='manual' instead.")
    } 
    
  }
  
  ### Read using manual pheno.txt file
  if(type=="manual"){
    if(is.data.frame(pheno_input)){
      pheno <- pheno_input
    }else{
      lines <- readLines(pheno_input, n=20)
      header <- which(grepl("^Sample_ID", lines))
      if(!length(header) == 1){
        stop("Error! Cannot find comma seperated header with first column named 'Sample_ID'")}
      pheno <- read.delim(pheno_input,  header=TRUE, sep=",")
    }
  }
  
  ## Order as counts
  if(!is.null(counts)){
    pheno_ord <- unlist(lapply(as.list(as.character(pheno$Sample_ID)), function(x){
      logi_colnam <- grepl(x, colnames(counts))
      if(sum(logi_colnam)>1){
        stop("Input pheno does not have unique Sample_IDs")
      }
      return(which(logi_colnam))
    }))
    if(any(duplicated(pheno_ord))){
      stop("Sample_IDs were not unique in pheno input.")
    }
    rownames(pheno) <- colnames(counts)[pheno_ord]
    pheno <- pheno[match(colnames(counts), rownames(pheno)),]
    stopifnot(identical(rownames(pheno), colnames(counts)))
    cat("Of", length(colnames(counts)), "sample names in counts,", sum(as.numeric(rownames(pheno) %in% colnames(counts))), "were found in pheno file path.\n")
    cat("\n")
    print(data.frame(pheno=as.character(pheno$Sample_ID), counts=colnames(counts)))
  }else{
    warning("\nNo count table was specified. Final Pheno will be unordered!\n")
  }
  
  ## Add progress report
  if(!is.null(progress_report)){
    cat("\nA progress report was specified, will attempt to match rownames...\n")
    
    # Fix sample names (illumina automatically exchanges "-" for "_"
    if(!sum(rownames(progress_report) %in% rownames(pheno)) == nrow(pheno)){
      rownames(progress_report) <- gsub("-", "_", rownames(progress_report))
    }
    prog_ord <- unlist(lapply(as.list(rownames(progress_report)), function(x){
      logi_colnam <- grepl(x, rownames(pheno))
      if(sum(logi_colnam)>1){
        stop("Input pheno does not have unique Sample_IDs")
      }
      return(which(logi_colnam))
    }))
    if(any(duplicated(prog_ord))){
      stop("Sample_IDs were not unique in pheno input.")
    }
    progress_report <- progress_report[match(rownames(pheno), rownames(progress_report)),]
    stopifnot(identical(rownames(pheno), rownames(progress_report)))
    pheno <- cbind(pheno, progress_report)
    cat("\n")
  }else{
    warning("\nNo progress report was specified. Will be missing in the final Pheno.\n")
  }
  cat("Done!\n")
  return(pheno)
}




