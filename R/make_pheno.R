#' Make Phenotype File
#'
#' \code{make_pheno} generates a Pheno data.frame object containing additional
#' metadata for unique samples.
#'
#' Given the path to a directory with a single csv file, the function will read
#' the file. Alternatively, a data.frame with sample information can be
#' provided. If provided with a count table, (see \code{\link{make_counts}}), it
#' will attempt to organize the the row names according to the column names in
#' the counts table.
#'
#' The function was originally written for Illumina SampleSheetUsed.csv files
#' downloaded from the BaseSpace server using the basemount software. It will,
#' however, be compatible with any comma seperated text file with any of its
#' first lines containing a header with the first column named 'Sample_ID'.
#'
#' The function will also attempt to read a 'Samples_stat.txt' (tab seperate
#' file) file also downloaded from Illumina BaseSpace. This file will add
#' original number of reads (past filter) generated in the sequencing. If
#' provided with a progress_report object (see \code{\link{make_counts}} the
#' script will attempt to add this data to the Pheno data.frame.
#'
#' @family PAC generation
#'
#' @seealso \url{https://github.com/junchaoshi/sports1.0} for download and
#'   documentation about running Sports. \url{https://github.com/Danis102} for
#'   updates on the current package.
#'
#' @param pheno Character vector with the path to a .csv file or a
#'   data.frame with a column named "Sample_ID".
#'
#' @param type Character indicating what type to file to expect. If
#'   type="manual", the function will attempt to read a file named 'pheno.csv',
#'   in which the a column has been named 'Sample_ID' containing the exact
#'   sample names matching basenames of the fastq sample files. If
#'   \emph{type="basespace"} the function will attempt to read and join
#'   information from 'SampleSheetUsed.csv' and 'Samples_stat.txt' in
#'   pheno that have been downloaded from Illumina Basespace. 
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
#' library(seqpac) 
#' 
#' ### First make counts 
#' 
#' input = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' counts  <- make_counts(input, threads=6, parse="default_neb", type="fastq",
#'                        trimming="seqpac", plot=TRUE, evidence=c(experiment=2, sample=1))
#'
#'
#' ### Then generate a phenotype table with make_pheno (herre using file names)
#'
#' #  Note:  'Sample_ID' column need to be the same IDs as colnames in the counts table.
#'
#' pheno <- as.data.frame(do.call("rbind", strsplit(list.files(input,pattern="*.fastq.gz"), "_|\\."))[,c(1,2,3,4)]) 
#' colnames(pheno) <- c("stage", "batch", "index", "sample") 
#' pheno$Sample_ID <- apply(pheno, 1, function(x){paste(x, collapse="_")}) 
#' pheno <- make_pheno(pheno=pheno, progress_report=counts$progress_report, counts=counts$counts)
#'
#'  
#' ### Lastly combine into PAC
#' 
#' # Note: a simple annotation table will be added automatically.
#' 
#' pac <- make_PAC(pheno=pheno, counts=counts, anno=NULL)
#'
#' @export
make_pheno<- function(pheno, type="manual", counts=NULL, progress_report=NULL){
  
  ### Read using basespace download files
  if(type=="basespace"){
    if(grepl("SampleSheetUsed.csv|SampleSheet.csv", pheno)){
      lines <- readLines(pheno, n=20)
      header <- which(grepl("\\<Sample_ID", lines))
      if(!length(header) == 1){
        stop("Error! Cannot find comma seperated header with first column named 'Sample_ID'")
      }
      pheno <- read.delim(pheno, skip= header-1,  header=TRUE, sep=",")
      cat("Illumina type SampleSheet.csv file was found.\n")
      ## Universal search for colnames:
      col_srch <- c("^Sample_ID$", "^Sample_Name$", "^SampleProject$|^Sample_Project$", "^Index$|^index$")
      col_nr <- unlist(lapply(col_srch, function(x){which(grepl(x, colnames(pheno)))})) 
      pheno <- pheno[, col_nr]
      ## Try to add Illumina stat
      # path_stat <- list.files(dirnames(pheno, pattern=".txt", full.names = TRUE)
      # if(length(path_stat) == 1){
      #               stat <- read.delim(paste0(pheno, "/Samples_stat.txt"), header=TRUE, sep="\t")
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
    if(is.data.frame(pheno)){
      header <- which(grepl("^Sample_ID|^sample_ID|^Sample_id|^sample_id", colnames(pheno)))
      if(!length(header) == 1){
        stop("Cannot find column named 'Sample_ID' \nor you have >1 columns named 'Sample_ID'")
      }
      colnames(pheno)[header] <- "Sample_ID"
      
    }else{
      lines <- readLines(pheno, n=20)
      header <- which(grepl("^Sample_ID|^sample_ID|^Sample_id|^sample_id", lines))
      if(!length(header) == 1){
        stop("Cannot find comma seperated header with first column \nnamed 'Sample_ID' or you have >1 columns named 'Sample_ID'")
        }
      head_1 <- stringr::str_count (lines[header], ",")
      row_1 <- stringr::str_count (lines[header+1], ",")
      if(head_1-row_1==0){
            pheno <- read.delim(pheno,  header=TRUE, sep=",")
      }else{
            pheno <- read.delim(pheno,  header=TRUE, sep="\t")
      }
    }
    pheno$Sample_ID <- as.character(gsub("-", "_", as.character(pheno$Sample_ID)))
    }
  
  
  ## Order as counts using grepl
  if(!is.null(counts)){
    # Search column names in pheno$Sample_ID
    typ <- "pheno"
    ord <- unlist(lapply(as.list(colnames(counts)), function(x){
      logi_colnam <- grepl(x, pheno$Sample_ID)
      if(sum(logi_colnam)>1){
        stop("Input pheno does not have unique Sample_IDs")
      }
      if(sum(logi_colnam) == 0){
        return(0)
      }else{
        return(which(logi_colnam))
      }
    }))
    # If n matches does not match n pheno do the other way around
    # Save what type you have done until later
    if(!length(ord[!ord==0])==nrow(pheno)){
      typ <- "counts"
      ord <- unlist(lapply(as.list(pheno$Sample_ID), function(x){
        logi_colnam <- grepl(x,  colnames(counts))
        if(sum(logi_colnam)>1){
          stop("Input pheno does not have unique Sample_IDs")
        }
        return(which(logi_colnam))
      }))
    }
    
    # Double check everything and report missing samples from typ=pheno 
    if(any(duplicated(ord[!ord==0]))){
      stop("Sample_IDs were not unique in pheno input.")
    }
    if(!length(ord[!ord==0])==nrow(pheno)){
      stop("\nNot all Sample_ID in pheno were available in counts column names.\nDouble check your Sample_ID column in pheno input.") 
    }
    
    
    # Reorder pheno according to counts sample names
    # Don't forget what type
    df <- as.data.frame(matrix(NA, nrow=ncol(counts), ncol=ncol(pheno)))
    rownames(df) <- colnames(counts)
    colnames(df) <- colnames(pheno)
    
    if(typ=="pheno"){
      for(i in 1:nrow(df)){
        if(!ord[i] == 0){
          df[i,] <- as.character(t(pheno[ord[i],]))
         }
      }
    }
    if(typ=="counts"){
      for(i in 1:nrow(df)){
       df[ord[i],] <- as.character(t(pheno[i,]))
      }
    }
    pheno <- df
      
    # Report outcome
    stopifnot(identical(rownames(pheno), colnames(counts)))
    logi_miss <- ord %in% 0
    if(any(logi_miss)){
      warning("Not all samples in counts were represented in pheno input.\n  These will have 'NA' in pheno.")
      }
    cat("Of", length(colnames(counts)), "sample names in counts,", sum(logi_miss), "were found in pheno file path.\n")
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




