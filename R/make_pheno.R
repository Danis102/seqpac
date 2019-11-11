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
#' @param pheno_file_path Character vector with the path to a .csv file.
#'
#' @param countTable Data.frame object with the same column names as in Sample_ID column of the .csv file.
#'
#' @param countTable Data.frame object with progress report.
#'
#' @return data.frame
#'
#' @examples
#' report <- make_pheno(<your_path_to_sports_output_directory>)
#' pheno <- make_pheno(pheno_file_path, countTable=countTable, progress_report=report)
#'
#' @export
make_pheno <- function(pheno_file_path, countTable=NULL, progress_report=NULL){
                            path <- list.files(pheno_file_path, pattern="*.csv", full.names = TRUE)
                                if(length(path) > 1){stop("Error! There are more than one .csv file to choose between in the specified folder.")}
                            lines <- readLines(path, n=20)
                            header <- which(grepl("\\<Sample_ID", lines))
                                if(!length(header) == 1){stop("Error! Cannot find comma seperated header with first column named 'Sample_ID'")}
                            pheno <- read.delim(paste0(pheno_file_path, "/SampleSheetUsed.csv"), skip= header-1,  header=TRUE, sep=",")

                                ## Illumina specific
                                if(header == 11){
                                      cat("Illumina type SampleSheetUsed.csv file was found.\n")
                                      colnames(pheno)[colnames(pheno) == "BASESPACE_ONLY_ORIGINAL_SAMPLE_NUMBER"] <- "Sample_number"
                                      pheno$Sample_number <- paste0("S", pheno$Sample_number)
                                      stopifnot(identical(as.character(pheno$Sample_ID), as.character(pheno$BASESPACE_ONLY_ORIGINAL_SAMPLE_ID)))
                                      stopifnot(identical(as.character(pheno$Sample_Name), as.character(pheno$BASESPACE_ONLY_ORIGINAL_SAMPLE_NAME)))
                                      pheno <- pheno[, c("Sample_ID", "Sample_Name", "SampleProject", "Sample_number", "Index")]
                                      pheno  <- cbind(data.frame(Sample = paste(pheno$Sample_Name, pheno$Sample_number, sep="_")), pheno)
                                      ## Try to add Illumina stat
                                      path_stat <- list.files(pheno_file_path, pattern=".txt", full.names = TRUE)
                                      if(length(path_stat) == 1){
                                                    stat <- read.delim(paste0(pheno_file_path, "/Samples_stat.txt"), header=TRUE, sep="\t")
                                                    stopifnot(any(as.character(pheno$Sample_ID) %in% as.character(stat$SampleId)))
                                                    pheno <- cbind(pheno, stat[match(as.character(pheno$Sample_ID), as.character(stat$SampleId)), !colnames(stat) %in% c("SampleId", "Name", "Index")])
                                          }

                                ## Other
                                }else{
                                      cat("A SampleSheetUsed.csv file was found but did not \n")
                                      cat("contain a header with 'Sample_ID' on the first column \n")
                                      cat("at row 11. This format is used for Illumina type \n")
                                      cat("SampleSheetUsed.csv. \n")
                                      cat("Cross-your fingers; I will attempt to read the file instead from the line where 'Sample_ID' was found...\n")
                                      pheno  <- cbind(data.frame(Sample = pheno[, "Sample_ID"]), pheno)
                                    }

                                rownames(pheno) <- pheno$Sample

                                ## Order as countTable
                                if(!is.null(countTable)){
                                    pheno <- pheno[match(colnames(countTable), rownames(pheno)),]
                                    stopifnot(identical(rownames(pheno), colnames(countTable)))
                                    cat("Of", length(colnames(countTable)), "sample names in countTable,", sum(as.numeric(rownames(pheno) %in% colnames(countTable))), "were found in pheno file path.\n")
                                    cat("\n")
                                    print(data.frame(pheno=rownames(pheno), countTable=colnames(countTable)))
                                }else{warning("No counTable was specified; Final Pheno will be unordered!\n")}

                                if(!is.null(progress_report)){
                                    cat("A progress report was specified, will attempt to match rownames...\n")
                                    progress_report <- progress_report[match(rownames(pheno), rownames(progress_report)),]
                                    stopifnot(identical(rownames(pheno), rownames(progress_report)))
                                    pheno <- cbind(pheno, progress_report)
                                    cat("\n")
                                }else{warning("No progress report was specified; will be missing the final Pheno.\n")}
                                cat("Done!")

                                return(pheno)
}




