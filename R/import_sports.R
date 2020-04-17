#' Imports Sports Output Files Into R
#'
#' \code{import_sports} uses parallel processing to read Sports result files and
#' store them in a list of dataframes.
#'
#' Given the path to Sports output directory, this function will import all
#' result output.txt files of a given Sports project and store them in a list of
#' data.frames. This is done while perserving the basenames of the original
#' fastq file (input file for Sports). Files will be imported
#'
#' @family PAC generation
#'
#' @seealso \url{https://github.com/junchaoshi/sports1.0} for download and
#'   documentation about running Sports. \url{https://github.com/Danis102} for
#'   updates on the current package.
#'
#' @param path A character vector with the path to Sport output directory.
#'
#' @param threads An integer vector stating the number of parallell jobs.
#'
#' @return
#'   The function will search the Sports output directory (path) for all files
#'   ending with _output.txt and import them into a list of data.frames.
#'
#'   Using threads > 1, files will be imported in parallell, one file per
#'   dedicated thread. Use \code{parallel::detectcores()} to see available
#'   threads on the computer. Do not exceed the number of available threads!
#'
#'
#' @examples
#'   sports_lst <- import_sports(<your_path_to_sports_output_directory>, threads=<number_of_jobs>)
#'   sports_lst <- import_sports(path=path, threads=12)
#'
#' @export
import_sports <- function(path, threads=1){
                        require(parallel)
                        require(pbmcapply)
                        ### Get paths and sample names for output.txt files
                        count_files <- list.files(path, pattern ="*_output.txt", full.names=TRUE, recursive=TRUE)
                        count_files_nams <- list.files(path, pattern ="*_output.txt", full.names=FALSE, recursive=TRUE)
                        count_files_nams <- do.call("rbind", strsplit(count_files_nams, "/"))[,1]
                        count_files_nams <- gsub("\\<1_|_merge.\\>", "", count_files_nams)
                        count_files_nams <- gsub("-", "_", count_files_nams)
                        ## Read Sport output files
                        cat("Now reading sports output files using ", threads, " parallel threaded processes of", paste0(detectCores(logical = FALSE)), "possible threads\n")
                        cat("This may still take some time (approx. 1 sample/thread/minute) \n")
                        cat("Started reading ", paste0(Sys.time()), "\n")
                        #files <- pbmclapply(count_files, mc.cores = threads,  function(x){y <- read.table(x, header=TRUE); return(y)}) # Adds progress par to mclapply. Didn't work. Forum indicates that it may work with newer versions of R
                        files <- mclapply(count_files, mc.cores = threads,  function(x){y <- read.table(x, header=TRUE); return(y)})
                        names(files) <- count_files_nams
                        cat("Finished reading ", paste0(Sys.time()), "\n")
                        return(files)
                       }
