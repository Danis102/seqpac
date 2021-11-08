
################################################################################
#  Generics PAC

#------------------------------------------------------------------------------
#' @rdname pheno
#' @return The Pheno data frame table from a PAC-object
#' @export
setGeneric("pheno", function(object){standardGeneric("pheno")})

#------------------------------------------------------------------------------
#' @rdname pheno
#' @export
setGeneric("pheno<-", function(value){standardGeneric("pheno<-")})

#------------------------------------------------------------------------------
#' @rdname anno
#' @return The Anno data frame table from a PAC-object
#' @export
setGeneric("anno", function(object){standardGeneric("anno")})

#------------------------------------------------------------------------------
#' @rdname anno
#' @export
setGeneric("anno<-", function(value){standardGeneric("anno<-")})

#------------------------------------------------------------------------------
#' @rdname counts
#' @return The Counts data frame table from a PAC-object
#' @export
setGeneric("counts", function(object){standardGeneric("counts")})

#------------------------------------------------------------------------------
#' @rdname counts
#' @export
setGeneric("counts<-", function(value){standardGeneric("counts<-")})

#------------------------------------------------------------------------------
#' @rdname norm
#' @return The list of normalized counts tables from a PAC-object
#' @export
setGeneric("norm", function(object){standardGeneric("norm")})

#------------------------------------------------------------------------------
#' @rdname norm
#' @export
setGeneric("norm<-", function(value){standardGeneric("norm<-")})

#------------------------------------------------------------------------------
#' @rdname summary
#' @return The list of summary tables from a PAC-object
#' @export
setGeneric("summary", function(object){standardGeneric("summary")})

#------------------------------------------------------------------------------
#' @rdname summary
#' @export
setGeneric("summary<-", function(value){standardGeneric("summary<-")})



################################################################################
#  Generics reanno

#------------------------------------------------------------------------------
#' @rdname overview
#' @return The overview table of a reanno-object as a data frame. 
#' @export
setGeneric("overview", function(x){standardGeneric("overview")})

#------------------------------------------------------------------------------
#' @rdname overview
#' @export
setGeneric("overview<-", function(value){standardGeneric("overview<-")})

#------------------------------------------------------------------------------
#' @rdname full
#' @return All the results that were imported into the reanno-object returned as
#'   a list of data frames.
#' @export
setGeneric("full", function(x){standardGeneric("full")})

#------------------------------------------------------------------------------
#' @rdname full
#' @export
setGeneric("full<-", function(value){standardGeneric("full<-")})



################################################################################
# S4 methods PAC

#------------------------------------------------------------------------------
pheno.PAC <- function(object){ object@Pheno }
#' sample information
#' @rdname pheno
#' @param pheno (PAC-object)
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("pheno", "PAC", pheno.PAC)

#------------------------------------------------------------------------------
anno.PAC <- function(object){ object@Anno }
#' annotation table
#' @rdname anno
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("anno", "PAC", anno.PAC)

#------------------------------------------------------------------------------
counts.PAC <- function(object){ object@Counts }
#' raw counts table
#' @rdname counts
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("counts", "PAC", counts.PAC)

#------------------------------------------------------------------------------
norm.PAC <- function(object){ object@norm }
#' list of normalized tables
#' @rdname norm
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("norm", "PAC", norm.PAC)

#------------------------------------------------------------------------------
summary.PAC <- function(object){ object@summary }
#' list of summary tables
#' @rdname summary
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("summary", "PAC", summary.PAC)

#------------------------------------------------------------------------------
names.PAC <- function(x){ names(as(x, "list")) }
#' @rdname names
#' @return The names of the content in a PAC-object. 
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("names", "PAC", names.PAC)

#------------------------------------------------------------------------------
rownames.PAC <- function(x){ rownames(x@Counts) }
#' sequence names
#' @rdname rownames
#' @return The sequences (row names) held by the Counts/Anno tables of a
#'   PAC-object.
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("rownames", "PAC", rownames.PAC)


#------------------------------------------------------------------------------
colnames.PAC <- function(x){ colnames(x@Counts) }
#' sample names
#' @rdname colnames
#' @return The sample names held by the Counts table (column names) and the
#'   Pheno table (row names) of a PAC-object.
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("colnames", "PAC", colnames.PAC)

#------------------------------------------------------------------------------
length.PAC <- function(x){ length(as(x, "list")) }
#' number of objects in PAC
#' @rdname length
#' @return The number of items in a PAC-object. 
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("length", "PAC", length.PAC)

#------------------------------------------------------------------------------
ncol.PAC <- function(x){ ncol(x@Counts) }
#' number of samples
#' @rdname ncol
#' @return The number of samples in a PAC-object. 
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("ncol", "PAC", ncol.PAC)

#------------------------------------------------------------------------------
nrow.PAC <- function(x){ nrow(x@Counts) }
#' number of sequences
#' @rdname nrow
#' @return The number of sequences in a PAC-object. 
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("nrow", "PAC", nrow.PAC)

#------------------------------------------------------------------------------
show.PAC <- function(object) {
  pac <- as(object, "list")
  cat("PAC object with: \n")
  cat("  ", nrow(pac$Pheno), "samples\n")
  cat("  ",nrow(pac$Anno), "sequences\n")
  avg <- round(mean(colSums(pac$Counts)), digits=0)
  mn <- round(min(colSums(pac$Counts)), digits=0)
  mx <- round(max(colSums(pac$Counts)), digits=0)
  mx_seq <- round(max(rowMeans(pac$Counts)), digits=0)
  mn_seq <- round(min(rowMeans(pac$Counts)), digits=0)
  cat(paste0("   mean total counts: ", avg, " (min:", mn, "/max:", mx, ")\n"))
  cat("   best sequence:", mx_seq, "mean counts\n")
  cat("   worst sequence:", mn_seq, "mean counts\n")
  if("norm" %in% names(object)){
    cat("normalized tables:", length(pac$norm),"\n")
    cat(names(pac$norm),"\n")
  }
  if("summary" %in% names(object)){
    cat("summarized tables:", length(pac$summary),"\n")
    cat(names(pac$summary),"\n")
  }
}
#' overview PAC object
#' @rdname show
#' @return Overview of a PAC-object.
#' @importFrom methods as new setMethod
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'
#' # extra fuctionality with s4 PAC-object:
#' pac_s4 <- as.PAC(pac)
#' pac_s4
#' names(pac_s4)
#' length(pac_s4)
#' nrow(pac_s4)
#' ncol(pac_s4)
#' rownames(pac_s4)
#' colnames(pac_s4)
#' pheno(pac_s4) 
#' head(anno(pac_s4))
#' head(counts(pac_s4))
#' head(norm(pac_s4)$cpm)
#' 
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage")) 
#' pac_s4                   
#' head(summary(pac_s4)$cpmMeans_stage)
#' 
#' @export
setMethod("show", signature(object="PAC"), show.PAC)



################################################################################
# S4 methods reanno

#------------------------------------------------------------------------------
overview.reanno <- function(x){ x@Overview }
#' overview table of reanno object
#' @rdname overview
#' @return Overview table of reanno-object.
#' @examples
#' ######################################################### 
#' ##### Create an reanno object
#' 
#' ##  First load a PAC- object
#' 
#'  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#'  pac$Anno <- pac$Anno[,1, drop = FALSE]
#'  
#'  
#' ##  Then specify paths to fasta references
#' # If you are having problem see the vignette small RNA guide for more info.
#'  
#'  trna_path <- system.file("extdata/trna", "tRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)  
#'  rrna_path <- system.file("extdata/rrna", "rRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)
#'  
#'  ref_paths <- list(trna= trna_path, rrna= rrna_path)
#'                                     
#' ##  Add output path of your choice.
#' # Here we use the R temporary folder depending on platform                                     
#'if(grepl("windows", .Platform$OS.type)){
#'  output <- paste0(tempdir(), "\\seqpac\\test")
#'}else{
#'  output <- paste0(tempdir(), "/seqpac/test")}
#' 
#' ## Make sure it is empty (otherwise you will be prompted for a question)
#' out_fls  <- list.files(output, recursive=TRUE)
#' closeAllConnections()
#' suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
#'
#' ##  Then map your PAC-object against the fasta references                                  
#'  map_reanno(pac, ref_paths=ref_paths, output_path=output,
#'                type="internal", mismatches=2,  import="biotype", 
#'                threads=2, keep_temp=FALSE)
#' 
#' ##  Then generate a reanno-object of the temporary bowtie-files 
#' reanno_object <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
#' 
#'## Accessing content and S4/S3 conversion:
#' names(reanno_object)
#' overview(reanno_object)
#' full(reanno_object)
#' rownames(reanno_object)
#' length(reanno_object)
#' nrow(reanno_object)
#' reanno_s3 <- as(reanno_object, "list")
#' reanno_s4 <- as.reanno(reanno_s3)
#' 
#' @export
setMethod("overview", "reanno", overview.reanno)

#------------------------------------------------------------------------------
full.reanno <- function(x){ x@Full_anno } 
#' full annotation list of reanno object
#' @rdname full
#' @return All the results that were imported into the reanno-object returned as
#'   a list of data frames.
#' @examples
#' ######################################################### 
#' ##### Create an reanno object
#' 
#' ##  First load a PAC- object
#' 
#'  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#'  pac$Anno <- pac$Anno[,1, drop = FALSE]
#'  
#'  
#' ##  Then specify paths to fasta references
#' # If you are having problem see the vignette small RNA guide for more info.
#'  
#'  trna_path <- system.file("extdata/trna", "tRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)  
#'  rrna_path <- system.file("extdata/rrna", "rRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)
#'  
#'  ref_paths <- list(trna= trna_path, rrna= rrna_path)
#'                                     
#' ##  Add output path of your choice.
#' # Here we use the R temporary folder depending on platform                                     
#'if(grepl("windows", .Platform$OS.type)){
#'  output <- paste0(tempdir(), "\\seqpac\\test")
#'}else{
#'  output <- paste0(tempdir(), "/seqpac/test")}
#' 
#' ## Make sure it is empty (otherwise you will be prompted for a question)
#' out_fls  <- list.files(output, recursive=TRUE)
#' closeAllConnections()
#' suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
#'
#' ##  Then map your PAC-object against the fasta references                                  
#'  map_reanno(pac, ref_paths=ref_paths, output_path=output,
#'                type="internal", mismatches=2,  import="biotype", 
#'                threads=2, keep_temp=FALSE)
#' 
#' ##  Then generate a reanno-object of the temporary bowtie-files 
#' reanno_object <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
#' 
#'## Accessing content and S4/S3 conversion:
#' names(reanno_object)
#' overview(reanno_object)
#' full(reanno_object)
#' rownames(reanno_object)
#' length(reanno_object)
#' nrow(reanno_object)
#' reanno_s3 <- as(reanno_object, "list")
#' reanno_s4 <- as.reanno(reanno_s3)  
#' @export
setMethod("full", "reanno", full.reanno)

#------------------------------------------------------------------------------
names.reanno <- function(x){ names(as(x, "list"))}
#' names of objects in reanno
#' @rdname names
#' @return Names of the items in the reanno-object.
#' @export
setMethod("names", "reanno", names.reanno)

#------------------------------------------------------------------------------
rownames.reanno <- function(x){ dplyr::pull(x@Overview[1])}
#' sequences in reanno
#' @rdname rownames
#' @return Sequences in the reanno-object.
#' @export
setMethod("rownames", "reanno", rownames.reanno)

#------------------------------------------------------------------------------
length.reanno <- function(x){ length(as(x, "list")) } 
#' length of reanno object
#' @rdname length
#' @return Number of items in the reanno-object.
#' @export
setMethod("length", "reanno", length.reanno)

#------------------------------------------------------------------------------
nrow.reanno <- function(x){ nrow(x@Overview) } 
#' number of sequences
#' @rdname nrow
#' @return Number of sequences in the reanno-object.
#' @export
setMethod("nrow", "reanno", nrow.reanno)
