
################################################################################
#  Generics PAC

#------------------------------------------------------------------------------
#' @rdname pheno
#' @export
setGeneric("pheno", function(object){standardGeneric("pheno")})

#------------------------------------------------------------------------------
#' @rdname pheno
#' @export
setGeneric("pheno<-", function(object){standardGeneric("pheno<-")})

#------------------------------------------------------------------------------
#' @rdname anno
#' @export
setGeneric("anno", function(object){standardGeneric("anno")})

#------------------------------------------------------------------------------
#' @rdname anno
#' @export
setGeneric("anno<-", function(object){standardGeneric("anno<-")})

#------------------------------------------------------------------------------
#' @rdname counts
#' @export
setGeneric("counts", function(object){standardGeneric("counts")})

#------------------------------------------------------------------------------
#' @rdname counts
#' @export
setGeneric("counts<-", function(object){standardGeneric("counts<-")})

#------------------------------------------------------------------------------
#' @rdname norm
#' @export
setGeneric("norm", function(object){standardGeneric("norm")})

#------------------------------------------------------------------------------
#' @rdname norm
#' @export
setGeneric("norm<-", function(object){standardGeneric("norm<-")})

#------------------------------------------------------------------------------
#' @rdname summary
#' @export
setGeneric("summary", function(object){standardGeneric("summary")})

#------------------------------------------------------------------------------
#' @rdname summary
#' @export
setGeneric("summary<-", function(object){standardGeneric("summary<-")})



################################################################################
#  Generics reanno

#------------------------------------------------------------------------------
#' @rdname overview
#' @export
setGeneric("overview", function(x){standardGeneric("overview")})

#------------------------------------------------------------------------------
#' @rdname overview
#' @export
setGeneric("overview<-", function(x){standardGeneric("overview<-")})

#------------------------------------------------------------------------------
#' @rdname full
#' @export
setGeneric("full", function(x){standardGeneric("full")})

#------------------------------------------------------------------------------
#' @rdname full
#' @export
setGeneric("full<-", function(x){standardGeneric("full<-")})



################################################################################
# S4 methods PAC

#------------------------------------------------------------------------------
pheno.PAC <- function(object){ object@Pheno }
#' sample information
#' @rdname pheno
#' @export
setMethod("pheno", "PAC", pheno.PAC)

#------------------------------------------------------------------------------
anno.PAC <- function(object){ object@Anno }
#' annotation table
#' @rdname anno
#' @export
setMethod("anno", "PAC", anno.PAC)

#------------------------------------------------------------------------------
counts.PAC <- function(object){ object@Counts }
#' raw counts table
#' @rdname counts
#' @export
setMethod("counts", "PAC", counts.PAC)

#------------------------------------------------------------------------------
norm.PAC <- function(object){ object@norm }
#' list of normalized tables
#' @rdname norm
#' @export
setMethod("norm", "PAC", norm.PAC)

#------------------------------------------------------------------------------
summary.PAC <- function(object){
 object <- as(object, "S3")
 tib_sum <- lapply(object$summary, function(y){tibble::as_tibble(y)})
 return(tib_sum)
}
#' list of summary tables
#' @rdname summary
#' @export
setMethod("summary", signature(object="PAC"), summary.PAC)

#------------------------------------------------------------------------------
names.PAC <- function(x){ names(as(x, "list")) }
#' @rdname names
#' @export
setMethod("names", "PAC", names.PAC)

#------------------------------------------------------------------------------
rownames.PAC <- function(x){ rownames(x@Counts) }
#' sequence names
#' @rdname names
#' @export
setMethod("rownames", "PAC", rownames.PAC)


#------------------------------------------------------------------------------
colnames.PAC <- function(x){ colnames(x@Counts) }
#' sample names
#' @rdname colnames
#' @export
setMethod("colnames", "PAC", colnames.PAC)

#------------------------------------------------------------------------------
length.PAC <- function(x){ length(as(x, "list")) }
#' number of objects in PAC
#' @rdname length
#' @export
setMethod("length", "PAC", length.PAC)

#------------------------------------------------------------------------------
ncol.PAC <- function(x){ ncol(x@Counts) }
#' number of samples
#' @rdname ncol
#' @export
setMethod("ncol", "PAC", ncol.PAC)

#------------------------------------------------------------------------------
nrow.PAC <- function(x){ nrow(x@Counts) }
#' number of sequences
#' @rdname nrow
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
  cat("   sequence with highest:", mx_seq, "mean counts\n")
  cat("    sequence with lowest:", mn_seq, "mean counts\n")
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
#' @export
setMethod("show", signature(object="PAC"), show.PAC)



################################################################################
# S4 methods reanno

#------------------------------------------------------------------------------
overview.reanno <- function(x){ x@Overview }
#' overview table of reanno object
#' @rdname overview
#' @export
setMethod("overview", "reanno", overview.reanno)

#------------------------------------------------------------------------------
full.reanno <- function(x){ x@Full_anno } 
#' full annotation list of reanno object
#' @rdname full
#' @export
setMethod("full", "reanno", full.reanno)

#------------------------------------------------------------------------------
names.reanno <- function(x){ names(as(x, "list"))}
#' names of objects in reanno
#' @rdname names
#' @export
setMethod("names", "reanno", names.reanno)

#------------------------------------------------------------------------------
rownames.reanno <- function(x){ dplyr::pull(x@Overview[1])}
#' sequence names in reanno
#' @rdname rownames
#' @export
setMethod("rownames", "reanno", rownames.reanno)

#------------------------------------------------------------------------------
length.reanno <- function(x){ length(as(x, "list")) } 
#' length of reanno object
#' @rdname length
#' @export
setMethod("length", "reanno", length.reanno)

#------------------------------------------------------------------------------
nrow.reanno <- function(x){ nrow(x@Overview) } 
#' number of sequences
#' @rdname nrow
#' @export
setMethod("nrow", "reanno", nrow.reanno)
