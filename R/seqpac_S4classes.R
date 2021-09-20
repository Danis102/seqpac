################################################################################
###-----------------------------------------------------------------------------
### This file contains all new S4 classes, generics and methods for seqpac. This
### involves:
### 1. The PAC object 
### 2. The reanno object 
### Generating these objects are controlled by (1) make_PAC and (2) make_reanno
### functions. From these functions, S3 versions of the objects can be
### generated using the 'output=' option. In addition S3 versions are available
### through coercion methods, eg. as(PAC,"S3"). The map object generated from
### the PAC_mapper function is a regular list.
################################################################################
###-----------------------------------------------------------------------------
### the S4 PAC object (make_PAC):
##
#' S4 object for PAC 
#' 
#' S4 version of the PAC object. 
#' 
#' Just like the S3 version, all subtables must have identical sequence (row)
#' names. To extract tables from S4 objects, use the getter commands or use @
#' (e.g. PAC@Pheno) just like $ for the S3 version (e.g. PAC$Pheno). You may
#' also use a coercion method, e.g. \code{PAC_S3 <- as(PAC_S4, "list")}
#' 
#' Holds up to 5 slots: \describe{
#'    \item{\strong{Pheno}}{data.frame > information about the samples}
#'    \item{\strong{Anno}}{data.frame > simplified sequences annotations}
#'    \item{\strong{Counts}}{data.frame > sequence counts over samples}
#'    \item{\emph{norm (optional)}}{list of data.frames > normalized counts}
#'    \item{\emph{summary (optional)}}{list of data.frames > summarized counts or normalized counts}
#'    }
#'    
#' @rdname PAC
#' @export
setClass(Class = "PAC",
          slots=c(Pheno = "data.frame", 
                  Anno= "data.frame",
                  Counts= "data.frame",
                  norm= "list",
                  summary="list",
                  reanno= "list"
                  ))

#-------------------------------
# Constructor PAC
#' @rdname PAC
#' @export
PAC <- function(Pheno, Anno, Counts, norm, summary, reanno){
  new("PAC", Pheno=Pheno, Anno=Anno, Counts=Counts,
      norm=norm, summary=summary, reanno=reanno)}

#-------------------------------
# Test validity
setValidity("PAC", function(object){
   # Is it empty
   if(is.null(length(object))){
     return("Empty object")
   }
   # Test 1
   pac <- as(object, "list")
   return(PAC_check(pac))
 })


#-------------------------------
# Coercion method PAC
setAs("PAC", "list",
      function(from){
        pac <- list(Pheno=from@Pheno,
                    Anno=from@Anno,
                    Counts=from@Counts)
        
        if(!is.null(from@norm[[1]])){
         pac <- c(pac, list(norm= from@norm))
        }
        if(!is.null(from@summary[[1]])){
         pac <- c(pac, list(summary= from@summary))
        }
        if(!is.null(from@reanno[[1]])){
         pac <- c(pac, list(reanno= from@reanno))
        }
        class(pac) <- c("PAC_S3","list")
        return(pac)
        })

###############################################
#' Converts an S3 PAC into a S4 PAC
#'
#' \code{as.PAC} Converts an S3 PAC object (list) to an S4 PAC object
#' 
#' Seqpac comes with two versions of the PAC object, either S4 or S3. The S3 PAC
#' is simply a list where each object can be received using the $-sign. The S4
#' version is a newer type of R object with predefined slots, that is received
#' using the @-sign. This function converts an S3 PAC object (list) to an S4 PAC
#' object. You can also use the S4 coercion method "as" to turn an S4 PAC into
#' an S3. See examples below.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param from S3 PAC-object containing at least a Pheno table with samples as
#'   row names, Anno table with sequences as row names and a Count table with
#'   raw counts (columns=samples, rows=sequences).
#' 
#' @return An S4 PAC object.
#'   
#' @examples
#' library(seqpac)
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'                   
#' # check type
#' 
#' class(pac)
#' isS4(pac)                        
#'  
#' # Turn S3 PAC object into a S4                                    
#' pac_s4 <- as.PAC(pac)
#'  
#' class(pac_s4)
#' isS4(pac_s4)   
#' 
#' # Turn S3 PAC object into a S4                                    
#' pac_s3 <- as(pac_s4, "list")
#' 
#' # Don't forget that in the summary and norm slots of the S4 PAC lies regular
#' # S3 lists. Thus, to receive these tables from an S4 PAC you need to combine
#' # both S4 and S3 receivers:
#' 
#' pac_s4 <- PAC_norm(pac_s4, norm = "cpm")
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", type = "means", 
#'                    pheno_target=list("stage"), merge_pac=TRUE)
#' 
#' pac_s4 
#' head(pac_s4@norm$cpm)
#' head(pac_s4@summary$cpmMeans_stage)
#' 
#' @export
#' 
as.PAC <- function(from){
       pac <- seqpac::PAC(Pheno=from$Pheno,
                   Anno=from$Anno,
                   Counts=from$Counts,
                   norm=list(NULL),
                   summary=list(NULL),
                   reanno=list(NULL)
                   )
       if("norm" %in% names(from)){
              pac@norm <- from$norm
        }
       if("summary" %in% names(from)){
              pac@summary <- from$summary
        }
       if("reanno" %in% names(from)){
              pac@reanno <- from$reanno
        }
       return(pac)
        }



#-------------------------------
# Converts an S3 PAC into a S4 PAC

# setAs("PAC_S3", "PAC",
#                 function(from){
#                      pac <- PAC(Pheno=from$Pheno,
#                                  Anno=from$Anno,
#                                  Counts=from$Counts,
#                                  norm=list(NULL),
#                                  summary=list(NULL),
#                                  reanno=list(NULL)
#                                  )
#                              if("norm" %in% names(from)){
#                                     pac@norm <- from$norm 
#                               }
#                              if("summary" %in% names(from)){
#                                     pac@summary <- from$summary 
#                               }
#                              if("reanno" %in% names(from)){
#                                     pac@reanno <- from$reanno 
#                               }
#                              return(pac)
#                               })


# as_PAC <- function(from){
#        pac <- PAC(Pheno=from$Pheno,
#                    Anno=from$Anno,
#                    Counts=from$Counts,
#                    norm=list(NULL),
#                    summary=list(NULL),
#                    reanno=list(NULL)
#                    )
#        if("norm" %in% names(from)){
#               pac@norm <- from$norm
#         }
#        if("summary" %in% names(from)){
#               pac@summary <- from$summary
#         }
#        if("reanno" %in% names(from)){
#               pac@reanno <- from$reanno
#         }
#        return(pac)
#         }

################################################################################
################################################################################
################################################################################
###-----------------------------------------------------------------------------
### the S4 reanno object (make_reanno):
#'
#' S4 object for reanno output 
#' 
#' Holds the imported information from mapping using the map_reanno function.
#' All information are held in tibble class (tbl/tbl_df) tables from the tibble
#' package. 
#' 
#' @return Contains two slots:  1. Overview: A table holding a summary of the
#' mapping. 2. Full_anno: Lists of tables holding the full imports per mismatch
#' cycle (mis0, mis1 etc) and reference (mi0-ref1, mis0-ref2, mis1-ref1,
#' mis1-ref2 etc).
#' 
#' @rdname reanno
#' 
#' @importClassesFrom tibble tbl_df
#' 
#' @export

setClass(Class = "reanno",
          slots=c(Overview = "tbl_df", 
                  Full_anno= "list"))

#-------------------------------
# Constructor reanno
#' @rdname reanno
#' @export
reanno <- function(Overview, Full_anno){
  new("reanno", Overview=Overview, Full_anno=Full_anno)}

#-------------------------------
# Test validity reanno

setValidity("reanno", function(object){
  # Is it empty
  if(is.null(length(object))){
    return("Empty object")
  }
  # Test 1
  seq_nam <- dplyr::pull(object@Overview[1])
  test1 <- lapply(object@Full_anno, function(y){
    lapply(y, function(z){identical(seq_nam, dplyr::pull(z[1]))})
  })
  logi_test1   <- unlist(test1)
  # Test 2
  test2 <- lapply(object@Full_anno, function(y){
    lapply(y, function(z){"tbl_df" %in%  class(z)})
  })
  logi_test2   <- unlist(test2)
  # Return
  if(any(!logi_test1)){ 
    return(cat("\nSequence names are not identical.",
               "\nPlease check @Full_anno tables:\n",
               names(logi_test1)[which(!logi_test1)]))
  }
  if(any(!logi_test2)){ 
    return(cat("\n@Full_anno is not a tibble.",
               "\nPlease check @Full_anno tables:\n",
               names(logi_test2)[which(!logi_test2)]))
  }
  TRUE
})

#-------------------------------
# Coercion method reanno
#
setAs("reanno", "list",
      function(from){
        rn <- list(Overview=from@Overview, 
                      Full_anno=from@Full_anno)
        class(rn) <- c("reanno_S3", "list")
        return(rn)
        })


###############################################
#' Converts an S3 PAC into a S4 PAC
#'
#' \code{as.reanno} Converts an S3 reanno object (list) to an S4 reanno object
#' 
#' Seqpac comes with two versions of the reanno object, either S4 or S3. The S3 
#' is simply a list where each object can be received using the $-sign. The S4
#' version is a newer type of R object with predefined slots, that is received
#' using the @-sign. This function converts an S3 reanno object (list) to an S4 
#' object. You can also use the S4 coercion method "as" to turn an S4 into
#' an S3. See examples below.
#' 
#' @family PAC reannotation
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param from S3 reanno object.
#' 
#' @return An S4 reanno object.
#'   
#' @examples
#' 
#'  library(seqpac)
#' 
#' ## load data
#'  load(system.file("extdata", "drosophila_sRNA_pac_filt.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'  pac = pac_cpm_filt
#'  
#'  ## Genome example
#'  # setup temporary output folder depending on windows or linux 
#'  
#'  if(grepl("windows", .Platform$OS.type)){
#'    output <- paste0(tempdir(), "\\seqpac\\test")
#'  }else{
#'    output <- paste0(tempdir(), "/seqpac/test")
#'  }
#'  
#'  Empty temp output folder
#'  out_fls  <- list.files(output, recursive=TRUE, full.names = TRUE)
#'  suppressWarnings(file.remove(out_fls))
#'  
#'  
#'  # Run reanno workflow loading test genome fasta
#'  mycoplasma_path <- paste0(getwd(), "/tests/testthat/data_for_tests/mycoplasma_genome")
#'  ref_paths <- list(genome1= list.files(mycoplasma_path, pattern=".fa", 
#'                                        full.names = TRUE),
#'                    genome2= list.files(mycoplasma_path, pattern=".fa", 
#'                                        full.names = TRUE))
#'  
#'  
#'
#'  map_reanno(PAC=pac, ref_paths=ref_paths, output_path=output,
#'                 type="internal", mismatches=3, import="genome", 
#'                 threads=8, keep_temp=FALSE)
#'
#'  reanno_genome <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE, 
#'                                 output="list")
#'                                 
#'                                 
#'  pac <- add_reanno(reanno_genome, type="genome", genome_max=10, 
#'                      mismatches=1, merge_pac=pac)
#' 
#' 
#' class(pac)
#' isS4(pac)                        
#'  
#' class(reanno_genome)
#' isS4(reanno_genome)
#' names(reanno_genome)   
#'  
#' # Turn S3 reanno object into a S4                                    
#' reanno_s4 <- as.reanno(reanno_genome)
#' class(reanno_s4)
#' isS4(reanno_s4) 
#'  
#' # Similar turn S3 PAC object into a S4                                    
#' pac_s4 <- as.PAC(pac)
#' class(pac_s4)
#' isS4(pac_s4)   
#' 
#' 
#' # Turn S3 reanno object back into a S4 using specific S4 coercion:                               
#' reanno_s3 <- as(reanno_s4, "list")
#' class(reanno_s3)
#' isS4(reanno_s3)
#' 
#' 
#' # Turn S3 PAC object into a S4                                    
#' pac_s3 <- as(pac_s4, "list")
#' 
#' # Don't forget that in the slots of S4 lies regular S3 objects. Thus,
#' # to receive these tables from an S4  you need to combine both S4 and S3
#' # receivers:
#' 
#' pac_s4 <- PAC_norm(pac_s4, norm = "cpm")
#' pac_s4 <- PAC_summary(pac_s4, norm = "cpm", type = "means", 
#'                    pheno_target=list("stage"), merge_pac=TRUE)
#' 
#' pac_s4 
#' head(pac_s4@norm$cpm)
#' head(pac_s4@summary$cpmMeans_stage)
#' 
#' @export
#' 
#' 
#'
#' 
as.reanno <- function(from){
       reanno <- seqpac::reanno(Overview=from$Overview,
                   Full_anno=from$Full_anno
                   )
       return(reanno)
        }



################################################################################
###-----------------------------------------------------------------------------




