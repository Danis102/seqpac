#' Annotate with GTF file
#'
#' This function will annotate a PAC object using input from a GTF/SAF file.
#'
#' Wrapper to Rsubreads::featureCounts to generate annotations for sequences in
#' PAC using GTF or SAF file input. Given a PAC object that have been annotated
#' with genome coordinates using \code{PAC_genome} (PAC$Anno should contain
#' "chr","start","end", and "strand" columns), this function will add additional
#' annotations provided by overlaps with genomic features like gene exons and
#' repeats. The function has been developed for Ensembl and repeatMasker GTF
#' files, but should work with any standard GTF or SAF file.
#' 
#' 
#' @family PAC analysis
#'   
#' @seealso \url{https://github.com/Danis102} for updates.
#'
#' @param PAC PAC-list object. The Anno object needs to contain genome
#'   coordinates in the following columns: "chr","start","end", and "strand".
#'   These coordinates are eaisly obtained by using \code{PAC_genome}
#'   
#' @param genome
#' 
#' @param data_custom
#' @param data_custom
#' @param graphs Logical whether graphs should be reported or not.
#'   
#' @return List of ggplot2 plots and the data used for generating the plots. Use
#'   ls.str() to explore each level.
#'   
#' @examples


input="/data/Data_analysis/Genomes/Drosophila/dm6/Ensembl/dm6_ensembl_release_101/gtf/Drosophila_melanogaster.gtf"
gtf_feature="exon"
gtf_attrbute="gene_id"
PAC_gtf <- function(PAC, input=NULL, type="gtf", gtf_feature="exon", gtf_attrbute="gene_id")
                            
  
  
                            if(type=="gtf"){
                              if(class(input)=="data.frame"){
                                new_anno <- input
                              }else{
                                gtf <- rtracklayer::import(input)
                                gtf <- gtf[gtf$type==gtf_feature,]
                                gtf <- as.data.frame(gtf)
                                new_anno <- as.data.frame(data.table::data.table(GeneID=gtf[,gtf_attrbute], Chr=gtf$seqnames, Start=gtf$start, End=gtf$end, Strand=gtf$strand))
                             }}
                           if(type=="saf"){
                                if(class(input)=="data.frame"){
                                  new_anno <- input
                                }else{
                                  new_anno <- read.delim(input, header=TRUE)
                                }
                          stopifnot(colnames(new_anno)==c("GeneID", "Chr", "Start", "End", "Strand")) 
                                  
                                  
                                  
                                  
                                  