context("PAC reannoation and advanced mapping")
library(seqpac)

invisible(capture.output(
  test_that("Testing reanno workflow ...", {  
        ################ 
        ### Setup
        
        load(system.file("extdata", "drosophila_sRNA_pac_filt.Rdata", package = "seqpac", mustWork = TRUE))
        pac = pac_cpm_filt
        
        ## Genome
        mycoplasma_path <- paste0(getwd(), "/data_for_tests/mycoplasma_genome")
        #mycoplasma_path <- paste0(getwd(), "/tests/testthat/data_for_tests/mycoplasma_genome")
        ref_paths <- list(genome1= list.files(mycoplasma_path, pattern=".fa", full.names = TRUE),
                          genome2= list.files(mycoplasma_path, pattern=".fa", full.names = TRUE))
        output <- paste0(tempdir(), "/seqpac/test")
        out_fls  <- list.files(output, recursive=TRUE)
        suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
        map_reanno(PAC=pac, ref_paths=ref_paths, output_path=output,
                           type="internal", mismatches=3, import="genome", threads=8, keep_temp=FALSE)
        reanno_genome <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
        pac <- add_reanno(reanno_genome, type="genome", genome_max=10, mismatches=1, merge_pac=pac)
        
        
        ## Biotype
        output <- paste0(tempdir(), "/seqpac/test")
        out_fls  <- list.files(output, recursive=TRUE)
        suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
        
        trna_path <- paste0(getwd(), "/data_for_tests/trna")
        rrna_path <- paste0(getwd(), "/data_for_tests/rrna")
        #trna_path <- paste0(getwd(), "/tests/testthat/data_for_tests/trna")
        #rrna_path <- paste0(getwd(), "/tests/testthat/data_for_tests/rrna")
        ref_paths <- list(trna= list.files(trna_path, pattern=".fa", full.names = TRUE),
                          rrna= list.files(rrna_path, pattern=".fa", full.names = TRUE))
        
        map_reanno(pac, ref_paths=ref_paths, output_path=output,
                    type="internal", mismatches=2,  import="biotype", threads=8, keep_temp=FALSE)
        
        reanno_biotype <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
        bio_search <- list(
                         rrna=c("5S", "5.8S", "12S", "16S", "18S", "28S", "pre_S45"),
                         trna =c("_tRNA", "mt:tRNA")
                        )
        pac <- add_reanno(reanno_biotype, bio_search=bio_search, type="biotype", bio_perfect=FALSE, mismatches = 2, merge_pac=pac)
        hierarchy <- list(rrna="rrna_",
                          trna="trna_"
                          )
        pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=0, bio_name="Biotypes_mis0", merge_pac=TRUE)
        pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=2, bio_name="Biotypes_mis2", merge_pac=TRUE)

################ 
### Tests 

# Test output from reannotation 
 
    expect_true(PAC_check(pac))
    
    test <- pac$Anno[!is.na(pac$Anno$mis0_genome1_genome),]
    expect_equal(dim(test), c(3,19))


# Test plotting after reannotation
 
    test <- try(PAC_stackbar(pac, anno_target=list("Biotypes_mis0"), summary="pheno", pheno_target=list("stage")))
    expect_equal(class(test), c("gg", "ggplot"))
    test <- try(PAC_stackbar(pac, anno_target=list("Biotypes_mis0"), summary="pheno", pheno_target=list("stage")))
    expect_equal(class(test), c("gg", "ggplot"))


# Test PAC_mapper and cov_plots

    ref <- paste0(getwd(), "/data_for_tests/trna_no_index/tRNA_copy.fa")
    #ref <- paste0(getwd(), "/tests/testthat/data_for_tests/trna_no_index/tRNA_copy.fa")
    
    output <- paste0(tempdir(), "/seqpac")
    out_fls  <- list.files(output, recursive=TRUE)
    suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
    
    
    pac <- PAC_summary(pac, norm = "cpm", type = "means", pheno_target=list("stage"), merge_pac = TRUE)
    map_object <- PAC_mapper(pac, ref=ref, N_up = "NNN", N_down = "NNN", mapper="reanno", mismatches=0, threads=8, report_string=TRUE)
    plts <- PAC_covplot(pac, map=map_object, summary_target= list("cpmMeans_stage"), map_target="Drosophila_melanogaster_tRNA-Ala-AGC-1-1")
    err <- try(plts[[1]])
    expect_true(class(map_object[[1]][[1]]) %in% "DNAStringSet")
    expect_equal(length(map_object), 317)
    expect_equal(class(plts[[1]]), c("gg", "ggplot"))
    expect_equal(class(err), c("gg", "ggplot"))
    
    map_object_rng <- map_rangetype(map_object, type="nucleotides", intervals=list(start=1:5, end=95:100))
    map_object_prc <- map_rangetype(map_object, type="percent", intervals=list(start=1:5, mid=45:50, end=95:100))
    
    expect_equal(stringr::str_count(paste0(map_object_rng[[1]][[1]]), pattern = "N"), 6)
    expect_equal(stringr::str_count(map_object_rng[[1]][[2]]$Align_string[1], pattern = "\\*"), 6)

    })))

  

  

