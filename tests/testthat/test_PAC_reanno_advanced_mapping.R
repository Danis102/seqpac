context("PAC reannoation and advanced mapping\n")
library(seqpac)


  test_that("Testing reanno workflow ...", {  
  ################ 
  ### Setup

  quiet <- function(x) { 
      sink(tempfile()) 
      on.exit(sink()) 
      invisible(force(x)) 
      }

  load(system.file("extdata", "drosophila_sRNA_pac_filt.Rdata", 
                   package = "seqpac", mustWork = TRUE))
  pac = pac_cpm_filt
  
  ## Genome
  if(grepl("windows", .Platform$OS.type)){
    output <- paste0(tempdir(), "\\seqpac\\test")
  }else{
    output <- paste0(tempdir(), "/seqpac/test")
  }
  # tests uses testhat folder for getwd()
  mycoplasma_path <- paste0(getwd(), "/data_for_tests/mycoplasma_genome") 
  #mycoplasma_path <- paste0(getwd(), "/tests/testthat/data_for_tests/mycoplasma_genome")
  ref_paths <- list(genome1= list.files(mycoplasma_path, pattern=".fa", 
                                        full.names = TRUE),
                    genome2= list.files(mycoplasma_path, pattern=".fa", 
                                        full.names = TRUE))
  
  out_fls  <- list.files(output, recursive=TRUE, full.names = TRUE)
  suppressWarnings(file.remove(out_fls))
  quiet(
      map_reanno(PAC=pac, ref_paths=ref_paths, output_path=output,
                 type="internal", mismatches=3, import="genome", 
                 threads=8, keep_temp=FALSE)
  )
  
  quiet( 
    reanno_genome <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
  )
  quiet(  
    pac <- add_reanno(reanno_genome, type="genome", genome_max=10, 
                      mismatches=1, merge_pac=pac)
  )

  
  ## Biotype
  #output <- paste0(tempdir(), "/seqpac/test")
  #out_fls  <- list.files(output, recursive=TRUE)
  #suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
  
  #doParallel::stopImplicitCluster()
  closeAllConnections()
  out_fls  <- list.files(output, recursive=TRUE, full.names = TRUE)
  if(grepl("windows", .Platform$OS.type)){
    out_fls <- gsub("/", "\\\\", out_fls)
  }
  suppressWarnings(file.remove(out_fls))
  
  trna_path <- paste0(getwd(), "/data_for_tests/trna")
  rrna_path <- paste0(getwd(), "/data_for_tests/rrna")
  #trna_path <- paste0(getwd(), "/tests/testthat/data_for_tests/trna")
  #rrna_path <- paste0(getwd(), "/tests/testthat/data_for_tests/rrna")
  ref_paths <- list(trna= list.files(trna_path, pattern=".fa", 
                                     full.names = TRUE),
                    rrna= list.files(rrna_path, pattern=".fa", 
                                     full.names = TRUE))
  quiet(  
    map_reanno(pac, ref_paths=ref_paths, output_path=output,
               type="internal", mismatches=2,  import="biotype", 
               threads=8, keep_temp=FALSE)
  )
  quiet(  
    reanno_biotype <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
  )
    bio_search <- list(
                rrna=c("5S", "5.8S", "12S", "16S", "18S", "28S", "pre_S45"),
                trna =c("_tRNA", "mt:tRNA")
                    )
    quiet(  
    pac <- add_reanno(reanno_biotype, bio_search=bio_search, 
                      type="biotype", bio_perfect=FALSE, 
                      mismatches = 2, merge_pac=pac)
    )
    hierarchy <- list(rrna="rrna_",
                      trna="trna_"
                      )
    quiet(  
    pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=0, 
                           bio_name="Biotypes_mis0", merge_pac=TRUE)
    )
    quiet(  
    pac <- simplify_reanno(input=pac, hierarchy=hierarchy, mismatches=2, 
                           bio_name="Biotypes_mis2", merge_pac=TRUE)
    )
    quiet(
    pac <- PAC_summary(pac, pheno_target=list("stage"), norm="cpm")
    )
################ 
### Tests 

# Test output from reannotation 
 
    expect_true(PAC_check(pac))
    
    test <- pac$Anno[!is.na(pac$Anno$mis0_genome1_genome),]
    expect_equal(dim(test), c(3,19))


# Test plotting after reannotation
    # stackbar
    quiet( 
      test <- try(PAC_stackbar(pac, anno_target=list("Biotypes_mis0"), 
                               summary="pheno", pheno_target=list("stage")))
    )
    expect_equal(class(test), c("gg", "ggplot"))
    # sizedist
    quiet( 
      test <- try(PAC_sizedist(pac, anno_target=list("Biotypes_mis0"), 
                               pheno_target=list("stage")))
    )
    expect_equal(class(test[[1]][[1]]), c("gg", "ggplot"))
    # nbias
    quiet( 
      test <- try(PAC_nbias(pac, anno_target=list("Biotypes_mis0"), 
                               summary_target=list("cpmMeans_stage")))
    )
    expect_equal(class(test[[1]][[1]]), c("gg", "ggplot"))
    expect_equal(length(test[[1]]), 3)
    quiet(  
      test <- try(PAC_nbias(pac, position = 10, range = c(20,30), 
                            anno_target=list("Biotypes_mis0"), 
                            summary_target=list("cpmMeans_stage")))
    )
    expect_equal(class(test[[1]][[1]]), c("gg", "ggplot"))
    
    rgn <- c(min(as.character(test[[2]][[1]]$length)), 
             max(as.character(test[[2]][[1]]$length)))
    expect_equal(rgn, c("20","30"))


# Test PAC_mapper, cov_plots
    ref <- paste0(getwd(), "/data_for_tests/trna_no_index/tRNA_copy.fa")
    ss <- paste0(getwd(), "/data_for_tests/trna_no_index/tRNA.ss")
    #ref <- paste0(getwd(), "/tests/testthat/data_for_tests/trna_no_index/tRNA_copy.fa")
    #ss <- paste0(getwd(), "/tests/testthat/data_for_tests/trna_no_index/tRNA.ss")
    
    closeAllConnections()
    if(grepl("windows", .Platform$OS.type)){
      output <- paste0(tempdir(), "\\seqpac")
    }else{
      output <- paste0(tempdir(), "/seqpac")
    } 
    out_fls  <- list.files(output, recursive=TRUE)
    suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
    
    quiet(
      map_object <- PAC_mapper(pac, ref=ref, N_up = "NNN", N_down = "NNN", 
                               mapper="reanno", mismatches=0, threads=8, 
                               report_string=TRUE)
    )
    quiet(
      plts <- PAC_covplot(pac, map=map_object, 
                          summary_target= list("cpmMeans_stage"), 
                          map_target="Drosophila_melanogaster_tRNA-Ala-AGC-1-1")
    )
    err <- try(plts[[1]])
    expect_true(class(map_object[[1]][[1]]) %in% "DNAStringSet")
    expect_equal(length(map_object), 290)
    expect_equal(class(plts[[1]]), c("gg", "ggplot"))
    expect_equal(class(err), c("gg", "ggplot"))
    
# Test map_rangetype    
    quiet( 
      map_object_rng <- map_rangetype(map_object, type="nucleotides", 
                                      intervals=list(start=1:5, end=95:100))
    )
    quiet( 
      map_object_prc <- map_rangetype(map_object, type="percent", 
                                      intervals=list(start=1:5, mid=45:50, 
                                                     end=95:100))
    )
    quiet( 
      map_object_ss <- map_rangetype(map_object, type="ss", 
                                ss=ss, min_loop_width=4)      
    )
    expect_equal(stringr::str_count(paste0(map_object_rng[[1]][[1]]), 
                                    pattern = "N"), 6)
    expect_equal(stringr::str_count(map_object_rng[[1]][[2]]$Align_string[1], 
                                    pattern = "\\*"), 6)
    expect_equal(stringr::str_count(map_object_ss[[1]][[3]][3], 
                                    pattern = ">"), 21)
    
# Test tRNA_class and PAC_trna    
    quiet(
        pac_trna <- tRNA_class(pac, map_object_ss, terminal=5)
    )
    quiet(
        trna_result <- PAC_trna(pac_trna, norm="cpm", filter = NULL,
            join = FALSE, top = 15, log2fc = TRUE,
            pheno_target = list("stage", c("Stage1", "Stage3")), 
            anno_target_1 = list("type"),
            anno_target_2 = list("class"))
    )
    expect_true(PAC_check(pac_trna))
    
    col_nams <- c("Size","class","decoder","acceptor","tRNA_ref","type")
    expect_equal(names(pac_trna$Anno), col_nams) 
    
    expect_equal(class(trna_result$plots$Log2FC_Anno_1), c("gg", "ggplot"))
    expect_equal(class(trna_result$plots$Percent_bars$Stage1), c("gg", "ggplot"))
    expect_equal(length(trna_result$plots$Log2FC_Anno_1$layers), as.integer(3))

  })


  

