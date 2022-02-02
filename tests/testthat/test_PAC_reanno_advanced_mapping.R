context("PAC reannoation and advanced mapping\n")
library(seqpac)


test_that("Testing reanno workflow ...", {

#----------------------------------------------
### Setup (quite function, create output folders, load pac)

  quiet <- function(x) { 
      sink(tempfile()) 
      on.exit(sink()) 
      invisible(force(x)) 
  }
  
  out_refs <- paste0(tempdir(), "/seqpac/test/refs")
  out_map <- gsub("/refs", "/map", out_refs)
  
  dir.create(out_refs, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_map, showWarnings = FALSE, recursive = TRUE)
  
  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                    package = "seqpac", mustWork = TRUE))
  
  pac$Anno <- pac$Anno[,1, drop=FALSE] 
  
#----------------------------------------------  
### Reannotate a genome
  # Move mycoplasma genome to temp and generate index  
  mycoplasma_path <- system.file("extdata/mycoplasma_genome", "mycoplasma.fa", 
                               package = "seqpac", mustWork = TRUE)
  
  file.copy(from=mycoplasma_path, to=out_refs, overwrite=TRUE, 
            recursive=FALSE, copy.mode=FALSE)
  
  mycoplasma_temp <- paste0(out_refs, "/mycoplasma.fa")
  
  Rbowtie::bowtie_build(mycoplasma_temp, 
                      outdir=gsub("mycoplasma.fa", "", mycoplasma_temp), 
                      prefix="mycoplasma", force=TRUE)
  
  # Run reanno workflow with genome mode
  ref_paths <- list(genome1= mycoplasma_temp,
                    genome2= mycoplasma_temp)
  
  quiet(
    map_reanno(PAC=pac, ref_paths=ref_paths, output_path=out_map,
                 type="internal", mismatches=3, import="genome", 
                 threads=2, keep_temp=TRUE, override = TRUE)
  )
  
  quiet( 
    reanno_genome <- make_reanno(out_map, PAC=pac, mis_fasta_check = TRUE, 
                                 output="list")
  )
  
  quiet(  
    pac <- add_reanno(reanno_genome, type="genome", genome_max=10, 
                      mismatches=1, merge_pac=pac)
  )
  
  # Clean up temp folder
  closeAllConnections()
  out_fls  <- list.files(out_map, recursive=TRUE, full.names = TRUE)
  suppressWarnings(file.remove(out_fls))
  
#----------------------------------------------  
## Biotype
  # Move tRNA and rRNA to temp and generate indexes   
  trna_path <- system.file("extdata/trna", "tRNA.fa", 
                           package = "seqpac", mustWork = TRUE)
  
  file.copy(from=trna_path, to=out_refs, overwrite=TRUE, 
            recursive=FALSE, copy.mode=FALSE)

  trna_temp <- paste0(out_refs, "/tRNA.fa")
  
  Rbowtie::bowtie_build(trna_temp, 
                        outdir=gsub("tRNA.fa", "", trna_temp), 
                        prefix="tRNA", force=TRUE)
  
  rrna_path <- system.file("extdata/rrna", "rRNA.fa", 
                           package = "seqpac", mustWork = TRUE)
  
  file.copy(from=rrna_path, to=out_refs, overwrite=TRUE, 
            recursive=FALSE, copy.mode=FALSE)
  
  rrna_temp <- paste0(out_refs, "/rRNA.fa")
  
  Rbowtie::bowtie_build(rrna_temp, 
                        outdir=gsub("rRNA.fa", "", rrna_temp), 
                        prefix="rRNA", force=TRUE)  
  
  
  ref_paths <- list(trna= trna_temp, rrna= rrna_temp)


  quiet(  
    map_reanno(pac, ref_paths=ref_paths, output_path=out_map,
               type="internal", mismatches=1,  import="biotype", 
               threads=2, keep_temp=FALSE, override = TRUE)
  )
  quiet(  
    reanno_biotype <- make_reanno(out_map, PAC=pac, mis_fasta_check = TRUE,
                                  output="list")
  )
    bio_search <- list(
                rrna=c("5S", "5.8S", "12S", "16S", "18S", "28S", "pre_S45"),
                trna =c("_tRNA", "mt:tRNA")
                    )
    quiet(  
    pac <- add_reanno(reanno_biotype, bio_search=bio_search, 
                      type="biotype", bio_perfect=FALSE, 
                      mismatches = 1, merge_pac=pac)
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
                           bio_name="Biotypes_mis1", merge_pac=TRUE)
    )
    quiet(
    pac <- PAC_summary(pac, pheno_target=list("stage"), norm="cpm")
    )
    
    # Clean up temp folder
    out_fls  <- list.files(gsub("/map", "", out_map), recursive=TRUE, 
                           full.names = TRUE)
    suppressWarnings(file.remove(out_fls))

#----------------------------------------------
### Tests reannotation workflow

# Test pac output

    expect_true(PAC_check(pac))
    test <- pac$Anno[!is.na(pac$Anno$mis0_genome1_genome),]
    expect_equal(dim(test), c(3,18))


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


#----------------------------------------------
### Test PAC_mapper, cov_plots

    quiet(
      pac  <- PAC_summary(pac, norm="cpm", pheno_target=list("stage"))
    )

    ref <- system.file("extdata/trna2", "tRNA2.fa",
                       package = "seqpac", mustWork = TRUE)

    ss <- system.file("extdata/trna2", "tRNA2.ss",
                       package = "seqpac", mustWork = TRUE)

    quiet(
      map_object <- PAC_mapper(pac, ref=ref, N_up = "NNN", N_down = "NNN",
                               mismatches=0, threads=2,
                               report_string=TRUE, override=TRUE)
    )
    quiet(
      plts <- PAC_covplot(pac, map=map_object,
                          summary_target= list("cpmMeans_stage"),
                          map_target="tRNA-Ala-AGC-1-1")
    )
    err <- try(plts[[1]])
    expect_true(class(map_object[[1]][[1]]) %in% "DNAStringSet")
    expect_equal(length(map_object), 290)
    expect_equal(class(plts[[1]]), c("gg", "ggplot"))
    expect_equal(class(err), c("gg", "ggplot"))

#----------------------------------------------
### Test map_rangetype
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

###----------------------------------------------
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


#----------------------------------------------
### Test PAC_gtf

  #First reload pac with original Anno
  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata",
                     package = "seqpac", mustWork = TRUE))

  # Create a gtf file
  anno <- pac$Anno
  anno <- anno[!grepl("Warning", anno$mis0_chromosomes_genome),]
  anno <- anno[!is.na(anno$mis0_chromosomes_genome),]
  coord <- anno$mis0_chromosomes_genome
  coord <- suppressWarnings(do.call("rbind", strsplit(coord, "\\|"))[,1])
  coord <- suppressWarnings(do.call("rbind", strsplit(coord, "\\;start=|;")))
  gr <- GenomicRanges::GRanges(seqnames=coord[,1],
                              IRanges::IRanges(as.numeric(coord[,2]),
                                               as.numeric(coord[,2])+anno$Size ),
                              strand=coord[,3])
  GenomicRanges::mcols(gr) <- data.frame(biotype=anno$Biotypes_mis3,
                                        bio_zero=as.character(anno$mis0_bio))
  spl <- sample(1:length(gr), round(length(gr)*0.3), replace=FALSE)
  gr1 <- gr[spl]
  gr2 <- gr[!1:length(gr) %in% spl]

  # Save artifical gtf in temp ref folder
  out1 <- paste0(out_refs, "/temp1.gtf")
  out2 <- paste0(out_refs, "/temp2.gtf")

  rtracklayer::export(gr1, out1, format="gtf")
  rtracklayer::export(gr2, out2, format="gtf")

  # Remove non-full coordinates from pac and make it smaller
  pac$Anno <-  pac$Anno [, grepl("chromosomes_genome|Size", colnames(pac$Anno))]

  rm_bad0 <- ifelse(grepl("Warning",  pac$Anno$mis0_chromosomes_genome), "rm", "keep")
  rm_bad1 <- ifelse(grepl("Warning",  pac$Anno$mis1_chromosomes_genome), "rm", "keep")
  rm_bad2 <- ifelse(grepl("Warning",  pac$Anno$mis2_chromosomes_genome), "rm", "keep")
  rm_bad3 <- ifelse(grepl("Warning",  pac$Anno$mis3_chromosomes_genome), "rm", "keep")
  pac$Anno$temp <- paste0(rm_bad0, rm_bad1, rm_bad2, rm_bad3)

  quiet(
    pac_test_large <- PAC_filter(pac, anno_target=list("temp", "keepkeepkeepkeep"))
  )
  quiet(
    pac_test_small <- PAC_filter(pac_test_large, size = c(20,22))
    )

  # Run PAC gtf in different settings
  gtf <- list(gtf1=out1, gtf2=out2)
  target <- list(gtf1=c("biotype","bio_zero"), gtf2=c("biotype","bio_zero") )

  quiet(
    simple_0 <- PAC_gtf(pac_test_large, mismatches=0, return="simplify",
                            gtf=gtf[1], target=target[1], threads=2)
    )
  quiet(
    simple_3 <- PAC_gtf(pac_test_small, mismatches=3, return="simplify",
                            gtf=gtf[1], target=target[1], threads=2)
    )
  quiet(
    full_1 <- PAC_gtf(pac_test_small, mismatches=1, return="full",
                      gtf=gtf[1], target=target[1], threads=2)
  )
  quiet(
    all_2 <- PAC_gtf(pac_test_small, mismatches=2, return="all",
                     gtf=gtf, target=target, threads=2)
  )

  # Generate error with wrong genome
  mycoplasma_path <- system.file("extdata/mycoplasma_genome", "mycoplasma.fa",
                                 package = "seqpac", mustWork = TRUE)

  file.copy(from=mycoplasma_path, to=out_refs, overwrite=TRUE,
            recursive=FALSE, copy.mode=FALSE)

  mycoplasma_temp <- paste0(out_refs, "/mycoplasma.fa")

  Rbowtie::bowtie_build(mycoplasma_temp,
                        outdir=gsub("mycoplasma.fa", "", mycoplasma_temp),
                        prefix="mycoplasma", force=TRUE)

  test_err  <- quiet(
    suppressWarnings(
      suppressMessages(try(PAC_gtf(pac_test_small, genome=mycoplasma_temp, mismatches=0,
                      return="simplify", gtf=gtf, target=target, threads=2),
                      silent = TRUE))
      )
    )

  # Tests for PAC_gtf output
  expect_true(sum(grepl("mis0|mis1|mis2|mis3", colnames(simple_3))) == 4)
  expect_true(sum(grepl("mis0", colnames(simple_0))) == 1)
  expect_true(grepl("Your chromosome names in", test_err))
  expect_true(sum(grepl("mis0|mis1|mis2|mis3", names(full_1))) > 800)
  expect_true(
    sum(colnames(full_1[[1]]) == c("seqid","start","end","strand"
                                   ,"biotype","bio_zero")) == 6)
  expect_true(sum(names(all_2) == c("simplify", "full")) == 2)
    
})

# Clean up temp
closeAllConnections()
fls_temp  <- list.files(tempdir(), recursive=TRUE, full.names = TRUE)
suppressWarnings(file.remove(fls_temp)) 
