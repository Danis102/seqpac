# global.R - Load libraries and source modular components

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(seqpac)
library(shinyjs)

# Source modules
source("modules/module_launch.R")
source("modules/module_load.R")
source("modules/module_filter.R")
source("modules/module_annotate.R")
source("modules/module_analyze.R")

# Enable file upload max size (e.g. 500MB for FASTQ files)
options(shiny.maxRequestSize = 500 * 1024^2)
