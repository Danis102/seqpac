
# 
# 
# 
# extract_file <- function(file_path, mode=NULL, search=NULL, cut= NULL, 
#                           remove_search=FALSE){
#   if(mode="fasta"){
#     search <- "^>"
#     remove_search=TRUE
#   }
#   if(mode="fasta_first"){
#     search <- "^>"
#     cut= " "
#     remove_search=TRUE
#   }
#   `%>%` <- dplyr::`%>%`
#   text <- readr::read_delim(file_path, col_names = FALSE,  delim="¤¤")
#   
#   if(is.null(cut)){
#     text <- text %>% 
#       dplyr::filter(grepl(search, X1)) %>%
#          dplyr::pull(X1)
#   }
#   if(!is.null(cut)){     
#     text <- text %>% 
#       dplyr::filter(grepl(search, X1)) %>%
#         tidyr::separate(col=X1, sep=cut, into=c("X1"), 
#                         remove=TRUE, extra="drop") %>%
#           dplyr::pull(X1)   
#   }
#   if(remove_search==TRUE){
#         text <-  gsub(pattern = search, replacement = "", text)
#   }
# }


  