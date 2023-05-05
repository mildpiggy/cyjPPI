

#' Enrich REACTOME complex
#'
#' Base on clusterProfiler::enricher, enrich REACTOME complex on giving protein list.
#'
#' @param proteinlist Vector, human uniprot protein ID
#'
#' @return
#'
#'
#' @export
#' A enrichResult instance
#' @examples
enrichComplex <- function(proteinlist){
  reactome_complex_gson = readRDS(system.file("extdata/ReactomeComplex","reactome_complex_pro_gson.rds",package = "cyjPPI"))
  enricher(proteinlist,gson = reactome_complex_gson, minGSSize = 3, maxGSSize = 100)
}

enrichComplex2 <- function(genelist){
  reactome_complex_gson2 = readRDS(system.file("extdata/ReactomeComplex","reactome_complex_gene_gson.rds",package = "cyjPPI"))
  enricher(genelist,gson = reactome_complex_gson2, minGSSize = 3, maxGSSize = 100)
}



