

#' Title
#'
#' @param PPIres
#' @param FC
#' @param layoutName
#'
#' @return
#' @export
#'
#' @examples
PPI2cyj <- function(PPIres, FC = NULL, layoutName){
  # print(class(PPIres))

  PPI_edges = PPIres
  if(all(c("from","to","combined_score") %in% colnames(PPI_edges))){
    PPI_edges <- PPI_edges %>% dplyr::rename(c(source = "from", target = "to",combined_score = "combined_score"))
  }else if(all(c("Gene1","Gene1","Score") %in% colnames(PPI_edges))){
    PPI_edges1111 <<- PPI_edges
    PPI_edges <- PPI_edges %>% dplyr::rename(c(source = "Gene1", target = "Gene2",combined_score = "Score"))
    PPI_edges2222 <<- PPI_edges
  }else{
    stop("PPIres should be the output data.frame from test_PPI or test_ReactomeFI")
  }

  PPI_edges$interaction <- T # here 'interaction' is a required column for dataFramesToJSON

  PPI_nodes <- data.frame(id = PPIres[,1:2] %>% unlist %>% unique())
  PPI_nodes$links_num = sapply(PPI_nodes[,1], function(x){
    c(PPIres[which(PPIres[,1] == x),2],
      PPIres[which(PPIres[,2] == x),1]
    ) %>% unique %>% length
  })
  # sig_tab <- get_signicant(DEPres,contrasts = "W4_vs_PBS",return_type = "table")
  # PPI_nodes$lfc = sig_tab[match(PPI_nodes$id,sig_tab$name),"W4_vs_PBS_ratio"]

  if(inherits(FC, "data.frame")){
    PPI_nodes$lfc = FC[match(PPI_nodes$id,FC$name), 2]
  }else if(is.vector(FC)){
    PPI_nodes$lfc = FC[match(PPI_nodes$id,names(FC))]
  }

  PPI_nodes$highlighted =
    ifelse(PPI_nodes$links_num > 3, T, F)
  # sample(c(T,F),nrow(PPI_nodes),replace = T)
  # rep(F,nrow(PPI_nodes))


  PPI_nodes111 <<- PPI_nodes
  graph.json <- dataFramesToJSON(PPI_edges, PPI_nodes)

  styles = cyjPPI::styles()
  cat("use style:", styles["Default"],"\n")
  cyjShiny(graph=graph.json,
           layoutName= layoutName,
           styleFile = styles["Default"]
  )
}


## styles
styles <- function(){
  style_path = system.file("extdata/styles",package = "cyjPPI")
  styles = c("Default" = file.path(style_path, "defaultStyle2.json"), ## 所设计的 style
             "Manual", # 自定义的 style, 基于 default 改造
             "Biological"= file.path(style_path,"biologicalStyle.js"),
             # "expermentStyle" = file.path(style_path,"expermentStyle.js"),
             "galFiltered" = file.path(style_path,"galFiltered-style.json"),
             "smallDemoStyle" = file.path(style_path,"smallDemoStyle.json"),
             "style01" = file.path(style_path,"style01.js"), # 等同于 biologicalStyle.js
             "styleCancer" = file.path(style_path,"styleCancer_21July2017.js"),
             "yeastGalactoseStyle" = file.path(style_path,"yeastGalactoseStyle.js")
  )
}


## The function to transform gene symbol to human homologous gene symbol.
# x, the gene symbol vector.
# species, Which species to transform including "Mouse","Rat","Canine"
Symbol_transform <- function(x,species){
  x = unique(x)
  all_reactome_gene_map_table = readRDS(system.file("extdata","homologousgene","reactome_gene_MouseRat_map.rds",package = "cyjPPI"))
  all_reactome_gene_map_table2 = all_reactome_gene_map_table[which(all_reactome_gene_map_table[,species] %in% x),c
                                                             ("Human","Human_ID",species,paste0(species,"_ID"))
  ]
  all_reactome_gene_map_table2 = all_reactome_gene_map_table2[!duplicated(all_reactome_gene_map_table2[,species]),]
  cat("Input ",length(x)," unique genes: ",length(all_reactome_gene_map_table2[,species]), " matched a homologous human gene in reactome FI data",
      length(all_reactome_gene_map_table2[,species])/length(x)*100,"% \n", sep="")
  return(all_reactome_gene_map_table2)
}


get_FI <- function(gene_symbol, species, reactome_FI_table){
  gene_symbol = na.omit(gene_symbol)
  if(species !="Human"){
    homogene_symbol = Symbol_transform(gene_symbol,species)
    gene_symbol <- homogene_symbol[,"Human"]
    if(length(gene_symbol) < 2) stop("Too few giving gene can be mapped to reactome human gene,please check your input")
  }

  sub_FI = reactome_FI_table %>% dplyr::filter( (Gene1 %in% gene_symbol) & (Gene2 %in% gene_symbol))
}


## the function to get ReactomeFI function, from DEP2::test_ppi
#' Title
#'
#' @param species
#' @param reactome_FI_data
#' @param score_cutoff A numeric less than 1, required lowest interaction scores. A score lowwer than 400 means the interaction is unconfident
#'
#' @return
#' @export
#'
#' @inheritParams DEP2::test_PPI
#'
#' @examples
#' data()
test_ReactomeFI <- function (x, contrasts = NULL, species = "Human",
                             reactome_FI_data = NULL, score_cutoff = 0.5)
{
  assertthat::assert_that(class(x) == "SummarizedExperiment" |
                            class(x) == "DEGdata" | class(x) == "character", is.null(contrasts) |
                            is.character(contrasts), is.character(species) && length(species) ==
                            1,
                          is.numeric(score_cutoff) | is.integer(score_cutoff),
                          is.null(reactome_FI_data)|is.data.frame(reactome_FI_data),
                          dplyr::between(score_cutoff,0,1))

  if(is.null(reactome_FI_data)){
    data(reactome_FI_table)
    reactome_FI_data = reactome_FI_table
    # reactome_FI_data <- readRDS(system.file("extdata/e"))
  }


  if (class(x) == "character") {
    gene = x
  }else if (class(x) == "DEGdata" | class(x) == "SummarizedExperiment") {
    sig = get_signicant(x, contrasts = contrasts, return_type = "subset")
    if (all(c("gene_name", "protein_ID") %in% colnames(rowData(x)))) {
      gene = rowData(x)$gene_name
    }else {
      gene = rownames(sig)
    }
  }
  if (length(gene) < 2) {
    stop("Giving candidate is too few, please check your input")
  }
  reactome_FI = get_FI(gene,species,reactome_FI_data)
  if(!is.null(reactome_FI)) reactome_FI = dplyr::filter(reactome_FI,Score >= score_cutoff)

  return(reactome_FI)
}


## Creat the PPI R6 shiny object
#' Creat a cytoscape.js shiny based on the
#'
#' @param PPI_res The result from test_ppi or test_reactome
#' @param Annotation_enrichment_res A enrichResult object (or its data.frame result) that stores enrichment analysis result
#' @param DEP_res The dep result from DEP2/DEP
#'
#' @return
#' A PPI.enrich.R6.Shiny R6 class object, containing the analysis results and shiny application
#'
#' @export
#'
#' @examples
creat_cyjFI = function(PPI_res, Annotation_enrichment_res, DEP_res){
  theapp.r6 = PPI.enrich.R6.Shiny$new(PPI_res = PPI_res,
                                      ORA_res = Annotation_enrichment_res,
                                      DEP_res = dep_pg)
}



get_lfc = function(dep, contrast = NULL){
  assertthat::assert_that(inherits(dep,"SummarizedExperiment"))

  if(is.null(contrast)){
    contrast = get_contrast(dep)[1]
  }
  if(!contrast %in% get_contrast(dep)){
    stop(contrast,"is not a contrast in input dep.")
  }
  lfc = get_df_wide(dep)[,paste0(contrast,"_diff")]
  names(lfc) = names(dep)
  return(lfc)
}

