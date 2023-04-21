PPI2cyj <- function(PPIres, L2FC_tab = NULL, layoutName){
  print(class(PPIres))
  PPI_nodes <- data.frame(id = PPIres[,1:2] %>% unlist %>% unique())
  PPI_nodes$links_num = sapply(PPI_nodes[,1], function(x){
    c(PPIres[which(PPIres[,1] == x),2],
      PPIres[which(PPIres[,2] == x),1]
    ) %>% unique %>% length
  })
  # sig_tab <- get_signicant(DEPres,contrasts = "W4_vs_PBS",return_type = "table")
  # PPI_nodes$lfc = sig_tab[match(PPI_nodes$id,sig_tab$name),"W4_vs_PBS_ratio"]

  if(inherits(L2FC_tab, "data.frame")){
    PPI_nodes$lfc = L2FC_tab[match(PPI_nodes$id,L2FC_tab$name), 2]
  }else if(is.vector(L2FC_tab)){
    PPI_nodes$lfc = L2FC_tab[match(PPI_nodes$id,names(L2FC_tab))]
  }

  PPI_nodes$highlighted =
    ifelse(PPI_nodes$links_num > 3, T, F)
  # sample(c(T,F),nrow(PPI_nodes),replace = T)
  # rep(F,nrow(PPI_nodes))

  PPI_edges = PPIres
  PPI_edges <- PPI_edges %>% dplyr::rename(c(source = "from", target = "to",combined_score = "combined_score"))
  PPI_edges$interaction <- T

  graph.json <- dataFramesToJSON(PPI_edges, PPI_nodes)

  cat("use style:", styles["Default"],"\n")
  cyjShiny(graph=graph.json,
           layoutName= layoutName,
           styleFile = styles["Default"]
  )
}


## styles
styles <- function(){
  style_path = system.file("inst/styles",package = "cyjPPI")
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
  all_reactome_gene_map_table = readRDS("all_reactome_gene_map.rds")
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

  sub_FI = reactome_FI_table %>% dplyr::filter(Gene1 %in% gene_symbol & Gene2 %in% gene_symbol)
}


## the function to get ReactomeFI function, from DEP2::test_ppi
test_ReactomeFI <- function (x, contrasts = NULL, species = "Human",
                             choose_scores = NULL, reactome_FI_data = reactome_FI_table, score_cutoff = 0.5)
{
  assertthat::assert_that(class(x) == "SummarizedExperiment" |
                            class(x) == "DEGdata" | class(x) == "character", is.null(contrasts) |
                            is.character(contrasts), is.character(species) && length(species) ==
                            1, is.null(choose_scores) | (is.character(choose_scores) &
                                                           all(choose_scores %in% c("combined_score", "neighborhood",
                                                                                    "fusion", "cooccurence", "coexpression", "experimental",
                                                                                    "database", "textmining"))),
                          is.numeric(score_cutoff) | is.integer(score_cutoff),
                          is.null(reactome_FI_data)|is.data.frame(reactome_FI_data),
                          dplyr::between(score_cutoff,0,1))

  if(is.null(reactome_FI_data)) reactome_FI_data <- readRDS("reactome_FI_table.rds")

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
  reactome_FI111 <<- reactome_FI
  if(!is.null(reactome_FI)) reactome_FI = dplyr::filter(reactome_FI,Score >= score_cutoff)

  return(reactome_FI)
}


## Creat the PPI R6 shiny object
creat.ppi.r6 = function(PPI_res, Annotation_enrichment_res, DEP_res){
  theapp.r6 = PPI.enrich.R6.Shiny$new(PPI_res = PPI_res,
                                      ORA_res = Annotation_enrichment_res,
                                      DEP_res = dep_pg)
}



