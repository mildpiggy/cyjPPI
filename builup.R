library(usethis)
library(devtools)

load_all()

reactome_complex_symbol_gson <- readRDS("inst/ReactomeComplex/reactome_complex_gene_gson.rds")
reactome_complex_symbol_gson@keytype
reactome_complex_uniprot_gson <- readRDS("inst/ReactomeComplex/reactome_complex_pro_gson.rds")
reactome_complex_uniprot_gson@keytype

reactome_complex_participants <- readRDS("inst/ReactomeComplex/participants_list2.rds")

use_data(reactome_complex_symbol_gson,compress = "xz")
use_data(reactome_complex_uniprot_gson,compress = "xz")
use_data(reactome_complex_participants,compress = "xz")

use_r("R6class")

load("../cyjShiny test/DEP2_results.rdata")

test_ReactomeFI(dep_pg)

reactome_FI_table <- readRDS("../cyjShiny test/reactome_FI_table.rds")

use_data(reactome_FI_table)


## use USP15 ip-ms resutl as example
library(DEP2)
PG <- read.csv("../USP15 ms_data/USP15 IP质谱/proteinGroups.txt",sep = "\t")
PG2 <- read.csv("../USP15 ms_data/USP15 IP质谱/proteinGroups2.txt",sep = "\t")
PG3 <- read.csv("../USP15 ms_data/USP15 IP质谱/proteinGroups3.txt",sep = "\t")

colnames(PG) %>% grep("10Dox|5Dox",.)
PG = PG[,-grep("10Dox|5Dox",colnames(PG))]

colnames(PG) = gsub("X50Dox|50Dox","Ip",colnames(PG)) %>% gsub("X0Dox|0Dox","Con",.)


identical(colnames(PG),colnames(PG2))
table(colnames(PG) %in% colnames(PG2))
(colnames(PG) %in% colnames(PG3))
colnames(PG)[!(colnames(PG) %in% colnames(PG3))]

uniq = make_unique(PG,"Gene.names", "Protein.IDs")
ecols = grep("LFQ.intensity.",colnames(uniq))
se = make_se_parse(uniq,ecols,"delim",sep = "_")
filt = filter_se(se,thr = 0, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
norm = normalize_vsn(filt)
imp = impute(norm,fun = "MinDet")
colData(imp)$condition
diff = test_diff(imp,type = "control",control = "Con",fdr.type = "BH")
dep = add_rejections(diff,alpha = 0.01,lfc = 2)
plot_heatmap(dep)

get_signicant(dep, return_type = "names") -> vector

# saveRDS(dep,"USP15_Flag.rds")

complex_result = enrichComplex2(vector)

dotplot(complex_result,showCategory = 30)

test_ReactomeFI(dep)

## 内源IP
PG <- read.csv("../USP15 ms_data/USP15 内源IP质谱/proteinGroups.txt",sep = "\t")
uniq = make_unique(PG,"Gene.names", "Protein.IDs")
colnames(PG)
ecols = grep("LFQ.intensity.",colnames(uniq))
se = make_se_parse(uniq,ecols,"delim",sep = "_")
filt = filter_se(se,thr = 0, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
norm = normalize_vsn(filt)
imp = impute(norm,fun = "MinDet")

colData(imp)$condition
diff = test_diff(imp,type = "control",control = "CON",fdr.type = "BH")
dep = add_rejections(diff,alpha = 0.01,lfc = 2)
plot_heatmap(dep)
get_signicant(dep, return_type = "names") -> vector

complex_result = enrichComplex2(vector)
as.data.frame(complex_result)
dotplot(complex_result,showCategory = 30)

test_ReactomeFI(dep) -> FI_res
test_PPI(dep) -> PPI_res

source("R/R6class.R")
FI_res
load_all()
complex_result2 = complex_result %>%
  clusterProfiler.dplyr::filter(qvalue<0.2,pvalue < 0.05)
the_r6 <-
  cyjPPI::PPI.enrich.R6.Shiny$new(PPI_res = FI_res,
                                  ORA_res = complex_result2,
                                  DEP_res = dep)
  # creat_cyjFI(PPI_res = PPI_res,DEP_res = dep, Annotation_enrichment_res = complex_result)
the_r6$run_app()
the_r6$server
shinyApp(ui = the_r6$ui,server = the_r6$server)
