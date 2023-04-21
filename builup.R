library(usethis)
library(devtools)

reactome_complex_symbol_gson <- readRDS("inst/ReactomeComplex/reactome_complex_gene_gson.rds")
reactome_complex_symbol_gson@keytype
reactome_complex_uniprot_gson <- readRDS("inst/ReactomeComplex/reactome_complex_pro_gson.rds")
reactome_complex_uniprot_gson@keytype

reactome_complex_participants <- readRDS("inst/ReactomeComplex/participants_list2.rds")

use_data(reactome_complex_symbol_gson,compress = "xz")
use_data(reactome_complex_uniprot_gson,compress = "xz")
use_data(reactome_complex_participants,compress = "xz")

use_r("R6class")
