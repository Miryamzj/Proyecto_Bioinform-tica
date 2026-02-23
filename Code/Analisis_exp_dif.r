# Proyecto de análisis de expresión diferencial utilizando el data set SRP117733 del repositorio recount3, que corresponde a un estudio sobre el síndrome de Klinefelter (KS) en la etapa prepuberal.

#Autor: Miryam Zamora Jimenez
#Fecha: 2026-02-20

# Cargamos la librería recount3 para acceder a los datos de expresión génica disponibles en el repositorio.
library("recount3") # Cargamos la librería recount3 para acceder a los datos de expresión génica disponibles en el repositorio.
library("limma") # Cargamos la librería limma para realizar el análisis de expresión diferencial utilizando el método de modelado lineal.
library("edgeR") # Cargamos la librería edgeR para realizar el análisis de expresión diferencial utilizando el método de modelado de conteos.
library("ggplot2") # Cargamos la librería ggplot2 para realizar visualizaciones de los resultados del análisis de expresión diferencial.
library("pheatmap") # Cargamos la librería pheatmap para realizar visualizaciones de mapas de calor de los resultados del análisis de expresión diferencial.
library("RColorBrewer") # Cargamos la librería RColorBrewer para utilizar paletas de colores en las visualizaciones de los resultados del análisis de expresión diferencial.

options(
    recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
)

## ---- 01_download --------------------------------------------------------
human_projects <- available_projects(organism = "human")

project_info <- subset(
    human_projects,
    project == "SRP117733" & project_type == "data_sources"
)
stopifnot(nrow(project_info) == 1)

rse_gene_SRP117733 <- create_rse(project_info)

## (OPCIONAL) guarda para no redescargar si se cae la sesión
dir.create("Data/rse", recursive = TRUE, showWarnings = FALSE)
saveRDS(rse_gene_SRP117733, "Data/rse/rse_gene_SRP117733_raw.rds")

## ---- 02_counts ----------------------------------------------------------
assay(rse_gene_SRP117733, "counts") <- compute_read_counts(rse_gene_SRP117733)

## ---- 03_define_groups ---------------------------------------------------
txt <- tolower(rse_gene_SRP117733$sra.sample_attributes)

stage <- ifelse(grepl("pre[- ]pubertal", txt), "prepuberal", NA_character_)
group <- ifelse(
    grepl("klinefelter", txt),
    "SK",
    ifelse(grepl("\\bnormal\\b", txt), "control", NA_character_)
)

print(table(stage, group, useNA = "ifany"))

keep <- (stage == "prepuberal") & (group %in% c("SK", "control"))
keep[is.na(keep)] <- FALSE

print(table(keep)) # debe dar 8 TRUE
print(table(group[keep])) # debe dar 4 control y 4 SK

rse_prepub <- rse_gene_SRP117733[, keep]
rse_prepub$group <- factor(group[keep], levels = c("control", "SK"))
print(table(rse_prepub$group))

## (OPCIONAL) guarda el subset
saveRDS(rse_prepub, "Data/rse/rse_prepub_SRP117733.rds")
