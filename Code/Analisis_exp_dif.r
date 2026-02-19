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

## 1) Cambiar URL de recount3 a AWS
options(
    recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
)

## 2) Cargar proyectos y elegir SRP117733
human_projects <- available_projects(organism = "human")

project_info <- subset(
    human_projects,
    project == "SRP117733" & project_type == "data_sources"
)

stopifnot(nrow(project_info) == 1)

## 3) Descargar RSE (genes)
rse_gene <- create_rse(project_info)

#limpiando los datos para quedarnos sólo con los prepuberal SK y los prepuberal control.

rse_gene_SRP117733$sra.sample_attributes[1:22]

txt <- tolower(rse_gene_SRP117733$sra.sample_attributes)

stage <- ifelse(grepl("pre-pubertal", txt), "prepuberal", NA)
group <- ifelse(
    grepl("klinefelter", txt),
    "SK",
    ifelse(grepl("normal", txt), "control", NA)
)

table(stage, group, useNA = "ifany")


## 4) Crear assay "counts" (cuentas por lectura) a partir de raw_counts
assay(rse_gene, "counts") <- compute_read_counts(rse_gene)
