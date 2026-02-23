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

# Fuente de datos: usar AWS (más estable/rápido que el servidor histórico).
options(
    recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
)

## ---- 01_download --------------------------------------------------------
# 1) Listar proyectos humanos disponibles y seleccionar SRP117733

human_projects <- available_projects(organism = "human")

project_info <- subset(
    human_projects,
    project == "SRP117733" & project_type == "data_sources"
)
stopifnot(nrow(project_info) == 1)

# 2) Descargar el objeto RangedSummarizedExperiment (genes) del proyecto
rse_gene_SRP117733 <- create_rse(project_info)

## (OPCIONAL) guardamos para no redescargar si se cae la sesión
dir.create("Data/rse", recursive = TRUE, showWarnings = FALSE)
saveRDS(rse_gene_SRP117733, "Data/rse/rse_gene_SRP117733_raw.rds")

## ---- 02_counts ----------------------------------------------------------
# recount3 entrega 'raw_counts' (conteos por base / nucleótido).
# compute_read_counts() convierte a conteos por lectura, más apropiados para DE.
assay(rse_gene_SRP117733, "counts") <- compute_read_counts(rse_gene_SRP117733)

## ---- 03_define_groups ---------------------------------------------------
# Definir variables biológicas a partir del texto en sra.sample_attributes.
# Esto permite filtrar SOLO las muestras prepuberales y clasificarlas como SK/control.
txt <- tolower(rse_gene_SRP117733$sra.sample_attributes)

# stage: identifica si la muestra es prepuberal (pre-pubertal / pre pubertal)
stage <- ifelse(grepl("pre[- ]pubertal", txt), "prepuberal", NA_character_)

# group: define condición de interés: SK si menciona klinefelter, control si menciona normal
group <- ifelse(
    grepl("klinefelter", txt),
    "SK",
    ifelse(grepl("\\bnormal\\b", txt), "control", NA_character_)
)

# Verificación: debe existir la celda prepuberal-control=4 y prepuberal-SK=4
print(table(stage, group, useNA = "ifany"))

# keep: vector lógico para quedarnos SOLO con (prepuberal) & (SK/control)
keep <- (stage == "prepuberal") & (group %in% c("SK", "control"))
keep[is.na(keep)] <- FALSE

print(table(keep)) # debe dar 8 TRUE
print(table(group[keep])) # debe dar 4 control y 4 SK

# Subset del RSE a las 8 muestras elegidas
rse_prepub <- rse_gene_SRP117733[, keep]

# Definir group como factor con 'control' como referencia (importante para el contraste)
rse_prepub$group <- factor(group[keep], levels = c("control", "SK"))
print(table(rse_prepub$group))

## (OPCIONAL) guardamos el subset para no tener que repetir este paso si se cae la sesión
saveRDS(rse_prepub, "Data/rse/rse_prepub_SRP117733.rds")


## ---- 04_normalization_voom ------------------------------------------------

# Construir objeto edgeR (DGEList) a partir de la matriz de conteos
dge <- DGEList(counts = assay(rse_prepub, "counts"))

keep_genes <- filterByExpr(dge, group = rse_prepub$group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

# Normalización TMM: corrige diferencias en tamaño de librería/composición
dge <- calcNormFactors(dge)
design <- model.matrix(~group, data = as.data.frame(colData(rse_prepub)))
# voom: transforma a logCPM + estima varianza (necesario para el modelo lineal limma)
v <- voom(dge, plot = FALSE)

## ---- 05_mds_plot ----------------------------------------------------------
# MDS: visualización exploratoria de similitud entre muestras (QC/estructura global)
mds <- plotMDS(v$E, plot = FALSE)

df_mds <- data.frame(
    sample = colnames(rse_prepub),
    group = rse_prepub$group,
    dim1 = mds$x,
    dim2 = mds$y
)
library(ggplot2)
ggplot(df_mds, aes(dim1, dim2, color = group, label = sample)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.6, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("control" = "blue", "SK" = "red")) +
    theme_bw(base_size = 14) +
    labs(
        title = "MDS Plot: Prepubertal SK vs Control",
        x = "Leading logFC dim 1",
        y = "Leading logFC dim 2",
        color = "Group"
    )

## ---- 06_differential_expression ------------------------------------------
# Modelo estadístico: efecto del grupo (SK vs control)
design <- model.matrix(~group, data = as.data.frame(colData(rse_prepub)))
colnames(design)

# Ajuste de modelo lineal y moderación Bayesiana
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extraer tabla completa de resultados para el coeficiente groupSK
de_results <- topTable(
    fit,
    coef = "groupSK",
    number = nrow(v$E),
    sort.by = "P"
)

# Resumen: cuántos genes con FDR < 0.05
table(de_results$adj.P.Val < 0.05)

## ---- 07_plots_basic ------------------------------------------------------
# MA plot y Volcano (rápidos para reporte)
plotMA(fit, coef = "groupSK")
volcanoplot(fit, coef = "groupSK", highlight = 10)

## ---- 08_heatmap_top50 ----------------------------------------------------
# Heatmap con top 50 genes por FDR
top50_ids <- rownames(de_results)[rank(de_results$adj.P.Val) <= 50]
expr_top50 <- v$E[top50_ids, ]

ann <- data.frame(Group = rse_prepub$group)
rownames(ann) <- colnames(expr_top50)

pheatmap(
    expr_top50,
    annotation_col = ann,
    show_colnames = FALSE,
    fontsize_row = 7
)

## ---- 09_save_results -----------------------------------------------------
# Guardar resultados a CSV para tu carpeta Data/
dir.create("Data/results", recursive = TRUE, showWarnings = FALSE)
write.csv(
    de_results,
    "Data/results/DE_limma_voom_prepuberal_SK_vs_control.csv",
    row.names = TRUE
)

# Guardar objetos clave (por reproducibilidad)
saveRDS(v, "Data/results/voom_object_v.rds")
saveRDS(fit, "Data/results/limma_fit.rds")

# definimos la significancia con FDR
sig <- de_results[de_results$adj.P.Val < 0.05, ]
nrow(sig)
""
" [1] 0
> sig
[1] logFC     AveExpr   t         P.Value   adj.P.Val B        
<0 rows> (or 0-length row.names)
"
""
