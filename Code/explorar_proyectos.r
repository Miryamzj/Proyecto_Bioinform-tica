library("recount3")
human_projects <- available_projects()
human_projects
# [1] "SRP000941" "SRP001007" "SRP001008" "SRP001009" "SRP001010" "SRP001011" "SRP001012" "SRP001013" "SRP001014" "SRP001015"

head(human_projects)

#exploramos en un table los dataSets disponibles para humanos

View(human_projects)

project_info <- subset(human_projects, project == "SRP117733")

project_info

rse_gene_SRP117733 <- create_rse(project_info)
rse_gene_SRP117733

# Consultando la información del proyecto en https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
# Se eligió el proyecto SRP117733, que corresponde a un estudio sobre el sindrome de Klinefelter (KS),
# que es una condición genética que afecta a los hombres y se caracteriza por la presencia de un cromosoma X adicional (XXY en lugar de XY).
# En este data Set se encuentran 22 muestras que corresponden a los siguientes grupos:
# Adultos SK: 3
#Adultos controles “matched”: 3
#Solo Sertoli: 4
#Espermatogénesis completa: 4
#Prepuberal SK: 4
#Prepuberal control: 4

#De los cuales sólo se tomarán los Prepuberal SK y los Prepuberal control, quedando un total de 8 muestras para el análisis, para poder comparar
# la expresión génica entre ambos grupos y así identificar posibles diferencias asociadas al síndrome de Klinefelter en la etapa prepuberal.
