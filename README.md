# Proyecto_Bioinform-tica
Proyecto de expresión diferencial, usando las herramientas de Rcount3 
# Análisis de Expresión Diferencial – SRP117733 (Síndrome de Klinefelter)

## Descripción del Proyecto

Este proyecto realiza un análisis de expresión diferencial utilizando datos de RNA-seq provenientes del repositorio **recount3**. El dataset corresponde al estudio **SRP117733**, el cual investiga patrones de expresión génica asociados al **síndrome de Klinefelter (SK)**.

El análisis se enfoca exclusivamente en muestras **prepuberales**, comparando:

- 4 muestras Control  
- 4 muestras con Síndrome de Klinefelter (SK)

El objetivo principal fue evaluar si existen diferencias transcriptómicas robustas entre ambos grupos en esta etapa del desarrollo.

---

## Pipeline de Análisis

El flujo de trabajo sigue un pipeline estándar para análisis de RNA-seq basado en **limma-voom**:

1. Descarga de datos desde recount3  
2. Conversión de conteos por base a conteos por lectura  
3. Filtrado de muestras (solo etapa prepuberal)  
4. Filtrado de genes con baja expresión (`filterByExpr`)  
5. Normalización TMM (edgeR)  
6. Modelado media-varianza con **voom**  
7. Ajuste de modelo lineal con **limma**  
8. Moderación Bayesiana de varianzas  
9. Corrección por múltiples pruebas mediante **Benjamini-Hochberg (FDR)**  

---

## Resultados Principales

- No se identificaron genes diferencialmente expresados con **FDR < 0.05**.
- Los análisis exploratorios (MDS, MA plot, Volcano plot y Heatmap) sugieren diferencias transcriptómicas sutiles entre grupos.
- La ausencia de genes significativos puede deberse a:
  - Tamaño muestral reducido (n = 4 por grupo)
  - Efectos biológicos leves en etapa prepuberal
  - Variabilidad interindividual

---

## Estructura del Repositorio
Code/  Scripts de análisis
Data/  Objetos RSE descargados y tablas de resultados
Plots/  Figuras generadas (MDS, MA, Volcano, Heatmap)
Reports/  Reporte en formato Markdown


---

## Figuras Generadas

- MDS plot (agrupamiento de muestras)
- MA plot (abundancia promedio vs log fold change)
- Volcano plot (logFC vs -log10 p-value)
- Heatmap de los 50 genes con menor FDR

---

## Software y Paquetes Utilizados

- R (≥ 4.5.1)
- recount3
- edgeR
- limma
- ggplot2
- pheatmap
- RColorBrewer

---

## Conclusión

Bajo las condiciones analizadas, no se detectaron diferencias transcriptómicas estadísticamente significativas entre muestras prepuberales con síndrome de Klinefelter y controles tras corrección por múltiples pruebas.

Este proyecto demuestra la implementación correcta de un pipeline reproducible de RNA-seq utilizando limma-voom y resalta la importancia del control del FDR y del tamaño muestral en estudios transcriptómicos.

---

## Autora

Miryam Zamora Jimenez  
2026