# Análisis de Expresión Diferencial
## SRP117733 – Muestras prepuberales (SK vs Control)

**Autor:** Miryam Zamora Jimenez  
**Fecha:** 2026-02-20  
**Dataset:** SRP117733 (recount3)  
**Organismo:** Homo sapiens  

---

## 1. Objetivo

Evaluar diferencias de expresión génica entre individuos con síndrome de Klinefelter (SK) y controles en etapa prepuberal, utilizando datos RNA-seq del proyecto SRP117733.

Se empleó un modelo limma-voom tras normalización TMM, con control del FDR para múltiples pruebas.

---

## 2. Diseño experimental

- 4 muestras prepuberales SK
- 4 muestras prepuberales control
- Total: 8 muestras

Se excluyeron muestras adultas para evitar efectos de confusión asociados a espermatogénesis y diferencias de composición celular.

---

## 3. Pipeline bioinformático

1. Descarga de datos mediante `recount3`
2. Conversión de raw_counts a counts por lectura
3. Subset a muestras prepuberales
4. Normalización TMM (edgeR)
5. Transformación voom (limma)
6. Ajuste de modelo lineal
7. Corrección por FDR (Benjamini-Hochberg)

donde:
- control = referencia
- groupSK = efecto SK − control

---

## 4. Resultados Exploratorios

### 4.1 MDS Plot

![MDS](/Plots/MDS_prepub.png)


**Interpretación:**

El análisis MDS no muestra una separación clara entre muestras SK y control en etapa prepuberal. Las muestras de ambos grupos se encuentran parcialmente mezcladas en el espacio bidimensional definido por las dos principales dimensiones de variación. Esto sugiere que el efecto transcriptómico global del síndrome de Klinefelter en etapa prepuberal es limitado en comparación con la variabilidad interindividual.

## 5. MA_plot

![MA](/Plots/MA_plot.png)

**Interpretación:**

El MA plot muestra la relación entre la abundancia promedio de expresión (AveExpr) y el cambio logarítmico de expresión (logFC) entre SK y control.

La mayoría de los genes se distribuyen simétricamente alrededor de logFC = 0, lo que indica ausencia de un desplazamiento global en la expresión génica entre grupos. 

Se observa mayor dispersión en genes de baja expresión, fenómeno esperado en datos RNA-seq debido a mayor variabilidad técnica y biológica en niveles bajos de conteo. En genes con mayor expresión, la variabilidad disminuye, lo que sugiere un comportamiento adecuado del modelo limma-voom.

En conjunto, el MA plot respalda el resultado estadístico de que no se detectan genes diferencialmente expresados con FDR < 0.05 en etapa prepuberal.

---

## 6. Heat_MAP

![Heat_map](/Plots/heat_map.png)

**Interpretación:**

Se generó un mapa de calor utilizando los 50 genes con menor valor ajustado (adj.P.Val). Aunque ninguno alcanzó significancia estadística (FDR < 0.05), estos representan los genes con mayor evidencia relativa de cambio entre grupos.

El heatmap muestra una tendencia de agrupamiento parcial por condición, donde las muestras control y SK tienden a organizarse en clústeres separados. Esto sugiere la existencia de diferencias transcriptómicas sutiles entre los grupos.

Sin embargo, la separación no es completamente marcada, lo que concuerda con la ausencia de genes diferencialmente expresados tras corrección por FDR. 

Estos resultados indican que, en etapa prepuberal, el efecto transcriptómico asociado a SK es leve o requiere mayor tamaño muestral para detectarse con poder estadístico adecuado.
---

## 7. Volcano Plot

![Volcano](/Plots/Volcano_plot.png)

**Interpretación:**

El volcano plot muestra la relación entre el cambio logarítmico de expresión (logFC) y la significancia estadística (-log10 P-value) para todos los genes analizados.

Se observa la distribución característica en forma de “V”, donde genes con mayor magnitud de cambio tienden a presentar valores de P más pequeños. Sin embargo, ninguno de los genes supera el umbral de significancia tras corrección por múltiples pruebas (FDR < 0.05).

La distribución simétrica alrededor de logFC = 0 sugiere ausencia de un desplazamiento global en la expresión génica entre muestras prepuberales con síndrome de Klinefelter y controles.

Estos resultados son consistentes con el MA plot y el heatmap, indicando que las diferencias transcriptómicas en esta etapa son sutiles o requieren mayor tamaño muestral para detectarse con suficiente poder estadístico.

---

## 8. Genes Diferencialmente Expresados

## 6. Genes diferencialmente expresados

La identificación de genes diferencialmente expresados se realizó utilizando el modelo lineal ajustado con limma sobre los datos transformados mediante voom. Para cada gen se estimó el cambio de expresión (logFC), el estadístico t moderado y el valor de p correspondiente.

Dado que se evaluaron miles de genes simultáneamente, se aplicó una corrección por múltiples pruebas utilizando el método de Benjamini-Hochberg para controlar la tasa de falsos descubrimientos (FDR). Se consideró como umbral de significancia un FDR < 0.05.

Tras aplicar este criterio, no se identificaron genes diferencialmente expresados entre muestras prepuberales con síndrome de Klinefelter y controles:

Número de genes con FDR < 0.05: 0

Aunque algunos genes presentaron valores de p relativamente bajos y cambios de expresión moderados (logFC elevados), estos no permanecieron significativos después de la corrección por múltiples pruebas.

Este resultado sugiere que, bajo las condiciones analizadas y con el tamaño muestral disponible (n = 4 por grupo), no se detectan alteraciones transcriptómicas robustas en etapa prepuberal asociadas al síndrome de Klinefelter.
La ausencia de genes significativos no implica ausencia de diferencias biológicas, sino que indica que las diferencias observadas no alcanzan evidencia estadística suficiente bajo los criterios establecidos.

## 9. Conclusiones

En este estudio se realizó un análisis de expresión diferencial utilizando datos del proyecto SRP117733 del repositorio recount3, enfocándose en muestras prepuberales con síndrome de Klinefelter (SK) y controles sanos.

El análisis incluyó:

- Transformación voom para modelar la relación media-varianza en datos de RNA-seq.
- Ajuste de modelos lineales mediante limma.
- Moderación Bayesiana de varianzas.
- Corrección por múltiples pruebas utilizando el método Benjamini-Hochberg (FDR).

Los análisis exploratorios (MDS, MA plot y heatmap) sugieren la presencia de diferencias transcriptómicas sutiles entre SK y controles. Sin embargo, tras la corrección por FDR, no se identificaron genes diferencialmente expresados con un umbral de significancia de 0.05.

La ausencia de genes significativos puede explicarse por:

1. Tamaño muestral reducido (n = 4 por grupo).
2. Efectos transcriptómicos leves en etapa prepuberal.
3. Variabilidad biológica interindividual.

En conjunto, los resultados indican que en etapa prepuberal el impacto transcriptómico del síndrome de Klinefelter es limitado o requiere mayor poder estadístico para ser detectado de forma robusta.

Futuros análisis podrían incluir:
- Mayor tamaño de muestra.
- Modelos ajustados por covariables adicionales.
- Análisis enfocados en cromosoma X.
- Enfoques de enriquecimiento funcional.

Este análisis demuestra la correcta implementación de un pipeline estándar de RNA-seq basado en limma-voom y destaca la importancia del control de múltiples pruebas en estudios de expresión génica.