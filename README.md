# RSeqXplorer

El análisis de expresión diferencial de genes es un proceso clave para poder identificar cambios en la expresión génica entre diferentes condiciones biológicas (por ejemplo entre controles vs casos o distintos tratamientos). Sin embargo, su correcta ejecución suele requerir de conocimientos básicos de programación, pudiendo suponer una barrera para parte de la comunidad científica. RSeqXplorer tiene como objetivo facilitar este tipo de análisis mediante una interfaz accesible e intuitiva.

--- 

## Pipeline principal 

El pipeline de análisis incluye los siguientes pasos:
1. **Carga de datos locales**
    - Matriz de conteos crudos
    - Metadatos del experimento 
    - Acepta tanto ficheros en formato `.csv` como `.txt` con distintos separadores
    - Soporte para datos de humano (*Homo sapiens*) y ratón (*Mus musculus*)

2. **Control de calidad**
    - Selección de la variable diseño y el grupo de referencia para los siguientes pasos
    - Filtrado de conteos bajos y normalización TMM mediante `edgeR`
    - Visualización pre y post filtrado de los datos mediante Boxplots y PCA comparativos

3. **Análisis de la expresión diferencial con `edgeR`**
    - Selección del grupo tratamiento
    - Selección de los parámetros del análisis
    - Modelado por GLM 
    - Visualización de los genes significativos (tabla y volcano plot)

4. **Análisis de enriquecimiento funcional con `clusterProfiler`**
    - Selección de parámetros para el análisis
    - Selección del método de enriquecimiento a realizar:
        - ORA (Over Representation Analysis)
        - GSEA (Gene Set Enrichment Analysis)
        - Tanto de términos GO como de rutas KEGG

---

## Getting started

### Paquetes necesarios
La aplicación utiliza los siguientes paquetes de R para su correcta ejecución. Para instalar automáticamente los faltantes se puede ejecutar el siguiente código:

```r
required_p <- c("shiny", "shinydashboard", "shinyWidgets", "shinycssloaders",
                "edgeR", "ggplot2", "DT", "PCAtools", "dplyr", 
                "EnhancedVolcano", "clusterProfiler", "AnnotationDbi",
                "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "stats")

installed_p <- rownames(installed.packages())
to_install <- setdiff(required_p, istalled_p)
```

### Iniciar la aplicación
Una vez instalados los paquetes, la aplicación se puede utilizar de forma local siguiendo los siguientes pasos:

1. Abrir RStudio
2. Abrir el fichero `app.R`
3. Hacer click en `Run App`

o 

```r
shiny::runApp("app.R")
```

---

## Ejemplo

En la carpeta [data](https://github.com/leyremp/shinyApp/tree/main/data) se adjuntan dos ficheros de ejemplo, uno de conteos (ex_counts.txt) y otro de metadatos (ex_metadata.txt) con los que se puede desarrollar el pipeline al completo. 

En la sección de **Carga de datos** hay que especificar los separadores para cada fichero, para conteos `Tabulador` y para los metadatos `Punto y coma`. Para comprobar si se han cargado correctamente, las dos tablas se pueden visualizar junto con un mensaje de comprobación, que si aparece en verde indica que están correctamente cargados, pero en caso contrario, el mensaje aparece en rojo con el error detectado. 

Con los datos cargados, se puede proceder con el resto del pipeline.

---

## Notas

- Los datos de entrada deben estar correctamente formateados para evitar errores durante el análisis.

- En la sección **Expresión diferencial** se debe seleccionar correctamente el tipo de identificador de los genes (ENSEMBL, SYMBOL, ENTREZ) para poder realizar el enriquecimiento funcional.

- El tiempo de análisis puede variar en función del tamaño del conjunto de los datos.

