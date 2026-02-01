#global config file

#librerías
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(edgeR)
library(ggplot2)
library(DT)
library(PCAtools)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(stats)


#aumentar el tamaño permitido
options(shiny.maxRequestSize=30*1024^2)

#funciones auxiliares
source("functions/ids.R")

#módulos
source("modules/m_data_upload.R")
source("modules/m_qc.R")
source("modules/m_deg.R")
source("modules/m_enrichment.R")

#metadatos app
a_title <- "RNA-Seq Pipeline" #no definitivo
version <- "1.0"
dev <- "Leyre Melià"
subj <- "Estudios in silico en biomedicina"
github <- "GitHub" #link cuando se tenga
year <- format(Sys.Date(), "%Y")
