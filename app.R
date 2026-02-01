#app.R v1
source("global.R")

#ui
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = a_title),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("Home", tabName = "t_home", icon = icon("home")), 
      menuItem("Carga de datos", tabName = "t_upload", icon = icon("upload")),
      menuItem("Control de calidad", tabName = "t_qc", icon = icon("chart-column")),
      menuItem("Expresión diferencial", tabName = "t_deg", icon = icon("dna")),
      menuItem("Enriquecimiento funcional", tabName = "t_enrich", icon = icon("project-diagram"))
    ),
    
    #info aplicación
    tags$style(HTML("
    .sidebar-footer { position: absolute; bottom: 0; width: 100%; padding: 10px; }
  ")),
    
    tags$div(
      class = "sidebar-footer",
      style = "padding: 10px 15px; font-size: 12px; color: #b8c7ce",
      tags$small(
        tags$div(
          tags$strong(a_title),
          paste0(" v ", version)
        ),
        tags$em(
          tags$div(
            "Desarrollado por:", tags$strong(dev)
          ),
          tags$div(
            paste0(subj, " © ", year)
          ),
          tags$div(
            paste0("GitHub: ", github) #cambiar a link cuando se tenga
          )
        )
      )
    )
  ), #sidebar
  
  #pestaña home
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "t_home",
        #encabezado de la app
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = T,
            title = NULL,
            h1(a_title, 
               style = "text-align: center; margin: 20px 0;"), #centrar el texto
            h4("Plataforma integral para el análisis de datos de RNA-Seq", #subtítulo
               style = "text-align: center; color: #555; margin: 20px 0;")
          )
        ),
        
        #info
        fluidRow(
          box(
            width = 6,
            status = "info",
            solidHeader = F,
            title = tagList(icon("info-circle"), "¿Para qué sirve RNA-Seq Pipeline?"),
            p("Esta aplicación realiza un análisis completo de datos de RNA-seq."),
            p("Desde la carga de datos hasta el enriquecimiento funcional de los genes
            diferencialmente expresados."),
            p("Está diseñada para el público investigador que quiera explorar de forma 
          sencilla los resultados del análisis.")
          ),
          
          box(
            width = 6,
            status = "info",
            solidHeader = F,
            title = tagList(icon("check-double"), "Características principales"),
            tags$ul(
              tags$li("Flexibilidad: Acepta experimentos de humano y ratón"),
              tags$li("Control de calidad: Filtrado de los conteos bajos y normalización mediante TMM"),
              tags$li("Análisis de la expresión diferencial: Modelado por GLM"),
              tags$li("Análisis de enriquecimiento funcional: GO (Gene Ontology) y KEGG pathways")
            )
          )
        ), #fluidrow
        
        #workflow
        fluidRow(
          box(
            width = 12,
            status = "primary",
            solidHeader = F,
            title = tagList(icon("route"), "Workflow"),
            
            fluidRow(
              column(6,
                     h4(icon("upload", style = "color: #D47A7A;"),
                        "Paso 1. Carga de datos", style = "font-weight: bold; color: #D47A7A;"),
                     p("Importa la matriz de conteos y el archivo de metadatos del experimento."),
                     br()
              ),
              column(6,
                     h4(icon("chart-column", style = "color: #C7915B;"),
                        "Paso 2. Control de calidad", style = "font-weight: bold; color: #C7915B;"),
                     p("Visualiza la distribución de los datos antes y después de filtrar y normalizar."),
                     br()
              )
            ),
            
            fluidRow(
              column(6,
                     h4(icon("dna", style = "color: #519657;"),
                        "Paso 3. Expresión diferecial", style = "font-weight: bold; color: #519657;"),
                     p("Identifica los genes diferencialmente expresados entre las condiciones experimentales."),
                     br()
              ),
              column(6,
                     h4(icon("project-diagram", style = "color: #4682B4;"),
                        "Paso 4. Enriquecimiento funcional", style = "font-weight: bold; color: #4682B4;"),
                     p("Identifica las rutas y procesos enriquecidos."),
                     br()
              )
            )
          )
        )
      ),#fin home
      
      #resto de pestañas
      tabItem(tabName = "t_upload", m_data_upload_ui("m_upload_1")),
      tabItem(tabName = "t_qc", m_qc_ui("m_qc_1")),
      tabItem(tabName = "t_deg", m_deg_ui("m_deg_1")),
      tabItem(tabName = "t_enrich", m_enrichment_ui("m_enrich_1"))
    )#tabitems
  )#dashboardbody
)

server <- function(input, output, session) {
 
  #cargar datos
  uploaded_d <- m_data_upload_server("m_upload_1")
  
  #qc (necesita uploaded_d, devuelve qc_res)
  qc_res <- m_qc_server("m_qc_1", i_data = uploaded_d)
  
  #deg
  deg_res <- m_deg_server("m_deg_1", i_data = qc_res)
  
  #enrichment
  enrich_res <- m_enrichment_server("m_enrich_1", i_data = deg_res)
}

shinyApp(ui = ui, server = server)