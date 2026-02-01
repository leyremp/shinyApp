#módulo de carga de datos (local)

#ui
m_data_upload_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
      fluidRow( #título de la sección
      column(12,
             h2(icon("upload", style = "color: #D47A7A;"), 
                span("Carga de datos", style = "color: #D47A7A; font-weight: bold;")),
             hr(style = "border-top: 1px solid #FFADAD; opacity: 0.5;")
      )
    ),
    
    fluidRow(
      #matriz de conteos
      box(
        title = tagList(icon("table", style = "color: #D47A7A;"),
                        "Matriz de conteos"),
        width = 6,
        solidHeader = F,
        status = "danger",
        
        #fichero de entrada
        fileInput(ns("f_counts"), "Seleccionar archivo (.csv, .txt):",
                  accept = c("text/csv", "text/comma-separated-values",
                             "text/plain", ".csv")),
        #formato fichero
        fluidRow(
          column(6, checkboxInput(ns("h_counts"), "Cabecera", T)),
          column(6, radioButtons(ns("sep_counts"), "Separador:",
                                 choices = c(Coma = ",", Tabulador = "\t",
                                             "Punto y coma" = ";"),
                                 selected = ",", inline = T))
        ),
        checkboxInput(ns("rownames_counts"), 
                      "Nombres de genes como primera columna", T)
      ),
      
      #metadatos
      box(
        title = tagList(icon("list-alt", style = "color: #D47A7A;"),
                        "Metadatos del experimento"),
        width = 6,
        solidHeader = F,
        status = "danger",
        
        #fichero de entrada
        fileInput(ns("f_meta"), "Seleccionar archivo (.csv, .txt)",
                  accept = c("text/csv", "text/comma-separated-values",
                             "text/plain", ".csv")),
        #formato fichero
        fluidRow(
          column(6, checkboxInput(ns("h_meta"), "Cabecera", T)),
          column(6, radioButtons(ns("sep_meta"), "Separador:",
                                 choices = c(Coma = ",", Tabulador = "\t",
                                             "Punto y coma" = ";"),
                                 selected = ",", inline = T))
    ),
    helpText(icon("info-circle"), 
             "Las filas de este fichero deben coincidir con las columnas de la matriz de conteos.")
    )
  ),
  
  #vista previa de los datos subidos
  fluidRow(
    box(
      width = 12,
      status = "danger",
      solidHeader = F,
      tabBox(
        title = span("Comprobación de los datos subidos", style = "color: #D47A7A; font-weight: bold;"),
        width = 12,
        id = ns("exp_tabs"),
        side = "right",
        
        #conteos
        tabPanel("Matriz de conteos", icon = icon("table"),
                 div(style = "overflow-x: scroll;", DTOutput(ns("exp_counts")))
        ),
        #metadatos
        tabPanel("Metadatos del experimento", icon = icon("list-alt"),
                 div(style = "overflow-x: scroll;", DTOutput(ns("exp_meta")))
        ),
        
        #checks
        tabPanel("Validación", icon = icon("check-circle"),
                 div(style = "padding: 30px; text-align: center;",
                     uiOutput(ns("validation_alert"))
                     
                 )
        )
      ) #tabbox
    )
    
   ) #fluidrow
 ) #tagList
}

#server
m_data_upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    #leer ficheros
    #conteos
    r_counts <- reactive({
      req(input$f_counts)
      #validar la lectura y devolver error si falla
      tryCatch({
        df <- read.csv(input$f_counts$datapath,
                       header = input$h_counts,
                       sep = input$sep_counts,
                       row.names = if(input$rownames_counts) 1 else NULL)
        return(df)
      }, error = function(e){
        showNotification(paste("Error en la lectura de la matriz de conteos:",
                               e$message), type = "error")
        return(NULL)
      })
    })
    
    #metadatos
    r_meta <- reactive({
      req(input$f_meta)
      #validar la lectura y devolver error si falla
      tryCatch({
        df <- read.csv(input$f_meta$datapath,
                       header = input$h_meta,
                       sep = input$sep_meta,
                       row.names = 1)
        return(df)
      }, error = function(e){
        showNotification(paste("Error en la lectura de los metadatos:",
                               e$message), type = "error")
        return(NULL)
      })
    })
    
    #validación
    val <- reactive({
      req(r_counts(), r_meta())
      counts <- r_counts()
      meta <- r_meta()
      
      #dimensiones
      if (ncol(counts) != nrow(meta)) {
        return(list(valid = F,
                    message = paste("Las dimensiones no coinciden:\n",
                                    "· Muestras en la matriz de conteos:", ncol(counts),
                                    "\n· Muestras en los metadatos:", nrow(meta))
        ))
      }
      
      #nombres muestras
      comunes <- intersect(colnames(counts), rownames(meta))
      if(length(comunes) != ncol(counts)) {
        return(list(valid = F,
                    message = paste("Los nombres de las muestras no coinciden:\n",
                                    "· Muestras en la matriz de conteos:", 
                                    paste(head(colnames(counts)), collapse = ", "), "\n",
                                    "\n· Muestras en los metadatos:", 
                                    paste(head(rownames(meta)), collapse = ", "))
        ))
      }
      
      return(list(valid = T,
                  message = "Datos correctos.")) #no se enseña
    })
      
    #outputs
    #tablas
    output$exp_counts <- DT::renderDT({
      req(r_counts())
      DT::datatable(head(r_counts(), 50), 
                    options = list(scrollX = T, dom = "t"))
    })  
    
    output$exp_meta <- DT::renderDT({
      req(r_meta())
      DT::datatable(head(r_meta(), 50), 
                    options = list(scrollX = T, dom = "t"))
    })
    
    #validación
    output$validation_alert <- renderUI({
      req(val())
      res <- val()
      if (res$valid) {
        div(class = "alert alert-success", style = "border-radius: 10px;",
            icon("check-circle"), "Datos cargados correctamente.")
      } else {
        div(class = "alert alert-danger",  style = "border-radius: 10px;",
            icon("exclamation-triangle"), "Error en la carga de datos:",
            br(), pre(res$message))
      }
    })
    
    #return si la carga es correcta
    data_checked <- reactive({
      req(val())
      if (val()$valid) {
        list(counts = r_counts(), meta = r_meta())
      } else NULL
    })
    return(data_checked)
  })
}