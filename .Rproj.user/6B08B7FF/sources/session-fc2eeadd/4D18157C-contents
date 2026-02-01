#módulo de control de calidad

#ui
m_qc_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow( #título de la sección
      column(12,
             h2(icon("chart-column", style = "color: #C7915B;"), 
                span("Control de calidad", style = "color: #C7915B; font-weight: bold;")),
             hr(style = "border-top: 1px solid #FFD6A5; opacity: 0.5;")
      )
    ),
    
    #config.(elegir columna grupos)
    fluidRow(
      box(title = tagList(icon("cog", style = "color: #C7915B;"),
                          "Configuración"),
          width = 12, 
          solidHeader = F,
          status = "warning",
          
          column(6, 
                 uiOutput(ns("ui_group_sel")), #porque se necesitan los datos del mod. de carga
                 uiOutput(ns("ui_group_ref")),
                 p(class = "text-muted", icon("info-circle"),
                   "La variable de diseño define los grupos y el grupo de referencia es el grupo con el que comparar (ej. control)"),
                 ),
          
          column(6,
                 br(), 
                 actionBttn(ns("run"),
                              "Ejecutar Control de Calidad y Normalización",
                            style = "jelly", color = "warning",
                            icon = icon("play")),
                 uiOutput(ns("genes_f"))
      )
     )
    ), #fluidRow
    
    #resultados (comp. raw y procesados)
    fluidRow(
      #boxplots
      box(title = span("Distribución de los datos (Boxplots)", style = "color: #C7915B; font-weight: bold;"),
          width = 12,
          status = "warning",
          column(6,
                 h4("Datos Crudos (raw)", style = "text-align: center;"),
                 withSpinner(plotOutput(ns("r_boxplot"), height = "500px"))
          ),
          
          column(6,
                 h4("Datos Normalizados", style = "text-align: center;"),
                 withSpinner(plotOutput(ns("n_boxplot"), height = "500px"))
          )
      )
    ), #fluidrow
    
    #PCAs
    fluidRow(
      #boxplots
      box(title = span("Análisis de Componentes Principales (PCA)", style = "color: #C7915B; font-weight: bold;"),
          width = 12,
          status = "warning",
          column(6,
                 h4("Datos Crudos (raw)", style = "text-align: center;"),
                 withSpinner(plotOutput(ns("r_pca"), height = "500px"))
          ),
          
          column(6,
                 h4("Datos Normalizados", style = "text-align: center;"),
                 withSpinner(plotOutput(ns("n_pca"), height = "500px"))
          )
      )
    )
  ) #taglist
}

#server
m_qc_server <- function(id, i_data) { #i_data: datos cargados del mod. de carga
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #seleccionar la columna con los grupos a comparar
    output$ui_group_sel <- renderUI({
      req(i_data())
      meta <- i_data()$meta
      
      selectInput(ns("group_sel"), "Seleccionar la variable de diseño:",
                  choices = colnames(meta),
                  selected = colnames(meta)[1]) #por defecto
    })
    
    #seleccionar grupo de referencia
    output$ui_group_ref <- renderUI({
      req(i_data(), input$group_sel)
      meta <- i_data()$meta
      
      lvls <- unique(as.character(meta[[input$group_sel]]))
      
      selectInput(ns("group_ref"), "Seleccionar el grupo de referencia:",
                  choices = lvls,
                  selected = lvls[1]) #por defecto
    })
    
    #procesamiento y normalización
    proc_results <- eventReactive(input$run, {
      req(i_data(), input$group_sel, input$group_ref)
      
      #datos a procesar
      counts <- i_data()$counts
      meta <- i_data()$meta
      v_group <- input$group_sel #variable que almacena el grupo
      ref <- input$group_ref #nivel de referencia
      
      #procesamiento
      withProgress(message = "Procesando los datos...", value = 0, {
        #1. crear el objeto DGE
        incProgress(0.25, detail = "Creando objeto DGEList...")
        dge <- DGEList(counts = counts, samples = meta)
        
        group_names <- as.factor(make.names(dge$samples[[v_group]])) #eliminar espacios en el nombre
        dge$samples$group <- group_names
        
        ref_name <- make.names(ref)
        tryCatch({
        dge$samples$group <- relevel(dge$samples$group, ref = ref_name)
        }, error = function(e) {
          showNotification("Error al asignar el grupo de referencia",
                           type = "error")
          return(NULL)
        })
        
        #2. filtrar y normalizar
        incProgress(0.5, detail = "Filtrando y normalizando los datos...")
        keep_genes <- filterByExpr(dge, group = dge$samples$group)
        f_dge <- dge[keep_genes, , keep.lib.sizes = F]
        n_dge <- calcNormFactors(f_dge, method = "TMM")
        
        #3. calcular cpm para las gráficas
        incProgress(0.75, detail = "Calculando CPM...")
        r_logcpm <- cpm(dge, log = T)
        n_logcpm <- cpm(n_dge, log = T)
        
        #datos que devolver
        list(dge = n_dge,
             r_lcpm = r_logcpm,
             n_lcpm = n_logcpm,
             meta = n_dge$samples,
             v_group = "group",
             ref_group = ref_name,
             genes = sum(keep_genes)
        )
      }) #withprogress

    })
    
    output$genes_f <- renderUI({
      res <- req(proc_results())
      div(style = "margin-top: 15px; font-size: 15px;",
          span(class = "label label-warning",
               paste("Genes restantes tras el filtrado:", res$genes)))
    })
    
    mycolors <- c("paleturquoise1", "lightsalmon1", "palegreen2",
                  "mistyrose1", "khaki1","lightblue2")
    #gráficos
    #raw
    output$r_boxplot <- renderPlot({
      res <- req(proc_results())
      col_map <- mycolors[as.numeric(res$meta$group)]
      boxplot(res$r_lcpm, ylab="Log2CPM", las=2,
              col = col_map)
    })
    
    output$r_pca <- renderPlot({
      res <- req(proc_results())
      r_pca <- PCAtools::pca(mat = res$r_lcpm, metadata = res$meta, scale = T)
      PCAtools::biplot(r_pca, colby = res$v_group,
                       colkey = mycolors, legendPosition = "bottom")
    })
    
    #norm
    output$n_boxplot <- renderPlot({
      res <- req(proc_results())
      col_map <- mycolors[as.numeric(res$meta$group)]
      boxplot(res$n_lcpm, ylab="Log2CPM", las=2,
              col = col_map)
    })
    
    output$n_pca <- renderPlot({
      res <- req(proc_results())
      n_pca <- PCAtools::pca(mat = res$n_lcpm, metadata = res$meta, scale = T)
      PCAtools::biplot(n_pca, colby = res$v_group,
                       colkey = mycolors, legendPosition = "bottom")
    })
    return(proc_results)
  })
}