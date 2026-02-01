#modulo de expresión diferencial (edgeR)

#ui
m_deg_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    fluidRow( #título de la sección
      column(12,
             h2(icon("dna", style = "color: #519657;"), 
                span("Expresión diferencial", style = "color: #519657; font-weight: bold;")),
             hr(style = "border-top: 1px solid #B2F2BB; opacity: 0.5;")
      )
    ),
    
    #resumen rápido
    fluidRow(
      valueBoxOutput(ns("upreg"), width = 4),
      valueBoxOutput(ns("downreg"), width = 4),
      valueBoxOutput(ns("total"), width = 4)
    ),
    
    #config. y resultados
    fluidRow(
      #config
      box(title = tagList(icon("cog", style = "color: #519657;"),
                          "Selección de parámetros:"),
          solidHeader = F,
          status = "success",
          width = 3,
          
          
          #selección del tipo de identificador de los genes
          selectInput(ns("id_type"),
                      "Tipo de identificador de los genes:",
                      choices = c("ENSEMBL ID" = "ENSEMBL",
                                  "ENTREZ ID" = "ENTREZ",
                                  "Gene Symbol" = "SYMBOL"),
                      selected = "ENSEMBL"),
          hr(),
          
          #selección de grupos para el contraste
          uiOutput(ns("ui_sel_contrast")),
          
          hr(),
          
          #selección de parámetros
          numericInput(ns("fdr_co"), "Valor umbral de FDR (p-valor ajustado):",
                       value = 0.05, #por defecto
                       min = 0,
                       max = 1,
                       step = 0.01),
          
          numericInput(ns("logfc_co"), "Valor umbral absoluto de logFC:",
                       value = 1, #por defecto
                       min = 0,
                       step = 0.5),
          hr(),
          actionBttn(ns("run"), "Ejecutar análisis", 
                     icon = icon("play"),
                     style = "jelly",
                     color = "success")),
    
      #resultados
      box(width = 9,
          status = "success",
          tabBox(width = 12,
                 title = span("Resultados del análisis de expresión diferencial mediante edgeR", 
                              style = "color: #519657; font-weight: bold;"),
                 id = ns("tab_res"),
                 tabPanel("Tabla de resultados",
                          withSpinner(DTOutput(ns("tt_table")))),
                 tabPanel("Volcano Plot",
                          withSpinner(plotOutput(ns("vplot"), height = "600px")))
                 ) #tabbox
          ) #box
      
    ) #fluidrow
  ) #taglist
}

#server
m_deg_server <- function(id, i_data) {
  moduleServer(id, function(input, output, session) {
    
    #grupos contraste
    output$ui_sel_contrast <- renderUI({
      req(i_data())
      ns <- session$ns
      
      lvls <- levels(i_data()$dge$samples$group)
      ref <- i_data()$ref_group
      treatment_groups <- setdiff(lvls, ref)
      tagList(
        selectInput(ns("treat_group"), "Grupo Tratamiento:",
                    choices = treatment_groups,
                    selected = treatment_groups[1]),
        p("Modelo de análisis: GLM"),
        p("Se compara el grupo tratamiento con el grupo de referencia:",
          span(class = "label label-success", ref))
      )
    })
    
    #deg
    res <- eventReactive(input$run, {
      req(i_data())
      dge <- i_data()$dge #datos post normalización del módulo qc
      ref <- i_data()$ref_group 
      
      validate(need(!is.null(dge), "Error: no se ha encontrado el objeto DGE."))
      
      #flujo normal de análisis
      withProgress(message = "Ejecutando análisis de expresión diferencial...", value = 0, {
        incProgress(0.2, detail = "Estimando la dispersión...")
        
        design <- model.matrix(~0 + group, data = dge$samples)
        colnames(design) <- levels(dge$samples$group)
        
        dge <- estimateDisp(dge, design = design)
        
        incProgress(0.4, detail = "Ajustando el modelo GLM...")
        fit <- glmFit(dge, design)
        
        incProgress(0.2, detail = "Calculando contrastes...")
        contrast <- paste0(input$treat_group, "-", ref)
        c <- makeContrasts(contrasts = contrast, levels = design)
        
        incProgress(0.3, detail = "Realizando LRT...")
        lrt <- glmLRT(fit, contrast = c)
        
        incProgress(0.1, detail = "Generando tabla de resultados...")
        tt <- topTags(lrt, n = Inf, adjust.method = "BH")$table
        tt$GeneID <- rownames(tt)
      })
      
      return(tt)
    })
    
    #filtrar por significancia
    proc_res <- reactive({
      req(res())
      df <- res()
      
      #definir columna significancia y regulación
      df$Regulation <- "No Sig"
      df$Regulation[df$FDR < input$fdr_co & df$logFC > input$logfc_co] <- "Up"
      df$Regulation[df$FDR < input$fdr_co & df$logFC < -input$logfc_co] <- "Down"
      
      #añadir el tipo de identificador para el enriquecimiento
      df$ID_type <- input$id_type
      
      return(df)
    })
    
    #visualización resultados
    output$upreg <- renderValueBox({
      req(proc_res())
      reg <- proc_res()$Regulation
      valueBox(sum(reg == "Up"), "Genes Up-regulated", icon = icon("arrow-up"),
               color = "red")
    })
    
    output$downreg <- renderValueBox({
      req(proc_res())
      reg <- proc_res()$Regulation
      valueBox(sum(reg == "Down"), "Genes Down-regulated", icon = icon("arrow-down"),
               color = "blue")
    })
    
    output$total <- renderValueBox({
      req(proc_res())
      reg <- proc_res()$Regulation
      valueBox(sum(reg != "No Sig"), "Genes significativos", icon = icon("list"),
               color = "purple")
    })
    
    output$tt_table <- renderDT({
      req(proc_res())
      tt <- proc_res()
      tt <- tt[, c("GeneID", setdiff(names(tt), "GeneID"))]
      datatable(tt,
                rownames = F, 
                extensions = c("Buttons"),
                options = list(
                  dom = "Bfrtip", #orden de las secciones botones, filter...
                  buttons = c("copy", "csv", "excel"),
                  scrollX = T,
                  pageLength = 20)) |> 
        formatSignif(columns = c("logFC", "logCPM", "LR", "PValue", "FDR"), digits = 4) |> 
        formatStyle("Regulation", target = "row", 
                    backgroundColor = styleEqual(c("Up", "Down"), c("#ffcccc", "#cceeff"))) #para que se detecten visualmente
    })
    
    output$vplot <- renderPlot({
      req(res(), input$treat_group)
      data <- res()
      c_title <- attr(data, "contrast_name")
     
      EnhancedVolcano(data, lab = data$GeneID,
                      x = "logFC", y = "FDR",
                      pCutoff = input$fdr_co,
                      FCcutoff = input$logfc_co,
                      title = c_title,
                      pointSize = 2.5,
                      labSize = 3.5,
                      legendPosition = "right",
                      col = c("azure4","springgreen3","indianred3","cornflowerblue"))
    })
    
    return(proc_res)
  })
}