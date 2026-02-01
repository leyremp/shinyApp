#modulo de enriquecimiento funcional (clusterProfiler)

#ui
m_enrichment_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    fluidRow( #título de la sección
      column(12,
             h2(icon("project-diagram", style = "color: #4682B4;"), 
                span("Enriquecimiento funcional", style = "color: #4682B4; font-weight: bold;")),
             hr(style = "border-top: 1px solid #A0C4FF; opacity: 0.5;")
      )
    ),
    
    #Config
    fluidRow(
      box(
        title = tagList(icon("cog", style = "color: #4682B4;"),
                        "Configuración"),
        solidHeader = F,
        status = "info",
        width = 3,
        
        tags$h5(
          span("Configuración Básica", style = "font-weight: bold;"),
          style = "border-bottom: 1px solid #e5e5e5; padding-bottom: 6px;"),
        
        #Selección de especie
        selectInput(ns("specie"),
                    "Especie:",
                    choices = c("Humano" = "human",
                                "Ratón" = "mouse"),
                    selected = "human"),
        
        tags$h5(
          span("Métodos de enriquecimiento", style = "font-weight: bold;"),
          style = "border-bottom: 1px solid #e5e5e5; padding-bottom: 6px;"),
        
        p("Selecciona uno o más métodos:"),
        
        checkboxGroupInput(ns("method"),
                           NULL,
                           choices = c(
                             "ORA - Gene Ontology (GO)" = "ora_go",
                             "ORA - KEGG Pathways" = "ora_kegg",
                             "GSEA - Gene Ontology (GO)" = "gsea_go",
                             "GSEA - KEGG Pathways" = "gsea_kegg"
                           )),
        
        tags$h5(
          span("Parámetros del análisis", style = "font-weight: bold;"),
          style = "border-bottom: 1px solid #e5e5e5; padding-bottom: 6px;"),              
        
        numericInput(ns("pvalue_co"),
                     "Valor umbral del p-valor:",
                     value = 0.05, #por defecto
                     min = 0,
                     max = 1, 
                     step = 0.01),
        
        numericInput(ns("minGSSize"),
                     "Valor mínimo del tamaño del grupo:",
                     value = 10, #por defecto
                     min = 5,
                     max = 100, 
                     step = 5),
        
        numericInput(ns("maxGSSize"),
                     "Valor máximo del tamaño del grupo:",
                     value = 100, #por defecto
                     min = 50,
                     max = 1000, 
                     step = 50),
        
        conditionalPanel(
          condition = "input.method && (input.method.includes('ora_go') || input.method.includes('gsea_go'))",
          ns = ns,
          
          tags$h5(
            span("Ontología GO", style = "font-weight: bold;"),
            style = "border-bottom: 1px solid #e5e5e5; padding-bottom: 6px;"),
          
          selectInput(ns("ont"),
                      "Selecciona una ontología:",
                      choices = c("Biological Process" = "BP",
                                  "Molecular Function" = "MF",
                                  "Cellular Component" = "CC"),
                      selected = "BP")
        ),
        
        
        hr(),
        actionBttn(ns("run"),
                   "Ejecutar análisis de enriquecimiento",
                   icon = icon("play"), style = "jelly",
                   color = "primary")
      ), #box
      
      #resultados
      box(
        width = 9,
        status = "info",
        solidHeader = F,
        title =  span("Resultados del análisis de enriquecimiento", style = "color: #4682B4; font-weight: bold;"),
        withSpinner(uiOutput(ns("res_vis")))
      ) #box
    ) #fluidrow
  )    
}

m_enrichment_server <- function(id, i_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #detectar IDs
    ids <- reactive({
      req(i_data())
      tt <- i_data()
      if ("ID_type" %in% colnames(tt)) {
        return(tt$ID_type[1])
      } else {
        return(detectIDs(tt$GeneID, specie = input$specie))
      }
    })
    
    #bases de datos/código kegg de organismo
    organism_db <- reactive({
      if (input$specie == "human") {
        return(org.Hs.eg.db)
      } else {
        return(org.Mm.eg.db)
      }
    })
    
    organism_kegg <- reactive({
      if (input$specie == "human") {
        return("hsa")
      } else {
        return("mmu")
      }
    })
    
    data <- eventReactive(input$run, {
      req(i_data(), ids())
      tt <- i_data()
      ori_id <- ids()
      
      validate(need(!is.null(tt), "Error. No hay datos de expresión diferencial."))
      validate(need("GeneID" %in% colnames(tt), "Error. No existe la columna GeneID"))
      
      #definir grupos de genes: universo, up y downregulados
      uni_genes <- unique(tt$GeneID)
      up_genes <- unique(tt$GeneID[tt$Regulation == "Up"])
      down_genes <- unique(tt$GeneID[tt$Regulation == "Down"])
      
      validate(need(length(up_genes) > 0 | length(down_genes) > 0,
                    "No hay genes diferencialmente expresados."))
      
      #convertir a ENSEMBL (si no lo son) y a ENTREZID
      withProgress(message = "Convirtiendo los identificadores...", value = 0.2, {
        if(ori_id != "ENSEMBL") {
          ensembl_ids <- convertIDs(uni_genes,
                                    ori_id = ori_id,
                                    new_id = "ENSEMBL",
                                    specie = input$specie)
            
          uni_ens <- ensembl_ids$ENSEMBL
          up_ens <- ensembl_ids$ENSEMBL[match(up_genes,ensembl_ids[[ori_id]])]
          up_ens <- unique(up_ens[!is.na(up_ens)])
          down_ens <- ensembl_ids$ENSEMBL[match(down_genes,ensembl_ids[[ori_id]])]
          down_ens <- unique(down_ens[!is.na(down_ens)])
        } else {
          uni_ens <- uni_genes
          up_ens <- up_genes
          down_ens <- down_genes
        }
        entrez_ids <- convertIDs(uni_genes,
                                 ori_id = ori_id,
                                 new_id = "ENTREZID",
                                 specie = input$specie)
        
        uni_entrez <- entrez_ids$ENTREZID
        up_entrez <- entrez_ids$ENTREZID[match(up_genes,entrez_ids[[ori_id]])]
        up_entrez <- unique(up_entrez[!is.na(up_entrez)])
        down_entrez <- entrez_ids$ENTREZID[match(down_genes,entrez_ids[[ori_id]])]
        down_entrez <- unique(down_entrez[!is.na(down_entrez)])
        
        #preparar los datos para GSEA GO
        geneList <- tt$logFC
        names(geneList) <- tt$GeneID

        if (ori_id != "ENSEMBL") {
          names(geneList) <- ensembl_ids$ENSEMBL[match(names(geneList), ensembl_ids[[ori_id]])]
        }
        geneList <- geneList[!is.na(names(geneList))]
        geneList <- sort(geneList, decreasing = T)
        
        #preparar los datos para GSEA KEGG
        geneList_kegg <- tt$logFC
        names(geneList_kegg) <- tt$GeneID
        if (ori_id != "ENTREZID") {
          names(geneList_kegg) <- entrez_ids$ENTREZID[match(names(geneList_kegg), 
                                                       entrez_ids[[ori_id]])]
        }
        geneList_kegg <- geneList_kegg[!is.na(names(geneList_kegg))]
        geneList_kegg <- sort(geneList_kegg, decreasing = T)
      })
      
      #devolver datos
      list(uni_ens = uni_ens,
           up_ens = up_ens,
           down_ens = down_ens,
           uni_entrez = uni_entrez,
           up_entrez = up_entrez,
           down_entrez = down_entrez,
           geneList = geneList,
           geneList_kegg = geneList_kegg)
    })
    
    #ORA - GO (si se selecciona)
    ora_go <- reactive({
      req(data())
      req("ora_go" %in% input$method)
      data <- data()
      
      withProgress(
        message = "Ejecutando Over-representation analysis con GO...",
        value = 0.2,
        {
          #ora para genes upregulados
          oraUp <- NULL
          if (length(data$up_ens) > 0) {
            oraUp <- tryCatch({
              enrichGO(gene = data$up_ens,
                       universe = data$uni_ens,
                       OrgDb = organism_db(),
                       keyType = "ENSEMBL",
                       ont = input$ont,
                       pAdjustMethod = "BH",
                       pvalueCutoff = input$pvalue_co,
                       minGSSize = input$minGSSize,
                       maxGSSize = input$maxGSSize)
            }, error = function(e) {
              showNotification(paste("Error en la ejecución de ORA GO (Up):",
                                     e$message), type = "error")
              NULL
            })
          }
            #ora para genes downregulados
            oraDown <- NULL
            if (length(data$down_ens) > 0) {
              oraDown <- tryCatch({
                enrichGO(gene = data$down_ens,
                         universe = data$uni_ens,
                         OrgDb = organism_db(),
                         keyType = "ENSEMBL",
                         ont = input$ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff = input$pvalue_co,
                         minGSSize = input$minGSSize,
                         maxGSSize = input$maxGSSize)
              }, error = function(e) {
                showNotification(paste("Error en la ejecución de ORA GO (Down):",
                                       e$message), type = "error")
                NULL
              })
          }
      }) #withprogress
      
      list(up = oraUp, down = oraDown)
    })
    
    #ORA - KEGG (si se selecciona)
    ora_kegg <- reactive({
      req(data())
      req("ora_kegg" %in% input$method)
      data <- data()
      
      withProgress(
        message = "Ejecutando Over-representation analysis con KEGG...",
        value = 0.2,
        {
          #ora para genes upregulados
          oraUpK <- NULL
          if (length(data$up_entrez) > 0) {
            oraUpK <- tryCatch({
              enrichKEGG(gene = data$up_entrez,
                         universe = data$uni_entrez,
                         organism = organism_kegg(),
                         keyType = "ncbi-geneid",
                         pvalueCutoff = input$pvalue_co)
            }, error = function(e) {
              showNotification(paste("Error en la ejecución de ORA KEGG (Up):",
                                     e$message), type = "error")
              NULL
            })
          }
          #ora para genes downregulados
          oraDownK <- NULL
          if (length(data$down_entrez) > 0) {
            oraDownK <- tryCatch({
              enrichKEGG(gene = data$down_entrez,
                         universe = data$uni_entrez,
                         organism = organism_kegg(),
                         keyType = "ncbi-geneid",
                         pvalueCutoff = input$pvalue_co)
            }, error = function(e) {
              showNotification(paste("Error en la ejecución de ORA KEGG (Down):",
                                     e$message), type = "error")
              NULL
            })
          }
        }) #withprogress
      
      list(up = oraUpK, down = oraDownK)
    })
    
    #GSEA - GO (si selecciona)
    gsea_go <- reactive({
      req(data())
      req("gsea_go" %in% input$method)
      data <- data()
      
      withProgress(
        message = "Ejecutando Gene Set Enrichment Analysis (GSEA) con GO...",
        value = 0.5,
        {
          gseaGO <- tryCatch({
            gseGO(geneList = data$geneList,
                  ont = input$ont,
                  OrgDb = organism_db(),
                  keyType = "ENSEMBL",
                  minGSSize = input$minGSSize,
                  maxGSSize = input$maxGSSize,
                  pAdjustMethod = "BH",
                  pvalueCutoff = input$pvalue_co,
                  verbose = T)
          }, error = function(e) {
            showNotification(paste("Error en GSEA (GO):",
                                   e$message), type = "error")
            NULL
          })
        })
      
      gseaGO
    })  
     
      #GSEA - KEGG (si se selecciona)
      gsea_kegg <- reactive({
        req(data())
        req("gsea_kegg" %in% input$method)
        data <- data()
        
        withProgress(
          message = "Ejecutando Gene Set Enrichment Analysis (GSEA) con KEGG...",
          value = 0.5,
          {
            gseaKEGG <- tryCatch({
              gseKEGG(geneList = data$geneList_kegg,
                      organism = organism_kegg(),
                      keyType = "ncbi-geneid",
                      pvalueCutoff = input$pvalue_co,
                      verbose = T)
             
            }, error = function(e) {
              showNotification(paste("Error en GSEA (KEGG):",
                                     e$message), type = "error")
              NULL
            })
          })
        gseaKEGG
      })
      
      #outputs
      output$res_vis <- renderUI({
        req(input$method)
        method <- input$method
        
        ont_name <- c("BP" = "Biological Process",
                      "MF" = "Molecular Function",
                      "CC" = "Cellular Component")
        ont_label <- ont_name[input$ont] #para los outputs
        
        tabs <- list()
        
        #ORA - GO
        if ("ora_go" %in% method) {
          tabs <- append(tabs, list(
            tabPanel(paste("Over Representation Analysis (ORA) - GO: ", ont_label),
            br(),
            tabsetPanel(
              tabPanel("Up-regulated",
                       h4("Análisis funcional: genes Up-regulados"),
                       plotOutput(ns("ora_go_up_p"), height = "500px"),
                       DTOutput(ns("ora_go_up_t"))),
              
              tabPanel("Down-regulated",
                       h4("Análisis funcional: genes Down-regulados"),
                       plotOutput(ns("ora_go_down_p"), height = "500px"),
                       DTOutput(ns("ora_go_down_t")))
              )
            )
          ))
        }
          
        #ORA - KEGG
        if ("ora_kegg" %in% method) {
          tabs <- append(tabs, list(
            tabPanel("Over Representation Analysis (ORA) - KEGG",
            br(),
            tabsetPanel(
              tabPanel("Up-regulated",
                       h4("Rutas: Genes Up-regulados"),
                       plotOutput(ns("ora_kegg_up_p"), height = "500px"),
                       DTOutput(ns("ora_kegg_up_t"))),
              
              tabPanel("Down-regulated",
                       h4("Rutas: Genes Down-regulados"),
                       plotOutput(ns("ora_kegg_down_p"), height = "500px"),
                       DTOutput(ns("ora_kegg_down_t")))
              )
            )
          ))
        } 
        
        #GSEA - GO
        if ("gsea_go" %in% method) {
          tabs <- append(tabs, list(
            tabPanel(paste("Gene Set Enrichment Analysis (GSEA) - GO: ", ont_label),
                     br(),
                     h4("Distribución de enriquecimiento"),
                     plotOutput(ns("gsea_go_p"), height = "500px"),
                     DTOutput(ns("gsea_go_t"))
            )
          ))
        }
        
        #GSEA -KEGG
        if ("gsea_kegg" %in% method) {
          tabs <- append(tabs, list(
            tabPanel("Gene Set Enrichment Analysis (GSEA) - KEGG",
                     br(),
                     h4("Distribución de enriquecimiento"),
                     plotOutput(ns("gsea_kegg_p"), height = "500px"),
                     DTOutput(ns("gsea_kegg_t"))
            )
          ))
        }
          
        if (length(tabs) > 0) {
          do.call(tabsetPanel, tabs) #enseñar las tabs seleccionadas
        }  else {
          p("Selecciona un método y pulsa 'Ejecutar análisis de enriquecimiento'.")
        }
      })
      
      #ORA - GO
      #gráficos
      output$ora_go_up_p <- renderPlot({
        req(ora_go())
        og <- ora_go()
        ont_name <- c("BP" = "Biological Process",
                      "MF" = "Molecular Function",
                      "CC" = "Cellular Component")
        ont_label <- ont_name[input$ont] #para el título
        
        if (!is.null(og$up) && nrow(og$up) > 0) {
          p <- barplot(og$up, showCategory = 10) + ggtitle(paste("GO Enrichment (Up-reg): ", ont_label)) +
            theme_minimal()
          print(p)
          } else {
          plot.new()
          text(0.5, 0.5, "No hay términos GO enriquecidos", cex = 1.2)
        }
      })
        
      output$ora_go_down_p <- renderPlot({
        req(ora_go())
        og <- ora_go()
        ont_name <- c("BP" = "Biological Process",
                      "MF" = "Molecular Function",
                      "CC" = "Cellular Component")
        ont_label <- ont_name[input$ont] #para el título
        
        if (!is.null(og$down) && nrow(og$down) > 0) {
          p <- barplot(og$down, showCategory = 10) + ggtitle(paste("GO Enrichment (Down-reg): ", ont_label)) +
            theme_minimal()
          print(p)
        } else {
          plot.new()
          text(0.5, 0.5, "No hay términos GO enriquecidos", cex = 1.2)
        }
      })
      
      #tablas
      output$ora_go_up_t <- renderDT({
        req(ora_go())
        og <- ora_go()
        if (!is.null(og$up)) {
          datatable(
            as.data.frame(og$up),            
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 10)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 4)
        }
      })
      
      output$ora_go_down_t <- renderDT({
        req(ora_go())
        og <- ora_go()
        if (!is.null(og$down)) {
          datatable(
            as.data.frame(og$down),
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 10)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 4)
        }
      })
      
      #ORA - KEGG
      
      #gráficos
      output$ora_kegg_up_p <- renderPlot({
        req(ora_kegg())
        ok <- ora_kegg()
        if (!is.null(ok$up) && nrow(ok$up) > 0) {
          p <- barplot(ok$up, showCategory = 10) + ggtitle("Rutas KEGG enriquecidas (Up-reguladas)") +
            theme_minimal()
         print(p)
        } else {
          plot.new()
          text(0.5, 0.5, "No hay Rutas KEGG enriquecidas", cex = 1.2)
        }
      })
      
      output$ora_kegg_down_p <- renderPlot({
        req(ora_kegg())
        ok <- ora_kegg()
        if (!is.null(ok$down) && nrow(ok$down) > 0) {
          p <- barplot(ok$down, showCategory = 10) + ggtitle("Rutas KEGG enriquecidas (Down-reguladas)") +
            theme_minimal()
          print(p)
        } else {
          plot.new()
          text(0.5, 0.5, "No hay Rutas KEGG enriquecidas", cex = 1.2)
        }
      })
      
      #tablas
      output$ora_kegg_up_t <- renderDT({
        req(ora_kegg())
        ok <- ora_kegg()
        if (!is.null(ok$up)) {
          datatable(
            as.data.frame(ok$up),
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 10)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 4)
        }
      })
      
      output$ora_kegg_down_t <- renderDT({
        req(ora_kegg())
        ok <- ora_kegg()
        if (!is.null(ok$down)) {
          datatable(
            as.data.frame(ok$down),
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 10)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 4)
        }
      })
      
      # GSEA - GO
      #gráfico
      output$gsea_go_p <- renderPlot({
        req(gsea_go())
        gg <- gsea_go()
        ont_name <- c("BP" = "Biological Process",
                      "MF" = "Molecular Function",
                      "CC" = "Cellular Component")
        ont_label <- ont_name[input$ont] #para el título
        
        if (!is.null(gg) && nrow(gg) > 0) {
          p <- dotplot(gg, showCategory = 10) + ggtitle(paste("GSEA GO: ", ont_label)) +
            theme_minimal()
          print(p)
        } else {
          plot.new()
          text(0.5, 0.5, "No hay términos GO enriquecidos", cex = 1.2)
        }
      })
      
      #tabla
      output$gsea_go_t <- renderDT({
        req(gsea_go())
        gg <- gsea_go()
        if (!is.null(gg)) {
          datatable(
            as.data.frame(gg),
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 20)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue", "NES"),
                         digits = 4)
        }
      })
      
      # GSEA - KEGG
      #gráfico
      output$gsea_kegg_p <- renderPlot({
        req(gsea_kegg())
        gk <- gsea_kegg()
        if (!is.null(gk) && nrow(gk) > 0) {
          p <- dotplot(gk, showCategory = 10) +
            theme_minimal()
          print(p)
        } else {
          plot.new()
          text(0.5, 0.5, "No hay Rutas KEGG enriquecidas", cex = 1.2)
        }
      })
      
      #tabla
      output$gsea_kegg_t <- renderDT({
        req(gsea_kegg())
        gk <- gsea_kegg()
        if (!is.null(gk)) {
          datatable(
            as.data.frame(gk),
            rownames = F,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip", #orden de las secciones botones, filter...
              buttons = c("copy", "csv", "excel"),
              scrollX = T,
              pageLength = 20)) |> 
            formatSignif(columns = c("pvalue", "p.adjust", "qvalue", "NES"), 
                         digits = 4)
        }
      })
      
      return(reactive({
        list(ora_go = ora_go(),
             ora_kegg = ora_kegg(),
             gsea_go = gsea_go(),
             gsea_kegg = gsea_kegg())
      }))
  })  
}