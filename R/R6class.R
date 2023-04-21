## Define R6 class----------------------------------------------------------------------------------------------------
PPI.enrich.R6.Shiny = R6Class(
  "PPI.enrich.R6.Shiny",

  # private values --------------------------------------------------------------------------------
  private = list(
    currentValue = NULL,
    cyj = NULL, # 储存第一时间的cyj对象
    network_json = NULL,
    session = NULL,
    layout = NULL, # 储存node position 信息
    style_arguments = NULL # 储存 manual style 的颜色信息
  ),

  # public values --------------------------------------------------------------------------------
  public = list(

    return_cyj = function(){
      return(private$cyj)
    }, ## 返回 cyj 对象

    return_argument = function(){
      return(private$style_arguments)
    },

    return_session = function(){
      return(private$session)
    },

    return_network_json = function(){
      return(private$network_json)
    },

    PPIres = NA,
    ORAres = NA,
    DEPres = NA,

    # the initialize function in R6
    initialize = function(PPI_res = NA, ORA_res = NA, DEP_res = NA){
      message("initializing Simple.R6.shiny app")
      private$currentValue <- 0
      self$PPIres = PPI_res
      self$ORAres = ORA_res
      self$DEPres = DEP_res
    },

    # the runapp function
    run_app = function(port = NULL){
      runApp(shinyApp(self$ui, self$server), port, launch.browser=TRUE)
    },

    # the ui function in public------------------------------------------------------------
    ui = function(){
      fluidPage(

        sidebarLayout(



          sidebarPanel(
            ## layout
            # uiOutput("selectnode_UI"),
            h4("Network-enrichment"),
            uiOutput("layoutname_ui"),

            br(),
            actionButton("saveLayout", "Store Layout"),
            actionButton("restoreLayout", "Restore Layout"),
            # actionButton("clearselection", "Deselect"),
            br(),br(),
            hr(),
            h4(strong("Highlight node")),
            actionButton("hightselection", "Highlight selection"),
            actionButton("unhightselection", "Unhighlight selection"),
            br(),br(),
            actionButton("invertselection", "Invert selection"),
            br(),hr(),
            h4(strong("Hide node")),
            actionButton("hideselection", "Hide selection"),
            actionButton("showall", "Show all nodes"),

            br(),hr(),
            h4(strong("Style options")),
            selectInput("visualStyleSelector", "Select Visual Style",
                        choices=styles),
            column(colourInput("high_color",label = "High color",value = "#E30ACF"), width = 4),
            column(colourInput("mid_color",label = "Mid color",value = "#F2F2F2"), width = 4),
            column(colourInput("low_color",label = "Low color",value = "#1BABBF"), width = 4),
            column(numericInput("col_limit",label = "Color outer limit",value = 4), width = 6),
            column(numericInput("mid_limit",label = "Color middle blank",value = 1), width = 6),
            actionButton(inputId = "reDraw","Redraw in style"),
            br(),
            actionButton("saveStyle", "Store manual style"),
            actionButton("restoreLayout", "Restore style"),
            # h6("Send random node 'lfc' attributes (visible only with Biological Style, mapped to color):"),
            # actionButton("randomNodeAttributes", "Send"),
            br(),
            hr(),
            actionButton("sfn", "Select First Neighbor"),
            actionButton("fit", "Fit Graph"),
            actionButton("fitSelected", "Fit Selected"),
            br(),br(),
            # h6("Try out png-saving capability, using the currently displayed network"),
            actionButton("savePNGbutton", "Generate PNG"),
            uiOutput("downPNG_ui"),
            br(),br(),
            actionButton("saveJSbutton", "Generate JSON"),
            uiOutput("downJS_ui"),
            width=3
            #style="margin-right:10px; padding-right:0px;"
          ),
          mainPanel(
            fluidRow(
              column(cyjShinyOutput('cyjShiny', height="800px"),width = 12)

            ),
            fluidRow(
              column(DTOutput("enriched_table"),width = 12)

            )
            #style="margin-left:0px; padding-left:0px;"
          )
        ) # sidebarLayout
      )}, # ui

    # shiny server function, initialize though the results data stored in S6 public ------------------------------------------------------------
    server = function(input, output, session,
                      PPIres = self$PPIres, ORAres = self$ORAres, DEPres = self$DEPres){


      ora_res2 <- ora_res %>%
        clusterProfiler.dplyr::filter(qvalue<0.2,pvalue < 0.05) %>%
        as.data.frame()

      observeEvent(input$fit, ignoreInit=TRUE, {
        fit(session, 80)
      })

      observeEvent(input$fitSelected,  ignoreInit=TRUE,{
        fitSelected(session, 100)
      })

      # save layout, transmitted to input$tbl.nodePositions
      observeEvent(input$saveLayout, ignoreInit=TRUE, {
        getNodePositions(session)
      })

      observeEvent(input$tbl.nodePositions, ignoreInit=TRUE, {
        printf("new posTable ready\n");
        tbl <- fromJSON(input$tbl.nodePositions) # layout is returned in tbl.nodePositions
        private$layout <- tbl  # store layout in
        # state$tbl.nodePositions <- tbl
        print(tbl)
      })


      # restore layout,
      observeEvent(input$restoreLayout, ignoreInit=TRUE, {
        tbl.pos <- private$layout
        # tbl.pos11 <<-tbl.pos
        if( (!is.null(tbl.pos)) && nrow(tbl.pos) > 0)
          setNodePositions(session, tbl.pos)
      })


      # save manual style
      observeEvent(input$saveStyle,ignoreInit = T,{

      })

      observeEvent(input$restoreStyle, ignoreInit=TRUE, {

      })

      # observeEvent(input$clearselection, ignoreInit=TRUE, {
      #   clearSelection(session)
      # })

      observeEvent(input$hideselection, ignoreInit=TRUE, {
        hideSelection(session)
      })

      # observeEvent(input$removeSelection, ignoreInit=TRUE, {
      #   session$sendCustomMessage(type = "removeSelection", message = list())
      # })

      # invert selection ----
      observeEvent(input$invertselection, ignoreInit=TRUE, {
        getSelectedNodes(session)
      })

      invert <- reactiveVal(0)
      observeEvent(input$selectedNodes, ignoreInit=TRUE, {
        if(invert() != input$invertselection & input$invertselection > 0){
          Sys.sleep(0.3)

          PPI_nodes <- data.frame(id = PPIres[,1:2] %>% unlist %>% unique())

          selected = input$selectedNodes
          unselected = PPI_nodes$id[which(!PPI_nodes$id %in% selected)]
          clearSelection(session)
          Sys.sleep(0.3)
          selectNodes(session, unselected)
          invert <- invert({
            input$invertselection
          })
        }

      })


      # hight selection -----
      hightlight <- reactiveVal(0)
      observeEvent(input$hightselection, ignoreInit=TRUE, {
        cat("awef")
        # session$sendCustomMessage("getSelected2", message = list())
        session$sendCustomMessage("getALL", message = list())
        getSelectedNodes(session)
      })

      observeEvent(input$selectedNodes, ignoreInit=TRUE, {
        if(hightlight() != input$hightselection & input$hightselection > 0){
          Sys.sleep(0.3)
          selected = input$selectedNodes
          selected111 <<- selected
          cat("hightlight nodes\n")

          # # get node attribution
          # all = input$allNodesdata
          # all111 <<- all
          # cat("awefaw nodes\n")
          # all2 <- matrix(all,nrow = names(all) %>% unique() %>% length()) %>% t %>% as.data.frame()
          # colnames(all2) = names(all) %>% unique()
          #
          # # Change highlighted attribution
          # highed = all2$highlighted %>% as.logical()
          # highed[which(all2$id %in% selected)] = TRUE

          cyjShiny::setNodeAttributes(session, attributeName="highlighted", nodes=selected, values=rep(T,length(selected)))

          hightlight <- hightlight({
            input$hightselection
          })
        }
      })

      # Unhightlight ----
      unhightlight <- reactiveVal(0)
      observeEvent(input$unhightselection, ignoreInit=TRUE, {
        cat("awef")
        # session$sendCustomMessage("getSelected2", message = list())
        session$sendCustomMessage("getALL", message = list())
        getSelectedNodes(session)
      })


      observeEvent(input$selectedNodes, ignoreInit=TRUE, {
        if(unhightlight() != input$unhightselection & input$unhightselection > 0){
          Sys.sleep(0.3)
          selected = input$selectedNodes
          # selected111 <<- selected
          cat("unhightlight nodes\n")

          # exit node attribution
          all = input$allNodesdata
          # all111 <<- all
          all2 <- matrix(all,nrow = names(all) %>% unique() %>% length()) %>% t %>% as.data.frame()
          colnames(all2) = names(all) %>% unique()
          # all22222<<-all2

          # Change highlighted attribution
          highed = all2$highlighted %>% as.logical()
          highed[which(all2$id %in% selected)] = FALSE
          cyjShiny::setNodeAttributes(session, attributeName="highlighted", nodes=all2$id, values=highed)

          unhightlight <- unhightlight({
            input$unhightselection
          })
        }
      })

      observeEvent(input$showall, ignoreInit=TRUE, {
        showAll(session)
      })

      observeEvent(input$sfn,  ignoreInit=TRUE,{
        selectFirstNeighbors(session)
      })

      ## change style
      init_style = reactiveVal(1)
      observeEvent(
        {
          input$reDraw
          input$visualStyleSelector
          Manual_style_jstext()
          input$Layoutname
        },
        ignoreInit=F,
        ignoreNULL = F,
        {
          newStyleFile <- input$visualStyleSelector
          newStyleFile222 <<- input$visualStyleSelector
          printf("newStyle: %s\n", newStyleFile)

          input_arguments <<- private$style_arguments
          if(newStyleFile == "Manual"){
            if(init_style() < 3
               &&
               ((!is.null(input_arguments)) && input_arguments$visualStyleSelector == "Manual")
            ){
              Sys.sleep(1)
              init_style = init_style(init_style() + 1)
            }

            Sys.sleep(0.5)
            cat("Use a manual style \n")
            # jsonText <- default_jsonText
            jsonText <- Manual_style_jstext()
            jsonText111 <<- jsonText
            message <- list(json = jsonText)
            session <- shiny::getDefaultReactiveDomain()
            session$sendCustomMessage("loadStyle", message)
          }else{
            loadStyleFile(newStyleFile)
          }
          # private$session = shiny::getDefaultReactiveDomain()
        }
      )

      # change layout, 根据 Layoutname dolayout
      init_layout <- reactiveVal({
        init_layout111 <<-
          private$layout
      })
      observeEvent(input$Layoutname,  ignoreInit=F,{
        if(input$Layoutname != ""){
          strategy <- input$Layoutname
          Sys.sleep(0.5)
          doLayout(session, strategy)
          print(paste0("Change layout to ",strategy))
          # later::later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
        }
        if( (!is.null(init_layout())) && nrow(init_layout()) > 0){
          Sys.sleep(1.5)
          cat("Restore the layout based on a saved position at app launch \n")
          print(head(init_layout111))
          setNodePositions(session, init_layout111)
          init_layout = init_layout(NULL)
        }
      })


      ## output cyj and enrichment table ----
      the_cjs <- reactive({
        # if(is.null(input$Layoutname)) {
        layoutName = "cose"
        # }else{layoutName = input$Layoutname}
        the_cjs <- PPI2cyj(PPIres = PPIres,
                           L2FC_tab = get_signicant(DEPres,contrasts = "W4_vs_PBS",
                                                    return_type = "table") %>% dplyr::select(name, W4_vs_PBS_ratio),
                           layoutName = layoutName) ## 这里可以如果layoutname 固定成 cose
      })


      # render outputs ----

      output$cyjShiny <- renderCyjShiny({
        private$cyj <- the_cjs() # save the initial cyj in private$cyj
        # private$session <- session

        the_cjs()
      })

      output$enriched_table <- renderDT({
        the_res <- as.data.frame(ora_res2)
        the_res[,c('pvalue', 'p.adjust','qvalue')] <- signif(the_res[,c('pvalue', 'p.adjust','qvalue')],3)
        the_res[,c('pvalue', 'p.adjust','qvalue')] <- format(the_res[,c('pvalue', 'p.adjust','qvalue')],scientific = T)
        DT::datatable(the_res,options = list(scrollX = TRUE))

      })

      observeEvent(input$enriched_table_rows_selected,ignoreInit = T,ignoreNULL = T,{
        print(input$enriched_table_rows_selected)
        selected_terms <- ora_res2[input$enriched_table_rows_selected,]
        print(selected_terms)
        selected_genes = selected_terms$geneID %>% sapply(.,strsplit,split = "/") %>% unlist %>% unique()
        print(selected_genes)
        clearSelection(session)
        Sys.sleep(0.5)
        selectNodes(session, selected_genes)
      })

      ## save as PNG ----
      observeEvent(input$savePNGbutton, ignoreInit=TRUE, {
        file.name <- tempfile(fileext=".png")
        # savePNGtoFile(session, file.name)
        session$sendCustomMessage(type = "savePNGtoFile", message = list())
      })

      PNG_binary <- reactive({
        if(!is.null(input$pngData)){
          cat("received pngData")
          theJSON <- input$pngData

          png.parsed <- fromJSON(theJSON)
          substr(png.parsed, 1, 30) # [1] "data:image/png;base64,iVBORw0K"
          nchar(png.parsed)  # [1] 768714
          png.parsed.headless <- substr(png.parsed, 23, nchar(png.parsed))  # chop off the uri header
          png.parsed.binary <- base64decode(png.parsed.headless)
        }else{NULL}
      })

      output$downPNG_ui <- renderUI({
        if(!is.null(PNG_binary())){
          downloadButton("downPNGbutton", "Download PNG")
        }
      })

      output$downPNGbutton <- downloadHandler(
        filename = function() {
          # session$sendCustomMessage(type = "savePNGtoFile11", message = list(file.name))
          # cat("awef\n")
          paste("foo", ".png", sep="")
        },
        content = function(file) {
          print("writing png to foo.png \n")

          conn = file(file, "wb")
          writeBin(PNG_binary(), conn)
          close(conn)
        }
      )

      Manual_style_jstext <- reactive({
        visualStyleSelector <- input$visualStyleSelector
        if(visualStyleSelector != "Manual") # use other style file
          return(11)

        cat("Generate a manual_style_jstext\n")
        obj = default_obj
        obj[obj$selector %>% grep("node\\[lfc",.,value=F),] -> temp

        # high_col <- c("#E30ACF")
        # mid_col <- c("#F2F2F2")
        # low_col <- c("#1BABBF")
        high_col = input$high_color
        mid_col = input$mid_color
        low_col = input$low_color

        col_lim = abs(input$col_limit)
        mid_value = abs(input$mid_limit)
        temp$css$`background-color`[1:2] = high_col
        temp$css$`background-color`[3] = paste0("mapData(lfc,",mid_value,",",col_lim,
                                                ",",mid_col,",",high_col,")")
        temp$css$`background-color`[4] = paste0("mapData(lfc,",-col_lim,",",-mid_value,
                                                ",",low_col,",",mid_col,")")
        temp$css$`background-color`[5:6] = low_col

        temp$selector[1] = paste0("node[lfc > ",col_lim,"]")
        temp$selector[2] = paste0("node[lfc = ",col_lim,"]")
        temp$selector[3] = paste0("node[lfc < ",col_lim,"][node[lfc >= ",mid_value,"]")
        temp$selector[4] = paste0("node[lfc <= ",-mid_value,"][node[lfc > ",-col_lim,"]")
        temp$selector[5] = paste0("node[lfc = ",-col_lim,"]")
        temp$selector[6] = paste0("node[lfc < ",-col_lim,"]")


        temp$selector[7] = paste0("node[lfc < ",mid_value,"][lfc > ",-mid_value,"]")
        temp$css$`background-color`[7] = mid_col

        obj[obj$selector %>% grep("node\\[lfc",.,value=F),] <- temp
        obj111 <<- obj
        default_jsonText <- as.character(toJSON(obj))
      })

      observeEvent(input$saveJSbutton, ignoreInit=TRUE, {
        print("generated a JSON data")
        session$sendCustomMessage(type = "saveJSONtoFile", message = list())
      })

      # observeEvent(input$JSONData, ignoreInit=TRUE, {
      #   print("received JSONData")
      #   theJSON11 <<- input$JSONData
      #   theJSON <- jsonlite::toJSON(theJSON11, pretty = TRUE, auto_unbox = TRUE)
      #   cat(theJSON , file = "foo.json", fill = FALSE, labels = NULL, append = FALSE)
      #
      # })

      JSON_data <- reactive({
        if(!is.null(input$JSONData)){
          # theJSON <- input$JSONData
          theJSON <- jsonlite::toJSON(input$JSONData, pretty = TRUE, auto_unbox = TRUE)
          private$network_json = theJSON
          return(theJSON)
        }else{NULL}
      })

      output$downJS_ui <- renderUI({
        if(!is.null(JSON_data())){
          cat("generated successfully. \n")
          downloadButton("downJSbutton", "Download network json")
        }
      })

      # * render options ----
      output$layoutname_ui <- renderUI({
        # cat("awgewfdw")
        # temp111 <<- private$style_arguments

        selectInput("Layoutname", "Select Layoutname:",
                    choices = c("preset", "cose", "cola", "circle",
                                "concentric", "breadthfirst", "grid",
                                "random", "dagre", "euler",
                                "fcose","springy","spread"),
                    selected = ifelse(is.null(private$style_arguments$Layoutname),
                                      "cose",
                                      private$style_arguments$Layoutname)
        )
      })

      # updata style options after initialize
      observeEvent(input$high_color,{
        # temp111 <<- input$high_color
        # temp222 <<-  private$style_arguments$high_color
        if(!is.null(private$style_arguments)){
          updateColourInput(session,"high_color", value = private$style_arguments$high_color)
          updateColourInput(session,"low_color", value = private$style_arguments$low_color)
          updateColourInput(session,"mid_color", value = private$style_arguments$mid_color)
          updateNumericInput(session,"col_limit", value = private$style_arguments$col_limit)
          updateNumericInput(session,"mid_limit", value = private$style_arguments$mid_limit)
          updateSelectInput(session,"visualStyleSelector",
                            choices = styles,
                            selected = private$style_arguments$visualStyleSelector)

          # ## 初始化style
          # newStyleFile <- input$visualStyleSelector
          # newStyleFile111 <<- newStyleFile
          # printf("newStyle: %s\n", newStyleFile)
          # if(newStyleFile == "Manual"){
          #   Sys.sleep(2)
          #   cat("First use a manual style \n")
          #   # jsonText <- default_jsonText
          #   jsonText <- Manual_style_jstext()
          #   jsonText111 <<- jsonText
          #   col_limit111 <<- input$col_limit
          #   message <- list(json = jsonText)
          #   session <- shiny::getDefaultReactiveDomain()
          #   session$sendCustomMessage("loadStyle", message)
          # }else{
          #   loadStyleFile(newStyleFile)
          # }
          # private$session = shiny::getDefaultReactiveDomain()

        }


      },ignoreNULL = F,ignoreInit = F,once = T)

      output$downJSbutton <- downloadHandler(
        filename = function() {
          # session$sendCustomMessage(type = "savePNGtoFile11", message = list(file.name))
          # cat("awef\n")
          paste("foo", ".json", sep="")
        },
        content = function(file) {
          print("writing JS to foo.json")

          conn = file(file, "wb")
          write(JSON_data(), conn)
          close(conn)
        }
      )

      # end session
      session$onSessionEnded(function() {
        isolate({
          private$style_arguments <- input
        })
      })

    } # server

  ) # public
) # class
