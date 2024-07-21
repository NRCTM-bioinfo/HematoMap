
#' @import shiny
ui <- fluidPage(
  navbarPage(
    "HematoMap v1.1.0",
    # tabPanel(
    #   "tab2", fluid = TRUE
    # ),
    navbarMenu(
      title = "scRNA",
      
      ## Input #####
      tabPanel(
        "#1 Input: Upload", 
        # tags$style('.form-group{margin-bottom: 0px;}'),
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(
              id = "data",
              tabPanel(
                title = "Seurat RDS",
                titlePanel(title = "Input: Seurat"),
                helpText("Upload Seurat RDS file"),
                tags$hr(),
                fileInput("scrna.rds", "data.rds",
                          multiple = FALSE,
                          accept = c(".rds")),
                tags$hr(),
              ),
              # tabPanel(
              #   title = "10x Raw",
              #   titlePanel(title = "Input: 10X"),
              #   helpText("Upload 10X standard files"),
              #   tags$hr(),
              #   fileInput("scrna.barcodes", "barcodes.tsv",
              #             multiple = FALSE,
              #             accept = c("text/tsv",
              #                        "text/comma-separated-values,text/plain",
              #                        ".tsv")),
              #   fileInput("scrna.genes", "genes.tsv",
              #             multiple = FALSE,
              #             accept = c("text/tsv",
              #                        "text/comma-separated-values,text/plain",
              #                        ".tsv")),
              #   fileInput("scrna.matrix", "matrix.mtx",
              #             multiple = FALSE,
              #             accept = c("text/comma-separated-values,text/plain",
              #                        ".mtx")),
              #   tags$hr(),
              # )
            ),
            actionButton(inputId = "upload.scrna", label = "Upload")
          ),
          mainPanel()
        )
      ),
      
      ## ThetaLike #####
      tabPanel(
        title = "#2 Plot: Theta-Like",
        sidebarLayout(
          sidebarPanel(
            titlePanel(title = "Theta-Like"),
            selectInput(inputId = "theta.like.plot.scell",
                        label = "plot.scell",
                        choices = c(""),
                        selected = "",
                        width = "300px"),
            sliderInput(inputId = "theta.like.point.size",
                        label = "point.size",
                        min = 0, max = 5,
                        value = 1, step = 0.5,
                        width = "300px"
                        ),
          ),
          mainPanel(
            plotOutput(outputId = "theta.like")
          )
        )
      ),
      
      
      ## CircleTree #####
      tabPanel(
        title = "#3 Plot: CircleTree",
        sidebarLayout(
          sidebarPanel(
            titlePanel(title = "CircleTree"),
            selectInput(inputId = "circle.tree.group.subc",
                        label = "group.subc",
                        choices = c("external", "internal", "reference"),
                        selected = "external",
                        width = "300px"
            ),
            selectInput(inputId = "circle.tree.top.subc",
                        label = "top.subc",
                        choices = as.character(2:10),
                        selected = "3",
                        width = "300px"
            ),
            selectInput(inputId = "circle.tree.color.mapping",
                        label = "color.mapping",
                        choices = c("cell.type", "cell.percentage"),
                        selected = "cell.type",
                        width = "300px"
            ),
            selectInput(inputId = "circle.tree.size.mode",
                        label = "size.mode",
                        choices = c("absolute.pct", "relative.pct"),
                        selected = "absolute.pct",
                        width = "300px"
            ),
            selectInput(inputId = "circle.tree.relative.compute",
                        label = "relative.compute",
                        choices = c("ratio", "difference"),
                        selected = "ratio",
                        width = "300px"
            ),
            # tags$hr(),
            # actionButton(inputId = "circle.tree", label = "Plot"),
            
            # textInput(inputId = "TimeFinderMax",
            #           label = "To:",
            #           value = "22.00",
            #           width = "100px"),
            # sliderInput(inputId = "RankOnTeam",
            #             label = "Select Swimmer Rank On Team",
            #             min = 1,
            #             max = 10,
            #             value = c(1,6),
            #             width = "300px")
          ),
          mainPanel(
            plotOutput(outputId = "circle.tree")
          )
        )
      ),
      
      
      ## ClusterTree #####
      tabPanel(
        title = "#4 Plot: ClusterTree",
        sidebarLayout(
          sidebarPanel(
            titlePanel(title = "ClusterTree"),
            selectInput(inputId = "cluster.tree.group.subc",
                        label = "group.subc",
                        choices = c("external", "internal", "reference"),
                        selected = "external",
                        width = "300px"
            ),
            selectInput(inputId = "cluster.tree.top.subc",
                        label = "top.subc",
                        choices = as.character(2:10),
                        selected = "3",
                        width = "300px"
            ),
            selectInput(inputId = "cluster.tree.color.mapping",
                        label = "color.mapping",
                        choices = c("cell.type", "cell.percentage"),
                        selected = "cell.type",
                        width = "300px"
            ),
            selectInput(inputId = "cluster.tree.size.mode",
                        label = "size.mode",
                        choices = c("absolute.pct", "relative.pct"),
                        selected = "absolute.pct",
                        width = "300px"
            ),
            selectInput(inputId = "cluster.tree.relative.compute",
                        label = "relative.compute",
                        choices = c("ratio", "difference"),
                        selected = "ratio",
                        width = "300px"
            ),
          ),
          mainPanel(
            plotOutput(outputId = "cluster.tree")
          )
        )
      ),
    ),
    navbarMenu(
      title = "Bulk",
      
      ## Input #####
      tabPanel(
        "#1 Input: Upload", 
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(
              id = "data",
              tabPanel(
                title = "Expression file",
                titlePanel(title = "Input: CSV"),
                helpText("Upload expression csv file"),
                tags$hr(),
                fileInput(inputId = "bulk.csv", "expr.csv",
                          multiple = FALSE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")),
                tags$hr(),
              ),
              tabPanel(
                title = "RDS file",
                titlePanel(title = "Input: RDS"),
                helpText("Upload RDS file"),
                tags$hr(),
                fileInput(inputId = "bulk.rds", "data.rds",
                          multiple = FALSE,
                          accept = c(".rds")),
                tags$hr(),
              ),
            ),
            actionButton(inputId = "upload.bulk", label = "Upload")
          ),
          mainPanel(
            # TODO text -> tips to steps
          )
        )
      ),
      
      ## LassoTree #####
      tabPanel(
        title = "#2 Plot: LassoTree",
        sidebarLayout(
          sidebarPanel(
            titlePanel(title = "LassoTree"),
            selectInput(inputId = "lasso.tree.color.mapping",
                        label = "color.mapping",
                        choices = c("value", "cell.type"),
                        selected = "value",
                        width = "300px"
            ),
          ),
          mainPanel(
            plotOutput(outputId = "lasso.tree")
          )
        )
      ),
      
      ## LassoLollipop #####
      tabPanel(
        title = "#3 Plot: LassoLollipop",
        sidebarLayout(
          sidebarPanel(
            titlePanel(title = "LassoLollipop"),
            sliderInput(inputId = "lasso.lollipop.point.size",
                        label = "point.size",
                        min = 1, max = 10,
                        value = 6,
                        width = "300px"),
          ),
          mainPanel(
            plotOutput(outputId = "lasso.lollipop", width = "100%")
          )
        )
      ),
      
    )
  )
)

#' @importFrom data.table fread
server <- function(input, output, session){
  message.time("Start App!")
  
  ## scRNA ######
  observeEvent(input$upload.scrna, {
    progress <- shiny::Progress$new()
    progress$set(message = "Start", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    updateProgress(value = 10, detail = "Read data")
    object <- NULL
    if(!is.null(input$scrna.rds)){
      object <- readRDS(input$scrna.rds$datapath)
      updateProgress(value = 20, detail = "Run Seurat Workflow")
      object <- runSeuratWorkflow(object, force.run = T)
      updateProgress(value = 40, detail = "Compute Subcluster")
      hemato.subc <- CreateSubclusterObject(object = object, n.core = 1)
      updateProgress(value = 80, detail = "Create Subcluster Object")
    }
    # hemato.subc <- readRDS("~/projects/scrna/tmp/hemato.subc.rds")
    
    # update - plot.scell
    updateSelectInput(session = session, inputId = "theta.like.plot.scell",
                      label = "plot.scell",
                      choices = rownames(hemato.subc@meta.data),
                      selected = rownames(hemato.subc@meta.data)[1])
    output$theta.like <- renderPlot({
      message.time("Plot Circle.Tree")
      plotThetaLike(hemato.subc=hemato.subc, plot.scell = input$theta.like.plot.scell, point.size = input$theta.like.point.size,
                    label.subc.size = 1.3, label.like.size = 4, text.size = 11, line.size = 0.5)
    }, height = 400, width = 400)
    
    output$circle.tree <- renderPlot({
      message.time("Plot Circle.Tree")
      plotCircleTree(hemato.subc=hemato.subc, group.subc=input$circle.tree.group.subc, top.subc = as.numeric(input$circle.tree.top.subc), 
                     color.mapping = input$circle.tree.color.mapping, size.mode = input$circle.tree.size.mode, 
                     relative.compute = input$circle.tree.relative.compute, 
                     label.size = 4.2, ribbon.size = 6, line.size = 0.5, point.size = 10, text.size = 11)
    }, height = 550, width = 500)
    
    output$cluster.tree <- renderPlot({
      message.time("Plot Cluster.Tree")
      plotClusterTree(hemato.subc=hemato.subc, group.subc=input$cluster.tree.group.subc, top.subc = as.numeric(input$cluster.tree.top.subc), 
                      color.mapping = input$cluster.tree.color.mapping, size.mode = input$cluster.tree.size.mode, 
                      relative.compute = input$cluster.tree.relative.compute, 
                      line.size = 0.5, ribbon.size = 6, point.size = 12, label.pop.size = 4.2, text.size = 11)
    }, height = 550, width = 500)
    
    
  })
  
  
  
  ## Bulk ######
  observeEvent(input$upload.bulk, {
    progress <- shiny::Progress$new()
    progress$set(message = "Start", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    updateProgress(value = 20, detail = "Read data")
    data.bulk <- NULL
    if(!is.null(input$bulk.rds)){
      data.bulk <- readRDS(input$bulk.rds$datapath)
      if(!is.data.frame(data.bulk) && !is.matrix(data.bulk)) stop("data should be data.frame or matrix")
    }else if(!is.null(input$bulk.csv)){
      data.bulk <- fread(input = input$bulk.csv, data.table = F)
    }
    if(!is.numeric(data.bulk[, 1])){
      rownames(data.bulk) <- data.bulk[,1]
      data.bulk <- data.bulk[,-1, drop = F]
    }
    
    updateProgress(value = 50, detail = "Compute Lasso Score")
    message.time("Compute Lasso Score")
    lasso.score <- computeLassoScore(as.matrix(data.bulk))
    updateProgress(value = 100, detail = "Complete")
    
    
    output$lasso.lollipop <- renderPlot({
      message.time("Plot Lasso.Lollipop")
      plotLassoLollipop(mat.score=lasso.score, point.size = input$lasso.lollipop.point.size, text.size = 12)
    }, height = 500, width = 300)
    
    output$lasso.tree <- renderPlot({
      message.time("Plot Lasso.Tree")
      plotLassoTree(mat.score=lasso.score, color.mapping = input$lasso.tree.color.mapping, 
                    line.size = 0.5, ribbon.size = 6, point.size = 3, label.pop.size = 4.2, 
                    text.size = 11, legend.barheight = unit(4, "mm"))
    }, height = 550, width = 500)
  })
  
}


#' Run Shiny-based Interactive User Interface
#' 
#' @param options Named options that should be passed to the runApp call (these 
#'   can be any of the following: "port", "launch.browser", "host", "quiet", 
#'   "display.mode" and "test.mode"). You can also specify width and height 
#'   parameters which provide a hint to the embedding environment about the 
#'   ideal height/width for the app.
#' 
#' @import shiny
#' @export
#' 
runApp <- function(options = list()){
  options(shiny.maxRequestSize=200*1024^2)
  shinyApp(ui = ui, server = server, options = list())
}

