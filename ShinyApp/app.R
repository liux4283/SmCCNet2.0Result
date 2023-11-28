library(shiny)
library(shinydashboard)
library(igraph)
plotMultiOmicsNetwork <- function(Abar, CorrMatrix, 
                                  ModuleIdx = 1, sub_type, CorrCutOff = 0.2,
                                  AddCorrSign = TRUE,
                                  ShowTypeLabel = TRUE,
                                  PlotTitle = "", NetLayout = "lgl",
                                  ShowNodes = TRUE,
                                  VertexLabelCex = 1, VertexSize = 1,
                                  EdgeColorIntensity = 1){
  
  p <- ncol(Abar)
  
  
  # Trim the module by CorrCutOff.
  M <- as.matrix(Abar)
  CorrSubMatrix <- as.matrix(CorrMatrix)
  M[which(abs(CorrSubMatrix) < CorrCutOff)] <- 0
  if(AddCorrSign){M <- M * sign(CorrSubMatrix)}
  M.node <- 1:nrow(Abar)
  newM.node <- M.node[which(apply(abs(M), 1, max) > 0)]
  
  # Get the actual names of the nodes from the row names of Abar
  node_names <- rownames(Abar)[newM.node]
  
  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(seq_len(p), ncol = 1)
    rownames(allidx) <- rownames(Abar)
    
    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                               diag = FALSE, mode = "undirected")
    
    # Vertex Color and Shape Definition
    k <- length(newM.node)
    unique_types <- unique(sub_type)
    type_colors <- c("darkviolet", "darkred", "midnightblue", "forestgreen", "deeppink")  # Darker color alternatives
    #type_colors <- c("purple", "dark orange", "blue", "green", "pink")  # Extend color list as needed
    type_shapes <- c("square", "circle", "triangle", "diamond", "rectangle")  # Extend shape list as needed
    vcol <- type_colors[match(sub_type[newM.node], unique_types)]
    vshape <- type_shapes[match(sub_type[newM.node], unique_types)]
    
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowTypeLabel){node_names <- rep(" ", length(newM.node))}
    
    # Define edge colors.
    ecol <- rep("gray80", igraph::ecount(net))
    ew <- abs(igraph::edge.attributes(net)$weight) * EdgeColorIntensity  # use the new argument here
    ecol[which(igraph::edge.attributes(net)$weight < 0)] <- grDevices::adjustcolor("blue", alpha.f = EdgeColorIntensity)  # adjust color as needed
    ecol[which(igraph::edge.attributes(net)$weight > 0)] <- grDevices::adjustcolor("red", alpha.f = EdgeColorIntensity)  # adjust color as needed
    # Define network layout.
    if(NetLayout == "circle"){
      l <- igraph::layout_in_circle(net)
    } else if(NetLayout == "sphere"){
      l <- igraph::layout_on_sphere(net)
    } else if(NetLayout == "fr"){
      l <- igraph::layout_with_fr(net)
    } else if(NetLayout == "lgl"){
      l <- igraph::layout_with_lgl(net)
    } else if(NetLayout == "kk"){
      l <- igraph::layout_with_kk(net)
    } else if(NetLayout == "random"){
      l <- igraph::layout_randomly(net)
    } else if(NetLayout == "grid"){
      l <- igraph::layout_on_grid(net)
    } else if(NetLayout == "drl"){
      l <- igraph::layout_with_drl(net)
    } else if(NetLayout == "star"){
      l <- igraph::layout_as_star(net)
    } else if(NetLayout == "tree"){
      l <- igraph::layout_as_tree(net)
    } else if(NetLayout == "bipartite"){
      l <- igraph::layout_as_bipartite(net)
    } else {
      stop("Unrecognized NetLayout input. Acceptable options are 'circle', 'sphere', 'fr', 'lgl', 'kk', 'random', 'grid', 'drl', 'star', 'tree', and 'bipartite'.")
    }
    
    
    # Plotting the network
    graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                   vertex.label.cex = VertexLabelCex, layout = l,
                   vertex.size = VertexSize, vertex.label = node_names, 
                   edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                   vertex.label.font = 2, main = PlotTitle)
    
  }
  
}


library(shiny)
library(shinydashboard)
library(igraph)




# Ensure the plotMultiOmicsNetwork function is loaded

ui <- dashboardPage(
  dashboardHeader(title = "SmCCNet Subnetwork Visualization"),
  dashboardSidebar(
    sidebarMenu(
      actionButton("info", "Info"),
      fileInput("file", "Upload .Rdata file",
                accept = c(".Rdata")
      ),
      sliderInput("CorrCutOff", "Correlation Cut-Off:", 
                  min = 0, max = 1, value = 0.2, step = 0.01
      ),
      selectInput("NetLayout", "Network Layout:",
                  choices = c("Large Graph Layout" = "lgl", 
                              "Circular Layout" = "circle", 
                              "Random Layout" = "random", 
                              "Grid Layout" = "grid", 
                              "Star Layout" = "star", 
                              "Tree Layout" = "tree")
      ),
      sliderInput("VertexLabelCex", "Vertex Label Size:", 
                  min = 0.5, max = 2, value = 1, step = 0.1
      ),
      sliderInput("VertexSize", "Vertex Size:", 
                  min = 0.5, max = 2, value = 1, step = 0.1
      ),
      sliderInput("EdgeColorIntensity", "Edge Intensity:", 
                  min = 1, max = 100, value = 1, step = 1
      ),
      actionButton("plot_network", "Plot Network"),
      downloadButton("download_plot", "Download Plot")
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
            .main-header .logo {
                font-size: 12px;
                overflow: visible;
                white-space: normal;
            }
        "))
    ),
    fluidRow(
      box(plotOutput("network_plot"), width = 12)
    )
  )
)

server <- function(input, output, session) {
  
  
  # Load example data at the start
  example_data <- local({
    # Provide the correct path to your .Rdata file
    load("example_subnetwork.Rdata") 
    list(M = M, correlation_sub = correlation_sub, sub_type = sub_type)
  })
  
  data <- reactiveVal(example_data) 
  
  observeEvent(input$file, {
    req(input$file)
    load(input$file$datapath)
    data(list(M = M, correlation_sub = correlation_sub, sub_type = sub_type))
  })
  
  observeEvent(input$plot_network, {
    req(data())
    output$network_plot <- renderPlot({
      plotMultiOmicsNetwork(data()$M, data()$correlation_sub, sub_type = data()$sub_type,
                            CorrCutOff = input$CorrCutOff,
                            NetLayout = input$NetLayout,
                            VertexLabelCex = input$VertexLabelCex,
                            VertexSize = input$VertexSize,
                            EdgeColorIntensity = input$EdgeColorIntensity)  # pass the new input here
    })
  })
  
  observeEvent(input$info, {
    showModal(
      modalDialog(
        title = "How to use this app",
        "This app allows you to visualize single/multi-omics networks with the subnetwork created and stored by SmCCNet. Here's how to use it:",
        tags$ul(
          tags$li("An example network can be plotted simply by clicking the 'Plot Network' without uploading anything."),
          tags$li("Upload a .Rdata file containing the variables M, correlation_sub, and sub_type, the .Rdata file should be in 'size_a_net_b.Rdata' format, where a is the pruned network size, and b is the network module index after hierarchical clustering."),
          tags$li("Adjust the Correlation Cut-Off slider to filter edges in the network based on their correlation values."),
          tags$li("Select a Network Layout from the dropdown menu to change the layout of the network visualization."),
          tags$li("Adjust the Vertex Label Size and Vertex Size sliders to change the size of the vertex labels and vertices, respectively."),
          tags$li("Adjust the Edge Intensity slider to change the color intensity and width of the edges."),
          tags$li("Click the 'Plot Network' button to generate the network visualization."),
          tags$li("Click the 'Download Plot' button to download the network visualization as a PDF.")
        ),
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("network_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(data())
      # Save the plot to the output file
      pdf(file)
      plotMultiOmicsNetwork(data()$M, data()$correlation_sub, sub_type = data()$sub_type,
                            CorrCutOff = input$CorrCutOff,
                            NetLayout = input$NetLayout,
                            VertexLabelCex = input$VertexLabelCex,
                            VertexSize = input$VertexSize,
                            EdgeColorIntensity = input$EdgeColorIntensity)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)
