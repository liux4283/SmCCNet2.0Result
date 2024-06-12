library(shiny)
library(shinydashboard)
library(igraph)
library(ggplot2)
library(grid)
library(tidyverse)
library(shiny)
library(shinydashboard)
library(igraph)
library(DT)
library(plotly)
library(reshape2)
library(visNetwork)
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
    type_shapes <- rep('circle',5)  # Extend shape list as needed
    vcol <- type_colors[match(sub_type[newM.node], unique_types)]
    vshape <- type_shapes[match(sub_type[newM.node], unique_types)]
    
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowTypeLabel){node_names <- rep(" ", length(newM.node))}
    
    # Define edge colors.
    ecol <- rep("gray80", igraph::ecount(net))
    ew <- abs(igraph::edge.attributes(net)$weight) * EdgeColorIntensity 
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



# Adding the UI
ui <- dashboardPage(
  dashboardHeader(title = tags$span("SmCCNet Subnetwork Visualization Tool", style = "font-size: 14px;")),
  dashboardSidebar(
    sidebarMenu(
      id = "sidebarMenu", # Add an ID to the sidebarMenu for tracking tab selections
      menuItem("Network Visualization", tabName = "network_viz", icon = icon("connectdevelop")),
      menuItem("PC Loading Visualization", tabName = "pc_loading_viz", icon = icon("bar-chart")),
      menuItem("Omics Correlation Data", tabName = "omics_corr_data", icon = icon("table")),
      menuItem("Correlation Heatmap", tabName = "corr_heatmap", icon = icon("th")),
      menuItem("PCA 3D Scatter Plot", tabName = "pca_scatter_plot", icon = icon("globe")),
      actionButton("info", "Info"),
      fileInput("file", "Upload .Rdata file", accept = c(".Rdata")),
      
      # ConditionalPanel for network visualization controls
      conditionalPanel(
        condition = "input.sidebarMenu === 'network_viz'",
        sliderInput("CorrCutOff", "Correlation Cut-Off:", min = 0, max = 1, value = 0.2, step = 0.01),
        selectInput("NetLayout", "Network Layout:",
                    choices = c("Circular Layout" = "circle", 
                                "Large Graph Layout" = "lgl", 
                                "Random Layout" = "random", 
                                "Grid Layout" = "grid", 
                                "Star Layout" = "star", 
                                "Tree Layout" = "tree")),
        sliderInput("VertexLabelCex", "Vertex Label Size:", min = 0.1, max = 10, value = 0.5, step = 0.1),
        sliderInput("VertexSize", "Vertex Size:", min = 0.1, max = 10, value = 0.5, step = 0.1),
        sliderInput("EdgeColorIntensity", "Edge Intensity:", min = 0.5, max = 100, value = 2, step = 0.5),
        downloadButton("download_plot", "Download Plot")
      )
    )
  ),
  dashboardBody(
    tabItems(
      # Network Visualization Tab
      tabItem(tabName = "network_viz",
              fluidRow(
                box(plotOutput("network_plot"), width = 12)
              )),
      # PC Loading Visualization Tab
      tabItem(tabName = "pc_loading_viz",
              selectInput("pc_choice", "Choose Principal Component:", choices = c("PC1", "PC2", "PC3")),
              plotOutput("pc_loading_plot"),
              downloadButton("download_pc_plot", "Download PC Plot")),
      tabItem(tabName = "omics_corr_data",
              DTOutput("omics_correlation_table"),
              downloadButton("download_table", "Download Table")  # The button for downloading the table
      ),
      tabItem(tabName = "corr_heatmap",
              plotOutput("heatmap_plot"),
              downloadButton("download_heatmap", "Download Heatmap")),
      tabItem(tabName = "pca_scatter_plot",
              plotlyOutput("pca_plot"))
    )
  )
)

# Adding the Server Logic
server <- function(input, output, session) {
  # Reactive value to store data
  data <- reactiveVal(list())
  
  # Load data when file is uploaded
  observeEvent(input$file, {
    req(input$file)
    load(input$file$datapath)
    row.names(omics_correlation_data) = NULL
    omics_correlation_data$type = sub_type
    omics_correlation_data$correlation = round(omics_correlation_data$correlation, 5)
    omics_correlation_data$p = round(omics_correlation_data$p, 5)
    omics_correlation_data = omics_correlation_data[,c(1,4,2,3)]
    colnames(omics_correlation_data) = c('name', 'type','correlation to phenotype', 'p-value')
    data(list(M = M, correlation_sub = correlation_sub, sub_type = sub_type, pc_loading = pc_loading, omics_correlation_data = omics_correlation_data, pca_x1_pc1 = pca_x1_pc1))
  })
  
  observeEvent(input$info, {
    showModal(modalDialog(
      title = "How to Use SmCCNet Visualization",
      "Welcome to the SmCCNet Visualization tool! Here's how to get started:",
      tags$ul(
        tags$li("First, ensure you have generated the '.Rdata' file containing -omics subnetwork information using the SmCCNet package."),
        tags$li("Network Visualization: Explore the single/multi-omics interactive network plots by adjusting the correlation cut-off, network layout, and visual properties."),
        tags$li("PC Loading Visualization: View the principal component (PC) loading plots. Select the PC of interest to see its loadings."),
        tags$li("Omics Correlation Data: Inspect the table of omics correlation data, including correlation to phenotype and p-values."),
        tags$li("Correlation Heatmap: Analyze the correlation heatmap of network features."),
        tags$li("PCA 3D Scatter Plot: Discover patterns in 3D space by viewing PCA results, color-coded by phenotype."),
        tags$li("Upload .Rdata file: Click 'Upload .Rdata file' to upload the data generated by SmCCNet. The file should include necessary variables for visualization."),
        tags$li("Download Plots: Downloadable plots are available for network visualization, PC loadings, and the correlation heatmap."),
        tags$li("For detailed instructions on using SmCCNet to generate the necessary '.Rdata' file, please refer to the SmCCNet package documentation.")
      ),
      size = "l",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Network plotting
  latestPlotFile <- tempfile(fileext = ".svg")
  output$network_plot <- renderPlot({
    req(data())
    if(!is.null(data()$M) && !is.null(data()$correlation_sub) && !is.null(data()$sub_type)){
      plotMultiOmicsNetwork(data()$M, data()$correlation_sub, sub_type = data()$sub_type,
                            CorrCutOff = input$CorrCutOff,
                            NetLayout = input$NetLayout,
                            VertexLabelCex = input$VertexLabelCex,
                            VertexSize = input$VertexSize,
                            EdgeColorIntensity = (input$EdgeColorIntensity))
      svg(file = latestPlotFile)  # Adjust size as needed
      plotMultiOmicsNetwork(data()$M, data()$correlation_sub, sub_type = data()$sub_type,
                            CorrCutOff = input$CorrCutOff,
                            NetLayout = input$NetLayout,
                            VertexLabelCex = input$VertexLabelCex,
                            VertexSize = input$VertexSize,
                            EdgeColorIntensity = (input$EdgeColorIntensity))
      dev.off()
    }
  })
  
  # PC Loading Visualization
  output$pc_loading_plot <- renderPlot({
    req(data()$pc_loading)
    pc_num <- match(input$pc_choice, c("PC1", "PC2", "PC3"))
    df <- data.frame(name = row.names(data()$pc_loading), loading = abs(data()$pc_loading[,pc_num]))
    ggplot(df, aes(x = reorder(name, loading), y = loading)) +
      geom_bar(stat = "identity", fill = "#4DAF4A") +
      coord_flip() +
      labs(x = "Features", y = "Loading", title = paste("Loadings for", input$pc_choice)) +
      theme_minimal() +
      theme(
        text = element_text(size = 14),  # Default text size for all text elements
        axis.title = element_text(size = 14),  # Font size for axis titles
        axis.text = element_text(size = 16),  # Font size for axis text (ticks)
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Font size and style for title
      )
  })
  
  output$omics_correlation_table <- renderDT({
    req(data()$omics_correlation_data)
    DT::datatable(data()$omics_correlation_data, options = list(pageLength = 10))
  }, server = FALSE)
  
  # correlation heatmap output
  output$heatmap_plot <- renderPlot({
    req(data()$correlation_sub)
    melted_cormat <- melt(data()$correlation_sub)
    ggplot(data = melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0, limits = c(-1, 1), space = "Lab", name="Correlation") +
      labs(title = "Correlation Heatmap of Network Features", x = "", y = "") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),  # Increase title size and center it
        axis.title = element_text(size = 14),  # Increase axis titles size (x and y)
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),  # Adjust x-axis text and increase size
        axis.text.y = element_text(angle = 0, size = 12),  # Adjust y-axis text size
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)   # Increase legend text size
      )
  })
  
  output$pca_plot <- renderPlotly({
    req(data()$pca_x1_pc1) # Make sure pca_x1_pc1 is loaded and available
    
    pca_data <- data()$pca_x1_pc1
    
    plot_ly(pca_data, x = ~pc1, y = ~pc2, z = ~pc3, type = 'scatter3d', mode = 'markers',
            marker = list(size = 5, opacity = 0.5, color = ~y, colorscale = 'Viridis'),
            color = ~y) %>%
      layout(title = "3D Scatter Plot of PCs (Color-coded by Phenotype)",
             scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')),
             legend = list(title = list(text = 'Phenotype')),
             autosize = FALSE,  # Disable autosizing
             width = 800,  # Set width to a larger value
             height = 600   # Set height to a larger value
      ) %>%
      colorbar(title = "Phenotype")
  })
  
  # Download plot handler remains unchanged
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("network-plot-", Sys.Date(), ".svg", sep="")
    },
    content = function(file) {
      file.copy(latestPlotFile, file)
    }
  )
  
  
  
  output$download_pc_plot <- downloadHandler(
    filename = function() {
      paste("PC_loading_plot_", input$pc_choice, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(data()$pc_loading)
      pc_num <- match(input$pc_choice, c("PC1", "PC2", "PC3"))
      df <- data.frame(name = row.names(data()$pc_loading), loading = abs(data()$pc_loading[,pc_num]))
      
      # Begin PDF plotting
      pdf(file, width = 8, height = 6)
      loading_bar <- ggplot(df, aes(x = reorder(name, loading), y = loading)) +
        geom_bar(stat = "identity", fill = "#4DAF4A") +
        coord_flip() +
        labs(x = "Features", y = "Loading", title = paste("Loadings for", input$pc_choice)) +
        theme_minimal() +
        theme(
          text = element_text(size = 14),  # Default text size for all text elements
          axis.title = element_text(size = 14),  # Font size for axis titles
          axis.text = element_text(size = 16),  # Font size for axis text (ticks)
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Font size and style for title
        )
      print(loading_bar)
      dev.off()  # Close the PDF device
    }
  )
  
  # For downloading the heatmap as PDF
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("correlation_heatmap", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(data()$correlation_sub)
      melted_cormat <- melt(data()$correlation_sub)
      # Begin PDF plotting
      pdf(file, width = 11, height = 9)
      heatmap_plot <- ggplot(data = melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0, limits = c(-1, 1), space = "Lab", name="Correlation") +
        labs(title = "Correlation Heatmap of Network Features", x = "", y = "") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 20),  # Increase title size and center it
          axis.title = element_text(size = 14),  # Increase axis titles size (x and y)
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),  # Adjust x-axis text and increase size
          axis.text.y = element_text(angle = 0, size = 16),  # Adjust y-axis text size
          legend.title = element_text(size = 14),  # Increase legend title size
          legend.text = element_text(size = 12)   # Increase legend text size
        )
      print(heatmap_plot)
      dev.off()  # Close the PDF device
    }
  )
  # download omics-phenotype correlation table
  output$download_table <- downloadHandler(
    filename = function() {
      paste("omics-correlation-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data()$omics_correlation_data)  # Ensure the data is available
      write.csv(data()$omics_correlation_data, file, row.names = FALSE)  # Write the table to the CSV file
    }
  )
  
}

# Launch the app
shinyApp(ui, server)

