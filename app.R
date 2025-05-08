library(shiny)
library(Seurat)
library(ggplot2)
library(ggnewscale)

ui <- fluidPage(
  titlePanel("Seurat UMAP Expression Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Dataset:",
                  choices = list(
                    "WT Early" = "data/wt_early_salivary_slim.rds",
                    "WT Late" = "data/wt_late_salivary_slim.rds"
                  )
      ),
      textInput("gene1", "Gene 1:", value = "Dr"),
      textInput("gene2", "Gene 2:", value = "trbl")
    ),
    mainPanel(
      plotOutput(
        "overlayPlot", height = "600px",
        brush = brushOpts(id = "plot_brush", resetOnNew = TRUE),
        dblclick = "plot_dblclick"
      ),
      plotOutput(
        "scatterPlot", height = "600px",
        brush = brushOpts(id = "scatter_brush", resetOnNew = TRUE),
        dblclick = "scatter_dblclick"
      ),
      verbatimTextOutput("stats")
    )
  )
)

server <- function(input, output, session) {
  seurat_obj <- reactive({
    req(input$dataset)
    readRDS(input$dataset)
  })
  wt_sg <- reactive({ subset(seurat_obj(), manual_celltypes == "Salivary Gland") })
  
  data_df <- reactive({
    obj <- wt_sg()
    coords <- Embeddings(obj, "umap")
    df <- as.data.frame(coords)
    df$gene1 <- FetchData(obj, vars = input$gene1)[,1]
    df$gene2 <- FetchData(obj, vars = input$gene2)[,1]
    df
  })
  
  stats <- reactive({
    df <- data_df(); n <- nrow(df)
    pct1 <- round(sum(df$gene1 > 0, na.rm = TRUE) / n * 100, 2)
    pct2 <- round(sum(df$gene2 > 0, na.rm = TRUE) / n * 100, 2)
    pctBoth <- round(sum(df$gene1 > 0 & df$gene2 > 0, na.rm = TRUE) / n * 100, 2)
    test <- cor.test(df$gene1, df$gene2, method = "pearson", use = "complete.obs")
    list(
      pct1 = pct1,
      pct2 = pct2,
      pctBoth = pctBoth,
      r = round(test$estimate[[1]], 2),
      p = signif(test$p.value, 3)
    )
  })
  
  # Reactive ranges for zooming UMAP plot
  ranges_overlay <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$plot_brush, {
    brush <- input$plot_brush
    ranges_overlay$x <- c(brush$xmin, brush$xmax)
    ranges_overlay$y <- c(brush$ymin, brush$ymax)
  })
  observeEvent(input$plot_dblclick, {
    ranges_overlay$x <- NULL
    ranges_overlay$y <- NULL
  })
  
  output$overlayPlot <- renderPlot({
    df <- data_df()
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = gene1), size = 5, alpha = 0.8) +
      scale_color_gradient(low = "gray", high = "red", name = input$gene1) +
      ggnewscale::new_scale_color() +
      geom_point(aes(color = gene2), size = 3, alpha = 0.6) +
      scale_color_gradient(low = "gray", high = "blue", name = input$gene2) +
      theme_classic() +
      labs(title = paste0("Overlay of ", input$gene1, " & ", input$gene2, " (Salivary Gland)"))
    if (!is.null(ranges_overlay$x)) {
      p <- p + coord_cartesian(xlim = ranges_overlay$x, ylim = ranges_overlay$y)
    }
    p
  })
  
  # Reactive ranges for zooming scatter plot
  ranges_scatter <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$scatter_brush, {
    brush <- input$scatter_brush
    ranges_scatter$x <- c(brush$xmin, brush$xmax)
    ranges_scatter$y <- c(brush$ymin, brush$ymax)
  })
  observeEvent(input$scatter_dblclick, {
    ranges_scatter$x <- NULL
    ranges_scatter$y <- NULL
  })
  
  output$scatterPlot <- renderPlot({
    df <- data_df()
    st <- stats()
    p2 <- ggplot(df, aes(x = gene1, y = gene2)) +
      geom_point(size = 4, alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      theme_classic() +
      labs(
        x = input$gene1,
        y = input$gene2,
        title = paste0("Correlation between ", input$gene1, " & ", input$gene2),
        subtitle = paste0("Pearson r = ", st$r, ", p = ", st$p)
      )
    if (!is.null(ranges_scatter$x)) {
      p2 <- p2 + coord_cartesian(xlim = ranges_scatter$x, ylim = ranges_scatter$y)
    }
    p2
  })
  
  output$stats <- renderPrint({
    st <- stats()
    cat(input$gene1, "> 0 in", st$pct1, "% cells\n")
    cat(input$gene2, "> 0 in", st$pct2, "% cells\n")
    cat("Both > 0 in", st$pctBoth, "% cells\n")
  })
}

shinyApp(ui, server)
