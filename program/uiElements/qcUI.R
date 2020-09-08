observe()
vln_ctrl_row_ui <- fluidRow(
    column(
        width = 3,
        selectizeInput("vlnFeatures",
                       label = "QC Violin Plot Features",
                       choices = colnames(obj()@meta.data),
                       selected = c("nCount_RNA"),
                       multiple = TRUE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("vlnSplit",
                       label = "QC Violin Plot Grouping Variable",
                       choices = colnames(obj()@meta.data),
                       selected = c("orig.ident"),
                       multiple = FALSE,
                       width = "100%"))
)
################################################
dimPlot_label_ui <- fluidRow(
    column(
        width = 6,
        align = "center",
        textOutput("DimPlotLabel"), style = "font-size:30px;font-style:bold"),
    column(
        width = 6,
        align = "center",
        textOutput("FeaturePlotLabel"), style = "font-size:30px;font-style:bold")
)

dimPlot_ctrl_row_ui <- fluidRow(
    column(
        width = 2,
        selectizeInput("dimPlotReduction",
                       label = "Dim Plot Reduction",
                       choices = names(obj()@reductions),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 2,
        selectizeInput("dimPlotGroup",
                       label = "Grouping Variable",
                       choices = colnames(obj()@meta.data),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 2,
        selectizeInput("dimPlotSplit",
                       label = "Splitting Variable",
                       choices = colnames(obj()@meta.data),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("featurePlotFeature",
                       label = "Gene to Plot",
                       choices = NULL,#rownames(obj()[["RNA"]]@counts),#Change to handle other slots later
                       selected = VariableFeatures(obj())[1],
                       multiple = TRUE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("featurePlotSlot",
                       label = "Normalization Method",
                       choices = names(obj()@assays),#Change to handle other slots later
                       selected = c("RNA"),
                       multiple = FALSE,
                       width = "100%"))
)

dimpPlot_row_ui <- fluidRow(
    column(
        width = 6,
        plotOutput("seuratDimPlot",
                   width = "100%",
                   height = "400px")),
    column(
        width = 6,
        plotOutput("seuratFeaturePlot",
                   width = "100%",
                   height = "400px"))
)

dimPlotMarkerTable_ui <- fluidRow(
    column(
        width = 12,
        dataTableOutput('markerTable'))
)
