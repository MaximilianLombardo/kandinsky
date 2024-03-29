#observe()
vln_ctrl_row_ui <- fluidRow(
    column(
        width = 3,
        selectizeInput("vlnFeatures",
                       label = "QC Violin Plot Features",
                       choices = NULL,#colnames(obj()@meta.data),
                       selected = NULL,#c("nCount_RNA"),
                       multiple = TRUE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("vlnSplit",
                       label = "QC Violin Plot Grouping Variable",
                       choices = NULL,#colnames(obj()@meta.data),
                       selected = NULL,#c("orig.ident"),
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
                       choices = NULL,#names(obj()@reductions),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 2,
        selectizeInput("dimPlotGroup",
                       label = "Grouping Variable",
                       choices = NULL,#colnames(obj()@meta.data),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 2,
        selectizeInput("dimPlotSplit",
                       label = "Splitting Variable",
                       choices = NULL,#colnames(obj()@meta.data),
                       multiple = FALSE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("featurePlotFeature",
                       label = "Gene to Plot",
                       choices = NULL,#rownames(obj()[["RNA"]]@counts),#Change to handle other slots later
                       selected = NULL,#VariableFeatures(obj())[1],
                       multiple = TRUE,
                       width = "100%")),
    column(
        width = 3,
        selectizeInput("featurePlotSlot",
                       label = "Normalization Method",
                       choices = NULL,#names(obj()@assays),#Change to handle other slots later
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
#############################################################
chord_ctrl_row_ui <- fluidRow(
    column(
        width = 3,
        selectizeInput("clusterNumber",
                       label = "Signalling Cluster",
                       choices = NULL,#levels(obj()@active.ident),
                       selected = c("01"),
                       multiple = FALSE,
                       width = "100%"))#,
    # column(
    #     width = 3,
    #     selectizeInput("vlnSplit",
    #                    label = "QC Violin Plot Grouping Variable",
    #                    choices = colnames(obj()@meta.data),
    #                    selected = c("orig.ident"),
    #                    multiple = FALSE,
    #                    width = "100%"))
)

chord_row_ui <- fluidRow(
    column(
        width = 12,
        plotOutput("chordDiagram", height = "800px"))
)

chord_table_ui <- fluidRow(
    column(
        width = 12,
        dataTableOutput("chordTable"))
)
####################################################3
umap_3d_ui <- fluidRow(
    column(
        width = 12,
        plotly::plotlyOutput("UMAP3D", height = "800px")
    )
)


###########################
#Hex Selection Plot UI Elements

hex_select_feature_control_ui <- fluidRow(
    column(width = 2),
    column(width = 2),
    column(
        width = 2, align = "center",
        checkboxInput("hexFeatureBool",
                      label = "Plot Features?",
                      value = FALSE))
)

hex_select_control_ui <- fluidRow(
    column(
        width = 2,
        selectizeInput("hexSelectReduction",
                       label = "Hex Selection Reduction",
                       choices = NULL,#names(obj()@reductions),
                       selected = NULL,#"umap",
                       multiple = FALSE,
                       width = "100%")),
    column(width = 2),
    column(
        width = 2,
        selectizeInput("hexSelectFeature",
                       label = "Hex Selection Feature",
                       choices = NULL,#rownames(obj()[["RNA"]]@scale.data),
                       selected = NULL,
                       multiple = FALSE,
                       width = "100%"))
)

hex_select_ui <- fluidRow(
    column(
        width = 6,
        plotly::plotlyOutput("hexSelectionPlot", height = "600px")
    )
)
