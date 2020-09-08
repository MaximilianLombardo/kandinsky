createCustomChordDiagram <- function(obj,
                                     cluster.number,
                                     links,
                                     background.color = "white"
                                     #thresh = 10
                                     ){
    ########################FUNCTION REQUIREMENTS###############
    require(circlize)
    require(Seurat)
    require(RColorBrewer)
    require(grDevices)
    require(d3heatmap)
    
    ########################ADD NODOSE ANNOTATIONS###############################
    # if(is.null(nodose.annotations)){
    #     
    #     nod.clusters <- mnod.cluster.nodes <- paste0(rep("MNOD", 26), formatC(0:25, digits = 1, format = "d", flag = "0"))
    #     nodose.annotations <- factor(x = character(length(nod.clusters)),
    #                                  levels = c("gut.innervating.nodose",
    #                                             "non.gut.innervating.nodose",
    #                                             "jugular"))
    #     names(nodose.annotations) <- nod.clusters
    #     
    #     
    #     nod.complete <- formatC(c(0:25), digits = 1, format = "d", flag = "0")
    #     # Annotate the nodose clusters that receive input from the gut
    #     gut.nods <- formatC(c(1,3,4,5,7,8,9,11,17,20), digits = 1, format = "d", flag = "0")#gut.nods <- formatC(c(1,3,4,8,9,11), digits = 1, format = "d", flag = "0")
    #     # Annotate Jugular clusters
    #     jugular <- formatC(c(14,15,22,23,24), digits = 1, format = "d", flag = "0")
    #     # Annotate non-gut innervating nodose clusters
    #     non.gut.nods <- setdiff(nod.complete, c(gut.nods, jugular))
    #     
    #     
    #     #create final annotations to be used with appropriate cluster names
    #     nodose.annotations[nod.clusters[nod.complete %in% gut.nods]] <- "gut.innervating.nodose"
    #     nodose.annotations[nod.clusters[nod.complete %in% non.gut.nods]] <- "non.gut.innervating.nodose"
    #     nodose.annotations[nod.clusters[nod.complete %in% jugular]] <- "jugular"
    #     
    # }else{
    #     nodose.annotations <- nodose.annotations
    # }
    
    
    ########################GET APPROPRIATE LINKS FROM COMPLETE LINKS AND CREATE APPROPRIATE SECTORS##################################
    links <- links[links$from.cluster %in% formatC(cluster.number, digits = 1, format = "d", flag = "0"),]
    links <- links[order(links$from.gene),]
    sectors <- unique(c(links$from.cluster, links$to.cluster))
    sectors <- factor(sectors, levels = sectors)#Preserve ordering when turning into a factor
    
    
    ########################SET CHART PARAMETERS AND INITIALIZE THE CHART/CIRCLE##########################
    ##Initializing parameters
    main.xlim <- 0.3#Size of the main signalling segment. parameterize?
    
    xrange <- c(0, seq(main.xlim, 1, by = (1 - main.xlim)/(length(sectors)-1)))
    xlim.data <- cbind(xrange[1:length(xrange)-1], xrange[2:length(xrange)])
    rownames(xlim.data) <- sectors
    
    small_gap <- 1
    big_gap <- 20
    
    if(background.color == "black"){
        par(bg = "black")
        border.color <- "black"
        dark.opacity <- 1
        light.opacity <- 1
    }else{
        par(bg = "white")
        border.color <- "white"
        dark.opacity <- 0.8
        light.opacity <- 0.4
    }
    
    
    circos.par(gap.after = c(big_gap, rep(small_gap, length(unique(as.character(links$to.cluster))) - 1), big_gap),
               start.degree = 226)
    
    circos.initialize(factors = sectors, xlim = xlim.data)
    
    ########################CREATE COLOR PALLETTES FOR SECTORS###########################
    #Create colors to code for different sectors -- NEWEWEWEW
    highlight.color <- adjustcolor("grey40",alpha.f = light.opacity)
    send.color <- adjustcolor("dodgerblue4",alpha.f = light.opacity)
    receive.color <- adjustcolor("dodgerblue4", alpha.f = light.opacity)
    sector.colors <- c(send.color, rep(receive.color, nrow(links)))
    
    ########################ADD HIGHLIGHT TRACKS TO GROUP SECTORS#######################
    #Add a grouping track to highlight
    circos.track(sectors, ylim = c(0, 0.2), track.height = 0.05, bg.border = border.color)#Changed to Black
    
    #Highlight the signalling (sending) cluster
    highlight.sector(sector.index = as.character(sectors[1]), track.index = 1, col = highlight.color, 
                     text = "signalling.cluster", text.vjust = -1, niceFacing = TRUE, facing = "bending.inside")
    
    #Highlight the receiving clusters
    highlight.sector(sector.index = sectors[-1], track.index = 1, col = highlight.color, 
                     text = "receiving.clusters", text.vjust = -1, niceFacing = TRUE, facing = "bending.inside")
    
    
    ########################@TODO ADD EXPRESSION LEVEL BAR /LINE CHART######################################################################################
    # ##Add a track which has line plots of the expression levels
    #Hormones for sending cluster and receptor expression levels for receiving cluster 
    
    #calculate average expression
    avg.exp <- AverageExpression(obj,
                                 features = unique(c(links$from.gene,links$to.gene)),
                                 assays = "RNA")#Make this a parameter later?
    avg.exp <- avg.exp$RNA
    colnames(avg.exp) <- formatC(as.integer(colnames(avg.exp)) + 1,
                                 digits = 1, format = "d", flag = "0")
    
    send.cluster.exp <- (avg.exp[order(unique(links$from.gene)), unique(links$from.cluster)])
    names(send.cluster.exp) <- unique(links$from.gene)#ensure we don't lose the naming of this vector
    receieve.cluster.exp <- (avg.exp[order(unique(links$to.gene)), unique(links$to.cluster)])
    
    
    
    ###
    precursor.links.colors <- list()
    precursors <- unique(as.character(links$from.gene))
    accent.pallete <- adjustcolor(brewer.pal(n = length(precursors), name = "PRGn"), alpha.f = dark.opacity)
    for(i in 1:length(precursors)){
        precursor.links.colors[[precursors[i]]] <- accent.pallete[i]
    }
    precursor.links.colors <- unlist(precursor.links.colors)
    
    # #Palettes for sending and receiving molecules
    set.seed(12345)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #thresh.receptor.genes <- sort(unique(rownames(which(receieve.cluster.exp > thresh, arr.ind = TRUE))))
    #receptor.pallette <- sample(col_vector, length(thresh.receptor.genes))
    receptor.colors <- sample(col_vector, length(unique(links$to.gene)))
    names(receptor.colors) <- unique(links$to.gene)
    
    circos.track(sectors, ylim = c(0, 0.2), track.height = 0.15, panel.fun = function(x,y){
        sector.index <- get.cell.meta.data("sector.index")
        xrange <- get.cell.meta.data("xlim")
        xmin <- xrange[1]
        xmax <- xrange[2]
        yrange <- get.cell.meta.data("ylim")
        ymin <- yrange[1]
        ymax <- yrange[2]
         
        if(grepl(unique(links$from.cluster), sector.index)){
            accent.pallete <- adjustcolor(brewer.pal(n = length(send.cluster.exp), name = "Accent"), alpha.f = dark.opacity)
            
            send.cluster.exp.normalized <-  send.cluster.exp/sum(send.cluster.exp)
            
            x.spot <- seq(from = xmin, to = xmax, length.out = length(send.cluster.exp.normalized))
            y.spot <- send.cluster.exp.normalized * ymax
            
            circos.lines(x = x.spot, y = y.spot,
                         col = adjustcolor(precursor.links.colors[names(send.cluster.exp.normalized)], alpha.f = dark.opacity),
                         type = "h", lwd = 10)
            
        }else{
            receive.sector.exp <- receieve.cluster.exp[, sector.index]
            names(receive.sector.exp) <- rownames(receieve.cluster.exp)
            #receive.sector.exp <- receive.sector.exp[receive.sector.exp > thresh]
            receive.sector.exp.normalized <-  receive.sector.exp/sum(receive.sector.exp)
            
            x.spot <- seq(from = xmin, to = xmax, length.out = length(receive.sector.exp.normalized))
            y.spot <- receive.sector.exp.normalized * ymax
            
            #Need to do something about those receptor colors
            circos.lines(x = x.spot, y = y.spot,
                         col = "black",#adjustcolor(receptor.colors[names(receive.sector.exp.normalized)], alpha.f = dark.opacity),
                         type = "h", lwd = 3)
            
        }
       
    }, bg.border = border.color, bg.col = adjustcolor("gray94", alpha.f = dark.opacity))
    
    ########################ADD SECTOR LABELS######
    circos.track(sectors, ylim = c(0, 0.2), track.height = 0.05, panel.fun = function(x,y){
        sector.index = get.cell.meta.data("sector.index")
        xcenter = get.cell.meta.data("xcenter")
        ycenter = get.cell.meta.data("ycenter")
        circos.text(xcenter, ycenter, sector.index, cex = 0.7)
    },bg.col = sector.colors ,bg.border = border.color)
    
    ########################ADD LINKS############################################################################################
    #Create a list with colors which represent different precursor genes
    precursor.links.colors <- list()
    precursors <- unique(as.character(links$from.gene))
    accent.pallete <- adjustcolor(brewer.pal(n = length(precursors), name = "PRGn"), alpha.f = dark.opacity)
    for(i in 1:length(precursors)){
        precursor.links.colors[[precursors[i]]] <- accent.pallete[i]
    }
    
    ##Link start and end points
    #Create an object which can dynamically keep track of a link's starting and ending points
    node.positions <- list()
    node.names <- c(unique(as.character(links$from.cluster)), unique(as.character(links$to.cluster)))
    
    #Initialized Values
    for(i in 1:length(node.names)){
        node.positions[[node.names[i]]] <- list()
        node.positions[[node.names[i]]]$link.start <- 0
        node.positions[[node.names[i]]]$link.end <- 0
    }
    
    #Add links with appropriate starting and ending points
    ###Implementation of the loop to correct for switching nodose nodes
    
    
    clusters <- unique(c(links$from.cluster,links$to.cluster))
    
    for(i in 1: nrow(links)){
        current.link <- links[i,]
        from.cluster <- current.link$from.cluster
        to.cluster <- current.link$to.cluster
        
        
        from.lim <- get.cell.meta.data("xlim", from.cluster, 1)
        to.lim <- get.cell.meta.data("xlim", to.cluster, 1)
        from.range <- get.cell.meta.data("xrange", from.cluster, 1)
        to.range <- get.cell.meta.data("xrange", to.cluster, 1)
        
        if(from.cluster %in% clusters){
            node.positions[[from.cluster]]$link.start <- node.positions[[from.cluster]]$link.start + from.lim["min.data"]
            node.positions[[from.cluster]]$link.end <- node.positions[[from.cluster]]$link.end + (node.positions[[from.cluster]]$link.start + (from.range/sum(links$from.cluster %in% from.cluster)))
        }else{
            node.positions[[from.cluster]]$link.start <- node.positions[[from.cluster]]$link.end
            node.positions[[from.cluster]]$link.end <- node.positions[[from.cluster]]$link.start + (from.range/sum(links$from.cluster %in% from.cluster))
        }
        
        if(to.cluster %in% clusters){
            node.positions[[to.cluster]]$link.start <- node.positions[[to.cluster]]$link.start + to.lim["min.data"]
            node.positions[[to.cluster]]$link.end <- node.positions[[to.cluster]]$link.end + (node.positions[[to.cluster]]$link.start + (to.range/sum(links$to.cluster %in% to.cluster)))
        }else{
            node.positions[[to.cluster]]$link.start <- node.positions[[to.cluster]]$link.end
            node.positions[[to.cluster]]$link.end <- node.positions[[to.cluster]]$link.start + (to.range/sum(links$to.cluster %in% to.cluster))
        }
        
        circos.link(from.cluster, c(node.positions[[from.cluster]]$link.start, node.positions[[from.cluster]]$link.end),
                    to.cluster, c(node.positions[[to.cluster]]$link.start, node.positions[[to.cluster]]$link.end),
                    col = precursor.links.colors[[as.character(current.link$from.gene)]])
        
        clusters <- clusters[clusters != from.cluster]
        clusters <- clusters[clusters != to.cluster]
        
    }
    
    ########################ADD LEGEND###############################################################
    #Add legend for Peptide precursors...
    legend("bottomleft",
           legend = names(precursor.links.colors),
           col = unlist(precursor.links.colors),
           pch = c(19),
           bty = "n",
           pt.cex = 2,
           cex = 1,#1.2
           text.col = "black",
           horiz = F ,
           inset = c(0.01, 0.01))
    
    #Add legend for Receptors precursors...
    # legend("topright",
    #        legend = names(receptor.pallette),
    #        col = receptor.pallette,
    #        pch = c(19),
    #        bty = "n",
    #        pt.cex = 2,
    #        cex = 0.8,
    #        text.col = "black",
    #        horiz = F ,
    #        inset = c(0, 0))
    
    #Add legend for Convertases...
    # legend("topleft",
    #        legend = names(convertase.pallete),
    #        col = convertase.pallete,
    #        pch = c(19),
    #        bty = "n",
    #        pt.cex = 2,
    #        cex = 1.2,
    #        text.col = "black",
    #        horiz = F ,
    #        inset = c(0.01, 0.01))
    
    
    circos.clear()
    
    
}
