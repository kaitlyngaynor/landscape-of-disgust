# random useful functions/work-arounds

#calc aic
aic <- function(U, n=NULL, p, aicc=FALSE){
	if(aicc==FALSE) aic <- U + (2*p)
	if(aicc==TRUE) aic <- U + (2*p) + ((2*p*(p+1))/(n-p-1))
	return(aic)
}

# extract legend
g_legend<-function(a.gplot){ 
 	tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  	legend <- tmp$grobs[[leg]] 
  	return(legend)
} 

# ticks on all axes of plot
ticks_all_sides <- function(plot){
	library(gtable)
	library(grid)

	gt2 <- ggplotGrob(plot)
	panel <-c(subset(gt2$layout, name=="panel", se=t:r))
	rn <- which(gt2$layout$name == "axis-b")
	axis.grob <- gt2$grobs[[rn]]
	axisb <- axis.grob$children[[2]]  # Two children - get the second
	xaxis = axisb$grobs[[1]]  # NOTE: tick marks first
	xaxis$y = xaxis$y - unit(0.25, "cm")  # Position them inside the panel
	gt2 <- gtable_add_rows(gt2, unit(0, "lines"), panel$t-1)
	gt2 <- gtable_add_grob(gt2, xaxis, l = panel$l, t = panel$t, r = panel$r, name = "ticks")
	panel <-c(subset(gt2$layout, name=="panel", se=t:r))
	rn <- which(gt2$layout$name == "axis-l")
	axis.grob <- gt2$grobs[[rn]]
	axisl <- axis.grob$children[[2]]  # Two children - get the second
	yaxis = axisl$grobs[[2]] # NOTE: tick marks second
	yaxis$x = yaxis$x - unit(0.25, "cm") # Position them inside the panel
	gt2 <- gtable_add_cols(gt2, unit(0, "lines"), panel$r)
	gt2 <- gtable_add_grob(gt2, yaxis, t = panel$t, l = panel$r+1, name = "ticks")
	gt2$layout[gt2$layout$name == "ticks", ]$clip = "off"
	return(gt2)
}

#convert matrix to multinet formalism
# net should be a named list of data.frames
# type is 'bipartite' or 'square'
# actor.attributes specifies attributes other than species name (should be a list of attributes where each attribute is a matrix, first column 'spp', column 2 attribute)
# include.weights is TRUE when matrix elements are used as interaction strengths, if FALSE all are 1
# interlayer.links should be a named data.frame with columns: spp, layer1, layer2, interaction weights
########## needs to be edited so that commas are output (potentially making it a one column matrix is best)
## added functionality to output data for muxViz
multinet_layout <- function(net, file, type='bipartite', actor.attributes=NULL,
    include.weights=TRUE, interlayer.links=NULL, intra.type="UNDIRECTED",
    inter.type='UNDIRECTED', muxviz=TRUE){

    lays <- length(net)
    ml <- matrix('', ncol=5, nrow=3+lays)
    ml[1,1] <- "#TYPE multilayer"
    ml[3,1] <- "#LAYERS" 

    if(type == 'bipartite'){
        if(class(net) == 'matrix') net <- lapply(net, as.data.frame)
        
        #intralayers ints   
        ml[(3+1):(3+lays),1] <- names(net)
        ml[(3+1):(3+lays),2] <- names(net)
        ml[(3+1):(3+lays),3] <- intra.type

        #interlayer ints
        for(i in 1:(lays-1)){
            ml <- rbind(ml, matrix(c(names(net)[i], names(net)[i+1], inter.type, '', ''), ncol=5))
        }

        # add actor attributes
        ml <- rbind(ml, matrix(rep('',5), ncol=5))

        ml <- rbind(ml, matrix(c('#ACTOR ATTRIBUTES', rep('',4)), ncol=5))
        ml <- rbind(ml, matrix(c('speciesname', 'STRING', rep('',3)), ncol=5))
        if(!is.null(actor.attributes)){
            for(i in 1:length(actor.attributes)){
                ml <- rbind(ml, matrix(c(names(actor.attributes[i]), 'STRING', rep('',3)), ncol=5)) 
            }
        }

        #add actors
        spp <- c()
        for(i in 1:length(net)){
            spp <- c(spp, rownames(net[[i]]), colnames(net[[i]]))
        }
        spp <- unique(spp)
        spp <- sapply(spp, function(x) gsub(' ', '_', x)); names(spp) <- NULL
        
        ml <- rbind(ml, matrix(rep('',5), ncol=5))
        ml <- rbind(ml, matrix(c('#ACTORS', rep('',4)), ncol=5))
        
        for(i in 1:length(spp)){
            ml <- rbind(ml, matrix(c(paste('spp',i,sep=''), spp[i], rep('',3)), ncol=5))   
        }

        if(!is.null(actor.attributes)){
            for(i in 1:length(actor.attributes)){
                value <- actor.attributes[[i]]
                for(x in 1:nrow(value)){
                    if(length(grep('_', value[x,1])) < 1) value[x,1] <- gsub(' ', '_', value[x,1])
                    ml[which(ml[,2] == value[x,1]),3] <- value[x,2]
                }
            }
        }

        # add edge attributes
        ml <- rbind(ml, matrix(rep('',5), ncol=5))
        ml <- rbind(ml, matrix(c('#EDGE ATTRIBUTES', rep('',4)), ncol=5))

        for(i in 1:lays){
            ml <- rbind(ml, matrix(c(names(net)[i], names(net)[i], 'weight', 'NUMERIC', ''), ncol=5))
        }   

        for(i in 1:(lays-1)){
            ml <- rbind(ml, matrix(c(names(net)[i], names(net)[i+1], 'weight', 'NUMERIC', ''), ncol=5))
        }

        # add edges
        ml <- rbind(ml, matrix(rep('',5), ncol=5))
        ml <- rbind(ml, matrix(c('#EDGES', rep('',4)), ncol=5))

        if(include.weights==TRUE){
            for(i in 1:length(net)){
                layer <- as.matrix(net[[i]])
                ints <- which(layer > 0, arr.ind=TRUE)
                for(x in 1:nrow(ints)){
                    ml <- rbind(ml, matrix(c(
                        ml[which(ml[,2] == gsub(' ', '_', rownames(layer)[ints[x,1]])),1],
                        names(net)[i],
                        ml[which(ml[,2] == gsub(' ', '_', colnames(layer)[ints[x,2]])),1],
                        names(net)[i],
                        layer[ints[x,1], ints[x,2]]
                        ), ncol=5))
                }  
            }

        }else{
            for(i in 1:length(net)){
                layer <- as.matrix(net[[i]])
                ints <- which(layer > 0, arr.ind=TRUE)
                for(x in 1:nrow(ints)){
                    ml <- rbind(ml, matrix(c(
                        ml[which(ml[,2] == gsub(' ', '_', rownames(layer)[ints[x,1]])),1],
                        names(net)[i],
                        ml[which(ml[,2] == gsub(' ', '_', colnames(layer)[ints[x,2]])),1],
                        names(net)[i],
                        1
                        ), ncol=5))
                }  
            } 
        }

        # add interlayer edges
        if(!is.null(interlayer.links)){
            for(i in 1:(lays-1)){
                for(x in 1:length(spp)){
                    if(spp[x] %in% c(rownames(net[[i]]), colnames(net)[[i]]) & spp[x] %in% c(rownames(net[[i+1]]), colnames(net)[[i+1]])){
                       ml <- rbind(ml, matrix(c(
                            spp[x],
                            names(net)[i],
                            spp[x],
                            names(net)[i+1],
                            interlayer.links$weight[which(interlayer.links$spp == spp[x] & interlayer.links$layer1 == names(net)[i] & interlayer.links$layer2 == names(net)[i+1])]), ncol=5)) 
                    }
                }
            }
        }else{
            for(i in 1:(lays-1)){
                for(x in 1:length(spp)){
                    if(gsub('_', ' ', spp[x]) %in% c(rownames(net[[i]]), colnames(net[[i]])) & gsub('_', ' ', spp[x]) %in% c(rownames(net[[i+1]]), colnames(net[[i+1]]))) {
                       ml <- rbind(ml, matrix(c(
                            ml[which(ml[,2] == spp[x]),1],
                            names(net)[i],
                            ml[which(ml[,2] == spp[x]),1],
                            names(net)[i+1],
                            1), ncol=5)) 
                    }
                }
            }
        }

        if(muxviz == TRUE){
            wfile <- paste(strsplit(file, '/')[[1]][1], strsplit(file, '/')[[1]][2], paste('edges', strsplit(file, '/')[[1]][3], sep='_'), sep='/')
            write.table(ml[(which(ml[,1] == '#EDGES')+1):nrow(ml),], file=wfile, quote=FALSE, sep=' ', col.names=FALSE, row.names=FALSE)
        
            wfile <- paste(strsplit(file, '/')[[1]][1], strsplit(file, '/')[[1]][2], paste('layout', strsplit(file, '/')[[1]][3], sep='_'), sep='/')
            nodes <- rbind(matrix(c('nodeID', 'nodeLabel'), ncol=2), ml[(which(ml[,1] == '#ACTORS')+1):(which(ml[,1] == '#EDGE ATTRIBUTES')-2),1:2])
            write.table(nodes, file=wfile, quote=FALSE, sep=' ', col.names=FALSE, row.names=FALSE)
        
            wfile <- paste(strsplit(file, '/')[[1]][1], strsplit(file, '/')[[1]][2], paste('layers', strsplit(file, '/')[[1]][3], sep='_'), sep='/')
            layers <- rbind(c('layerID', 'layerLabel'), cbind(c(1:lays),names(net)))
            write.table(layers, file=wfile, quote=FALSE, sep=' ', col.names=FALSE, row.names=FALSE)
        }

        # write the table mate!
        if(muxviz == FALSE) write.table(ml, file=file, sep=',', col.names=FALSE, row.names=FALSE, quote=FALSE)

        }else if(type == 'square'){
            if(class(net) == 'matrix') net <- as.data.frame(net)

            stop('Only supported for bipartite nets right now, sorry dude.')
        }
}

#pull out assigned edge weights from a multilayer
edge_weight <- function(ml){
    library(multinet)

    g <- c() # place to put weights
    e <- edges.idx.ml(ml) # edges in network based on ID
    n <- nodes.ml(ml)

    for(i in 1:nrow(e)){
        edge <- cbind(n[e[i,'from'],], n[e[i,'to'],])
        g <- c(g, get.values.ml(d, 'weight', edge=edge))
    }

    return(g)
}

# plot is then multinet_plot(ml, edge.width=as.numeric(edge_weight(ml)))

#plotting multilayer network with real edge control
multinet_plot <- function (x, layout = NULL, grid = NULL, mai = 0.1, vertex.shape = 16, 
    vertex.cex = 1, vertex.color = NULL, vertex.labels = NULL, 
    vertex.labels.pos = 3, vertex.labels.offset = 0.5, vertex.labels.cex = 1, 
    edge.type = 1, edge.width = 1, edge.color = 1, edge.arrow.length = 0.1, 
    edge.arrow.angle = 20, com = NULL, com.cex = 1, ...){

	library(multinet)
    num.cols = num.layers.ml(x)
    num.rows = 1
    if (!is.null(grid)) {
        if (!length(grid) == 2) 
            stop("argument grid must have two elements")
        num.rows = grid[1]
        num.cols = grid[2]
    }
    if (is.null(layout)) {
        layout <- layout.multiforce.ml(x)
    }
    x_coord <- function(xyz_coord) {
        xyz_coord$x + xyz_coord$z%%num.cols * width
    }
    y_coord <- function(xyz_coord) {
        xyz_coord$y + (num.rows - 1 - xyz_coord$z%/%num.cols) * 
            height
    }
    x.min = min(layout$x)
    y.min = min(layout$y)
    x.max = max(layout$x)
    y.max = max(layout$y)
    width = x.max - x.min + mai * (x.max - x.min)
    x.min = x.min - mai/2 * (x.max - x.min)
    height = y.max - y.min + mai * (y.max - y.min)
    y.min = y.min - mai/2 * (y.max - y.min)
    plot(NA, type = "n", xlim = c(x.min, x.min + width * num.cols), 
        ylim = c(y.min, y.min + height * num.rows), xaxt = "n", 
        yaxt = "n", bty = "n", xlab = "", ylab = "")
    segments((0:num.cols * width) + x.min, y.min, (0:num.cols * 
        width) + x.min, y.min + height * num.rows)
    segments(x.min, (0:num.rows * height) + y.min, x.min + width * 
        num.cols, (0:num.rows * height) + y.min)
    if (!is.null(com) && nrow(com) > 0) {
        num.com <- max(com$cid) + 1
        palette = rainbow(num.com, alpha = 0.5)
        draw.areas <- function(d) {
            xc <- x_coord(layout[d$aid, ])
            yc <- y_coord(layout[d$aid, ])
            os <- par()$cxy * com.cex
            xc <- c(xc + os[1]/2, xc + os[1]/2, xc - os[1]/2, 
                xc - os[1]/2)
            yc <- c(yc + os[2]/2, yc - os[2]/2, yc + os[2]/2, 
                yc - os[2]/2)
            extreme.points = chull(xc, yc)
            xspline(xc[extreme.points], yc[extreme.points], open = F, 
                shape = 1, border = NA, col = palette[d$cid + 
                  1])
        }
        c.list <- get.community.list.ml(com, x)
        lapply(c.list, draw.areas)
    }
    e <- edges.idx.ml(x)
    if(length(edge.width) == nrow(e)) e$weight <- edge.width #make sure they are the right order
    draw_edge <- function(d) {
        if (d["dir"] == 0) {
            segments(x_coord(layout[as.numeric(d["from"]), ]), y_coord(layout[as.numeric(d["from"]), 
                ]), x_coord(layout[as.numeric(d["to"]), ]), y_coord(layout[as.numeric(d["to"]), 
                ]), lty = edge.type, lwd = as.numeric(d["weight"]), col = edge.color)
        }
        if (d["dir"] == 1) {
            arrows(x_coord(layout[as.numeric(d["from"]), ]), y_coord(layout[as.numeric(d["from"]), 
                ]), x_coord(layout[as.numeric(d["to"]), ]), y_coord(layout[as.numeric(d["to"]), 
                ]), length = edge.arrow.length, angle = edge.arrow.angle, 
                lty = edge.type, lwd = as.numeric(d["weight"]), col = edge.color)
        }
    }
    apply(e, 1, draw_edge)
    if (is.null(vertex.color)) 
        vertex.color = layout$z + 1
    points(x_coord(layout), y_coord(layout), pch = vertex.shape, 
        col = vertex.color, cex = vertex.cex)
    if (is.null(vertex.labels)) 
        vertex.labels = layout$actor
    text(x_coord(layout), y_coord(layout), labels = vertex.labels, 
        pos = vertex.labels.pos, offset = vertex.labels.offset, 
        cex = vertex.labels.cex)
}

# extracts coordinates from kml files
fuck_you_kml <- function(file){
    f <- readLines(file)
    names=grep('SimpleData name=', f)
    points=grep('<Point><coordinates>', f)
    df <- matrix(NA, nrow=length(names), ncol=3)
    colnames(df) <- c('ID', 'Latitude', 'Longitude')
    for(i in 1:length(names)){
        df[i,1] <- strsplit(substr(f[names[i]], 25, 35), '<')[[1]][1]
        df[i,2] <- as.numeric(strsplit(strsplit(f[points[i]], ',')[[1]][2], '<')[[1]][1])
        df[i,3] <- as.numeric(strsplit(strsplit(f[points[i]], ',')[[1]][1], '>')[[1]][3])
    }
    df <- as.data.frame(df)
    return(df)
}



