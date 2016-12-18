HEATMAPS <- function(db, 
                     db.na, 
                     cols = NULL,
                     na.col,
                     go.term,
                     name,
                     p.value,
                     line.size = LS,
                     Dendro = TRUE,
                     p.value.cutoff = 0.05){
  
  
  dist <- dist(db)
  dd.col <- as.dendrogram(hclust(dist))
  col.ord <- order.dendrogram(dd.col)
  
  DB2 <- DB.na[col.ord,]
  
  m.DB <- melt(DB2)
  
  ddata_y <- dendro_data(dd.col, 
                         labels = dd.col$labels,
                         type = 'rectangle')
  
  theme_none <- theme(
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_text(colour=NA),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = NA, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(colour = NA),
    plot.margin = unit(c(1,1,1,1), 'lines')
  )
  
  if('white' %in% cols){
    resc <- rescale(x = 0,
                    to = c(0,1),
                    from = range(m.DB$value, na.rm = TRUE))
    
    wh <- which(cols == 'white')
    
    n.bef <- wh - 1
    
    n.aft <- length(cols) - wh
    
    seq.bef <- seq(from = 0, to = resc, length.out = n.bef)
    seq.aft <- seq(from = resc, to = 1, length.out = n.aft + 1)[-1]
    
    vls <- c(seq.bef, resc, seq.aft)
  } else {
    vls <- NULL
  }
  
  
  if(check.for.significance == 'y'){
    ggtit <- paste('log',logbase,'(fold-change). ', 
                   na.col, ' indicates non-significant fold change (p-value > ',
                   p.value.cutoff, ')','\n', sep = '')
  } else {
    ggtit <- paste('log', logbase, '(fold-change).', sep = '')
  }
  
  hmap <- ggplot(m.DB, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = value), colour = "black", size = line.size) +
    scale_x_discrete('', expand = c(0,0)) + 
    scale_y_discrete('', expand = c(0,0)) +
    scale_fill_gradientn(name = '',
                         colours = cols,
                         values = vls, 
                         na.value = na.col) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.margin = unit(c(1,1,1,5), 'lines')) 
  
  panel_height <- unit(0.75,"npc") - sum(ggplotGrob(hmap)[["heights"]][-3]) - unit(0.5,"line")
  
  hmap <- hmap + 
    theme(# legend.justification = 'left',
      legend.margin = unit(1, 'lines'),    
      legend.position = 'left',
      legend.key.width = unit(2, 'lines'),
      legend.direction = 'vertical',
      legend.title = element_blank(
        #             size = 18, 
        #                                       face = "italic",
        #                                       angle = 45,
        #                                       hjust = 1
      )
      ) + 
    guides(fill= guide_colorbar(barheight=panel_height))
  
  gbuildHmap <- ggplot_build(ggdraw(switch_axis_position(hmap, 'y')))
  
  
  dendro <- ggplot(segment(ddata_y)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none +
    ## ggtitle(paste('log',logbase,'(fold-change)\n', sep = '')) + 
    scale_x_continuous('', expand = c(1/(2*nrow(m.DB)/ncol(m.DB)),0)) +
    scale_y_continuous(breaks = seq(from = min(segment(ddata_y)$y), 
                                    to = max(segment(ddata_y)$y), 
                                    length.out = ncol(db)),
                       labels = colnames(db))
  
  gHmap <- ggplot_gtable(gbuildHmap) 
  gDendro <- ggplot_gtable(ggplot_build(dendro))
  
  gHmap$vp <- viewport(0.7, 0.75, x = 0.35, y = 0.425)
  gDendro$vp <- viewport(0.3, 0.75, x = 0.8, y = 0.425)
  
  Title <- textGrob(label = paste(paste(name, collapse = ' '),
                                  go.term,
                                  paste('p-value:', p.value, sep = ' '),
                                  sep = '\n'),
                    just = 'center',
                    gp = gpar(fontsize = 24, fontface = 'bold'),
                    vp = viewport(0.5, 0.2, x = 0.5, y = 0.9))
  subtitle <- textGrob(label = ggtit,
                       just = 'center',
                       gp = gpar(fontsize = 16, fontface = 'bold'),
                       vp = viewport(0.5, 0.05, x = 0.5, y = 0.825))
  
  grid.draw(gHmap)
  grid.draw(gDendro)
  grid.draw(Title)
  grid.draw(subtitle)
  
}

## Function that will let us choose a folder for heatmaps
choose.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}

## Function to do the log-to-fold transformation
logFC <- function(db,mutants,WT,logBase=2,pseudo=1) {
  if (length(WT) !=1 ) {
    stop('WT must refer to a single gene/column')
  }
  if (is.numeric(logBase)==FALSE) {
    stop('logBase must be a numeric value')
  }
  
  dbb <- (db[,mutants]+pseudo)/(db[,WT]+pseudo)
  dbb <- log(dbb,logBase)
  
  names(dbb) <- rownames(db)
  
  return(dbb)
}
