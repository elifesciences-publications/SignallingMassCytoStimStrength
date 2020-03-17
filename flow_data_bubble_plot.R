## function to take tables of %+ and MFI from flow/cytof data and plot in bubble plot
## author: Arianne Richard

## color.bar function adapted from 
## https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  ## makes color bar for bubble plot
  ## lut is output of colorRampPalette
  ## min and max are min and max values for color bar
  ## nticks is how many steps in color bar
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, labels=round(ticks, 2), las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

color.gradient <- function(x, minmax, colors=viridis(5), colsteps=100) {
  ## adapted from 
  ## https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  ## generates a color gradient across values in x
  ## minmax is a 2-value vector specifying the minimum and maximum values in the color scale
  suppressWarnings(
    if(!is.na(minmax)){
      return(colorRampPalette(colors) (colsteps) [findInterval(x, seq(minmax[1], minmax[2], 
                                                                      length.out=colsteps), all.inside=TRUE)])
    } else {
      return(colorRampPalette(colors) (colsteps) [findInterval(x, seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE), 
                                                                      length.out=colsteps), all.inside=TRUE)])
    }
  )
}

bubble_plot <- function(percent_table, MFI_table, scale_MFI=TRUE, plot_height=10, plot_width=10, title=NA, color_scale=NA) {
  # percent table is a data.frame of %+ for each marker in each sample - columns are markers, rows are samples
  # MFI table is a data.frame of MFI for the positive fraction of each marker in each sample - columns are markers, rows are samples
  # MFI scale refers to scaling the MFIs by marker
  # color_scale is a vector of length 2 delimiting the minimum and maximum of the color scale
  
  ## check same markers in each table and reorder if needed
  stopifnot(colnames(percent_table) %in% colnames(MFI_table))
  stopifnot(colnames(MFI_table) %in% colnames(percent_table))
  stopifnot(rownames(percent_table) %in% rownames(MFI_table))
  stopifnot(rownames(MFI_table) %in% rownames(percent_table))
  
  ## reorder the tables so the columns and rows in same order
  percent_table <- percent_table[rownames(MFI_table),colnames(MFI_table)]
  
  ## figure out the maximum size a circle could be
  ## this may need to be adjusted with more markers added or a bigger plot made
  max_height <- (plot_height-1.5)/nrow(percent_table)
  max_width <- (plot_width-1.5)/ncol(percent_table)
  max_circle <- min(max_height, max_width)/3.2
  
  ## log transform the MFI table
  MFI_table <- log10(MFI_table)
  
  ## scale the MFIs if desired - this centers the MFIs for each marker
  ## and gives them a standard deviation of 1
  ## Scaling means that you cannot compare colors between makers but all 
  ## measurements of the same marker are colored according to their relative values
  if(scale_MFI){
    MFI_table <- apply(MFI_table, 2, function(x){scale(x, center=TRUE, scale=TRUE)})
    print(MFI_table)
  }
  
  require(viridis)
  MFIs <- unlist(MFI_table)
  ## testing color gradient function
  print(color.gradient(MFIs, minmax=color_scale))
  
  ## plot circles so that area corresponds to percent of positive cells
  ## and color corresponds to the MFI of the positive fraction
  if(!is.na(title)){
    pdf(title, height=plot_height, width=plot_width)
    layout(matrix(c(1,2,3), nrow=1, ncol=3), widths=c(5,1.25, 1.25))
  }
  par(las=2, cex.axis=0.8, cex.lab=0.9, mgp=c(3.5,0.5,0), mar=c(5,5,2,1)+0.1)
  symbols(x=rep(c(1:ncol(percent_table)), each=nrow(percent_table)), 
          y=rep(c(1:nrow(percent_table)), times=ncol(percent_table)), 
          circles=unlist(sqrt(percent_table+1)), 
          fg=color.gradient(MFIs, minmax=color_scale),
          bg=color.gradient(MFIs, minmax=color_scale), 
          inches=max_circle,
          xaxt="n",
          yaxt="n",
          xlab='Signalling proteins',
          ylab='Samples')
  axis(1, at=c(1:ncol(percent_table)), labels=colnames(percent_table))
  axis(2, at=c(1:nrow(percent_table)), labels=rownames(percent_table))
  
  ## add color bar
  par(las=2, cex.axis=0.8, cex.lab=0.8, mgp=c(3.5,0.5,0), mar=c(5,2,3,1)+0.1, cex.main=0.9, font.main=1)
  suppressWarnings(if(!is.na(color_scale)){
    color.bar(colorRampPalette(viridis(5))(100), 
              color_scale[1], color_scale[2], 
              title=paste('scaled MFI of', 'positive fraction', sep='\n')) 
  } else {
    color.bar(colorRampPalette(viridis(5))(100), 
              min(MFIs, na.rm=TRUE), max(MFIs, na.rm=TRUE), 
              title=paste('scaled MFI of', 'positive fraction', sep='\n'))
  }
  )
  
  ## add circle key
  par(las=2, cex.axis=0.8, cex.lab=0.8, mgp=c(3.5,0.5,0), mar=c(5,2,3,1)+0.1, cex.main=0.9, font.main=1, bty='n')
  symbols(x=rep(1, times=nrow(percent_table)), 
          y=c(1:nrow(percent_table)), 
          circles=sqrt(seq(from=1, to=100, length.out=nrow(percent_table))), 
          fg='black',
          inches=max_circle,
          xaxt="n",
          yaxt="n",
          xlab='',
          ylab='',
          main=paste('percent of cells', 'in positive fraction', sep='\n'))
  axis(2, at=c(1:nrow(percent_table)), labels=round(seq(from=1, to=100, length.out=nrow(percent_table)), 2))
  
  if(!is.na(title)){
    dev.off()
  }
}

add.spaces <- function(data_table, rs=NA){
  ## this function adds spaces in table that will be used to generate the bubble plot
  ## data_table is the data that will be used for the plot
  ## rs is the rows AFTER which a space is desired
  suppressWarnings(if(is.na(rs)){
    return(data_table)
  }else{
    new_names <- vector('character', length(rs))
    for(i in 1:length(new_names)){
      new_names[i] <- paste(replicate(i, " "), collapse = "")
    }
    empty_df <- data.frame(matrix(rep(NA, times=ncol(data_table)*length(rs)), nrow=length(rs)))
    rownames(empty_df) <- new_names
    colnames(empty_df) <- colnames(data_table)
    for(j in 1:length(rs)){
      k <- rs[j] + j - 1
      data_table <- rbind(data_table[1:k,],empty_df[j,],data_table[-(1:k),])
    }
    return(data_table)
  })
}

### applying function to data ####

## uncomment below and load data from appropriate files
# tm <- read.table(<MFI table>)
# tp <- read.table(<percent positive table>)

tp <- add.spaces(tp, rs=c(2,4))
tm <- add.spaces(tm, rs=c(2,4))

## run the bubble plot function
bubble_plot(tp, tm, scale_MFI=TRUE, plot_height=5, plot_width=5.5, title='my bubble plot', color_scale=c(-2,2))

