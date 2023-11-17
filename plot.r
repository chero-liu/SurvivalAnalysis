# 定义一个保存图像的函数

plot_pic = function(x,pname,w=3.5,h=3.5) {
  png(paste0(pname,'.png'), width=w, height=h, units='in', res=400)
  print(x) ; dev.off()
  pdf(paste0(pname,'.pdf'), width=w, height=h,onefile = FALSE)
  print(x) ; dev.off()
}

subset_group <- function(data, geneA, geneB, a, b, c, d) {
  subset1 <- data[data[[geneA]] %in% a & data[[geneB]] %in% b, ]
  subset2 <- data[data[[geneA]] %in% c & data[[geneB]] %in% d, ]
  
  result <- rbind(subset1, subset2)
  
  return(result)
}


process_cutoff <- function(data, geneA, geneB, cutoff, outdir, prefix) {
  if (cutoff == 'best') {
    res.cut <- surv_cutpoint(data, time = 'times', event = 'patient.vital_status', variables = matching_genes)
    best_cutoff_dir <- paste0(outdir, '/bestcutoff/')
    dir.create(best_cutoff_dir, showWarnings = FALSE)
    
    cutpoint_file <- paste0(best_cutoff_dir, geneA, '_', geneB, '_optimal_cutpoint_of_variables.txt')
    write.table(res.cut$cutpoint, cutpoint_file)
    
    plotA <- plot(res.cut, geneA, palette = "npg")
    plotB <- plot(res.cut, geneB, palette = "npg")
    
    plotA_file <- paste0(best_cutoff_dir, geneA, "_", cutoff)
    plotB_file <- paste0(best_cutoff_dir, geneB, "_", cutoff)
    
    plot_pic(plotA, plotA_file, w = 5, h = 5)
    plot_pic(plotB, plotB_file, w = 5, h = 5)
    
    res.cat <- surv_categorize(res.cut)
  } else if (cutoff == 'median') {
    data[, geneA] <- ifelse(data[, geneA] <= median(data[, geneA]), 'low', 'high')
    data[, geneB] <- ifelse(data[, geneB] <= median(data[, geneB]), 'low', 'high')
    res.cat <- data
  } else {
    cutoff1 <- as.numeric(strsplit(cutoff, "_")[[1]][1])
    cutoff2 <- as.numeric(strsplit(cutoff, "_")[[1]][2])
    
    cut_value1_1 <- quantile(data[, geneA], cutoff1)
    cut_value2_1 <- quantile(data[, geneA], 1 - cutoff1)
    data[, geneA] <- ifelse(data[, geneA] <= cut_value1_1, 'low',
                                ifelse(data[, geneA] > cut_value2_1, 'high', 'NA'))
    
    data <- data[data[, geneA] %in% c('low', 'high'), ]
    
    cut_value1_2 <- quantile(data[, geneB], cutoff2)
    cut_value2_2 <- quantile(data[, geneB], 1 - cutoff2)
    data[, geneB] <- ifelse(data[, geneB] <= cut_value1_2, 'low',
                                ifelse(data[, geneB] > cut_value2_2, 'high', 'NA'))
    
    data <- data[data[, geneB] %in% c('low', 'high'), ]
    res.cat <- data
  }
  
  return(res.cat)
}

plot_ggsurvplot <- function(
    fit,                     # The survival fit object to be plotted
    xlab,                    # Label for the x-axis
    HR,                      # Hazard Ratio to annotate on the plot
    outdir = "./",           # Directory to save the output plots
    title = "",              # Title of the plot
    size = 0.6,              # Size of the plotted lines
    pval_size = 5,           # Font size for the p-value annotation
    pval_coord = c(0, 0.23), # Coordinates for the p-value annotation
    palette = c("#CD2626","#104E8B"), # Color palette for the plot
    font_main = 1,           # Main font size
    font_x = 12,             # Font size for the x-axis label
    font_y = 12,             # Font size for the y-axis label
    legend_title = "",       # Title for the legend
    legend_coords = c(0.7, 0.9), # Coordinates for the legend
    font_legend = 11,        # Font size for the legend
    risk_table = TRUE,       # Whether to include a risk table
    risk_table_title = "No. at risk:", # Title for the risk table
    tables_height = 0.15,    # Height for the risk table
    risk_table_fontsize = 3  # Font size for the risk table
) {
  
  pic <- ggsurvplot(
    fit,
    title = title,
    ggtheme = theme_survminer(),
    size = size,
    pval = TRUE, 
    pval.size = pval_size, 
    pval.coord = pval_coord,
    palette = palette,
    surv.median.line = "hv",
    conf.int = FALSE,    
    font.main = font_main,
    xlab = xlab,
    ylab = 'Overall Survival',
    font.x = font_x,
    font.y = font_y,
    legend.title = legend_title,
    legend = legend_coords,
    font.legend = font_legend,
    risk.table = risk_table,
    tables.theme = theme_cleantable(),
    risk.table.title = risk_table_title,
    risk.table.y.text = FALSE, 
    risk.table.y.text.col = TRUE, 
    tables.height = tables_height,
    risk.table.fontsize = risk_table_fontsize
  )
  
  plot_pic(pic, outdir, w = 5, h = 5)
  
  pic$plot <- pic$plot + ggplot2::annotate("text", x = max(fit$time) / 4, 
                                           y = 0.15, size = 5, label = paste("HR =", HR))
  
  plot_pic(pic, paste0(outdir,'_HR'), w = 5, h = 5)
}