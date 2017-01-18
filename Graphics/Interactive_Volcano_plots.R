# dependencies:
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))

# path to the DiffBind table with genes
input_file <- "/Users/steve/Bioinformatics/DiffBind_scripts_reports/DiffBind_Volcano_Plot_report/input/diff_bind.Treatment4-ChIPSeq-vs-Control-ChIPSeq.p100.csv"

# read the file into a dataframe
diff_df <- read.delim(file = input_file,header = TRUE,sep = ',')

# check some attributes of the data
colnames(diff_df)

##  [1] "seqnames"                   "start"                     
##  [3] "end"                        "width"                     
##  [5] "strand"                     "Conc"                      
##  [7] "Conc_Treatment4.ChIPSeq"    "Conc_Control.ChIPSeq"      
##  [9] "Fold"                       "p.value"                   
## [11] "FDR"                        "Sample1.Treatment4.ChIPSeq"
## [13] "Sample2.Treatment4.ChIPSeq" "Sample3.Treatment4.ChIPSeq"
## [15] "Sample7.Control.ChIPSeq"    "Sample8.Control.ChIPSeq"   
## [17] "Sample9.Control.ChIPSeq"    "feature"                   
## [19] "external_gene_name"         "gene_biotype"              
## [21] "start_position"             "end_position"              
## [23] "insideFeature"              "distancetoFeature"         
## [25] "shortestDistance"           "fromOverlappingOrNearest"

dim(diff_df)

## [1] 3192   26

# keep only the fields needed for the plot
# FDR = false discovery rate = adjusted p value = significance 
diff_df <- diff_df[c("external_gene_name", "Fold", "FDR")]

# preview the dataset; data required for the plot
head(diff_df)

##   external_gene_name  Fold      FDR
## 1      RP11-431K24.1 -4.13 1.04e-05
## 2              UBE4B  2.42 8.71e-06
## 3             UBIAD1  4.27 5.50e-06
## 4             UBIAD1  1.89 2.16e-04
## 5             UBIAD1  2.74 8.59e-05
## 6             UBIAD1  3.42 1.39e-08

# add a grouping column; default value is "not significant"
diff_df["group"] <- "NotSignificant"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['Fold']) < 1.5 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df[which(diff_df['FDR'] > 0.05 & abs(diff_df['Fold']) > 1.5 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['Fold']) > 1.5 ),"group"] <- "Significant&FoldChange"


# Find and label the top peaks..
top_peaks <- diff_df[with(diff_df, order(Fold, FDR)),][1:5,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-Fold, FDR)),][1:5,])


# Add gene labels to the plot
# Single Gene Annotation example
# m <- diff_df[with(diff_df, order(Fold, FDR)),][1,]
# a <- list(
#   x = m[["Fold"]],
#   y = -log10(m[["FDR"]]),
#   text = m[["external_gene_name"]],
#   xref = "x",
#   yref = "y",
#   showarrow = TRUE,
#   arrowhead = 7,
#   ax = 20,
#   ay = -40
# )

# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks))) {
  m <- top_peaks[i, ]
  a[[i]] <- list(
    x = m[["Fold"]],
    y = -log10(m[["FDR"]]),
    text = m[["external_gene_name"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}


# make the Plot.ly plot
p <- plot_ly(data = diff_df, x = Fold, y = -log10(FDR), text = external_gene_name, mode = "markers", color = group) %>% 
  layout(title ="Volcano Plot") %>%
  layout(annotations = a)
p

# to save plot to a HTML file:
htmlwidgets::saveWidget(as.widget(p), "graph.html")

# System Information
sessionInfo()


## R version 3.3.1 (2016-06-21)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.6 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] plotly_3.6.0  ggplot2_2.1.0
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7        knitr_1.14         magrittr_1.5      
##  [4] munsell_0.4.3      colorspace_1.2-6   R6_2.1.3          
##  [7] stringr_1.1.0      httr_1.2.1         plyr_1.8.4        
## [10] tools_3.3.1        grid_3.3.1         gtable_0.2.0      
## [13] htmltools_0.3.5    assertthat_0.1     yaml_2.1.13       
## [16] digest_0.6.10      tibble_1.2         gridExtra_2.2.1   
## [19] RColorBrewer_1.1-2 formatR_1.4        tidyr_0.6.0       
## [22] viridis_0.3.4      base64enc_0.1-3    htmlwidgets_0.7   
## [25] evaluate_0.9       rmarkdown_1.0      stringi_1.1.1     
## [28] scales_0.4.0       jsonlite_1.1

