library(dplyr)
library(circlize)
library(ComplexHeatmap)
# Load necessary package
install.packages("gplots")
library(gplots)

# Load the data
data <- read.table("/Users/hy/Documents/GitHub/thesis/graph/raw_data/heatmap_wordsize4_myReadout.txt", header = TRUE, sep = "\t", row.names = 1)

matrix <- as.matrix(data)

col_fun = colorRamp2(c(0, 20, 40), c("white", "#78BFFF", "#0085FF"))
heat <- Heatmap(matrix,
                heatmap_legend_param = list(
                  title = "Bits score"),
                col = col_fun,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))
draw(heat)

