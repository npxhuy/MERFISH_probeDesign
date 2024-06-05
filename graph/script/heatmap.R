library(dplyr)
library(circlize)
library(ComplexHeatmap)
# Load necessary package
install.packages("gplots")
library(gplots)

# Load the data
data <- read.table("/Users/hy/Documents/GitHub/thesis/graph/raw_data/heatmap_wordsize4_myReadout.txt", header = TRUE, sep = "\t", row.names = 1)

data2 <- read.table("/Users/hy/Documents/GitHub/thesis/graph/raw_data/heatmap_wordsize4_exampleReadout.txt", header = TRUE, sep = "\t", row.names = 1)

matrix <- as.matrix(data)
matrix2 <- as.matrix(data2)
col_fun = colorRamp2(c(0, 20, 40), c("white", "#4B47FF",'yellow'))
heat <- Heatmap(matrix,
                heatmap_legend_param = list(
                  title = "Bits score"),
                col = col_fun,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))
draw(heat)

heat2 <- Heatmap(matrix2,
                heatmap_legend_param = list(
                  title = "Bits score"),
                col = col_fun,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))
draw(heat2)
