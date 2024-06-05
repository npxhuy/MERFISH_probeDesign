# Load necessary package
library(ggplot2)
number_of_probe <- read.table("/Users/hy/Documents/GitHub/thesis/graph/raw_data/histogram_no_of_probe.txt", 
                              header = FALSE, sep = "\t")
colnames(number_of_probe) <- c("Species", "Numbers of probes")

# Create a bar plot of the number of probes
ggplot(number_of_probe, aes(x = Species, y = `Numbers of probes`)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "Bar Chart", x = "Species", y = "Number of Probes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
