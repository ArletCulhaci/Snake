#loading the needed libraries.
library(ggplot2)

#getting the input file with absolute path and the output file with absolute path from the cli arguments.
args <- commandArgs()
input <- args[6]
output <- args[7]

#reading the data, save as data_raw.
data_raw <- read.delim(input, sep=";", header = TRUE)

#plot the remaining reads of all entries(indivduals), the lowest will be lighter of color.
a <- ggplot(data_raw, aes(x=data_raw$Filename,y=data_raw$Total, fill=data_raw$Total)) + geom_bar(stat="identity") +
  labs(x="Individuals", y="# initial reads") +
  ggtitle("# reads per sample") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10), , axis.text=element_text(size=12))
a+scale_fill_gradient(low="lightblue", high="darkblue")

#save the image as pdf
pdf(output, width=15, height=15)
plot(a)
dev.off()
