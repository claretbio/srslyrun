library(stats)
library(ggplot2)


args <- commandArgs(trailingOnly=TRUE)

infile <- args[1]
libname <- args[2]
outfile <- args[3]


histo <- read.table(infile, sep="\t")
histo$V1 <- abs(histo$V1)
cutoff <- quantile(histo$V1, seq(0,1,0.005))[200]
lengths <- as.data.frame(histo[which(histo$V1 <= cutoff),])
histo_table <- as.data.frame(table(lengths))
libsum <- sum(histo_table$Freq)
histo_table$perc <- histo_table$Freq / libsum
histo_table$perc <- as.numeric(histo_table$perc)
histo_table$lengths <- as.numeric(histo_table$lengths)
plt <- ggplot(histo_table) + aes(x=lengths, y=perc) + geom_line() + geom_point() + scale_x_continuous("Length (bp)") + labs(title="Mapped template molecule lengths", subtitle = libname, caption='Plot brought to you by your friends at Claret Bio ;)') + ylab('Percent library at length') 

ggsave(plt, file=outfile, device="pdf", width=7, height=7)

