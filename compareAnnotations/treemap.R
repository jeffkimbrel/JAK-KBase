library("tidyverse")
library("treemapify")
library("jakR")

df = read_delim("/Users/kimbrel1/Science/scripts/JAK-KBase/compareAnnotations/temp.txt", delim = "\t")

df$GROUP = gsub("\\[", "", df$GROUP)
df$GROUP = gsub("\\]", "", df$GROUP)
df$GROUP = gsub(" ", "", df$GROUP)

ggplot(df, aes(area = COUNT, fill = COUNT, label = GROUP)) +
  geom_treemap(color = "black") +
  geom_treemap_text(color = "black", place = "centre", grow = TRUE) +
  scale_fill_gradientn(colors=c("white","navyblue","red")) +
  #scale_fill_manual(values = tol$rainbow(15)) +
  theme(legend.position = "bottom") +
  ggtitle("PDIF272563.7 Genes")
