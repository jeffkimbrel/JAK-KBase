library("tidyverse")
library("jakR")

ggplot_theme <- function(base_size = 10, base_family = "", keySize = 0.3){

  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, color = "black"),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.key = element_rect(color = NA, fill = NA),
      legend.key.size = unit(keySize, "cm"),
      legend.text = element_text(size = 10),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, margin = ggplot2::margin(1,0,1,0, "mm")),
      strip.text.y = element_text(angle = 0, size = 10)
    )
}

filename = '/Users/kimbrel1/Dropbox/scripts/JAK-KBase/compareAnnotations/csc.txt'

df = read_delim(filename, delim = "\t")

df.g = df %>%
  transform(DESCRIPTION = factor(DESCRIPTION, levels = DESCRIPTION)) %>%
  gather("type", "metric", 2:5) %>%
  transform(type = factor(type, levels = rev(c("BUFFER", "OVERLAP", "ADDED", "TOTAL"))))

df.g.bar = df.g %>% filter(type != "TOTAL")
df.g.line = df.g %>% filter(type == "TOTAL")

#pdf(file = "~/Desktop/PDIF272563.modelseed.pdf", pointsize = 9, height = 7, width = 6)
ggplot(df.g.bar, aes(x = DESCRIPTION, y = metric, fill = type, color = type)) +
  ggplot_theme() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values  = c("BUFFER" = NA, "OVERLAP" = "#EA6F21", "ADDED" = "#6EBD75", "TOTAL" = "black")) +
  scale_color_manual(values = c("BUFFER" = NA, "OVERLAP" = "black",   "ADDED" = "black", "TOTAL" = "black")) +
  geom_text(position = "stack", vjust = -0.5, aes(label = ifelse(type == "ADDED",   metric, NA))) +
  geom_text(position = "stack", vjust = 1.5,  aes(label = ifelse(type == "OVERLAP", metric, NA))) +
  theme(legend.position = "none") +
  #geom_line(data = df.g.line, aes(group = type)) +
  xlab("Annotation Event") +
  ylab("Count") +
  ggtitle("PDIF272563 ModelSeed RXNs")
#dev.off()
