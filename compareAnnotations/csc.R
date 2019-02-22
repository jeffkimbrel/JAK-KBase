#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("tidyverse")

ggplot_theme <- function(base_size = 15, base_family = "", keySize = 0.3){

  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, color = "black"),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 15, margin = ggplot2::margin(1,0,1,0, "mm")),
      strip.text.y = element_text(angle = 0, size = 15)
    )
}

filename = args[1]

df = read_delim(filename, delim = "\t")

df.g = df %>%
  transform(DESCRIPTION = factor(DESCRIPTION, levels = DESCRIPTION)) %>%
  gather("type", "metric", 2:5) %>%
  transform(type = factor(type, levels = rev(c("BUFFER", "OVERLAP", "ADDED", "TOTAL"))))

df.g.bar = df.g %>% filter(type != "TOTAL")
df.g.line = df.g %>% filter(type == "TOTAL")

#png(file = args[2], height = 700, width = 600)
pdf(file = args[2], height = 10, width = 7, pointsize = 9)
ggplot(df.g.bar, aes(x = DESCRIPTION, y = metric, fill = type, color = type)) +
  ggplot_theme() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values  = c("BUFFER" = NA, "OVERLAP" = "#888888", "ADDED" = "#6EBD75", "TOTAL" = "black")) +
  scale_color_manual(values = c("BUFFER" = NA, "OVERLAP" = "black",   "ADDED" = "black", "TOTAL" = "black")) +
  geom_text(position = "stack", size = 8, vjust = -0.25, aes(label = ifelse(type == "ADDED",   metric, NA))) +
  geom_text(position = "stack", size = 8, vjust = 1.25,  aes(label = ifelse(type == "OVERLAP", metric, NA))) +
  theme(legend.position = "none") +
  #geom_line(data = df.g.line, aes(group = type)) +
  xlab("Annotation Event") +
  ylab("Count") +
  ggtitle(args[3])
dev.off()
