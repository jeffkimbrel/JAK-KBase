library('tidyverse')
library("jakR")
library("vegan")

df = read_delim(file = '/Users/kimbrel1/Science/scripts/JAK-KBase/compareAnnotations/temp.txt', delim = "\t")

df = as.data.frame(df)
rownames(df) = df$RXN
df$RXN = NULL
df$Prokka = NULL

df.dist = vegdist(t(df), method = "jaccard")
NMDS = metaMDS(df.dist)
stressplot(NMDS)

NMDSpoints = as.data.frame(NMDS$points)
NMDSpoints$ONTOLOGY = rownames(NMDSpoints)
#NMDSpoints = left_join(NMDSpoints, groups, by = "genome")
#colnames(NMDSpoints) = c("MDS1", "MDS2", "Genome", "Cluster")
#NMDSpoints$Cluster = as.factor(NMDSpoints$Cluster)

ggplot(NMDSpoints, aes(x = MDS1, y = MDS2, fill = ONTOLOGY)) +
  jak_theme() +
  geom_point(size = 5, pch = 21, alpha = 0.9) +
  scale_fill_manual(values = tol$bright(6)) +
  ggrepel::geom_text_repel(aes(label = ONTOLOGY)) +
  theme(legend.position = "none")
