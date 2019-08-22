##
## 10. GSEA analysis
##
##
##

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(fgsea)) BiocManager::install("fgsea")
if(!require(qusage)) BiocManager::install("qusage")

library(fgsea)
library(qusage)

# Create ranks
gseaDat <- resultAB
gsea.ranks <- gseaDat[ ,"log2FoldChange"]
names(gsea.ranks) <- row.names(gseaDat)

barplot(sort(gsea.ranks, decreasing=T))

# Load pathways
#pathway <- read.gmt("D:/LYM/Projects/Data/h.all.v7.0.symbols.gmt")
pathway <- read.gmt("D:/LYM/Projects/Data/c2.all.v7.0.symbols.gmt")

# Conduct analysis
gsea.result <- fgsea(pathway, gsea.ranks, minSize=15, maxSize=500, nperm=1000)

head(gsea.result[order(padj, -abs(NES)),], n=10) # top 10 results

# Save tableplot
png(filename="tableplot_curated.png", width=800, height=12000)
topUp <- gsea.result %>% filter(ES > 0) %>% top_n(n=10, wt=-padj)
topDown <- gsea.result %>% filter(ES < 0) %>% top_n(n=10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>%  arrange(-ES)
plotGseaTable(pathway[topPathways$pathway], 
              gsea.ranks[complete.cases(gsea.ranks)], 
              gsea.result, 
              gseaParam = 0.5)
dev.off()
