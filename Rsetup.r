#setup script for phsophoprot
packs = c("ggplot2", "xlsx", "gridExtra", "ggfortify", "VennDiagram", "amap", "pheatmap")
bcpacks = c()

for (pack in packs) {
  if(!(requireNamespace(pack, quietly=TRUE))) {
    install.packages(pack)
  }
}

for (bcpack in bcpacks) {
  if(!(requireNamespace(bcpack, quietly=TRUE))) {
    library(BiocManager)
    BiocManager::install(bcpack)
  }
}