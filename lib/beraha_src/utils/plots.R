library(reticulate)
source_python("utils/python_utils.py")

barplot.clusters <- function(rhoChain, correctCluster, savefile=NULL) {
  clusters = countClusters(rhoChain)
  cols = c(rep("grey", ncol(clusters)))
  heights = numeric(ncol(clusters))
  for (i in 1:ncol(clusters)) {
    heights[i] = clusters[[i]]
    if (areEqualClusters(names(clusters)[i], correctCluster)) {
      cols[i] = "blue"
    }
  }

  if (! is.null(savefile)) {
    png(savefile)
    par(mar = c(7, 4, 2, 2) + 0.2)
    barplot(heights / sum(heights), names.arg=names(clusters), las=3, col=cols)
    dev.off()
  } else {
    par(mar = c(7, 4, 2, 2) + 0.2)
    barplot(heights / sum(heights), names.arg=names(clusters), las=3, col=cols)
  }
}
