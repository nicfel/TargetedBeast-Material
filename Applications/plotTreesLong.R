library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)
library("colorblindr")
library("coda")

# Clear workspace
rm(list=ls())

burnin = 0.2
totOffset = 0.35
# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

targeted="#d95f02"
default="#1b9e77"

colors = c("Targeted"="#e41a1c", Default="#377eb8", Intervals="#984ea3")
additional_colors = c("Parsimony Only"="#ff7f00", "other1"="#4daf4a", "other2"="#f781bf", "other3"="#a65628","other4"="#ffff33", "colorgrey"="#999999")

# get all log files in out that contain default
logFiles <- list.files(path = "out_long", pattern = "\\.trees", full.names = TRUE)

data_convergence = data.frame()
posterior = data.frame()
posterior_combined = data.frame()
time_to_convergence = data.frame()

# remove any files that contain Default
logFiles <- logFiles[!grepl("Default", logFiles)]
logFiles <- logFiles[!grepl("Intervals", logFiles)]
# remove sars from the log files
logFiles <- logFiles[!grepl("sars", logFiles, ignore.case = TRUE)]

# loop over the files
# for (logFile in logFiles) {
#   print(logFile)
#   # make a mcc tree
#   system(paste0("\\/Applications\\/BEAST\\ 2.7.7\\/bin\\/logcombiner -resample 1000000 -burnin 10 -log ",
#          logFile, " -o ", gsub(".trees", ".combined", logFile)))
#   system(paste0("\\/Applications\\/BEAST\\ 2.7.7\\/bin\\/treeannotator -burnin 0 -height keep ", gsub(".trees", ".combined", logFile),
#                 " ", gsub(".trees", ".tree", logFile)))
# }

# plot all tree using ggtree
require(ggtree)
require(treeio)
mcc_trees = list.files(path = "out_long", pattern = "\\.tree$", full.names = TRUE)

x_axis_labels = c("A) H3N2 long term", "B) H3N2 short term", "C) Influenza B long term","D) MPXV  USA", "E) WNV USA")

c = 1
plots = list()
for (mcc_tree in mcc_trees) {
  print(mcc_tree)
  tree <- read.beast(mcc_tree)
  # get the name of the tree
  tree_name <- gsub(".tree", "", basename(mcc_tree))
  # get teh group number f
  # split all tree#phylo$tiplabel on | and take the second group
  dates = lapply(strsplit(tree@phylo$tip.label, "\\|"), function(x) x[length(x)])
  mrsi = max(as.Date(unlist(dates)))
  
  # plot the tree
  p <- ggtree(tree, mrsd = mrsi, size=0.2) +
    geom_tippoint(size=1.1) +
    geom_tippoint(size=0.5, color="grey") +
    ggtitle(x_axis_labels[c]) +
    theme_tree2(panel.grid.x.major = element_line(color = "grey", size = 0.5, linetype = "dashed"),
                panel.grid.x.minor = element_blank(),
                panel.grid.x.major = element_line(color = "grey", size = 0.5, linetype = "dashed"),
                panel.grid.y.minor = element_blank())
  plot(p)
    
  c=c+1
  plots[[tree_name]] <- p
}

# arrange the plots in a grid
grid.arrange(grobs = plots, ncol = 2)
# save the plots
ggsave("../Figures/trees_long.pdf", arrangeGrob(grobs = plots, ncol = 2), width = 12, height = 12)

