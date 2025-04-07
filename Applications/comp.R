library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)
library("colorblindr")

# Clear workspace
rm(list=ls())

burnin = 0.1

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

targeted="#d95f02"
default="#1b9e77"


clade_file = read.table("./out3/clade_comp.tsv", header = F)

# find the biggest differences between V2 and V3
clade_file$diff = abs(clade_file$V2 - clade_file$V3)
# get the top 5 clades
top5 = clade_file[order(clade_file$diff, decreasing = T),][1:10,]
