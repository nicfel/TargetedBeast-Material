library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)
library("colorblindr")
library(coda)


# Clear workspace
rm(list=ls())


# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# 
# # open h3n2_recent_Targeted_5000_rep2.log in out4
# out4 <- read.table("./out4/h3n2_recent_Targeted_1000_rep2.log", header=TRUE, sep="\t")
# # remove 5 % burnin
# burnin = 0.05
# out4 <- out4[round(nrow(out4)*burnin):nrow(out4),]
# 
# # cluster the treeLikelihood values in two groups, plot them against the sample number
# out4$treeLikelihood <- as.numeric(out4$treeLikelihood)
# out4$sample <- 1:nrow(out4)

# get every line starting with # in ~/Documents/beast27/TargetedBeast/examples/mpxv_all_Targeted_733_rep0.clades
a = readLines("~/Documents/beast27/TargetedBeast/examples/mpxv_all_Targeted_733_rep0.clades")
# remove any line that does not start with #
a = a[str_detect(a, "^#")]
# remove the # from the beginning of each line
a = str_replace(a, "^#", "")
# get the NR to ID mapping by splitting every line on tab
b = strsplit(a, split="\t")

clade_isolate = c()
for (i in 1:length(b)){
  clade_isolate = c(clade_isolate, b[[i]][2])
}


t = read.table("~/Documents/beast27/TargetedBeast/examples/mpxv_all_Targeted_733_rep0.clades", sep="\t", header=F, skip=1)
log = read.table("~/Documents/beast27/TargetedBeast/examples/mpxv_all_Targeted_733_rep0.log", sep="\t", header=T )
burnin=0.1
l_before = nrow(t)
#skip the first 10 %
t = t[round(nrow(t)*burnin):4:nrow(t),2:ncol(t)]
log = log[round(nrow(log)*burnin):4:nrow(log),]

# get how many samples were skipped
n_skipped = l_before - nrow(t)

likelihood = log[,3]
prior = log[,4]
# flatten all element in t into a vector
talt = c()
for (i in 1:nrow(t)){
  talt = c(talt, t[i,])
}
# remove any element of talt if it has more than 5 ,
# talt = talt[str_count(talt, ",")<=5]

# get the frequency of all unique clades in talt
unique_clades = unique(talt)
uc = unique_clades
# get the freuqency of each clade
c=1
for (i in unique_clades){
   ind = which(talt==i)
   freq = length(ind)/nrow(t)
   talt = talt[-ind]
   if (freq>0.8 || freq <0.2){
     print(paste(length(ind), length(uc), length(talt)))
     uc = uc[-c]
   }else{
     c=c+1
   }
}

# for each clade, get the ess based on presence absence in t for each iteration
candidate_clades=c()
for (j in 1:length(uc)){
  clade = uc[j]
  ess = c()
  for (i in 1:nrow(t)){
    ess = c(ess, sum(t[i,2:ncol(t)] == clade))
  }
  if (mean(ess)==1){next}
  
  # calculate the ess_val
  ess_val = effectiveSize(ess)
  if (ess_val==0){ds}
  if (ess_val>10){next}
    
    
  candidate_clades = c(candidate_clades, clade)
  print(paste(j, length(uc), clade, ess_val))
  # plot the likelihood and log and color by on or off clade
  plot_dat = data.frame(x=1:nrow(t)*4-3, likelihood=likelihood, prior=prior, clade=ess)
  p1=ggplot(plot_dat, aes(x=x+n_skipped, y=likelihood, color=clade)) + 
    geom_line() + ggtitle(paste("Clade: ", clade, "ESS: ", ess_val))
  p2=ggplot(plot_dat, aes(x=x+n_skipped, y=prior, color=clade)) +
    geom_line() + ggtitle(paste("Clade: ", clade, "ESS: ", ess_val))
  
  # make a list with all the elements in b based on the clade numbers
  clade_nr = strsplit(clade[[1]], ",")[[1]]
  clade_nr = as.numeric(clade_nr)+1
  sample_names = clade_isolate[clade_nr]
  # plot the same names as a horizontal table
  p3 = ggplot() + geom_text(data=data.frame(x=1, y=seq(1,length(sample_names)), label=sample_names), aes(x=x, y=y, label=label)) + theme_void()
  
  p=grid.arrange(p1, p2, p3, ncol=1)
  plot(p)
  ggsave(plot=p, paste("./Figures/clade_", j, "_ess_", ess_val, ".pdf", sep=""), width=10, height=10)

}

