library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)
library("colorblindr")
library(treeio)
library(ggtree)
library(coda)

# Clear workspace
rm(list=ls())

burnin = 0.1

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

targeted="#d95f02"
default="#1b9e77"

# combine the trees files and run
#/Applications/BEAST\ 2.7.7/bin/logcombiner -resample 500000 -burnin 10 -log h3n2_recent_Targeted_dta_5000_rep*.trees -o Targeted.trees
# system("\\/Applications\\/BEAST\\ 2.7.7\\/bin\\/logcombiner -resample 500000 -burnin 10 -log out_dta/h3n2_recent_Targeted_dta_10000_rep0.trees -o out_dta/Targeted.trees")
# /Applications/BEAST\ 2.7.7/bin/treeannotator -burnin 0 -height keep Targeted.trees Targeted.tree 
# system("\\/Applications\\/BEAST\\ 2.7.7\\/bin\\/treeannotator -burnin 0 -lowMem true -height keep out_dta/Targeted.trees out_dta/Targeted.tree")
# system("\\/Applications\\/BEAST\\ 2.7.7\\/bin\\/logcombiner -burnin 10 -log out_dta/h3n2_recent_Targeted_dta_5000_rep*.log -o out_dta/Targeted.log")

# read in tree
tree = read.beast("out_dta/Targeted.tree")

# make a new node and tip trait called max_loc
tree@data$max_loc = NA
tree@data$max_prob = NA
# for each of the nodes, loop over and
for (i in seq(1, length(tree@data$max_loc), 1)) {
  # get the node
  probs = as.numeric(tree@data$l.set.prob[i][[1]])
  locs = tree@data$l.set[i][[1]]
  # get the index of the max prob
  max_index = which.max(probs)
  tree@data$max_loc[i] = locs[max_index]
  tree@data$max_prob[i] = probs[max_index]
  if (is.na(tree@data$max_loc[i])) {
    dsa
  }
}

# for each tip, get the date
mrsd = as.Date("1000-01-01")
for (i in tree@phylo$tip.label){
  tmp = strsplit(i, split="\\|")[[1]]
  mrsd = max(mrsd, as.Date(tmp[length(tmp)]))
}

# ColorBrewer’s Dark2 palette
region_colors <- c(
  Central      = "#1b9e77",  # teal/green
  East         = "#d95f02",  # burnt orange
  West         = "#7570b3",  # purple
  Canada           = "#e7298a",  # magenta/pink
  Mexico           = "#66a61e",  # green
  Central_America  = "#e6ab02",  # golden
  Caribbean        = "#a6761d"   # brown
)

color_labels = c("Central USA", "Eastern USA", "Western USA", "Canada", "Mexico", "Central America", "Caribbean")

p <- ggtree(tree, aes(color = max_loc, alpha = max_prob),
            size = 0.2, mrsd = mrsd) +
  geom_tippoint(colour = "black", size = 0.35) +
  geom_tippoint(aes(color = max_loc), size = 0.2) +
  theme_tree2() +
  scale_color_manual(name = "location with\nhighest posterior",
                     values = region_colors, labels = color_labels) +
  scale_alpha_continuous(name = "posterior\nsupport for location") +
  # put the legend into the figure on the top left
  theme(
    legend.position      = c(0, 1),   # (x, y) in npc units: 0 = left, 1 = top
    legend.justification = c(0, 1),   # anchor the legend’s own top-left corner
    legend.box.just      = "left",    # keep multiple guides left-aligned
    ## remove the white panels
    legend.background    = element_blank(),   # outer box
    legend.key           = element_blank(),   # keys behind symbols
    legend.box.background = element_blank()   # box around grouped legends
  ) +
  
  ## 2. Use two columns for every guide that appears
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 3)),
    alpha = guide_legend(ncol = 2)
  )+
  ylim(1, 10000)

# add rectangle layers (they will be last for now)
p <- p +
  annotate("rect", xmin = 2020.5, xmax = 2021, ymin = -Inf, ymax = Inf,
           fill = "grey95", alpha = 1) +
  annotate("rect", xmin = 2021.5, xmax = 2022, ymin = -Inf, ymax = Inf,
           fill = "grey95", alpha = 1) +
  annotate("rect", xmin = 2022.5, xmax = 2023, ymin = -Inf, ymax = Inf,
           fill = "grey95", alpha = 1)+
  annotate("rect", xmin = 2023.5, xmax = 2024, ymin = -Inf, ymax = Inf,
           fill = "grey95", alpha = 1)+
  annotate("rect", xmin = 2024.5, xmax = 2025, ymin = -Inf, ymax = Inf,
           fill = "grey95", alpha = 1)


# move the three newest layers (rectangles) to the bottom
p$layers <- append(p$layers[ (length(p$layers)-4) : length(p$layers) ],
                   p$layers[ 1 : (length(p$layers)-5) ])

plot(p)

ggsave("dta.pdf", p, width = 9, height = 6)


# read in the log file
log_large = read.table("out_dta/h3n2_recent_Targeted_dta_10000_rep0.log", header=TRUE, sep="\t")
# plot the posterior, likelihood and prior vs. samples
# calculate the posterior ESS
ess = effectiveSize(log_large$posterior)
# calculate the likelihood ESS
ess_likelihood = effectiveSize(log_large$likelihood)
# calculate the prior ESS
ess_prior = effectiveSize(log_large$prior)
values = data.frame(
  sample = log_large$Sample,
  value = log_large$posterior,
  quantity = paste0("Posterior (ESS=", round(ess, 0), ")")
)
values = rbind(values, data.frame(
  sample = log_large$Sample,
  value = log_large$likelihood,
  quantity = paste0("Likelihood (ESS=", round(ess_likelihood, 0), ")")
))
values = rbind(values, data.frame(
  sample = log_large$Sample,
  value = log_large$prior,
  quantity = paste0("Prior (ESS=", round(ess_prior, 0), ")")
))

# reorder to posterior likelohood prior
values$quantity = factor(values$quantity, levels=c(
  paste0("Posterior (ESS=", round(ess, 0), ")"),
  paste0("Likelihood (ESS=", round(ess_likelihood, 0), ")"),
  paste0("Prior (ESS=", round(ess_prior, 0), ")")
))

p2 = ggplot(values, aes(x=sample/10^6, y=value)) +
  geom_line(size=0.1) +
  theme_minimal()+
  facet_wrap(~quantity, ncol=1, scales="free_y")+
  ylab("") + xlab("Million MCMC samples")
plot(p2)
ggsave("dta_traces.pdf", p2, width = 7, height = 4)


# make a combined plot next to each other where the tree takes 2/3 of the space
# and p2 takes 1/3
p3 = grid.arrange(p, p2, ncol=2, widths=c(3/4, 1/4))
  
ggsave("dta.pdf", p3, width = 12, height = 7)



log_alt = read.table("out_dta/h3n2_recent_Targeted_dta_1000_rep0.log", header=TRUE, sep="\t")
# remove the first 10 % as burnin
log_alt = log_alt[round(nrow(log_alt)*burnin):nrow(log_alt),]
red_data = data.frame()
# loop over all labels in log
for (l in labels(log_large)[[2]]){
  # if starts with geoSubstModelLogger.trait.relGeoRate traitClockRate.trait
  if (startsWith(l, "geoSubstModelLogger.trait.relGeoRate") || startsWith(l, "traitClockRate.trait")){
    # get the values for this label
    value = HPDinterval(as.mcmc(log(log_large[, l])))
    values_alt = HPDinterval(as.mcmc(log(log_alt[,l])))
    
    width = value[2] - value[1]
    width_alt = values_alt[2] - values_alt[1]
    # compute the reduction in the width of the interval
    reduction = (width/width_alt)
    # get the mean of the value in log
    mean_val = mean(log_large[, l])
    
    # get the name of the label
    name = strsplit(l, "\\.")[[1]][3]
    # add to data frame
    red_data = rbind(red_data, data.frame(label=name, value=reduction, mean=mean_val))
  }
}

red_data$point = "relative rate"
red_data[is.na(red_data$label),"point"] = "trait clock rate"
  
p = ggplot(red_data, aes(x=mean, y=value, color=point, size=point)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, color="black") +
  scale_y_log10() +
  ylab("Reduction in width of 95% HPD interval") +
  xlab("Parameter value") +
  scale_color_manual(values=c("black", "#DC143C"), name="Parameter type") +
  scale_size_manual(values=c("trait clock rate"=3, "relative rate"=1), guide=FALSE) +
  theme_minimal()
plot(p)
ggsave("dta_reduction.pdf", p, width = 7, height = 5)
