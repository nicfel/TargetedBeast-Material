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


ESS_Target = 100

effectiveSize <- function(x, max_lag = NULL) {
  N <- length(x)
  x_demeaned <- x - mean(x)
  if (is.null(max_lag)) {
    max_lag <- min(N - 1, 10000)  # Usually you wouldn't need lags beyond a few thousand
  }
  acf_vals <- acf(x_demeaned, lag.max = max_lag, plot = FALSE,
                  type = "correlation", demean = FALSE)$acf
  acf_lags <- acf_vals[-1]  # remove the lag 0 element

  # Sum autocorrelations while they remain positive
  pos_sum <- 0
  for (r in acf_lags) {
    if (r <= 0) break
    pos_sum <- pos_sum + r
  }
  ess <- N / (1 + 2 * pos_sum)
  return(ess)
}


# get all log files in out that contain default
logFiles <- list.files(path = "out_long", pattern = "\\.log", full.names = TRUE)

data_convergence = data.frame()
posterior = data.frame()
posterior_combined = data.frame()
time_to_convergence = data.frame()
# leafs = c("infB"=1978, "h3n2"=5000, "wnv"=2664, "mpxv"=1579)

# remove any files that contain Default
logFiles <- logFiles[!grepl("Default", logFiles)]
# remove sars from the log files
logFiles <- logFiles[!grepl("sars", logFiles, ignore.case = TRUE)]


# get all potential datasets based on the names of the log files
datasets = c()
for (logFile in logFiles) {
  # remove out/, .log and replace _Default_ with _ a
  tmp <- str_replace(logFile, "out_long/", "")
  tmp <- str_replace(tmp, "\\.log", "")
  # split on _
  name = str_split(tmp, "_")[[1]]
  datasets = c(datasets, paste(name[[1]], name[[2]], "\n", name[[4]], "samples"))
}  

datasets = unique(datasets)

# loop over the files
for (logFile in logFiles) {
  print(logFile)
  # remove out/, .log and replace _Default_ with _ a
  tmp <- str_replace(logFile, "out_long/", "")
  tmp <- str_replace(tmp, "\\.log", "")
  # split on _
  name = str_split(tmp, "_")[[1]]
  # read the log file
  log <- read.table(logFile, header = TRUE, sep="\t",fill = TRUE)
  # remove any line if there is a NA anywhere
  if (any(is.na(log))){
    ind = which(apply(log, 1, function(x) any(is.na(x))))
    log = log[-ind,]
    # overwrite the log file
    write.table(log, logFile, sep="\t", row.names=FALSE)
  }
  
  # read in the corresponding out file line by line and look for any lines with seconds in them
  outLine = readLines(gsub(".log", ".out", logFile))
  # only keep lines that contain seconds
  outLine = outLine[grepl("seconds", outLine)]
  # split on \\s and convert the fourth group to numeric
  tmp = str_split(outLine, "\\s+")
  total_seconds = 0
  for (i in seq(1, length(tmp), 1)){
    if (length(tmp[[i]])>4){
      total_seconds = total_seconds + as.numeric(tmp[[i]][4])
    }
  }

  # read in the default file for comparison, if is targeted, use intervals, if is intervals, use default
  if (grepl("Targeted", name[[3]])){
    companme = "Targeted"
  } else {
    companme = "Default"
  }
  burnincomp = "Targeted"
  
  log.default <- read.table(str_replace(logFile, name[[3]], companme), header = TRUE, sep="\t", fill = TRUE)
  log.burnincomp <- read.table(str_replace(logFile, name[[3]], burnincomp), header = TRUE, sep="\t", fill = TRUE)
  
  ## get the corresponding out file and total computation time
  outLine.default = readLines(gsub(".log", ".out", str_replace(logFile, name[[3]], companme)))
  # only keep lines that contain seconds
  outLine.default = outLine.default[grepl("seconds", outLine.default)]
  # split on \\s and convert the fourth group to numeric
  tmp.default = str_split(outLine.default, "\\s+")
  total_seconds.default = 0
  for (i in seq(1, length(tmp.default), 1)){
    if (length(tmp.default[[i]])>4){
      total_seconds.default = total_seconds.default + as.numeric(tmp.default[[i]][4])
    }
  }
  
  if (any(is.na(log.default))){
    ind = which(apply(log.default, 1, function(x) any(is.na(x))))
    log.default = log.default[-ind,]
  }
  
  # calculate the ESS value for the posterior, likelihood and prior
  for (header in c("posterior", "likelihood", "prior")) {
    subset_log = log[round(nrow(log)*burnin):nrow(log),]
    subset_log_default = log.default[round(nrow(log.default)*burnin):nrow(log.default),]

    ess <- effectiveSize(subset_log[,header])
    ess.default <- effectiveSize(subset_log_default[,header])

    low_ESS = ess<ESS_Target
    low_ESS_default = ess.default<ESS_Target

    samples = max(subset_log[, "Sample"])-min(subset_log[, "Sample"])
    samples.default = max(subset_log_default[, "Sample"])-min(subset_log_default[, "Sample"])
    
    ratio = (ess/samples) / (ess.default/samples.default)

    # make a conservative estimate of the 95% HPD in log[, header] using only
    # the last 50% of samples
    limits_upper <- quantile(log.burnincomp[round(nrow(log.burnincomp)/2):nrow(log.burnincomp),header], 0.975)
    limits_lower <- quantile(log.burnincomp[round(nrow(log.burnincomp)/2):nrow(log.burnincomp),header], 0.025)

    # get the number of the first sample larger than the mean
    ind = which(log[,header] > limits_lower & log[,header] < limits_upper)
    ind.default = which(log.default[,header] > limits_lower & log.default[,header] < limits_upper)
    # find the first ind[i] for which ind[i]+1=ind[i+1]
    first_sample = NA
    for (i in seq(1, length(ind)-10)){
      is_first = T
      for (j in seq(1, 10)){
        if (ind[i]+j != ind[i+j]){
          is_first = F
          break
        }
      }
      if (is_first){
        first_sample = ind[i]
        break
      }
    }
    for (i in seq(1, length(ind.default)-10)){
      is_first = T
      for (j in seq(1, 10)){
        if (ind.default[i]+j != ind.default[i+j]){
          is_first = F
          break
        }
      }
      if (is_first){
        first_sample_default = ind.default[i]
        break
      }
    }

    if (is.na(first_sample_default)){
      first_sample_default = nrow(log.default)
    }

    dataset_name = paste(name[[1]],name[[2]], "\n",  name[[4]], "samples")
    
    x_val = which(datasets == dataset_name)
    
    if (grepl("Targeted", name[[3]])){
      offset = totOffset
    } else {
      offset = 0
    }
    
    comparison_method = "Default"
    plot=T
    if (grepl("Targeted", name[[3]])){
      comparison_method = "Intervals"
      plot=F
    }
    
    # compute the minimum time to convergence for estimated burnin periods where
    # the ESS is above 100
    print("start time to ess calc")
    if (!low_ESS){
      vals = log[,header]
      k=NA
      for (j in seq(first_sample, length(vals)-ESS_Target)){
        # calculate the ESS
        ess_val = effectiveSize(vals[j:length(vals)])
        if (ess_val>ESS_Target){
          k = length(vals)
          while (ess_val>ESS_Target){
            # find the first sample where the ESS is less than 100
            momentum = (k - k*ESS_Target/ess_val)/4
            momentum = as.numeric(ceiling(momentum))
            k = round(k-momentum)
            ess_val = effectiveSize(vals[j:k])
          }
          # increase k again until the ESS is above 100
          for (i in seq(k, length(vals)-ESS_Target)){
            ess_val = effectiveSize(vals[j:i])
            if (ess_val>ESS_Target){
              k = i
              break
            }
          }
          samples_to_converges = log[k, "Sample"]
          break
        }
      }
      # calculate the time to convergence
      # get the ratio of sample
      print("a")
      ratio = samples_to_converges/max(log[, "Sample"])
    }
    if (is.na(k) || low_ESS){
      ratio = 1
    }
    time_to_convergence = rbind(time_to_convergence,
                                 data.frame( x = x_val+offset,
                                            dataset = dataset_name,
                                            rep = name[[5]],
                                            samples = samples_to_converges,
                                            header = header,
                                            low_ESS = low_ESS,
                                            method=name[[3]],
                                            time = total_seconds*ratio))
    print("b")
    if (plot){
      if (!low_ESS_default){
        # do the same for default
        vals = log.default[,header]
        k=NA
        print("b1")
        for (j in seq(first_sample_default, length(vals)-ESS_Target)){
          # calculate the ESS
          ess_val = effectiveSize(vals[j:length(vals)])

          if (ess_val>ESS_Target){
            k = length(vals)
            while (ess_val>ESS_Target){
              # find the first sample where the ESS is less than 100
              momentum = (k - k*ESS_Target/ess_val)/4
              momentum = ceiling(momentum)
              k = round(k-momentum)
              ess_val = effectiveSize(vals[j:k])
            }
            print("b2")
            # increase k again until the ESS is above 100
            for (i in seq(k, length(vals)-ESS_Target)){
              ess_val = effectiveSize(vals[j:i])
              if (ess_val>ESS_Target){
                k = i
                break
              }
            }
            samples_to_converges.default = log.default[k, "Sample"]
            break
          }
        }
        # calculate the time to convergence
        # get the ratio of sample
        print("a")
        ratio = samples_to_converges.default/max(log.default[, "Sample"])
      }
      if (is.na(k) || low_ESS_default){
        ratio = 1
      }
      time_to_convergence = rbind(time_to_convergence,
                                   data.frame( x = x_val-totOffset,
                                              dataset = dataset_name,
                                              rep = name[[5]],
                                              samples = samples_to_converges.default,
                                              header = header,
                                              low_ESS = low_ESS_default,
                                              method=comparison_method,
                                              time = total_seconds.default*ratio))
    }
    
    
    data_convergence <- rbind(data_convergence, 
                              data.frame( x = x_val+offset,
                                        dataset = dataset_name, 
                                         rep = name[[5]],
                                         name = "individual",
                                         samples = samples,
                                         header = header, 
                                         method=name[[3]],
                                         ess = ess,
                                         plot=T,
                                         time = total_seconds*(1-burnin),
                                         burnin = log[first_sample, "Sample"]))
    data_convergence <- rbind(data_convergence, 
                              data.frame( x = x_val-totOffset,
                                          dataset = dataset_name, 
                                         rep = name[[5]],
                                         name = "individual",
                                         samples = samples,
                                         header = header, 
                                         method=comparison_method,
                                         ess = ess.default,
                                         plot=plot,
                                         time = total_seconds.default*(1-burnin),
                                         burnin = log.default[first_sample_default, "Sample"]))
    
    
    posterior <- rbind(posterior, 
                       data.frame(dataset = dataset_name, 
                                  rep = name[[5]],
                                  samples = log[round(nrow(log)*0.01):nrow(log), "Sample"],
                                  header = header, 
                                  method=name[[3]],
                                  value=log[round(nrow(log)*0.01):nrow(log),header]))
    posterior <- rbind(posterior,
                       data.frame(dataset = dataset_name, 
                                  rep = name[[5]],
                                  samples = log.default[round(nrow(log.default)*0.01):nrow(log.default), "Sample"],
                                  header = header, 
                                  method=comparison_method,
                                  value=log.default[round(nrow(log.default)*0.01):nrow(log.default),header]))
    
    posterior_combined  = rbind(posterior_combined, 
                                data.frame(dataset = dataset_name, 
                                           samples = log[round(nrow(log)*burnin):nrow(log), "Sample"],
                                           rep = name[[5]],
                                           header = header, 
                                           method=name[[3]],
                                           time = total_seconds*(1-burnin),
                                           value=log[round(nrow(log)*burnin):nrow(log),header]))
    if (plot){
      posterior_combined  = rbind(posterior_combined,
                                  data.frame(dataset = dataset_name, 
                                             samples = log.default[round(nrow(log.default)*burnin):nrow(log.default), "Sample"],
                                             rep = name[[5]],
                                             header = header, 
                                             method=comparison_method,
                                             time = total_seconds.default*(1-burnin),
                                             value=log.default[round(nrow(log.default)*burnin):nrow(log.default),header]))
    }
  }
  
}



methods = c("Default",  "Intervals", "Targeted")
offset = c(-totOffset, 0, totOffset)
for (dataset in unique(posterior_combined$dataset)){
  for (header in c("posterior", "likelihood", "prior")){
    
    cc=1
    for (c in methods){
      total_samples = 0
      total_time = 0
      for (r in unique(posterior_combined$rep)){
        # get the values for the samples
        vals = posterior_combined[posterior_combined$dataset == dataset & 
                                  posterior_combined$method == c &
                                  posterior_combined$header == header &
                                  posterior_combined$rep == r, c("samples", "time")]
        total_samples = total_samples + max(vals[[1]]) - min(vals[[1]])
        total_time = total_time + vals[[2]][1]
      }

      vals = posterior_combined[posterior_combined$dataset == dataset & 
                                  posterior_combined$method == c &
                                  posterior_combined$header == header, "value"]
      
      x_val = which(datasets == dataset)
      ess = effectiveSize(vals)
      
      
      data_convergence <- rbind(data_convergence, 
                                data.frame( x = x_val+offset[cc],
                                            dataset = dataset, 
                                            rep = "combined",
                                            name = "combined",
                                            samples = total_samples,
                                            header = header, 
                                            method=c,
                                            ess = ess,
                                            plot=T,
                                            time = total_time,
                                            burnin = NA))
      cc = cc+1
      
    }
  }
}

# make lines between any combination of rep in Default and rep in Targeted for 
# rep0 rep1 and rep2 for the ess values
lines = data.frame()
improvement = data.frame() # keep track of the mean improvement

# define the comparisons
comp_from = c("Default", "Intervals", "Default")
comp_to = c("Intervals", "Targeted", "Targeted")
off_from = c(-totOffset, 0, -totOffset)
off_to = c(0, totOffset, totOffset)
plot = c(T, T, F)
v_off = c(0.8, 0.8, 1)
text_addition = c("","","combined:\n")

for (dataset in unique(posterior_combined$dataset)){
  for (header in c("posterior", "likelihood", "prior")){
    val = c()
    val.height = c()
    val.burnin = c()
    val.burnin.height = c()
    for (c in seq(1,length(comp_to))){
      cf = comp_from[c]
      ct = comp_to[c]
      
      from = data_convergence[data_convergence$dataset == dataset & 
                          data_convergence$header == header &
                          data_convergence$name == "combined" &
                          data_convergence$method == cf, ]
      to = data_convergence[data_convergence$dataset == dataset & 
                          data_convergence$header == header &
                          data_convergence$name == "combined" &
                          data_convergence$method == ct, ]
      ess_from = from$ess/from$time*60*60
      ess_to = to$ess/to$time*60*60

      x_val = which(datasets == dataset)
      
      lines = rbind(lines, data.frame(xfrom = x_val+off_from[c], yfrom = ess_from, 
                              xto = x_val+off_to[c], yto = ess_to,
                              dataset = dataset, header = header, 
                              text_addition = text_addition[c],
                              method = cf,
                              plot=plot[c]))
      
      # calculate the average improvement
      # if (c==3){
        improvement = rbind(improvement, data.frame(x=x_val +(off_from[c] +off_to[c])/2, y=mean(c(ess_to, ess_from))*v_off[c], 
                                                    dataset = dataset, header = header, 
                                                    text_addition = text_addition[c],
                                                    improvement = ess_to/ess_from, burnin = ess_to/ess_from,
                                                    method = "Default"))
      # }
    }
  }
}

# reorder facets for header
data_convergence$header <- factor(data_convergence$header, levels = c("posterior", "likelihood", "prior"))
data_convergence$name = paste(data_convergence$name, "runs")
data_convergence$name <- factor(data_convergence$name, levels = c("combined runs", "individual runs"))

improvement$header <- factor(improvement$header, levels = c("posterior", "likelihood", "prior"))
lines$header <- factor(lines$header, levels = c("posterior", "likelihood", "prior"))

# make an additional dataframe that points to the improvement, but only in the posterior for H3N2 HA
explanation = improvement[improvement$header == "posterior" & improvement$dataset == "h3n2 HA \n 1000 samples",]
# improvement = improvement[-which(improvement$header == "posterior" & improvement$dataset == "h3n2 HA \n 1000 samples"),]

# change the names of the datasets to 
x_axis_labels = c("H3N2\nlong term", "H3N2\nshort term", "Influenza B\nlong term","MPXV\n USA", "SARS-CoV-2\nGlobal","WNV\nUSA")
x_axis_labels = c("H3N2\nlong term", "H3N2\nshort term", "Influenza B\nlong term","MPXV\n USA", "WNV\nUSA")

# plot the ESS values
p = ggplot() + 
  annotate("rect",xmin  = 0.5, xmax = 1.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 2.5, xmax = 3.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 4.5, xmax = 5.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  geom_point(data=data_convergence[data_convergence$plot & data_convergence$name=="combined runs",], aes(x=x, y=ess/time*60*60, fill=method, size=name), color="black", shape = 21) +
  geom_point(data=data_convergence[data_convergence$plot & data_convergence$name=="individual runs",], aes(x=x, y=ess/time*60*60, fill=method, size=name), color="black", shape = 21) +
  labs(x = "", y = "ESS per hour") +
  facet_grid(header~.) +
  scale_y_log10() +
  theme_minimal() +
  scale_fill_manual(name="operator setup", values=colors) +
  scale_color_manual(name="operator setup",values=colors) +
  geom_segment(data=lines[lines$plot,], aes(x=xfrom, y=yfrom, xend=xto, yend=yto), color="black", size=0.05) +
  scale_x_continuous(breaks = seq(1, length(datasets)), labels=x_axis_labels) +
  scale_size_manual(name="", values=c(3, 1)) +
  geom_label(data = improvement[improvement$text_addition=="",],
             aes(x = x, y = y, label = paste0(round(improvement, 1))),
             size = 3,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  geom_label(data = improvement[improvement$text_addition!="",],
             aes(x = x, y = 100, label = paste0(round(improvement, 1))),
             size = 3,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  coord_cartesian(ylim = c(0.25, 500))+
  theme(panel.grid.minor.x = element_blank())+
  # make the major lines wider and grey
  theme(panel.grid.major.x = element_blank())
plot(p)
ggsave("../Figures/convergence_longerruns.pdf", p, width = 8, height = 4)



# calculate the average reduction in burnin using comp_form and comp_to
burnin_improvement = data.frame()
for (dataset in unique(posterior_combined$dataset)){
  for (header in c("posterior", "likelihood", "prior")){
    val = c()
    for (c in seq(1,length(comp_to))){
      cf = comp_from[c]
      ct = comp_to[c]
      
      from = data_convergence[data_convergence$dataset == dataset & 
                          data_convergence$header == header &
                          data_convergence$name == "individual runs" &
                          data_convergence$method == cf & data_convergence$plot, "burnin"]
      to = data_convergence[data_convergence$dataset == dataset & 
                          data_convergence$header == header &
                          data_convergence$name == "individual runs" &
                          data_convergence$method == ct & data_convergence$plot, "burnin"]
      
      avg = c()
      # calculate the ratio for any pair
      for (f in from){
        for (t in to){
          avg = c(avg, log(f/t))
        }
      }
      
      x_val = which(datasets == dataset)

      burnin_improvement = rbind(burnin_improvement, 
                                  data.frame(x=x_val +(off_from[c] +off_to[c])/2, y=mean(c(from, to))*v_off[c], 
                                             dataset = dataset, header = header, 
                                             improvement = exp(mean(avg)),
                                             plot=plot[c],
                                             method = cf))
    }
  }
}

burnin_improvement$header <- factor(burnin_improvement$header, levels = c("posterior", "likelihood", "prior"))


# make the same figure for the burnin
p = ggplot(data_convergence[data_convergence$plot,], aes(x=x, y=burnin, fill=method)) + 
  annotate("rect",xmin  = 0.5, xmax = 1.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 2.5, xmax = 3.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 4.5, xmax = 5.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  geom_point(aes(color=method), shape = 21) +
  labs(x = "", y = "Number of samples for burnin phase") +
  facet_grid(header~.) +
  scale_y_log10() +
  theme_minimal() +
  scale_fill_manual(values=colors, name="operator setup") +
  scale_color_manual(values=colors, name="operator setup") +
  scale_x_continuous(breaks = seq(1, length(datasets)), labels=x_axis_labels)+
  geom_label(data=burnin_improvement[burnin_improvement$plot,], aes(x=x, y=y, label = paste0(signif(round(improvement,1), digits=2))),size = 3,
              # make bold
              fontface = "bold",
              color = "black",
              fill = NA,    # White background
              label.size = NA,   # Remove label border if desired
              label.padding = unit(0., "lines"),  # Make box as tight as possible
              vjust=0) +
  geom_label(data=burnin_improvement[!burnin_improvement$plot,], aes(x=x, y=10^4, label = paste0(signif(improvement, digits=2))),size = 3,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  theme(panel.grid.minor.x = element_blank())+
  # make the major lines wider and grey
  theme(panel.grid.major.x = element_blank())


plot(p)
ggsave("../Figures/burnin_longerruns.pdf", p, width = 8, height = 4)





# calculate the average reduction in burnin using comp_form and comp_to
time_improvement = data.frame()
for (dataset in unique(time_to_convergence$dataset)){
  for (header in c("posterior", "likelihood", "prior")){
    val = c()
    for (c in seq(1,length(comp_to))){
      cf = comp_from[c]
      ct = comp_to[c]
      
      from = time_to_convergence[time_to_convergence$dataset == dataset & 
                                   time_to_convergence$header == header &
                                   time_to_convergence$method == cf,]
      to = time_to_convergence[time_to_convergence$dataset == dataset & 
                              time_to_convergence$header == header &
                              time_to_convergence$method == ct,]
      avg = c()
      # calculate the ratio for any pair
      for (f in from$time){
        for (t in to$time){
          avg = c(avg, log(f/t))
        }
      }
      x_val = which(datasets == dataset)
      
      if (any(from$low_ESS) || any(to$low_ESS)){
        # if either is low ESS, set the improvement to 1
        uncertain = "+"
      }
      else{
        uncertain = ""
      }
      
      time_improvement = rbind(time_improvement, 
                                 data.frame(x=x_val +(off_from[c] +off_to[c])/2, y=mean(c(from$time, to$time))*v_off[c], 
                                            dataset = dataset, header = header, 
                                            improvement = exp(mean(avg)),
                                            uncertain = uncertain,
                                            plot=plot[c],
                                            method = cf))
    }
  }
}

time_improvement$header <- factor(time_improvement$header, levels = c("posterior", "likelihood", "prior"))
time_to_convergence$header <- factor(time_to_convergence$header, levels = c("posterior", "likelihood", "prior"))


# plot the time to convergence
p = ggplot() + 
  annotate("rect",xmin  = 0.5, xmax = 1.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 2.5, xmax = 3.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 4.5, xmax = 5.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  
  geom_point(data=time_to_convergence, aes(x=x, y=time/60/60, fill=method, color=method, shape=low_ESS)) +
  labs(x = "", y = "Hours to ESS > 100") +
  facet_grid(header~.) +
  scale_shape_manual(values=c(21, 24), name="ESS below threshold") +
  scale_y_log10() +
  theme_minimal() +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  scale_x_continuous(breaks = seq(1, length(datasets)), labels=datasets) +
  geom_label(data=time_improvement[time_improvement$plot,], aes(x=x, y=y/60/60, label = paste0(signif(round(improvement,1), digits=2), uncertain)),size = 2,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  geom_label(data=time_improvement[!time_improvement$plot,], aes(x=x, y=1, label = paste0(signif(improvement, digits=2), uncertain)),size = 2,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  # remove the minor vertical lines
  theme(panel.grid.minor.x = element_blank())+
  # make the major lines wider and grey
  theme(panel.grid.major.x = element_blank())
  
  # theme(panel.grid.major.x = element_line(size=5, color="grey95"))

plot(p)

ggsave("../Figures/time_to_convergence_longerruns.pdf", p, width = 8, height = 4)




# plot the time to convergence
p = ggplot() + 
  annotate("rect",xmin  = 0.5, xmax = 1.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 2.5, xmax = 3.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  annotate("rect",xmin  = 4.5, xmax = 5.5, ymin  = 0, ymax =  Inf, fill  = "grey70", alpha = .3) +
  labs(x = "", y = "Hours to posterior ESS > 100") +
  # facet_grid(header~.) +
  scale_shape_manual(values=c(21, 24), name="ESS below threshold") +
  geom_point(data=time_to_convergence[time_to_convergence$header=="posterior",], aes(x=x, y=time/60/60, fill=method, color=method, shape=low_ESS)) +
  scale_y_log10() +
  theme_minimal() +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  scale_x_continuous(breaks = seq(1, length(datasets)), labels=x_axis_labels) +
  geom_label(data=time_improvement[time_improvement$plot & time_improvement$header=="posterior",], aes(x=x, y=y/60/60, label = paste0(signif(round(improvement,1), digits=2), uncertain)),size = 3,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  geom_label(data=time_improvement[!time_improvement$plot & time_improvement$header=="posterior",], aes(x=x, y=1, label = paste0(signif(improvement, digits=2), uncertain)),size = 3,
             # make bold
             fontface = "bold",
             color = "black",
             fill = NA,    # White background
             label.size = NA,   # Remove label border if desired
             label.padding = unit(0., "lines"),  # Make box as tight as possible
             vjust=0) +
  # remove the minor vertical lines
  theme(panel.grid.minor.x = element_blank())+
  # make the major lines wider and grey
  theme(panel.grid.major.x = element_blank())
plot(p)

ggsave("../Figures/posterior_time_to_convergence_longerruns.pdf", p, width = 8, height = 3)



# 
# 
# posterior$header <- factor(posterior$header, levels = c("posterior", "likelihood", "prior"))
# # plot the convergence of the posterior
# p = ggplot(posterior, aes(x=samples, y=-value, color = method )) + 
#   geom_line(aes(group=interaction(method,rep)), size=0.1) +
#   labs(x = "Samples", y = "Traces after 1% burnin") +
#   facet_wrap(header~dataset, scales = "free_y") +
#   theme_minimal() +
#   scale_color_OkabeIto() +
#   # scale_x_log10() +
#   # make log scale onthe y axis, but for negative values
#   scale_y_log10()+
#   # flip the y axis
#   scale_y_reverse() 
# plot(p)
# 
# ggsave("traces_longerruns.png", p, width = 8, height = 7)
# p = ggplot(posterior[posterior$method=="Targeted" & posterior$dataset=="h3n2 5000 samples",], aes(x=samples, y=-value, color = method )) + 
#   geom_line(aes(group=interaction(method,rep)), size=0.1) +
#   labs(x = "Samples", y = "Traces after 1% burnin") +
#   facet_wrap(header~rep, scales = "free_y") +
#   theme_minimal() +
#   scale_color_manual(values=colors) +
#   # scale_x_log10() +
#   # make log scale onthe y axis, but for negative values
#   scale_y_log10()+
#   # flip the y axis
#   scale_y_reverse() 
# plot(p)
# 
# 
# posterior_combined$header <- factor(posterior_combined$header, levels = c("posterior", "likelihood", "prior"))
# 
# # make a density plot of the combined posterior
# p = ggplot(posterior_combined, aes(x=value, color = method, fill=method)) + 
#   geom_density(size=1, alpha=0.5) +
#   labs(x = "Posterior density", y = "Density") +
#   facet_wrap(header~dataset, scales = "free") +
#   theme_minimal() +
#   scale_color_OkabeIto()+
#   scale_fill_OkabeIto()
# plot(p)
# 
# # compute the ks distance for dfferent sample numbers between 
# ks_dist = data.frame()
# for (dataset in unique(posterior_combined$dataset)){
#   for (header in c("posterior", "likelihood", "prior")){
#     # get the values for the samples
#     values = posterior_combined[posterior_combined$dataset == dataset & 
#                                 posterior_combined$header == header &
#                                   posterior_combined$method == "Targeted"  , "value"]
#     values_default = posterior_combined[posterior_combined$dataset == dataset & 
#                                         posterior_combined$header == header & 
#                                         posterior_combined$method == "Default", "value"]
#     
#     min_vals = min(length(values), length(values_default))
#     x_vals = round(exp(seq(log(0.1*min_vals), log(min_vals), length.out=100)))
#     # calculate the ks_distance between values and values_default
#     for (x in x_vals){
#       ks = ks.test(values[1:x], values_default[1:x])
#       ks_dist = rbind(ks_dist, data.frame(dataset = dataset, header = header, x = x, ks = ks$statistic))
#     }
#   }
# }
# 
# # make a density plot of the combined posterior
# p = ggplot(ks_dist, aes(x=x, y=ks)) + 
#   geom_line(size=1) +
#   labs(x = "Posterior density", y = "Density") +
#   facet_wrap(header~dataset) +
#   scale_y_log10() +
#   theme_minimal()
# plot(p)
# 
# 
