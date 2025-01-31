library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)

# Clear workspace
rm(list=ls())

burnin = 0.0911

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

targeted="#d95f02"
default="#1b9e77"


ess_tracer_style <- function(x, max_lag = NULL) {
  # x: a numeric vector representing a single MCMC chain
  # max_lag: optional maximum lag to consider. If NULL, use full length minus 1.
  
  N <- length(x)
  # Subtract the mean (demean)
  x_demeaned <- x - mean(x)
  
  # If not specified, set max_lag to something sensible (e.g., N-1 or some fraction of N)
  if (is.null(max_lag)) {
    max_lag <- min(N - 1, 10000)  # Usually you wouldn't need lags beyond a few thousand
  }
  
  # We can use the 'acf' function to get autocorrelations up to max_lag
  # 'type = "correlation"' returns the usual ACF,
  # 'plot = FALSE' to suppress plots,
  # 'demean = FALSE' because we've already demeaned x.
  acf_vals <- acf(x_demeaned, lag.max = max_lag, plot = FALSE, 
                  type = "correlation", demean = FALSE)$acf
  
  # 'acf_vals' is an array of length max_lag + 1, with acf_vals[1] = lag 0 correlation = 1
  # We want acf_vals for lag = 1..max_lag:
  acf_lags <- acf_vals[-1]  # remove the lag 0 element
  
  # Sum autocorrelations while they remain positive
  pos_sum <- 0
  for (r in acf_lags) {
    if (r <= 0) break
    pos_sum <- pos_sum + r
  }
  
  # Calculate the ESS from the sum of positive autocorrelations
  ess <- N / (1 + 2 * pos_sum)
  return(ess)
}


# get all log files in out that contain default
logFiles <- list.files(path = "out", pattern = ".log", full.names = TRUE)
logFiles <- logFiles[str_detect(logFiles, "Default")]

data_convergence = data.frame()

leafs = c("infB"=1978, "h3n2"=5000)

# loop over the files
for (logFile in logFiles) {
  # remove out/, .log and replace _Default_ with _ a
  name <- str_replace(logFile, "out/", "")
  name <- str_replace(name, ".log", "")
  name <- str_replace(name, "Default_", "")
  # split on _
  name <- str_split(name, "_")[[1]]
  # read the log file
  log.default <- read.table(logFile, header = TRUE, sep="\t")
  # remove the first 10 %
  log.default <- log.default[round(nrow(log.default)*burnin):nrow(log.default),]
  # read in the corresponding targeted file
  log.targeted <- read.table(gsub("Default", "Targeted", logFile), header = TRUE, sep="\t")
  # remove the first 10 %
  log.targeted <- log.targeted[round(nrow(log.targeted)*burnin):nrow(log.targeted),]
  
  
  # get the amount of samples in the log file
  default.samples = max(log.default$Sample) - min(log.default$Sample)
  targeted.samples = max(log.targeted$Sample) - min(log.targeted$Sample)
  precentage = name[[3]]
  # multiply the percentage by the leafs
  totalSamples = leafs[name[[1]]] * as.numeric(precentage)
  
  # calculate the ESS value for the posterior, likelihood and prior
  for (header in c("posterior", "likelihood", "prior")) {
    ess.default <- ess_tracer_style(log.default[,header])
    ess.targeted <- ess_tracer_style(log.targeted[,header])
    
    # make a text for the actual ESS values plotted with one digit
    text = paste0("D: ", round(ess.default, 1), "\nT: ", round(ess.targeted, 1))
    
    # add the ESS values to the data frame
    data_convergence <- rbind(data_convergence, 
                              data.frame(dataset = name[[1]], 
                                         samples = totalSamples,
                                         text = text,
                                         header = header, 
                                         ratio = (ess.targeted/targeted.samples)/(ess.default/default.samples)))
  }
}

# reorder facets for header
data_convergence$header <- factor(data_convergence$header, levels = c("posterior", "likelihood", "prior"))
# plot the ESS values
p = ggplot(data_convergence[data_convergence$header=="posterior",], aes(x=samples, y=ratio, color = dataset)) + 
  geom_point() +
  facet_grid(header~.) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=c(default, targeted)) +
  labs(x = "number of samples", y = "Ratio of ESS\nper MCMC step using\ntargeted over\nuntargeted\nMCMC moves", fill = "Type") +
  # geom_text(aes(label = text), vjust = 0) +
  theme_minimal()
  # theme(legend.position = "none") +
  # scale_y_log10(limit = c(1, 30))
plot(p)
ggsave("convergence.png", p, width = 10, height = 5, units = "cm")

