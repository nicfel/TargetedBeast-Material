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
logFiles <- list.files(path = "out", pattern = "\\.log", full.names = TRUE)

data_convergence = data.frame()

# leafs = c("infB"=1978, "h3n2"=5000, "wnv"=2664, "mpxv"=1579)

# remove any files that contain Default
logFiles <- logFiles[!grepl("Default", logFiles)]

# loop over the files
for (logFile in logFiles) {
  # remove out/, .log and replace _Default_ with _ a
  name <- str_replace(logFile, "out/", "")
  name <- str_replace(name, "\\.log", "")
  # split on _
  name <- str_split(name, "_")[[1]]
  # read the log file
  log <- read.table(logFile, header = TRUE, sep="\t")
  # remove the first 10 %
  log <- log[round(nrow(log)*burnin):nrow(log),]
  
  # read in the default file for comparison
  log.default <- read.table(str_replace(logFile, name[[3]], "Default"), header = TRUE, sep="\t")
  log.default <- log.default[round(nrow(log.default)*burnin):nrow(log.default),]

  # get the amount of samples in the log file
  samples = max(log$Sample) - min(log$Sample)
  if (samples<10000){
    next
  }
  sample.default = max(log.default$Sample) - min(log.default$Sample)
  if (sample.default<10000){
    next
  }
    
  precentage = name[[4]]
  # multiply the percentage by the leafs
  totalSamples = as.numeric(precentage)
  
  # calculate the ESS value for the posterior, likelihood and prior
  for (header in c("posterior", "likelihood", "prior")) {
    ess <- ess_tracer_style(log[,header])
    ess.default <- ess_tracer_style(log.default[,header])
    
    if (ess<10 || ess.default<10){
      next
    }
    
    ratio = (ess/samples) / (ess.default/sample.default)

    # make a text for the actual ESS values plotted with one digit

    # add the ESS values to the data frame
    data_convergence <- rbind(data_convergence, 
                              data.frame(dataset = name[[1]], 
                                         samples = samples,
                                         leaves = totalSamples,
                                         header = header, 
                                         method=name[[3]],
                                         ess = ratio))
  }
}

# reorder facets for header
data_convergence$header <- factor(data_convergence$header, levels = c("posterior", "likelihood", "prior"))
# plot the ESS values
p = ggplot(data_convergence, aes(x=leaves, y=ess, color = method, fill=method)) + 
  geom_point(aes(shape=dataset), size=2) +
  facet_grid(header~.) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_color_manual(values=c(default, targeted)) +
  labs(x = "number of samples", y = "ESS per sample", fill = "Type") +
  # geom_text(aes(label = text), vjust = 0) +
  facet_grid(header~.) +
  scale_y_log10() +
  scale_x_log10()+
  geom_smooth(se=T, method="lm") +
  theme_minimal() +
  scale_color_manual(values=c("Targeted" = "#56B4E9", "Intervals" = "#009E73")) +
  scale_fill_manual(values=c("Targeted" = "#56B4E9", "Intervals" = "#009E73")) +
  coord_cartesian(ylim = c(1, 20))
  # theme(legend.position = "none") +
plot(p)
ggsave("convergence.png", p, width = 12, height = 7, units = "cm")

