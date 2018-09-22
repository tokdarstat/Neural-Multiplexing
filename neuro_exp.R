source("dynamic_neural_model.R")

synth.data <- rsynth(ntrials = c(15, 20, 20), pr.flat = 1, intervals = list(c(0,.1), c(.45,.55), c(.9,1)), wts = c(1/3, 1/3, 1/3), span = c(.1,.9), period = c(500, 1500))

spike.counts <- list()
breaks <- seq(0, 1e3, 25)
spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)
spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
spike.counts$bin.width <- diff(breaks)[1]

fit.post <- dynamic.model.fit(spike.counts, plot = TRUE, lengthScale = c(75, 125, 200, 300))
plot(fit.post, true.alphas = synth.data$alphas)
plot(fit.post, true.alphas = synth.data$alphas, tilt = TRUE)
print(summary(fit.post, tilt = TRUE))
