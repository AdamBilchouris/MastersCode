sunspotsDat <- read.csv("SN_y_tot_V2.0.csv", sep=';')[, 1:2]
str(sunspotsDat)
colnames(sunspotsDat) <- c("year", "sunspots")
str(sunspotsDat)
plot(sunspotsDat$year, sunspotsDat$sunspots, type='l', xlab="Year", ylab="Yearly Mean Sunspots")

# periodogram
sunspotsPer <- spectrum(sunspotsDat$sunspots)
plot(sunspotsPer$freq, sunspotsPer$spec, xlab='Frequency', ylab='Spectrum', main='Periodogram of Yearly Mean Sunspots', type='l')
1 / (sunspotsPer$freq[which.max(sunspotsPer$spec)])
