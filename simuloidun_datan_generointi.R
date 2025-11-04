library(dplyr)
library(TraMineR)
library(ggseqplot)

# alkuperäinen aineisto ositettu mallinnuksen perusteella neljään osaan (klusteriin)
seqdata_mhmm_K4 <- readRDS("seqdata_mhmm_K4.rds")

# lasketaan jokaisesta klusterista keskimääräiset kuukausittaiset 
# siirtymätodennäköisyydet havaintokategorioiden välillä
markov_mallit <- lapply(
  seqdata_mhmm_K4, \(x) {
    list(pi = prop.table(table(x$y[, 1])[1:6]), A = seqtrate(x$y, time.varying = TRUE))
  }
)
# funktio joka simuloi yllä saatujen todennäköisyyksien perusteella uusia havaintoja
simuloi_sekvenssi <- function(malli) {
  pi <- malli$pi
  A <- malli$A
  n <- dim(A)[3] + 1
  m <- nrow(A)
  y <- numeric(n)
  y[1] <- sample(m, size = 1, prob = pi)
  for (i in 2:n) {
    y[i] <- sample(m, size = 1, prob = A[y[i - 1], , i - 1])
  }
  y
}

# nimet ja värit
obs <- alphabet(seqdata_mhmm_K4[[1]]$y)
cols <- cpal(seqdata_mhmm_K4[[1]]$y)
# simuloi 4x500 sekvenssiä
set.seed(1)
simuloidut_sekvenssit <- lapply(markov_mallit, \(x) {
  y <- replicate(500, simuloi_sekvenssi(x))
  seqdef(t(y), states = obs, cpal = cols)
})

# tallennetaan simuloitu aineisto
saveRDS(simuloidut_sekvenssit, file = "simuloidut_sekvenssit.rds")
