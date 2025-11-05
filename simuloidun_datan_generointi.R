library(ggseqplot)
library(dplyr)
library(seqHMM)
library(TraMineR)
library(forcats)

d <- readRDS(file = "sla_data.rds") |> 
  mutate(
    id = as.factor(id),
    aad_std = (aad - 40) / 10,
    status = fct_recode(
      status,
      Töissä = "working",
      Sairauslomalla = "sick leave",
      Muu = "other"
    ),
    lääkitys = fct_recode(
      medication,
      "Ei lääkitystä" = "no", Lääkitys = "yes"
    ),
    y = interaction(status, lääkitys, sep = " / ", lex.order = TRUE)
  ) |> select(-medication) |> 
  filter(time == 0)

fit <- readRDS("fit_mhmm_K4.rds")
mpc <- most_probable_cluster(fit)
d$cluster <- mpc

# kovariaattien klusterikohtaiset (marginaali)jakaumat:

(p_history <- prop.table(table(d$cluster, d$history), margin = 1))
(p_comorbidity <- prop.table(table(d$cluster, d$comorbidity), margin = 1))
(p_female <- prop.table(table(d$cluster, d$female), margin = 1))
(p_education <- prop.table(table(d$cluster, d$education), margin = 1))
(p_aad_std <- prop.table(table(d$cluster, d$aad_std), margin = 1))

history <- factor(colnames(p_history), levels = levels(d$history))
comorbidity <- 0:1
female <- factor(colnames(p_female), levels = levels(d$female))
education <- factor(colnames(p_education), levels = levels(d$education))

set.seed(1)
dx <- bind_rows(
  lapply(1:4, \(i) {
    data.frame(
      id = (i - 1) * 500 +  1:500, 
      history = sample(history, size = 500, replace = TRUE, prob = p_history[i, ]),
      comorbidity = sample(comorbidity, size = 500, replace = TRUE, prob = p_comorbidity[i, ]),
      female = sample(female, size = 500, replace = TRUE, prob = p_female[i, ]),
      education = sample(education, size = 500, replace = TRUE, prob = p_education[i, ]),
      aad_std = sample(seq(-1,1, by = 0.1), size = 500, replace = TRUE, prob = p_aad_std[i, ])
    )
  }
  ),
  .id = "cluster"
)

# alkuperäinen aineisto ositettu mallinnuksen perusteella neljään osaan (klusteriin)
seqdata_mhmm_K4 <- readRDS("W:/JouniH/sla/seqdata_mhmm_K4.rds")

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

# nimet
obs <- alphabet(seqdata_mhmm_K4[[1]]$y)
# simuloi 4x500 sekvenssiä
dy <- bind_rows(
  lapply(1:4, \(i) {
    y <- replicate(500, simuloi_sekvenssi(markov_mallit[[i]]))
    data.frame(
      time = rep(0:59, each = 500),
      id = (i - 1) * 500 +  1:500, 
      y = c(t(y))
    )
  }), 
  .id = "cluster"
)
dy <- dy |> mutate(y = obs[y]) |> 
  tidyr::separate_wider_delim(y, " / ", names = c("status", "lääkitys")) |> 
  mutate(
    status = factor(status, levels = c("Töissä", "Sairauslomalla", "Muu")),
    lääkitys = factor(lääkitys, levels = c("Ei lääkitystä", "Lääkitys"))
  )
# yhdistetään simuloidut kovariaatit ja sekvenssit
dsim <- left_join(dy, dx, by = c("id", "cluster"))
# tallennetaan simuloitu aineisto
saveRDS(dsim, file = "simuloidut_sekvenssit.rds")
