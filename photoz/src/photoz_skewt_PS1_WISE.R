#------------------------------------------------------------
# Purpose: calculate the photo-z of quasars based on Pan-STARRS(PS1) and WISE photometry
# Version: July 11, 2017
# Author: Qian Yang, qianyang.astro@gmail.com
#------------------------------------------------------------

library(sn)
library(pracma)
library(data.table)
library(bit64)

file <- "./../../data/qso_ps1_wise.csv"
data <-  fread(file)
index <- data$FR>0 & data$FEG>0 & data$FER>0 & data$FEI>0 & data$FEZ>0 & data$FEY>0 & data$FEW1>0 & data$FEW2>0
dats <- data[index]

ra <- dats$RA
dec <- dats$DEC
zs <- dats$REDSHIFT
fg <- dats$FG
fr <- dats$FR
fi <- dats$FI
fz <- dats$FZ
fy <- dats$FY
fw1 <- dats$FW1
fw2 <- dats$FW2

feg <- dats$FEG
fer <- dats$FER
fei <- dats$FEI
fez <- dats$FEZ
fey <- dats$FEY
few1 <- dats$FEW1
few2 <- dats$FEW2

ebv <- dats$EBV

Arr <- c(3.172, 2.271, 1.682, 1.322, 1.087)
mg <- dats$GPSFMAG - ebv * Arr[1]
mr <- dats$RPSFMAG - ebv * Arr[2]
mi <- dats$IPSFMAG - ebv * Arr[3]

x <- mg - mi
a0 <- -0.01808
a1 <- 0.13595
a2 <- 0.01941
a3 <- -0.00183
dats$mg <- mg - (a0 + a1 * x + a2 * x^2 + a3 * x^3)
mg <- dats$mg
a0 <- -0.01836
a1 <- -0.03577
a2 <- 0.02612
a3 <- -0.00558
dats$mr <- mr - (a0 + a1 * x + a2 * x^2 + a3 * x^3)
mr <- dats$mr

nall <- length(dats$mr)

fgr <- fg / fr
fir <- fi / fr
fzr <- fz / fr
fyr <- fy / fr
fw1r <- fw1 / fr
fw21 <- fw2 / fw1

fegr <- (feg^2*fr^2 + fer^2 * fg^2)/(fr^4)
feir <- (fei^2*fr^2 + fer^2 * fi^2)/(fr^4)
fezr <- (fez^2*fr^2 + fer^2 * fz^2)/(fr^4)
feyr <- (fey^2*fr^2 + fer^2 * fy^2)/(fr^4)
few1r <- (few1^2*fr^2 + fer^2 * fw1^2)/(fr^4)
few21 <- (few2^2*fw1^2 + few1^2 * fw2^2)/(fw1^4)


fe_gi <- (fg*fi*fer^2)/(fr^4)
fe_gz <- (fg*fz*fer^2)/(fr^4)
fe_gy <- (fg*fy*fer^2)/(fr^4)
fe_gw1 <- (fg*fw1*fer^2)/(fr^4)

fe_iz <- (fi*fz*fer^2)/(fr^4)
fe_iy <- (fi*fy*fer^2)/(fr^4)
fe_iw1 <- (fi*fw1*fer^2)/(fr^4)

fe_zy <- (fz*fy*fer^2)/(fr^4)
fe_zw1 <- (fz*fw1*fer^2)/(fr^4)

fe_yw1 <- (fy*fw1*fer^2)/(fr^4)

fe_w12 <- (-fw2*few1^2)/(fr * fw1^2)

fit <- readRDS("./../model/qso_skewt_ps1_wise.rds")
redshift <- fit$redshift
numm <- length(redshift)
redarr <- fit$redarr
magarr <- fit$magarr
nmg <- length(magarr)
has <- fit$has
idarr <- seq(1, nmg, 1)

zprior <- readRDS("./../model/zprior_g.rds")
rsarr <- zprior$redshift
marr <- zprior$marr
nmg_prior <- length(marr)

photoz <- numeric(nall)
pdftotalq <- numeric(nall)
pdftotals <- numeric(nall)
photoz1 <- numeric(nall)
photoz2 <- numeric(nall)
prob <- numeric(nall)
photoz_other <- numeric(nall)
photoz1_other <- numeric(nall)
photoz2_other <- numeric(nall)
prob_other <- numeric(nall)
photoz_peak <- numeric(nall)
log_qso <- numeric(nall) + 99.0

for (i in 1:nall){
  idd <- which.min(abs(mg[i] - marr))
  if (mg[i] < 0){idd <- nmg_prior}
  zpr <- zprior[[idd]]

  color <- c(fgr[i], fir[i], fzr[i], fyr[i], fw1r[i], fw21[i])
  error <- diag(c(fegr[i], feir[i], fezr[i], feyr[i], few1r[i], few21[i]))
  error[1, 2] <- fe_gi[i]
  error[1, 3] <- fe_gz[i]
  error[1, 4] <- fe_gy[i]
  error[1, 5] <- fe_gw1[i]

  error[2, 3] <- fe_iz[i]
  error[2, 4] <- fe_iy[i]
  error[2, 5] <- fe_iw1[i]

  error[3, 4] <- fe_zy[i]
  error[3, 5] <- fe_zw1[i]

  error[4, 5] <- fe_yw1[i]

  error[5, 6] <- fe_w12[i]
  error[2, 1] <- error[1, 2]
  error[3, 1] <- error[1, 3]
  error[4, 1] <- error[1, 4]
  error[5, 1] <- error[1, 5]
  error[3, 2] <- error[2, 3]
  error[4, 2] <- error[2, 4]
  error[5, 2] <- error[2, 5]
  error[4, 3] <- error[3, 4]
  error[5, 3] <- error[3, 5]
  error[5, 4] <- error[4, 5]
  error[6, 5] <- error[5, 6]

  pdf <- double(numm+1)
  for (jn in 1:numm) {
    if (sum(has[jn, ])>0){
      id_one <- idarr[has[jn, ] == 1]
      if (length(id_one)>1){
          mag_one <- mr[i]
        id_one <- id_one[which.min(abs(magarr[id_one] - mag_one))]
      }
      id <- (jn-1)*nmg + id_one
      dpp <- fit[[id]]$dp.complete
      dpp$Omega <- dpp$Omega + error
      pdf[jn] = dmst(color, dp = dpp) * zpr[jn]
    }
  }

  pdf[is.na(pdf)] <- 0
  pdftotalq[i] <- sum(pdf)
  # max
  id <- which.max(pdf)
  photoz_peak[i] <- redshift[id]
  # ---- finding peaks in PDF -----
  peaks <- findpeaks(pdf, npeaks=2, sortstr=TRUE)
  np <- length(peaks[,1])
  if (np >0){
    if (np > 1){
      peak1 <- redshift[peaks[1, 2]]
      peak2 <- redshift[peaks[2, 2]]
      pdf1 <- sum(pdf[peaks[1, 3]: peaks[1, 4]])
      pdf2 <- sum(pdf[peaks[2, 3]: peaks[2, 4]])
      pdfarr <- c(pdf1, pdf2)
      pdf_peak <- max(pdfarr)
      prob[i] <- pdf_peak/pdftotalq[i]
      idpeak <- which.max(pdfarr)
      photoz[i] <- redshift[peaks[idpeak, 2]]
      photoz1[i] <- redshift[peaks[idpeak, 3]]
      photoz2[i] <- redshift[peaks[idpeak, 4]]
      #
      pdf_other <- min(pdfarr)
      idother <- which.min(pdfarr)
      photoz_other[i] <- redshift[peaks[idother, 2]]
      photoz1_other[i] <- redshift[peaks[idother, 3]]
      photoz2_other[i] <- redshift[peaks[idother, 4]]
      prob_other[i] <- pdf_other/pdftotalq[i]
    } else {
      photoz[i] <- redshift[peaks[2]]
      prob[i] <- 1.0
      photoz1[i] <- redshift[peaks[3]]
      photoz2[i] <- redshift[peaks[4]]
    }
  }
  if (pdftotalq[i]>0){
    log_qso[i] <- log10(pdftotalq[i])
  }
}
res <- cbind(ra, dec, zs, mg, mr, photoz, pdftotalq, photoz1, photoz2, prob, photoz_other, photoz1_other, photoz2_other, prob_other, log_qso)
write.table(res, file = "./../result/photoz_skewt_PS1_WISE.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
