#------------------------------------------------------------
# Purpose: calculate the photo-z of quasars based on SDSS photometry
# Version: July 11, 2017
# Author: Qian Yang, qianyang.astro@gmail.com
#------------------------------------------------------------

library(sn)
library(pracma)
library(data.table)
library(bit64)

file <- "./../../data/qso_sdss.csv"
data <-  fread(file)
index <- data$FR>0 & data$FEU>0 & data$FEG>0 & data$FER>0 & data$FEI>0 & data$FEZ>0
dats <- data[index]
exr <- c(5.155, 3.793, 2.751, 2.086, 1.479)
exr <- exr/exr[1]

ra <- dats$RA
dec <- dats$DEC
zs <- dats$REDSHIFT
fu <- dats$FU
fg <- dats$FG
fr <- dats$FR
fi <- dats$FI
fz <- dats$FZ
feu <- dats$FEU
feg <- dats$FEG
fer <- dats$FER
fei <- dats$FEI
fez <- dats$FEZ
dats$mg <- dats$PSFMAG_G - exr[2] * dats$EXTINCTION_U
dats$mr <- dats$PSFMAG_R - exr[3] * dats$EXTINCTION_U
mg <- dats$mg
mr <- dats$mr
nall <- length(dats$mr)

fur <- fu / fr
fgr <- fg / fr
fir <- fi / fr
fzr <- fz / fr
feur <- (feu^2*fr^2 + fer^2 * fu^2)/(fr^4)
fegr <- (feg^2*fr^2 + fer^2 * fg^2)/(fr^4)
feir <- (fei^2*fr^2 + fer^2 * fi^2)/(fr^4)
fezr <- (fez^2*fr^2 + fer^2 * fz^2)/(fr^4)

fe_ug <- (fu*fg*fer^2)/(fr^4)
fe_ui <- (fu*fi*fer^2)/(fr^4)
fe_uz <- (fu*fz*fer^2)/(fr^4)
fe_gi <- (fg*fi*fer^2)/(fr^4)
fe_gz <- (fg*fz*fer^2)/(fr^4)
fe_iz <- (fi*fz*fer^2)/(fr^4)

fit <- readRDS("./../model/qso_skewt_sdss.rds")
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
  zpr <- zprior[[idd]]

  color <- c(fur[i], fgr[i], fir[i], fzr[i])
  error <- diag(c(feur[i], fegr[i], feir[i], fezr[i]))
  error[1, 2] <- fe_ug[i]
  error[1, 3] <- fe_ui[i]
  error[1, 4] <- fe_uz[i]
  error[2, 3] <- fe_gi[i]
  error[2, 4] <- fe_gz[i]
  error[3, 4] <- fe_iz[i]
  error[2, 1] <- error[1, 2]
  error[3, 1] <- error[1, 3]
  error[4, 1] <- error[1, 4]
  error[3, 2] <- error[2, 3]
  error[4, 2] <- error[2, 4]
  error[4, 3] <- error[3, 4]

  pdf <- double(numm+1)
  for (jn in 1:numm) {
    if (sum(has[jn, ])>0){
      id_one <- idarr[has[jn, ] == 1]
      if (length(id_one)>1){
        mag_one <- mr[i]
        id_one <- id_one[which.min(abs(magarr[id_one] - mag_one))]
      }
      id <- (jn-1)*nmg + id_one
      dpp <- fit[[id]]$dp
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
write.table(res, file = "./../result/photoz_SDSS.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
