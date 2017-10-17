#------------------------------------------------------------
# Purpose: calculate the photo-z of quasars based on DECaLS and WISE photometry
# Version: July 11, 2017
# Author: Qian Yang, qianyang.astro@gmail.com
#------------------------------------------------------------

library(sn)
library(pracma)
library(data.table)
library(bit64)

file <- "./../../data/qso_decals_wise.csv"
data <-  fread(file)
index <- data$DEPTH_G>0 & data$DEPTH_R>0 & data$DEPTH_Z>0 & data$FR>0
dats <- data[index]

ra <- dats$RA
dec <- dats$DEC
zs <- dats$REDSHIFT
fg <- dats$FG
fr <- dats$FR
fz <- dats$FZ
fw1 <- dats$FW1
fw2 <- dats$FW2
fw3 <- dats$FW3
fw4 <- dats$FW4

feg <- dats$FEG
fer <- dats$FER
fez <- dats$FEZ
few1 <- dats$FEW1
few2 <- dats$FEW2
few3 <- dats$FEW3
few4 <- dats$FEW4

mg <- dats$MG
mr <- dats$MR

nall <- length(dats$MR)

fgr <- fg / fr
fzr <- fz / fr
fw1r <- fw1 / fr
fw21 <- fw2 / fw1
fw31 <- fw3 / fw1
fw41 <- fw4 / fw1
fegr <- (feg^2*fr^2 + fer^2 * fg^2)/(fr^4)
fezr <- (fez^2*fr^2 + fer^2 * fz^2)/(fr^4)
few1r <- (few1^2*fr^2 + fer^2 * fw1^2)/(fr^4)
few21 <- (few2^2*fw1^2 + few1^2 * fw2^2)/(fw1^4)
few31 <- (few3^2*fw1^2 + few1^2 * fw3^2)/(fw1^4)
few41 <- (few4^2*fw1^2 + few1^2 * fw4^2)/(fw1^4)

fe_gz <- (fg*fz*fer^2)/(fr^4)
fe_gw1 <- (fg*fw1*fer^2)/(fr^4)

fe_zw1 <- (fz*fw1*fer^2)/(fr^4)

fe_w12 <- (-few1^2*fw2)/(fr*fw1^2)
fe_w13 <- (-few1^2*fw3)/(fr*fw1^2)
fe_w14 <- (-few1^2*fw4)/(fr*fw1^2)

fe_w23 <- (fw2*fw3*few1^2)/(fw1^4)
fe_w24 <- (fw2*fw4*few1^2)/(fw1^4)

fe_w34 <- (fw3*fw4*few1^2)/(fw1^4)

fit <- readRDS("./../model/qso_skewt_decals_wise.rds")
redshift <- fit$redshift
numm <- length(redshift)
redarr <- fit$redarr
magarr <- fit$magarr
nmg <- length(magarr)
has <- fit$has
hasall <- fit$hasall

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

  color <- c(fgr[i], fzr[i], fw1r[i], fw21[i])
  error <- diag(c(fegr[i], fezr[i], few1r[i], few21[i]))
  error[1, 2] <- fe_gz[i]
  error[1, 3] <- fe_gw1[i]
  error[2, 3] <- fe_zw1[i]
  error[3, 4] <- fe_w12[i]
  error[2, 1] <- error[1, 2]
  error[3, 1] <- error[1, 3]
  error[3, 2] <- error[2, 3]
  error[4, 3] <- error[3, 4]

  pdf <- double(numm+1)
  for (jn in 1:105) {
    if (hasall[jn]>0){
      if (sum(has[jn, ])>0){
        id_one <- idarr[has[jn, ] == 1]
        if (length(id_one)>1){
            mag_one <- mr[i]
          id_one <- id_one[which.min(abs(magarr[id_one] - mag_one))]
        }
        id <- (jn-1)*(nmg + 1) + id_one
      } else{
        # print("hasall")
        id <- jn * (nmg + 1)
      }
      dpp <- fit[[id]]$dp.complete
      dpp$Omega <- dpp$Omega + error
      pdf[jn] <- dmst(color, dp = dpp) * zpr[jn]
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
write.table(res, file = "./../result/photoz_skewt_DECaLS_WISE.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
