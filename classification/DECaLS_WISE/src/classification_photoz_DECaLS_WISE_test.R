# PURPOSE: classify quasars and calculate the photo-z of quasars.
# AUTHOR: Qian Yang
# E-MAIL: qianyang.astro@gmail.com
# VERSION: July 11, 2017
# CRITERIA (quasar candidate): log_qso>-2 & (log_qso - log_star>4.4|log_star>90) & (log_qso-log_galaxy>5.4|log_galaxy>90)

library(sn)
library(mvtnorm)
library(pracma)
library(data.table)
library(bit64)

file <- "./../../../data/decals_S82_test.csv"
data <-  fread(file)
index <- data$BRICK_PRIMARY & data$DEPTH_G>0 & data$DEPTH_R>0 & data$DEPTH_Z>0 & data$FR>0 & data$FEW1>0 & data$FEW2>0

dats <- data[index]

ra <- dats$RA
dec <- dats$DEC
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

fg_ap <- dats$FG_AP
fr_ap <- dats$FR_AP
fz_ap <- dats$FZ_AP
fw1_ap <- dats$FW1_AP
fw2_ap <- dats$FW2_AP

feg_ap <- dats$FEG_AP
fer_ap <- dats$FER_AP
fez_ap <- dats$FEZ_AP
few1_ap <- dats$FEW1_AP
few2_ap <- dats$FEW2_AP

mg <- dats$MG
mr <- dats$MR
mz <- dats$MZ

mg_ap <- dats$MG_AP
mr_ap <- dats$MR_AP
mz_ap <- dats$MZ_AP

nall <- length(dats$MR)

# w2/w1 (QSO variability when using w2/r)
fgr <- fg / fr
fzr <- fz / fr
fw1r <- fw1 / fr
fw21 <- fw2 / fw1
fegr <- (feg^2*fr^2 + fer^2 * fg^2)/(fr^4)
fezr <- (fez^2*fr^2 + fer^2 * fz^2)/(fr^4)
few1r <- (few1^2*fr^2 + fer^2 * fw1^2)/(fr^4)
few21 <- (few2^2*fw1^2 + few1^2 * fw2^2)/(fw1^4)
fe_gz <- (fg*fz*fer^2)/(fr^4)
fe_gw1 <- (fg*fw1*fer^2)/(fr^4)
fe_zw1 <- (fz*fw1*fer^2)/(fr^4)
fe_w12 <- (-few1^2*fw2)/(fr*fw1^2)

# w2/r (smaller r data uncertainty)
fw2r <- fw2 / fr
few2r <- (few2^2*fr^2 + fer^2 * fw2^2)/(fr^4)
fe_gw2 <- (fg*fw2*fer^2)/(fr^4)
fe_zw2 <- (fz*fw2*fer^2)/(fr^4)
fe_w12_2 <- (fw1*fw2*fer^2)/(fr^4)

# ap
fgr_ap <- fg_ap/fr_ap
fzr_ap <- fz_ap/fr_ap
fw1r_ap <- fw1_ap/fr_ap
fw2r_ap <- fw2_ap/fr_ap
fegr_ap <- (feg_ap^2*fr_ap^2 + fer_ap^2 * fg_ap^2)/(fr_ap^4)
fezr_ap <- (fez_ap^2*fr_ap^2 + fer_ap^2 * fz_ap^2)/(fr_ap^4)
few1r_ap <- (few1_ap^2*fr_ap^2 + fer_ap^2 * fw1_ap^2)/(fr_ap^4)
few2r_ap <- (few2_ap^2*fr_ap^2 + fer_ap^2 * fw2_ap^2)/(fr_ap^4)

# ------------------------------------------------------------------
# model QSO
fitq <- readRDS("./../model/qso_skewt_decals_wise.rds")
redshift <- fitq$redshift
numm <- length(redshift)
redarr <- fitq$redarr
magarr <- fitq$magarr
nmg <- length(magarr)
has <- fitq$has
hasall <- fitq$hasall

# QSO priority from QLF (NPD et al.)
idarr <- seq(1, nmg, 1)
zpriorq <- readRDS("./../model/zprior_g.rds")
marrq <- zpriorq$marr
nmg_prior <- length(marrq)

# star colors and star counts from Besancon model and Galaxia code
fits <- readRDS("./../model/star_decam_wise.rds")
marrs <- fits$marr

# galaxy colors from templates in Brown et al.
fitg <- readRDS("./../model/galaxy_decam_wise.rds")
types <- fitg$types
names <- fitg$names
ng <- length(names)

# galaxy priority from BPZ and Size distribution (Shen et al.)
priorg <- readRDS("./../model/galaxy_zprior_decam.rds")
marrg <- priorg$marr
redshiftg <- priorg$redshift
nzg <- length(redshiftg)
# ------------------------------------------------------------------
photoz <- numeric(nall)
pdftotalq <- numeric(nall)
pdftotals <- numeric(nall)
# probability (ratio in PDF of this peak within photo-z -p1 +p2)
# photo-z error
photoz1 <- numeric(nall)
photoz2 <- numeric(nall)
prob <- numeric(nall)
photoz_other <- numeric(nall)
photoz1_other <- numeric(nall)
photoz2_other <- numeric(nall)
prob_other <- numeric(nall)
photoz_peak <- numeric(nall)

pdftotals <- numeric(nall)

pdftotalg <- numeric(nall)
type <- array("", dim = nall)
photoz_galaxy <- numeric(nall)
log_qso <- numeric(nall) + 99.0
log_star <- numeric(nall) + 99.0
log_galaxy <- numeric(nall) + 99.0

for (i in 1:nall){
  # ========================= star =============================
  idd <- which.min(abs(mr_ap[i] - marrs))
  color <- c(fgr_ap[i], fzr_ap[i], fw1r_ap[i], fw2r_ap[i])
  error <- diag(c(fegr_ap[i], fezr_ap[i], few1r_ap[i], few2r_ap[i]))
  ps <- dmvnorm(fits[[idd]], color, error, log = FALSE)
  pdftotals[i] <- sum(ps)/20.0
  if (pdftotals[i]>1e-323){
    log_star[i] <- log10(pdftotals[i])
  }
  # ========================= galaxy =============================
  idd <- which.min(abs(mr[i] - marrg))
  frac <- 0.0
  pp <- numeric(ng)
  zz <- numeric(ng)
  for (j in 1:ng){
    if (types[j] == "E"){
      frac <- 0.2 # 5 E galaxy templates
    }
    if (types[j] == "Spiral"){
      frac <- 0.1 # 10 Spiral galaxy templates
    }
    if (types[j] == "Irr"){
      frac <- 1.0/3.0 # 3 Irr galaxy templates
    }
    zpr <- priorg[[j]][idd,]
    ft <- fitg[[j]]
    cs <- cbind(ft$fgr, ft$fzr, ft$fw1r, ft$fw2r)
    pt <- dmvnorm(cs, color, error, log = FALSE)
    p <- pt * zpr * frac
    ind <- which.max(p)
    zz[j] <- redshiftg[ind]
    pp[j] <- sum(p)
  }
  ind <- which.max(pp)
  type[i] <- types[ind]
  photoz_galaxy[i] <- zz[ind]
  pdftotalg[i] <- sum(pp)
  if (pdftotalg[i]>1e-323){
    log_galaxy[i] <- log10(pdftotalg[i])
  }
  # ================= QSO =====================
  idd <- which.min(abs(mg[i] - marrq))
  if (mg[i] < 0){idd <- nmg_prior}
  zpr <- zpriorq[[idd]]

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
  for (jn in 1:numm) {
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
      dpp <- fitq[[id]]$dp.complete
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

res <- cbind(ra, dec, mg, mr, mz, photoz, pdftotalq, photoz1, photoz2, prob, photoz_other, photoz1_other, photoz2_other, prob_other, pdftotals, type, photoz_galaxy, pdftotalg, log_qso, log_star, log_galaxy)

write.table(res, file = "./../result/photoz_qso_decals_wise_test.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
