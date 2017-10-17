;+
; Name: flux_decals_wise
; Purpose: Calculate DECaLS (DR3) and unWISE flux
; Input: tractor fits file
; Output: output fits file including
;    ra, dec, shape (DCHISQ[0, *] - DCHISQ[1, *]),
;    ng, nr, nz, nw1, nw2, nw3, nw4,
;    BRICK_PRIMARY, depth_g, depth_r, depth_z,
;    mg_ap, mr_ap, mz_ap, (magnitude without galactic extinction correction)
;    mw1_ap, mw2_ap, mw3_ap, mw4_ap,
;    mg, mr, mz, (magnitude with galactic extinction correction)
;    mw1, mw2, mw3, mw4
;    meg, mer, mez, mew1, mew2, mew3, mew4 (magnitude and magnitude error)
;    fg_ap, fr_ap, fz_ap, (flux without galactic extinction correction)
;    feg_ap, fer_ap, fez_ap, (flux error without galactic extinction correction)
;    fw1_ap, fw2_ap, fw3_ap, fw4_ap,
;    few1_ap, few2_ap, few3_ap, fw4_ap
;    fg, fr, fz, (flux with galactic extinction correction)
;    feg, fer, fez, (flux error with galactic extinction correction)
;    fw1, fw2, fw3, fw4, few1, few2, few3, few4
; Version:
;  July 11, 2017
; Author: Qian Yang, qianyang.astro@gmail.com
;-
;------------------------------------------------------------

pro flux_decals_wise
file = 'tractor-0001m002.fits' ;tractor file in DECaLS DR3
data = mrdfits(file, 1)
DECAM_NOBS = data.DECAM_NOBS

ra = data.ra
dec = data.dec
BRICK_PRIMARY = data.BRICK_PRIMARY
depth = data.decam_depth
DECAM_FLUX = data.DECAM_FLUX
DECAM_FLUX_IVAR = data.DECAM_FLUX_IVAR
DECAM_MW_TRANSMISSION = data.DECAM_MW_TRANSMISSION
WISE_FLUX = data.WISE_FLUX
WISE_FLUX_IVAR = data.WISE_FLUX_IVAR
WISE_MW_TRANSMISSION = data.WISE_MW_TRANSMISSION
DCHISQ = data.DCHISQ
DECAM_NOBS = data.DECAM_NOBS
WISE_NOBS = data.WISE_NOBS
f0 = 10.d^(-(48.6+22.5)/2.5) ;erg/cm2/Hz
;
shape = DCHISQ[0, *] - DCHISQ[1, *]
ng = DECAM_NOBS[1, *]
nr = DECAM_NOBS[2, *]
nz = DECAM_NOBS[4, *]
fg = DECAM_FLUX[1, *]
fr = DECAM_FLUX[2, *]
fz = DECAM_FLUX[4, *]
depth_g = depth[1, *]
depth_r = depth[2, *]
depth_z = depth[4, *]
tg = DECAM_MW_TRANSMISSION[1, *]
tr = DECAM_MW_TRANSMISSION[2, *]
tz = DECAM_MW_TRANSMISSION[4, *]
ivar2_g = DECAM_FLUX_IVAR[1, *]
ivar2_r = DECAM_FLUX_IVAR[2, *]
ivar2_z = DECAM_FLUX_IVAR[4, *]
;
mg_ap = 22.5 - 2.5 * alog10(fg)
mr_ap = 22.5 - 2.5 * alog10(fr)
mz_ap = 22.5 - 2.5 * alog10(fz)
feg = 1.d/sqrt(ivar2_g)
fer = 1.d/sqrt(ivar2_r)
fez = 1.d/sqrt(ivar2_z)
meg = (2.5 / alog(10.d)) * (feg/fg)
mer = (2.5 / alog(10.d)) * (fer/fr)
mez = (2.5 / alog(10.d)) * (fez/fz)
fg_ap = fg * f0
fr_ap = fr * f0
fz_ap = fz * f0
feg_ap = feg * f0
fer_ap = fer * f0
fez_ap = fez * f0
fg = fg/tg
fr = fr/tr
fz = fz/tz
mg = 22.5 - 2.5 * alog10(fg)
mr = 22.5 - 2.5 * alog10(fr)
mz = 22.5 - 2.5 * alog10(fz)
fg = fg * f0
fr = fr * f0
fz = fz * f0
feg = feg/tg * f0
fer = fer/tr * f0
fez = fez/tz * f0
ind = where((feg le 0) or (depth_g le 0), nt)
if (nt gt 0) then begin
  feg_ap[ind] = -1.0
  feg[ind] = -1.0
  mg[ind] = -1.0
  meg[ind] = -1.0
  mg_ap[ind] = -1.0
endif
ind = where((fer le 0) or (depth_r le 0), nt)
if (nt gt 0) then begin
  fer_ap[ind] = -1.0
  fer[ind] = -1.0
  mr[ind] = -1.0
  mer[ind] = -1.0
  mr_ap[ind] = -1.0
endif
ind = where((fez le 0) or (depth_z le 0), nt)
if (nt gt 0) then begin
  fez_ap[ind] = -1.0
  fez[ind] = -1.0
  mz[ind] = -1.0
  mez[ind] = -1.0
  mz_ap[ind] = -1.0
endif
;
ind = where((fg le 0), nt)
if (nt gt 0) then begin
  mg[ind] = -1.0
  meg[ind] = -1.0
  mg_ap[ind] = -1.0
endif
ind = where((fr le 0), nt)
if (nt gt 0) then begin
  mr[ind] = -1.0
  mer[ind] = -1.0
  mr_ap[ind] = -1.0
endif
ind = where((fz le 0), nt)
if (nt gt 0) then begin
  mz[ind] = -1.0
  mez[ind] = -1.0
  mz_ap[ind] = -1.0
endif
;
nw1 = WISE_NOBS[0, *]
nw2 = WISE_NOBS[1, *]
nw3 = WISE_NOBS[2, *]
nw4 = WISE_NOBS[3, *]
fw1 = WISE_FLUX[0, *]
fw2 = WISE_FLUX[1, *]
fw3 = WISE_FLUX[2, *]
fw4 = WISE_FLUX[3, *]
tw1 = WISE_MW_TRANSMISSION[0, *]
tw2 = WISE_MW_TRANSMISSION[1, *]
tw3 = WISE_MW_TRANSMISSION[2, *]
tw4 = WISE_MW_TRANSMISSION[3, *]
ivar2_w1 = WISE_FLUX_IVAR[0, *]
ivar2_w2 = WISE_FLUX_IVAR[1, *]
ivar2_w3 = WISE_FLUX_IVAR[2, *]
ivar2_w4 = WISE_FLUX_IVAR[3, *]
few1 = 1.d/sqrt(ivar2_w1)
few2 = 1.d/sqrt(ivar2_w2)
few3 = 1.d/sqrt(ivar2_w3)
few4 = 1.d/sqrt(ivar2_w4)
mew1 = (2.5 / alog(10.d)) * (few1/fw1)
mew2 = (2.5 / alog(10.d)) * (few2/fw2)
mew3 = (2.5 / alog(10.d)) * (few3/fw3)
mew4 = (2.5 / alog(10.d)) * (few4/fw4)
mw1_ap = 22.5 - 2.5 * alog10(fw1)
mw2_ap = 22.5 - 2.5 * alog10(fw2)
mw3_ap = 22.5 - 2.5 * alog10(fw3)
mw4_ap = 22.5 - 2.5 * alog10(fw4)
fw1_ap = fw1 * f0
fw2_ap = fw2 * f0
fw3_ap = fw3 * f0
fw4_ap = fw4 * f0
few1_ap = few1 * f0
few2_ap = few2 * f0
few3_ap = few3 * f0
few4_ap = few4 * f0
fw1 = fw1/tw1
fw2 = fw2/tw2
fw3 = fw3/tw3
fw4 = fw4/tw4
mw1 = 22.5 - 2.5 * alog10(fw1)
mw2 = 22.5 - 2.5 * alog10(fw2)
mw3 = 22.5 - 2.5 * alog10(fw3)
mw4 = 22.5 - 2.5 * alog10(fw4)
fw1 = fw1 * f0
fw2 = fw2 * f0
fw3 = fw3 * f0
fw4 = fw4 * f0
few1 = few1/tw1 * f0
few2 = few2/tw2 * f0
few3 = few3/tw3 * f0
few4 = few4/tw4 * f0
ind = where(ivar2_w1 le 0, nt)
if (nt gt 0) then begin
  few1_ap[ind] = -1.0
  few1[ind] = -1.0
  mw1[ind] = -1.0
  mew1[ind] = -1.0
  mw1_ap[ind] = -1.0
endif
ind = where(ivar2_w2 le 0, nt)
if (nt gt 0) then begin
  few2_ap[ind] = -1.0
  few2[ind] = -1.0
  mw2[ind] = -1.0
  mew2[ind] = -1.0
  mw2_ap[ind] = -1.0
endif
ind = where(ivar2_w3 le 0, nt)
if (nt gt 0) then begin
  few3_ap[ind] = -1.0
  few3[ind] = -1.0
  mw3[ind] = -1.0
  mew3[ind] = -1.0
  mw3_ap[ind] = -1.0
endif
ind = where(ivar2_w4 le 0, nt)
if (nt gt 0) then begin
  few4_ap[ind] = -1.0
  few4[ind] = -1.0
  mw4[ind] = -1.0
  mew4[ind] = -1.0
  mw4_ap[ind] = -1.0
endif
;
ind = where(fw1 le 0, nt)
if (nt gt 0) then begin
  mw1[ind] = -1.0
  mew1[ind] = -1.0
  mw1_ap[ind] = -1.0
endif
ind = where(fw2 le 0, nt)
if (nt gt 0) then begin
  mw2[ind] = -1.0
  mew2[ind] = -1.0
  mw2_ap[ind] = -1.0
endif
ind = where(fw3 le 0, nt)
if (nt gt 0) then begin
  mw3[ind] = -1.0
  mew3[ind] = -1.0
  mw3_ap[ind] = -1.0
endif
ind = where(fw4 le 0, nt)
if (nt gt 0) then begin
  mw4[ind] = -1.0
  mew4[ind] = -1.0
  mw4_ap[ind] = -1.0
endif

nall = n_elements(ra)

result = replicate({ra:0.d, dec:0.d, shape:0.0, $
  ng:0, nr:0, nz:0, nw1:0, nw2:0, nw3:0, nw4:0, $
  BRICK_PRIMARY:'', depth_g:0.0, depth_r:0.0, depth_z:0.0, $
  mg_ap:0.0, mr_ap:0.0, mz_ap:0.0, $
  mw1_ap:0.0, mw2_ap:0.0, mw3_ap:0.0, mw4_ap:0.0, $
  mg:0.0, mr:0.0, mz:0.0, mw1:0.0, mw2:0.0, mw3:0.0, mw4:0.0, $
  meg:0.0, mer:0.0, mez:0.0, mew1:0.0, mew2:0.0, mew3:0.0, mew4:0.0, $
  fg_ap:0.d, fr_ap:0.d, fz_ap:0.d, feg_ap:0.d, fer_ap:0.d, fez_ap:0.d, $
  fw1_ap:0.d, fw2_ap:0.d, fw3_ap:0.d, fw4_ap:0.d, $
  few1_ap:0.d, few2_ap:0.d, few3_ap:0.d, few4_ap:0.d, $
  fg:0.d, fr:0.d, fz:0.d, feg:0.d, fer:0.d, fez:0.d, $
  fw1:0.d, fw2:0.d, fw3:0.d, fw4:0.d, $
  few1:0.d, few2:0.d, few3:0.d, few4:0.d $
}, nall)

result.ra = ra
result.dec = dec
result.shape = transpose(shape)
result.ng = transpose(ng)
result.nr = transpose(nr)
result.nz = transpose(nz)
result.nw1 = transpose(nw1)
result.nw2 = transpose(nw2)
result.nw3 = transpose(nw3)
result.nw4 = transpose(nw4)
result.BRICK_PRIMARY = BRICK_PRIMARY
result.depth_g = transpose(depth_g)
result.depth_r = transpose(depth_r)
result.depth_z = transpose(depth_z)
result.mg_ap = transpose(mg_ap)
result.mr_ap = transpose(mr_ap)
result.mz_ap = transpose(mz_ap)
result.mw1_ap = transpose(mw1_ap)
result.mw2_ap = transpose(mw2_ap)
result.mw3_ap = transpose(mw3_ap)
result.mw4_ap = transpose(mw4_ap)
result.mg = transpose(mg)
result.mr = transpose(mr)
result.mz = transpose(mz)
result.mw1 = transpose(mw1)
result.mw2 = transpose(mw2)
result.mw3 = transpose(mw3)
result.mw4 = transpose(mw4)
result.meg = transpose(meg)
result.mer = transpose(mer)
result.mez = transpose(mez)
result.mew1 = transpose(mew1)
result.mew2 = transpose(mew2)
result.mew3 = transpose(mew3)
result.mew4 = transpose(mew4)

result.fg_ap = transpose(fg_ap)
result.fr_ap = transpose(fr_ap)
result.fz_ap = transpose(fz_ap)
result.feg_ap = transpose(feg_ap)
result.fer_ap = transpose(fer_ap)
result.fez_ap = transpose(fez_ap)
result.fw1_ap = transpose(fw1_ap)
result.fw2_ap = transpose(fw2_ap)
result.fw3_ap = transpose(fw3_ap)
result.fw4_ap = transpose(fw4_ap)
result.few1_ap = transpose(few1_ap)
result.few2_ap = transpose(few2_ap)
result.few3_ap = transpose(few3_ap)
result.few4_ap = transpose(few4_ap)

result.fg = transpose(fg)
result.fr = transpose(fr)
result.fz = transpose(fz)
result.feg = transpose(feg)
result.fer = transpose(fer)
result.fez = transpose(fez)
result.fw1 = transpose(fw1)
result.fw2 = transpose(fw2)
result.fw3 = transpose(fw3)
result.fw4 = transpose(fw4)
result.few1 = transpose(few1)
result.few2 = transpose(few2)
result.few3 = transpose(few3)
result.few4 = transpose(few4)
mwrfits, result, 'flux_tractor-0001m002.fits'
end
