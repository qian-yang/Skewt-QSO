;+
; Name: get_ps1_flux
; Purpose: Convert Pan-STARRS(PS1) magnitude to flux
; Input: input fits file including GPSFMAG, RPSFMAG, IPSFMAG, ZPSFMAG, YPSFMAG, GPSFMAGERR, RPSFMAGERR, IPSFMAGERR, ZPSFMAGERR, YPSFMAGERR, EBV
; Output: output fits file including fg, fr, fi, fz, fy, feu, feg, fer, fei, fez, fey
; Version:
;  July 11, 2017
; Author: Qian Yang, qianyang.astro@gmail.com
;-
;------------------------------------------------------------

pro get_ps1_flux
file = 'ps1_ebv.fits'
data = mrdfits(file, 1)
g = data.GPSFMAG
r = data.RPSFMAG
i = data.IPSFMAG
z = data.ZPSFMAG
y = data.YPSFMAG
eg = data.GPSFMAGERR
er = data.RPSFMAGERR
ei = data.IPSFMAGERR
ez = data.ZPSFMAGERR
ey = data.YPSFMAGERR
ebv = data.EBV
ps1_mag2flux, g, r, i, z, y, eg, er, ei, ez, ey, ebv, fg = fg, fr = fr, fi = fi, fz = fz, fy = fy, feg = feg, fer = fer, fei = fei, fez = fez, fey = fey
struct_add_field, data, 'fg', fg
struct_add_field, data, 'fr', fr
struct_add_field, data, 'fi', fi
struct_add_field, data, 'fz', fz
struct_add_field, data, 'fy', fy
struct_add_field, data, 'feg', feg
struct_add_field, data, 'fer', fer
struct_add_field, data, 'fei', fei
struct_add_field, data, 'fez', fez
struct_add_field, data, 'fey', fey
mwrfits, data, 'ps1_ebv_flux_ap.fits'
end
