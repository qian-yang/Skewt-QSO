;+
; Name: get_sdss_flux
; Purpose: Convert SDSS magnitude to flux
; Input: input fits file including psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z, psfmagerr_u, psfmagerr_g, psfmagerr_r, psfmagerr_i, psfmagerr_z, extinction_u
; Output: output fits file including fu, fg, fr, fi, fz, feu, feg, fer, fei, fez
; Version:
;  July 11, 2017
; Author: Qian Yang, qianyang.astro@gmail.com
;-
;------------------------------------------------------------

pro get_sdss_flux
file = 'sdss_flag_ch.fits'
data = mrdfits(file, 1)
u = data.psfmag_u
g = data.psfmag_g
r = data.psfmag_r
i = data.psfmag_i
z = data.psfmag_z
eu = data.psfmagerr_u
eg = data.psfmagerr_g
er = data.psfmagerr_r
ei = data.psfmagerr_i
ez = data.psfmagerr_z
au = data.extinction_u
sdss_mag2flux, u, g, r, i, z, eu, eg, er, ei, ez, au, fu = fu, fg = fg, fr = fr, fi = fi, fz = fz, feu = feu, feg = feg, fer = fer, fei = fei, fez = fez
struct_add_field, data, 'fu', fu
struct_add_field, data, 'fg', fg
struct_add_field, data, 'fr', fr
struct_add_field, data, 'fi', fi
struct_add_field, data, 'fz', fz
struct_add_field, data, 'feu', feu
struct_add_field, data, 'feg', feg
struct_add_field, data, 'fer', fer
struct_add_field, data, 'fei', fei
struct_add_field, data, 'fez', fez
mwrfits, data, 'sdss_flag_ch_flux.fits'
end
