;+
; Name: sdss_mag2flux
; Purpose: Correct galactic extionction, convert to AB magnitude, and calculate flux for SDSS photometry
; Input: u, g, r, i, z, eu, eg, er, ei, ez, au
; Output: fu, fg, fr, fi, fz, feg, fer, fei, fez
; Version:
;  July 11, 2017
; Author: Qian Yang, qianyang.astro@gmail.com
;-
;------------------------------------------------------------

function mag2flux, mag
	;erg*s-1*cm-2*Hz-1
	flux = (10.0d^(mag/(-2.5d))) * (3.631*(1e-20))
	nv = -1.0
	index = where(mag le 0, n)
	if (n gt 0) then flux[index] = nv
	return, flux
end

function mag2fluxerr, mag, magerr, flux
	fluxerr = flux * 0.4 * alog(10.d) * magerr
	nv = -1.0
	index = where(magerr le 0, n)
	if (n gt 0) then fluxerr[index] = nv
	return, fluxerr
end

pro sdss_mag2flux, u, g, r, i, z, eu, eg, er, ei, ez, au, fu = fu, fg = fg, fr = fr, fi = fi, fz = fz, feu = feu, feg = feg, fer = fer, fei = fei, fez = fez
c = 2.99792458*1e18 ;A/s
exr = [5.155, 3.793, 2.751, 2.086, 1.479]
exr = exr/exr[0]
u = u - exr[0] * au  - 0.04
g = g - exr[1] * au
r = r - exr[2] * au
i = i - exr[3] * au
z = z - exr[4] * au + 0.02
;
fu = mag2flux(u)
fg = mag2flux(g)
fr = mag2flux(r)
fi = mag2flux(i)
fz = mag2flux(z)
feu = mag2fluxerr(u, eu, fu)
feg = mag2fluxerr(g, eg, fg)
fer = mag2fluxerr(r, er, fr)
fei = mag2fluxerr(i, ei, fi)
fez = mag2fluxerr(z, ez, fz)
end
