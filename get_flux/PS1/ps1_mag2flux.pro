;+
; Name: ps1_mag2flux
; Purpose: Correct galactic extionction, and calculate flux for PS1 photometry
; Input: g, r, i, z, y, eu, eg, er, ei, ez, ey, ebv
; Output: fg, fr, fi, fz, fy, feg, fer, fei, fez, fey
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

pro ps1_mag2flux, g, r, i, z, y, eg, er, ei, ez, ey, ebv, fg = fg, fr = fr, fi = fi, fz = fz, fy = fy, feg = feg, fer = fer, fei = fei, fez = fez, fey = fey
c = 2.99792458*1e18 ;A/s
Arr = [3.172, 2.271, 1.682, 1.322, 1.087]
g = g - Arr[0] * ebv
r = r - Arr[1] * ebv
i = i - Arr[2] * ebv
z = z - Arr[3] * ebv
y = y - Arr[4] * ebv
;
fg = mag2flux(g)
fr = mag2flux(r)
fi = mag2flux(i)
fz = mag2flux(z)
fy = mag2flux(y)
feg = mag2fluxerr(g, eg, fg)
fer = mag2fluxerr(r, er, fr)
fei = mag2fluxerr(i, ei, fi)
fez = mag2fluxerr(z, ez, fz)
fey = mag2fluxerr(y, ey, fy)
end
