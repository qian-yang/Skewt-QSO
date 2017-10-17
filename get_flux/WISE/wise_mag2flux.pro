;+
; Name: wise_mag2flux
; Purpose: Convert to AB magnitude, and calculate flux for WISE photometry
; Input: w1, w3, w3, w4, w1e, w2e, w3e, w4e
; Output: fw1, fw2, fw3, fw4, few1, few2, few3, few4
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

pro wise_mag2flux, w1, w2, w3, w4, w1e, w2e, w3e, w4e, fw1 = fw1, fw2 = fw2, fw3 = fw3, fw4 = fw4, few1 = few1, few2 = few2, few3 = few3, few4 = few4
c = 2.99792458*1e18 ;A/s
w1 = w1 + 2.699
w2 = w2 + 3.339
w3 = w3 + 5.174
w4 = w4 + 6.620
;
fw1 = mag2flux(w1)
fw2 = mag2flux(w2)
fw3 = mag2flux(w3)
fw4 = mag2flux(w4)
few1 = mag2fluxerr(w1, w1e, fw1)
few2 = mag2fluxerr(w2, w2e, fw2)
few3 = mag2fluxerr(w3, w3e, fw3)
few4 = mag2fluxerr(w4, w4e, fw4)
end
