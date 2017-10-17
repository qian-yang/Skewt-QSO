;+
; Name: get_wise_flux
; Purpose: Convert WISE magnitude to flux
; Input: input fits file including w1mpro, w2mpro, w3mpro, w4mpro, w1sigmpro, w2sigmpro, w3sigmpro, w4sigmpro
; Output: output fits file including fw1, fw2, fw3, fw4, few1, few2, few3, few4
; Version:
;  July 11, 2017
; Author: Qian Yang, qianyang.astro@gmail.com
;-
;------------------------------------------------------------


function nanreplace,array,value
if ~keyword_set(value) then begin
value=-1.0
endif
index=WHERE(FINITE(array,/nan))
if (min(index) ge 0) then begin
array[index]=value
endif
return,array
end

pro get_wise_flux
file = 'allwise.fits'
data = mrdfits(file, 1)
w1 = nanreplace(data.w1mpro)
w2 = nanreplace(data.w2mpro)
w3 = nanreplace(data.w3mpro)
w4 = nanreplace(data.w4mpro)
w1e = nanreplace(data.w1sigmpro)
w2e = nanreplace(data.w2sigmpro)
w3e = nanreplace(data.w3sigmpro)
w4e = nanreplace(data.w4sigmpro)
wise_mag2flux, w1, w2, w3, w4, w1e, w2e, w3e, w4e, fw1 = fw1, fw2 = fw2, fw3 = fw3, fw4 = fw4, few1 = few1, few2 = few2, few3 = few3, few4 = few4
struct_add_field, data, 'fw1', fw1
struct_add_field, data, 'fw2', fw2
struct_add_field, data, 'fw3', fw3
struct_add_field, data, 'fw4', fw4
struct_add_field, data, 'few1', few1
struct_add_field, data, 'few2', few2
struct_add_field, data, 'few3', few3
struct_add_field, data, 'few4', few4
mwrfits, data, 'wise_flux.fits'
end
