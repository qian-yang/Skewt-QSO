#------------------------------------------------------------
# Purpose: Get EBV value
# Version: July 11, 2017
# Author: Qian Yang, qianyang.astro@gmail.com
#------------------------------------------------------------

def get_ebv(ra, dec):
    qresult = query(ra, dec, coordsys='equ', mode = 'sfd')
    ebv = qresult['EBV_SFD']
    return(ebv)

import numpy as np
import pyfits
import math
from query_argonaut import *
filename = 'ps1.fits'
orig_table = pyfits.open(filename)[1].data
ra = orig_table['RA']
dec= orig_table['DEC']
ra = ra.tolist()
dec= dec.tolist()
ebv = get_ebv(ra, dec)
#
orig_cols = orig_table.columns
new_cols = pyfits.ColDefs([pyfits.Column(name='EBV', format='D', array = ebv)])
hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto('ps1_ebv.fits')
