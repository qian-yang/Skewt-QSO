# Skewt-QSO
## Quasar Photometric Redshifts and Candidate Selection
* The photometric redshift (photo-z) estimation code: photoz
* The candidate selection code: classification

* The photo-z and classification code is written in R language (https://cran.r-project.org), and depends on R packages, including sn, mvtnorm, pracma, data.table and bit64.
* To install a R package, for example sn, please run `install.packages("sn")` in R.
* To run a .R file, please run `source("filename.R")` in R.

## Flux
* The flux in the input files are calculated via IDL codes in the get_flux directory.

## Photo-z results
* The SDSS, SDSS-WISE, PS1-WISE, and DECaLS-WISE photo-z results are in files in the photoz/result directory.
* SDSS: photoz_SDSS.csv.zip
* SDSS-WISE: photoz_skewt_SDSS_WISE.csv.zip
* PS1-WISE: photoz_skewt_PS1_WISE.csv.zip
* DECaLS-WISE: photoz_skewt_DECaLS_WISE.csv.zip

## Photometric quasars in the SDSS Stripe 82 region
* In region 340< R.A.<45 deg and −1.25<Decl.<1.25 deg (roughly b < −50 deg).
Photometry used: DECaLS (DR3) and unWISE photometry.
File: catalog/photoz_qso_S82.fits

* Any question please fell free to contact Qian (qianyang.astro@gmail.com).
