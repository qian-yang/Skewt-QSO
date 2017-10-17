# Convert PS1 magnitude to flux
* get EBV values using Python code get_ps1_ext.py, which depends on code query_argonaut.py from https://github.com/gregreen/map3d.

* convert PS1 magnitude to flux using IDL code get_ps1_flux.pro, which depends on ps1_mag2flux.pro.

* Magnitude and magnitude errors from PS1 casjob (http://mastweb.stsci.edu/ps1casjobs).
* Table: StackObjectThin in PanSTARRS_DR1
* Magnitude: GPSFMAG, RPSFMAG, IPSFMAG, ZPSFMAG, YPSFMAG
* Magnitude error: GPSFMAGERR, RPSFMAGERR, IPSFMAGERR, ZPSFMAGERR, YPSFMAGERR
