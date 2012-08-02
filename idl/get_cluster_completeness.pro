;+
; NAME:
;       GET_CLUSTER_COMPLETENESS
;
; PURPOSE:
;       Get the completeness function in the stellar model space
;       given the mag completeness functions.  
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       A = GET_CLUSTER_COMPLETENESS()
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       smodname: string that idenifies the stellar model IDL save file
;                 [e.g., smodname='av025_rv05']
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Started     : Karl Gordon (MPIA, 18 Mar 2012)
;       23 Jul 2012 : Cleaned and added full documentation added (KDG)
;-

function get_cluster_completeness,smodname=smodname

if (not keyword_set(smodname)) then smodname = 'av025_rv05'
restore,'fit_sed_band_seds_'+smodname+'.sav'

; currently assuming a f814w step function cutoff
size_cube = size(band_seds.band_grid_seds)

comp_func = replicate(0.0,size_cube[2],size_cube[3],size_cube[4],size_cube[5])

f814w_vals = reform(band_seds.band_grid_seds[3,*,*,*,*],size_cube[2],size_cube[3],size_cube[4],size_cube[5])

; convert to mags
es_points = [3.2596e-18,1.34e-18,1.820979E-19,7.033186E-20,1.5233e-20,1.9106e-20]
vegamag = [22.65,23.46,26.16,25.52,26.07,24.70]
waves = [275.,336.,475.,814.,1100.,1600.]*1e-3
zero_points = es_points/10^(-0.4*vegamag)

f814w_mags = f814w_vals*0.0 + 100.
indxs = where(f814w_vals GT 0.0,n_indxs)
f814w_mags[indxs] = -2.5*alog10(f814w_vals[indxs]/zero_points[3])

; now add construct the completeness function
indxs = where(f814w_mags LT 24.0,n_indxs)
comp_func[indxs] = 1.0
; code for linear completeness function between two mags
;indxs = where((f814w_mags GE 23.0) AND (f814w_mags LT 24.0),n_indxs)
;comp_func[indxs] = -1.*(f814w_mags[indxs] - 24.)

;fits_write,'test_f814w_cfunc.fits',comp_func

return, comp_func

end
