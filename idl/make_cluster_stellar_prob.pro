;+
; NAME:
;       MAKE_CLUSTER_STELLAR_PROB
;
; PURPOSE:
;       Create the grid of stellar probabilities for a set of cluster
;       parameters; alpha (IMF slope) and age.  This creates an
;       output FITS file that is used by fit_cluster_prob.pro.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       MAKE_CLUSTER_STELLAR_PROB
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       run_tag : string to use to identify this run in output filenames
;       cluster_alpha_range : range of cluster alpha [IMF slope]
;       cluster_alpha_delta : delta of cluster alpha [IMF slope]
;       cluster_age_range : range of cluster ages [in log10 Myears units]
;       cluster_age_delta : delta of cluster ages [in log10 Myears units]
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
; 	Started     : Karl Gordon (23 Jul 2012)
;-

pro make_cluster_stellar_prob,smodname=smodname,run_tag=run_tag, $
                     cluster_age_range=cluster_age_range,cluster_age_delta=cluster_age_delta, $
                     cluster_alpha_range=cluster_alpha_range,cluster_alpha_delta=cluster_alpha_delta

if (not keyword_set(run_tag)) then run_tag = 'make_cluster'

; age in units of log10 Myears
if (not keyword_set(cluster_age_range)) then cluster_age_range = alog10([1,100.0])
if (not keyword_set(cluster_age_delta)) then cluster_age_delta = 0.25

if (not keyword_set(cluster_alpha_range)) then cluster_alpha_range = [0.5,3.0]
if (not keyword_set(cluster_alpha_delta)) then cluster_alpha_delta = 0.25

; create the cluster grid (cluster model parameters)
cluster_age_npts = fix((cluster_age_range[1] - cluster_age_range[0])/cluster_age_delta) + 1
cluster_age_vals = cluster_age_range[0] + findgen(cluster_age_npts)*cluster_age_delta
cluster_age_vals = 10^cluster_age_vals

cluster_alpha_npts = fix((cluster_alpha_range[1] - cluster_alpha_range[0])/cluster_alpha_delta) + 1
cluster_alpha_vals = cluster_alpha_range[0] + findgen(cluster_alpha_npts)*cluster_alpha_delta

; get the grid of models
if (not keyword_set(smodname)) then smodname = 'av025_rv05'
restore,'fit_sed_grid_seds_'+smodname+'.sav'

; get the temp and grav values
temp_vals = max(grid_seds[*,*].logt,dimension=2)
grav_vals = max(grid_seds[*,*].logg,dimension=1)

; get the points in the temp and grav vectors that have models
gtindxs = where(temp_vals GT 0.)
temp_vals = temp_vals[gtindxs]
n_gravs = n_elements(grav_vals)
ggindxs = where(grav_vals[1:n_gravs-1] GT 0.)
grav_vals = [grav_vals[0],(grav_vals[1:n_gravs-1])[ggindxs]]

; setup the stellar probability grid
;   2D for each cluster alpha and age
size_grid = size(grid_seds[*,*])
stellar_prob = fltarr(size_grid[1],size_grid[2],cluster_alpha_npts,cluster_age_npts)
; temp 2D version for ease of creation
tstellar_prob = fltarr(size_grid[1],size_grid[2])

; get the stellar evolutionary tracks
get_stellar_tracks,logl,logt,logg,radius,mass,bmass,gtag,age,no_inter=0

size_st_grid = size(mass)

; age in years from Myrs
tcluster_age = cluster_age_vals*1e6

; loop over cluster parameters
for x = 0,(cluster_alpha_npts-1) do begin
    for y = 0,(cluster_age_npts-1) do begin

        print,'cluster alpha, age = ', cluster_alpha_vals[x], cluster_age_vals[y]

        ; make sure the temp version starts fully initialized to zero
        tstellar_prob[*,*] = 0.0

        ; loop over masses
        for j = 1,(size_st_grid[2]-2) do begin
            indxs = where(finite(age[*,j]),n_indxs)
            if (min(age[indxs,j]) GT tcluster_age[y]) then begin
                cur_logt = logt[indxs[0],j]
                cur_logg = logg[indxs[0],j]
            endif else if (max(age[indxs,j]) LT tcluster_age[y]) then begin
                cur_logt = !values.f_nan
                cur_logg = !values.f_nan
            endif else begin
                cur_logt = interpol(logt[indxs,j],age[indxs,j],tcluster_age[y])
                cur_logg = interpol(logg[indxs,j],age[indxs,j],tcluster_age[y])
            endelse

            if (finite(cur_logt)) then begin

                dist = (grid_seds[*].logt - cur_logt)^2 + (grid_seds[*].logg - cur_logg)^2
                sindxs = sort(dist)
        
                dbmass = 0.5*(bmass[0,j-1] - bmass[0,j+1])
                tstellar_prob[sindxs[0]] += dbmass*(bmass[0,j]^(-1.0*cluster_alpha_vals[x]))

            endif
        endfor

        ; normalize the tstellar_prob to have a sum/area of 1
        tot_stellar_prob = total(tstellar_prob)
        tstellar_prob /= tot_stellar_prob

        ; now store for the full output
        stellar_prob[*,*,x,y] = tstellar_prob

    endfor
endfor

fxhmake,header,stellar_prob,/initialize
sxaddpar,header,'ALPHMIN',cluster_alpha_range[0],' Min cluster alpha'
sxaddpar,header,'ALPHMAX',cluster_alpha_range[1],' Max cluster alpha'
sxaddpar,header,'ALPHDELT',cluster_alpha_delta,' delta cluster alpha'
sxaddpar,header,'AGEMIN',cluster_age_range[0],' Min cluster age'
sxaddpar,header,'AGEMAX',cluster_age_range[1],' Max cluster age'
sxaddpar,header,'AGEDELT',cluster_age_delta,' delta cluster age'
fits_write,run_tag+'_stellar_prob_'+smodname+'.fits',stellar_prob,header

end
