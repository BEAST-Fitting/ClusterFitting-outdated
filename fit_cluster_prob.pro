;+
; NAME:
;       FIT_CLUSTER_PROB
;
; PURPOSE:
;       Fit a model of a stellar cluster to a set of cluster plus
;       nearby field stars using Hierarchical Bayesian model.  The
;       cluster model has a stellar (age, IMF slope), interstellar
;       dust (log-normal M31 dust), and membership
;       (based on integrated flux cluster fitting).
;
;       The star model+data information is contained in nD
;       likelyhood/probability distributions that are calculated from
;       the fit_seds.pro code (Bayesian stellar+insterstellar
;       fitting). 
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       FIT_CLUSTER_PROB
;
; INPUTS:
;       star_filebase : filebase for the individual star nD
;                       probability distributions (search string)
;
; KEYWORD PARAMETERS:
;       run_tag : string to use to identify this run in output filenames
;       field_eps : field population epsilon [default is 0]
;       silent : surpress all screen output (except error messages)
;       debug : settings for debugging (fewer files, etc.)
;
;     The following parameters all have defaults if not set.  Run code
;     w/o any parameters and they will be printed to the screen.
;       cluster_age_range : range of cluster ages [in log10 Myears units]
;       cluster_age_delta : delta of cluster ages [in log10 Myears units]
;       cluster_alpha_range : range of cluster alpha [IMF slope]
;       cluster_alpha_delta : delta of cluster alpha [IMF slope]
;       cluster_av_range : range of cluster A(V) assuming a log-normal
;                          distribution
;       cluster_av_delta : delta of cluster A(V)
;       cluster_av_sigma_range : range of width of cluster A(V)
;                                assuming a log-normal distribution
;       cluster_av_sigma_range : delta of width of cluster A(V)
;       cluster_rv_range : range of cluster R(V) assuming delta function
;       cluster_rv_delta : range of cluster R(V)
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
; 	Started     : Karl Gordon (NYU, Camp Hogg, 12 Mar 2012)
;       27 Mar 2012 : added start of documentation (KDG)
;       16 Jul 2012 : Fleshed out to the full 5D cluster space (MPIA, KDG)
;       23 Jul 2012 : Cleaned and added full documentation added (KDG)
;-

pro fit_cluster_prob,star_filebase,run_tag=run_tag,field_eps=field_eps,debug=debug, $
                     cluster_age_range=cluster_age_range,cluster_age_delta=cluster_age_delta, $
                     cluster_alpha_range=cluster_alpha_range,cluster_alpha_delta=cluster_alpha_delta, $
                     cluster_av_range=cluster_av_range,cluster_av_delta=cluster_av_delta, $
                     cluster_av_sigma_range=cluster_av_sigma_range,cluster_av_sigma_delta=cluster_av_sigma_delta, $
                     cluster_rv_range=cluster_rv_range,cluster_rv_delta=cluster_rv_delta, $
                     silent=silent
                     
; setup the range of cluster parameters (theta) to search over
; default are present if they are not set

; age in units of log10 Myears
if (not keyword_set(cluster_age_range)) then cluster_age_range = alog10([1,100.0])
if (not keyword_set(cluster_age_delta)) then cluster_age_delta = 0.25

if (not keyword_set(cluster_alpha_range)) then cluster_alpha_range = [0.5,3.0]
if (not keyword_set(cluster_alpha_delta)) then cluster_alpha_delta = 0.25

if (not keyword_set(cluster_av_range)) then cluster_av_range = [0.0,3.0]
if (not keyword_set(cluster_av_delta)) then cluster_av_delta = 0.25

if (not keyword_set(cluster_av_sigma_range)) then cluster_av_sigma_range = [0.1,0.5]
if (not keyword_set(cluster_av_sigma_delta)) then cluster_av_sigma_delta = 0.1

if (not keyword_set(cluster_rv_range)) then cluster_rv_range = [2.5,5.0]
if (not keyword_set(cluster_rv_delta)) then cluster_rv_delta = 0.5

; check we got at least 1 parameters
if (N_params() LT 1) then begin
    print,"Syntax - fit_cluster_prob,star_filebase,run_tag='',field_eps=0.0,/debug,/silent"
    print,'Other Keywords: 
    print,'          cluster_alpha_range, cluster_alpha_delta,'
    print,'          cluster_age_range, cluster_age_delta,'
    print,'          cluster_av_range, cluster_av_delta,'
    print,'          cluster_av_sigma_range, cluster_av_sigma_delta,'
    print,'          cluster_rv_range, cluster_rv_delta'
    print,'      range = 2 element vector, delta = scaler'
    print,''
    print,'default cluster parameters are:'
    print,'Alpha [IMF slope]: ' + strtrim(string(cluster_alpha_range[0],format="(F10.2)"),2) + $
          ' to ' + strtrim(string(cluster_alpha_range[1],format="(F10.2)"),2) + ' in steps of ' + $
          strtrim(string(cluster_alpha_delta,format="(F10.2)"),2)
    print,'Age [log10(Myrs)]: ' + strtrim(string(cluster_age_range[0],format="(F10.2)"),2) + $
          ' to ' + strtrim(string(cluster_age_range[1],format="(F10.2)"),2) + ' in steps of ' + $
          strtrim(string(cluster_age_delta,format="(F10.2)"),2)
    print,'       A(V) [mag]: ' + strtrim(string(cluster_av_range[0],format="(F10.2)"),2) + $
          ' to ' + strtrim(string(cluster_av_range[1],format="(F10.2)"),2) + ' in steps of ' + $
          strtrim(string(cluster_av_delta,format="(F10.2)"),2)
    print,' sigma A(V) [mag]: ' + strtrim(string(cluster_av_sigma_range[0],format="(F10.2)"),2) + $
          ' to ' + strtrim(string(cluster_av_sigma_range[1],format="(F10.2)"),2) + ' in steps of ' + $
          strtrim(string(cluster_av_sigma_delta,format="(F10.2)"),2)
    print,'       R(V) [mag]: ' + strtrim(string(cluster_rv_range[0],format="(F10.2)"),2) + $
          ' to ' + strtrim(string(cluster_rv_range[1],format="(F10.2)"),2) + ' in steps of ' + $
          strtrim(string(cluster_rv_delta,format="(F10.2)"),2)
    return
endif

if (not keyword_set(run_tag)) then run_tag = 'fit_cluster_prob'

; temporary? keyword for field_eps
if (not keyword_set(field_eps)) then field_eps = 0.0

; get the files with the nD probability distributions for each
; potential cluster star
files = file_search(star_filebase,count=n_files)

; check that we have some stellar files
if (n_files LE 0) then begin
    print,'No stellar files match the input filebase = ' + star_filebase
    return
endif

; debugging on fewer files
if (keyword_set(debug)) then n_files = min([50,n_files])

; go through the list of stars and read in the nD PDF for each one
if (not keyword_set(silent)) then print,'# of individual star files = ' + strtrim(string(n_files),2)

star_ra = dblarr(n_files)
star_dec = dblarr(n_files)
for i = 0,(n_files-1) do begin
    print,'reading in star file = ' + files[i]
    fits_open,files[i],ifcb
    fits_read,ifcb,sed,main_header,exten_no=1

    ; get (ra,dec) for cluster membership calculation
    star_ra[i] = fxpar(main_header,'RA')
    star_dec[i] = fxpar(main_header,'DEC')

    if (i EQ 0) then begin
        ; get the A(V) and R(V) vectors (used later)
        fits_read,ifcb,av_info,av_header,extname="AV_PROB"
        star_av_vals = av_info[*,0]
        fits_read,ifcb,rv_info,rv_header,extname="RV_PROB"
        star_rv_vals = rv_info[*,0]
    endif

    ; get the nd PDF
    fits_read,ifcb,oneprob,prob_header,extname="FULL_PROB"

    fits_close,ifcb
    
    if (i EQ 0) then begin
        ; initialize the full n_files x nD prob that
        ; will store all the stars nD probabilities
        size_oneprob = size(oneprob)
        fullprob = dblarr(size_oneprob[1],size_oneprob[2],size_oneprob[3],size_oneprob[4],n_files)
    endif

    ; save in the array (normalized)
    ;   - may want to revisit this normalization at a later date
    ;   - can we use the overall level of the likelihood function?
    fullprob[*,*,*,*,i] = oneprob/total(oneprob)
endfor

; initialize the model cluster PDF
p_gamma_theta_cluster = dblarr(size_oneprob[1],size_oneprob[2],size_oneprob[3],size_oneprob[4])

; read in the p(gamma|theta) for the field stars
;   created with create_p_gamma_theta_field.pro
;fits_read,cluster_p_gamma_theta_field,p_gamma_theta_field,pgtf_header
;  - right now it is a flat PDF
p_gamma_theta_field = p_gamma_theta_cluster*0.0 + 1.0


; get the membership probabilities for each star, 1D vector by star
;p_membership = pmem(star_ra, star_dec, cluster_rc, cluster_c, cluster_ra, cluster_dec)
;mean_ra = mean(star_ra)
;mean_dec = mean(star_dec)
;radius = sqrt((star_ra - mean_ra)^2 + (star_dec - mean_dec)^2)
;cluster_radius = max(radius)

;p_membership = radius/cluster_radius

; create the theta grid (cluster model parameters)

cluster_age_npts = fix((cluster_age_range[1] - cluster_age_range[0])/cluster_age_delta) + 1
cluster_age_vals = cluster_age_range[0] + findgen(cluster_age_npts)*cluster_age_delta
cluster_age_vals = 10^cluster_age_vals

cluster_alpha_npts = fix((cluster_alpha_range[1] - cluster_alpha_range[0])/cluster_alpha_delta) + 1
cluster_alpha_vals = cluster_alpha_range[0] + findgen(cluster_alpha_npts)*cluster_alpha_delta

cluster_av_npts = fix((cluster_av_range[1] - cluster_av_range[0])/cluster_av_delta) + 1
cluster_av_vals = cluster_av_range[0] + findgen(cluster_av_npts)*cluster_av_delta

cluster_av_sigma_npts = fix((cluster_av_sigma_range[1] - cluster_av_sigma_range[0])/cluster_av_sigma_delta) + 1
cluster_av_sigma_vals = cluster_av_sigma_range[0] + findgen(cluster_av_sigma_npts)*cluster_av_sigma_delta

cluster_rv_npts = fix((cluster_rv_range[1] - cluster_rv_range[0])/cluster_rv_delta) + 1
cluster_rv_vals = cluster_rv_range[0] + findgen(cluster_rv_npts)*cluster_rv_delta

theta_prob = replicate(0.d0,cluster_alpha_npts,cluster_age_npts, $
                       cluster_av_npts,cluster_av_sigma_npts,cluster_rv_npts)

; get the completeness function
; not a function of position at this point 
;    (probably will be in the future)
; this is currently in a highly developmental state
;  - will need to figure out how we are specificying this 
;    (currently hard codeds)
completeness_function = get_cluster_completeness()

; multiply by both the field by the completeness function *and* normalize to 1
p_gamma_theta_field *= completeness_function
p_gamma_theta_field /= total(p_gamma_theta_field)

; loop over the theta values of interest and generate the nD PDF for
; the cluster paramters

for i = 0,(cluster_alpha_npts-1) do begin
    for j = 0,(cluster_age_npts-1) do begin

        ; get the stellar portion of p(gamma|theta)
        ;   2D in log(teff) & log(g)
        p_stellar = get_cluster_stellar_prob(cluster_alpha_vals[i],cluster_age_vals[j],grid_seds=grid_seds)

        if (keyword_set(debug)) then begin
            fits_write,'test_stellar_prob_alpha_' + strtrim(string(cluster_alpha_vals[i],format='(F10.2)'),2) + $
                       '_age_' + strtrim(string(cluster_age_vals[j],format='(F10.2)'),2) + '.fits', $
                       p_stellar
        endif

        for k = 0,(cluster_av_npts-1) do begin
            for l = 0,(cluster_av_sigma_npts-1) do begin

                ; get the lognormal screen distribution
                ;   1D in single star A(V) space
                if (cluster_av_vals[k] GT 0) then begin
                    p_av = prob_dust_lognormal(star_av_vals,cluster_av_vals[k],cluster_av_sigma_vals[l])
                endif else begin ; special case if A(V) is set to 0
                    p_av = star_av_vals*0.
                    p_av[0] = 1.0
                endelse
                p_av /= total(p_av)

                for m = 0,(cluster_rv_npts-1) do begin
                    if (not keyword_set(silent)) then begin
                        print,'calculating theta [alpha, age, A(V), delta A(V), R(V)] = ', $
                              cluster_alpha_vals[i], cluster_age_vals[j], $
                              cluster_av_vals[k], cluster_av_sigma_vals[l], cluster_rv_vals[m]
                    endif

                    ; create the full 4D p(gamma|theta) for the cluster stars
                    ;  - currently does not include R(V)
                    for x = 0,(size_oneprob[3]-1) do begin
                        for y = 0,(size_oneprob[4]-1) do begin
                            p_gamma_theta_cluster[*,*,x,y] = p_stellar*p_av[x]
                        endfor
                    endfor
                                            
                    ; multiply by the completeness function 
                    ; *and* normalize to 1
                    p_gamma_theta_cluster *= completeness_function
                    p_gamma_theta_cluster /= total(p_gamma_theta_cluster)

                    ; may move inside following for loop
                    ; if completeness is position dependent
                    p_gamma_theta = (1.0 - field_eps)*p_gamma_theta_cluster + field_eps*p_gamma_theta_field

                    for c = 0,(n_files-1) do begin

                        ; compute the total probability for this star
                        tprob_full = fullprob[*,*,*,*,c]*p_gamma_theta
                        tprob = total(tprob_full)
                        
                        ; code to penalize non-fitting stars when there is no field in the model
                        if ((field_eps LE 0.0) AND (tprob LT 1e-20)) then tprob = 1e-20

                        theta_prob[i,j,k,l,m] += alog(tprob)

                    endfor

                    if (not keyword_set(silent)) then print,'theta_prob = ', theta_prob[i,j,k,l,m]

                endfor
            endfor
        endfor
    endfor
endfor

; find and renomalize by the maximum probability
max_tp = max(theta_prob)
if (not keyword_set(silent)) then print,'max theta_prob = ', max_tp
theta_prob -= max_tp

fxhmake,header,theta_prob,/initialize
sxaddpar,header,'MAX_PROB',max_tp,' Maximum probability used to normalize (subtract)'
fits_write,run_tag+'_theta_prob.fits',theta_prob,header

end
