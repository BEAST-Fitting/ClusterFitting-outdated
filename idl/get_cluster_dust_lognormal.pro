;+
; NAME:
;       GET_CLUSTER_DUST_LOGNORMAL
;
; PURPOSE:
;       A log-normal dust distribution is computed.  Nominal usage is
;       for the Hierachical Bayesian cluster fitting for PHAT.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       GET_CLUSTER_DUST_LOGNORMAL, x, mu, sig
;
; INPUTS:
;       x: values to compute the log-normal distribution on.
;       mu: mean value of the log-normal
;       sig: sigma or dispersion of the log-normal
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       returns the log-normal distribution sampled on the x values
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Started     : Morgan/Dan (NYU, Camp Hogg, Mar 2012)
;       16 Jul 2012 : modified name (Karl Gordon, MPIA)
;       23 Jul 2012 : modified again (KDG)
;-

function logNormal, x, mu, sig
;+
; Returns the log Normal distribution on x space
;
; INPUTS:
;     x   -  where to compute the distribution
;     mu  -  average (linear-scale)
;     sig -  dispersion (linear scale))
;
; OUTPUTS:
;     f(r) - logNormal distribution
;-

indxs = where(x GT 0,n_indxs)
vals = x*0.
vals[indxs] = (1.D/(x[indxs] * sqrt(2*!PI) * sig )) * exp( - 0.5 * (alog(x[indxs]/mu)/sig)^2 )

return, vals
end ; function logNormal

function get_cluster_dust_lognormal, xAv, mu, sig, AvMW = AvMW
;+
; Compute the likelihood of the dust
; We consider a lognrmal distribution of the dust and we account for the Milky
; Way foreground dust as well.
;
; INPUTS:
;    xAv    - values of Av at wich the likelihood need to be computed
;    mu     - average value expected
;    sig    - expeted variance
; KEYWORDS:
;    AvMW   - Milky Way amplitude of dust (default = 0.0 )
; OUTPUTS:
;  f(xAv,..)- likelihood values over xAv
;-
  if (n_elements(AvMW) EQ 0) then AvMW = 0.
  return, AvMW + logNormal( xAv, mu, sig )
end
