import numpy as np
import pylab as plt
import pyfits
import pdb
import sys


def lognorm(x,mu,sig):
    indxs = np.where(x > 0)
    vals = x*0.
    vals[indxs] = (1/(x[indxs]*np.sqrt(2*np.pi)*sig)) * np.exp(-0.5 * (np.log(x[indxs]/mu)/sig)**2)
    return vals

def prob_dust( xAv, mu, sig, AvMW=0.0):
    return AvMW + lognorm(xAv, mu, sig)

#def step_completeness():

nfiles=1 # change to correct number of files when I scale this up


#read in single star .fits file
hdulist = pyfits.open('mock4_randfield_jul12_b1_s10_av025_rv05.fits.gz')
scidata = hdulist[1].data
star_av_vals = hdulist['AV_PROB'].data[0,:]
star_rv_vals = hdulist['RV_PROB'].data[0,:]

#get the set of stellar probabilities for this set of stellar cluster parameters
hdulist1 = pyfits.open('../data/cluster_test_stellar_prob_av025_rv05.fits')
full_stellar_prob = np.transpose(hdulist1['PRIMARY'].data)


#assign and normalize full ND probability
fullprob = np.transpose(hdulist['FULL_PROB'].data)
fullprob /= np.sum(fullprob)

#only valid for sinle star -- redo for multiple stars
p_gamma_theta_cluster = fullprob*0.

#create grid of cluster points

age_grid = np.arange(6., 8.01, 0.25) # log age from 6 to 8 with 0.25 steps
alpha_grid = np.arange(0.5, 3.1, 0.25) # alpha from 0.5 to 3 with 0.25 steps
av_grid = np.arange(0.0, 3.1, 0.25) # av range from 0 to 3 with 0.25 steps
av_sig_grid = np.arange(0.1, 0.51, 0.1) # av sigrange 0.1 to 0.5 with 0.1 steps
rv_grid = np.arange(2.5, 5.1, 0.5) #rv range 2.5 to 5.0 with steps of 0.5





field_eps = 0.0 #epilson value for field stars

# FUNCTION TO GET COMPLETENESS
# completeness_function = get_cluster_completeness()


# PDF for field, set flat by default
p_gamma_theta_field = p_gamma_theta_cluster*0.0 + 1.0
p_gamma_theta_field *= 1. # need to add in compneteness here
p_gamma_theta_field /= np.sum(p_gamma_theta_field)


#loop over theta values of interest and generate ND PDF for cluster parameters

# XXXX make inline for loops once it runs XXXX


#loop over alpha


theta_prob = np.zeros((len(alpha_grid), len(age_grid), len(av_grid), len(av_sig_grid), len(rv_grid)))

for i in range(len(alpha_grid)):
    for j in range(len(age_grid)):
        # get stellar portion of p(gamma_theta), 2D in log(teff), log(g)
        p_stellar = full_stellar_prob[:,:,i,j]
        for k in range(len(av_grid)):
            for l in range(len(av_sig_grid)):
                if (av_grid[k] > 0.):
                    p_av = prob_dust(star_av_vals, av_grid[k], av_sig_grid[l])
                else:
                    p_av = star_av_vals*0.
                    p_av[0] = 1.
                p_av /= np.sum(p_av)
                for m in range(len(rv_grid)):
                    #Does nothing for now

                    for x in range(len(fullprob[0,0,:,0])):
                        for y in range(len(fullprob[0,0,0,:])):
                            p_gamma_theta_cluster[:,:,x,y] = p_stellar*p_av[x]
                                  
               
                    completeness_function = 1.
                    field_eps = 0.
                    p_gamma_theta_cluster *= completeness_function
                    p_gamma_theta_cluster /= p_gamma_theta_cluster.sum()
                    p_gamma_theta = (1. - field_eps)*p_gamma_theta_cluster + field_eps*p_gamma_theta_field
                    # compute total probability for star 'c'
                    #for c in range(len(nfiles)):
                    tprob_full = fullprob[:,:,:,:]*p_gamma_theta
                    tprob = np.sum(tprob_full)
                    
                    if ((field_eps <= 0.) & (tprob < 1e-20)):
                        tprob = 1e-20
                    theta_prob[i,j,k,l, m] += np.log(tprob)
                    print alpha_grid[i], age_grid[j], av_grid[k], av_sig_grid[l], rv_grid[m], theta_prob[i,j,k,l,m]





                                     
                                     
                                     
                                   
                     
                                     
                    
