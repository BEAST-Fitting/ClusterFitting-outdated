import numpy as np
import pylab as plt
import pyfits
import pdb
import sys
import emcee
import idlsave
import time
import itertools as it


def lognorm(x,mu,sig):
    #lognormal distribution
    vals = np.zeros_like(x)
    if (mu > 0.):
        vals = (1/(x*np.sqrt(2*np.pi)*sig)) * np.exp(-0.5 * (np.log(x/mu)/sig)**2)
    elif (mu == 0):
        #special case of Av = 0.
        vals[0,:] = 1. 
    #normalize in loop -- is loop necessary?
    for i in range(len(vals[0,:])):
        vals[:,i] /= vals[:,i].sum()
    return vals

def prob_dust(xAv, mu, sig, AvMW=0.0):
    #lognormal dust with no MW foreground
    return AvMW + lognorm(xAv, mu, sig)

def get_completeness():
    #completenes set to 1 everywhere for now
    return 1.

def stellar_grid():
    #read in stellar parameters
    star_pars = idlsave.read('fit_sed_grid_seds_av025_rv05.sav')
    #get logg, teff grid
    grav_vals = star_pars.grid_seds['logg']
    temp_vals = star_pars.grid_seds['logt'] 
    #find non-zero elements (i.e., have models)
    gtindxs = np.where(temp_vals > 0.)
    temp_vals = temp_vals[gtindxs]
    n_gravs = len(grav_vals)
    ggindx = np.where(grav_vals[1:n_gravs-1] > 0.)
    grav_vals = [grav_vals[0], (grav_vals[1:n_gravs-1])[ggindxs]]

    #get_isochrones,logl,logt,logg,radius,mass,bmass,gtag,age,no_inter=0
    

def star_prob(alpha, age):
    #look up closest alpha, age points on grid
    walpha = np.abs((alpha_grid-alpha)).argmin()
    wage = np.abs((age_grid-age)).argmin()
    return full_stellar_prob[:,:,walpha, wage], walpha, wage

#read in list of stellar .fits files
dirname = 'mock4_cl/'
filename = 'temp'
filelist = np.loadtxt(dirname+filename, dtype="string")
nfiles=len(filelist) 
#nfiles=1

#get the set of stellar probabilities for this set of stellar cluster parameters 
full_stellar_prob = np.transpose(pyfits.getdata('cluster_test_stellar_prob_av025_rv05.fits', extname='PRIMARY'))

#create grid of cluster points that match Karl's outputs
age_grid = 10**np.arange(6., 8.01, 0.25) # log age from 6 to 8 with 0.25 steps
alpha_grid = np.arange(0.5, 3.1, 0.25) # alpha from 0.5 to 3 with 0.25 steps
av_grid = np.arange(0.0, 3.1, 0.25) # av range from 0 to 3 with 0.25 steps
av_sig_grid = np.arange(0.1, 0.51, 0.1) # av sigrange 0.1 to 0.5 with 0.1 steps
rv_grid = np.arange(2.5, 5.1, 0.5) #rv range 2.5 to 5.0 with steps of 0.5
theta_prob = np.zeros((len(alpha_grid), len(age_grid), len(av_grid), len(av_sig_grid), len(rv_grid)))

#set field contribution to zero for now
field_eps = 0.0

#initialize likelihood
theta_prob = 0.
tprob_full = 0.

#reminder to un-hardcode expected grid sizes
star_av_vals = np.zeros((41, nfiles))
star_rv_vals = np.zeros((10, nfiles))
star_bmass_vals = np.zeros((50, nfiles))
star_bmass_prob = np.zeros((50, nfiles))
fullprob = np.zeros((76, 51, 41, 10, nfiles))
p_gamma_theta_cluster = np.zeros((76, 51, 41, 10, nfiles))
p_gamma_theta_field = np.ones((76, 51, 41, 10, nfiles))

#apply completeness and noramlize field population
p_gamma_theta_field *= get_completeness()
p_gamma_theta_field /= np.sum(p_gamma_theta_field)

#vector to store most probable stellar masses from Karl's fits
masses = np.zeros(nfiles)




#read in each star's .fits.gz ND probability
for c in range(nfiles):
    print 'Reading in Star Number ', c, 'of ', nfiles
    
    f = pyfits.open(dirname+filelist[c])

    star_av_vals[:,c] = np.array(f["AV_PROB"].data, dtype=float)[0,:] + 0.05
    star_rv_vals[:,c] = np.array(f["RV_PROB"].data, dtype=float)[0,:]

    fullprob[:,:,:,:,c] = np.array(f["FULL_PROB"].data, dtype=float).T
    fullprob[:,:,:,:,c] /= fullprob[:,:,:,:,c].sum()

    star_bmass_vals[:,c], star_bmass_prob[:,c] = np.array(f["BMASS_PROB"].data, dtype=float)
    f.close()
    
    masses[c] = star_bmass_vals[:,c][star_bmass_prob[:,c] == star_bmass_prob[:,c].max()]


def lnprob(theta, filelist=filelist, theta_prob=theta_prob, p_gamma_theta_cluster=p_gamma_theta_cluster, p_gamma_theta_field=p_gamma_theta_field, tprob_full=tprob_full):
    theta = [0.5, 1.0e6, 0.0, 0.1]
    #set priors on acceptable ranges for alpha, age, av, av_sig, etc
    priorcheck = [(theta[0] >= alpha_grid.min()), (theta[0] <= alpha_grid.max()), (theta[1] >= age_grid.min()), (theta[1] <= age_grid.max()), (theta[2] >= av_grid.min()), (theta[2] <= av_grid.max()), (theta[3] >= av_sig_grid.min()), (theta[3] <= av_sig_grid.max())]
    
    #check to make sure prior conditions are met, if not set prob to -infinity
    if not (False in priorcheck):
        STRT = time.time()
        #look up stellar probability on grid for alpha, age from emcee
        p_stellar, walpha, wage = star_prob(theta[0], theta[1])
        #get dust probability
        p_av = prob_dust(star_av_vals, theta[2], theta[3])
        #combine them to get PDF for group of stars -- seems like a less ineffiicent way to do this

        # new_p_gamma_theta_cluster = p_stellar[:,:,None,None] * p_av[None,None,:,:]  # USE ME!!!
        p_gamma_theta_cluster = p_stellar[:,:,None,None] * p_av[None,None,:,:]
        #apply completeness function
        p_gamma_theta_cluster *= get_completeness() 

        #normalize combined P(gamma|theta) -- correct? efficient?
        s = time.time()
        #for c in range(nfiles):
        #    p_gamma_theta_cluster[:,:,:,c] /= np.sum(p_gamma_theta_cluster[:,:,:,c])

        p_gamma_theta_cluster /= np.sum(p_gamma_theta_cluster)

        #Ignoring field populations for the time being -- seems to be a really slow part of the process otherwise
        p_gamma_theta =  p_gamma_theta_cluster #(1. - field_eps)*p_gamma_theta_cluster + field_eps*p_gamma_theta_field
        
        '''
        pdb.set_trace()
        #compute log probability for all stars in cluster -- correct? efficient?
        tprob=0
        s = time.time()
        print "tprob_full1", time.time() - s
        s = time.time()
        #for i in range(nfiles):
        #    tprob_full = fullprob[:,:,:,:,i]*p_gamma_theta
        #    tprob += np.log(np.sum(tprob_full))
        #pdb.set_trace()
        #tprob_full = np.log(np.sum(tprob_full))
        '''

        #This sum does what it's supposed to do...BUT...
        #Current difference with Karl's code is the 'catch' feature...
        tprob_list = np.sum(np.sum(np.sum(np.sum(fullprob[:,:,:,:,:]*p_gamma_theta[:,:,:,None,:], axis=0), axis=0), axis=0),axis=0)

        # the following line is an unjustified HACK in Karl's code
        # When the field model exists, this line should be disabled.
        tprob_list = np.clip(tprob_list, 1e-20, np.infty)

        # what's called `tprob_full` here is called `alog(tprob)` in Karl's code
        tprob_full = np.sum(np.log(tprob_list))
        print 'tprob_list', tprob_list.shape, tprob_list.min(), tprob_list.max(), 'tprob_full', tprob_full.shape

        #trying to replicate what Karl's code was doing
        #He added a catch to exclude(?) really bad fits

        #tprob = np.sum(tprob_full)
            #pdb.set_trace()
            #print c, theta, tprob, p_gamma_theta_cluster.sum()
        #theta_prob += np.log(tprob)
            #print theta, tprob, theta_prob
        #if ((field_eps <= 0.) & (tprob < 1e-50)): 
        #    theta_prob_full = -np.infty
        #else:        
        #    theta_prob += np.log(tprob)
        #pdb.set_trace()
    # if priors not True, then tprob_full should be -infinity
    else:
        tprob_full = -np.infty
    
    print theta[0], theta[1], theta[2], theta[3], tprob_full
    #pdb.set_trace()
    return tprob_full


#for c in range(nfiles):
#        tprob_full = fullprob[:,:,:,:]*p_gamma_theta
#        tprob = np.sum(tprob_full)
#        if ((field_eps <= 0.) & (tprob < 1e-50)): tprob = 1e-50
#        theta_prob += np.log(tprob)
#    print theta, theta_prob
#    return theta_prob



#initialization for emcee
#small number of steps, threads for laptop usage

nwalkers = 8
ndim = 4
nburn = 50
nsteps = 50
nthreads = 1


initial = [np.array([np.random.uniform(alpha_grid.min(), alpha_grid.max()), np.random.uniform(age_grid.min(), age_grid.max()), np.random.uniform(av_grid.min(), av_grid.max()), np.random.uniform(av_sig_grid.min(), av_sig_grid.max())]) for i in xrange(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=nthreads)
start = time.time()
pos,prob,state = sampler.run_mcmc(initial, nburn)
sampler.reset()
sampler.run_mcmc(np.array(pos),nsteps, rstate0=state)
duration = time.time()-start
print 'Sampling done in', duration, 's'




                                     
                                     
                                     
                                   
                     
                                     
                    
