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
    indxs = np.where(x > 0)
    vals = x*0.
    vals[indxs] = (1/(x[indxs]*np.sqrt(2*np.pi)*sig)) * np.exp(-0.5 * (np.log(x[indxs]/mu)/sig)**2)
    return vals

def prob_dust(xAv, mu, sig, AvMW=0.0):
    return AvMW + lognorm(xAv, mu, sig)

def get_completeness():
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
    walpha = np.abs((alpha_grid-alpha)).argmin()
    wage = np.abs((age_grid-age)).argmin()
    return full_stellar_prob[:,:,walpha, wage], walpha, wage


filelist = np.loadtxt('mock4_cl/temp', dtype="string")
nfiles=len(filelist) # change to correct number of files when I scale this up
#nfiles=1

#get the set of stellar probabilities for this set of stellar cluster parameters
full_stellar_prob = np.transpose(pyfits.getdata('cluster_test_stellar_prob_av025_rv05.fits', extname='PRIMARY'))

#create grid of cluster points

age_grid = np.arange(6., 8.01, 0.25) # log age from 6 to 8 with 0.25 steps
alpha_grid = np.arange(0.5, 3.1, 0.25) # alpha from 0.5 to 3 with 0.25 steps
av_grid = np.arange(0.0, 3.1, 0.25) # av range from 0 to 3 with 0.25 steps
av_sig_grid = np.arange(0.1, 0.51, 0.1) # av sigrange 0.1 to 0.5 with 0.1 steps
rv_grid = np.arange(2.5, 5.1, 0.5) #rv range 2.5 to 5.0 with steps of 0.5
theta_prob = np.zeros((len(alpha_grid), len(age_grid), len(av_grid), len(av_sig_grid), len(rv_grid)))


field_eps = 0.0 #epilson value for field stars
theta_prob = 0.

#first = pyfits.getdata('mock4_randfield_jul12_b1_s10_av025_rv05.fits.gz', extname='AV_PROB')[0,:]

star_av_vals = np.zeros((41, nfiles))
star_rv_vals = np.zeros((10, nfiles))
fullprob = np.zeros((76, 51, 41, 10, nfiles))



for c in range(nfiles):
    print 'Reading in Star Number ', c, 'of ', nfiles
    star_av_vals[:,c] = pyfits.getdata('mock4_cl/'+filelist[c], extname='AV_PROB')[1,:]
    star_rv_vals[c] = pyfits.getdata('mock4_cl/'+filelist[c], extname='RV_PROB')[1,:]
    fullprob[:,:,:,:,c] = np.transpose(pyfits.getdata('mock4_cl/'+filelist[c], extname='FULL_PROB'))
    fullprob[:,:,:,:,c] /= fullprob[:,:,:,:,c].sum()


def lnprob(theta, filelist=filelist, theta_prob=theta_prob):
    for c in range(nfiles):
        
        
    #initizlie model cluster PDF
        p_gamma_theta_cluster = fullprob[:,:,:,:,c]*0
        #print 'p_gamma_theta_cluster = fullprob*0'
    #PDF for field, set flat by default
        p_gamma_theta_field = p_gamma_theta_cluster*0.0 + 1.0
        #print 'p_gamma_theta_field = p_gamma_theta_cluster*0.0 + 1.0'
        p_gamma_theta_field *= get_completeness()
        #print 'p_gamma_theta_field *= get_completeness()'
        p_gamma_theta_field /= np.sum(p_gamma_theta_field)
        #print 'p_gamma_theta_field /= np.sum(p_gamma_theta_field)'

        p_stellar, walpha, wage = star_prob(theta[0], theta[1])
        #print 'p_stellar done'
        if (theta[2] > 0.):
            p_av = prob_dust(star_av_vals[:,c], theta[2], theta[3])
            p_av /= p_av.sum()
            
        elif (theta[2] == 0):
            p_av = star_av_vals[:,c]*0.
            p_av = 1
            p_av /= p_av.sum()        
        elif (theta[2] < 0.):
            p_av = star_av_vals[:,c]*0. + 1e-50
        #p_gamma_theta_cluster1 = p_gamma_theta_cluster  
        #p_gamma_theta_cluster1[:,:,walpha,wage] = p_stellar*p_av[walpha]
        #for x, y in it.product(np.arange(len(fullprob[0,0,:,0])), np.arange(len(fullprob[0,0,0,:]))):
        #    p_gamma_theta_cluster[:,:,x,y] = p_stellar*p_av[x]
        #print (p_gamma_theta_cluster1 -p_gamma_theta_cluster).sum()
        #pdb.set_trace()
        
            #print x, y
        
        #pdb.set_trace()
        for x in range(len(fullprob[0,0,:,0])):
            for y in range(len(fullprob[0,0,0,:])):
                p_gamma_theta_cluster[:,:,x,y] = p_stellar*p_av[x]
        #pdb.set_trace()
        p_gamma_theta_cluster *= get_completeness()
        p_gamma_theta_cluster /= np.sum(p_gamma_theta_cluster)
        p_gamma_theta = (1. - field_eps)*p_gamma_theta_cluster + field_eps*p_gamma_theta_field
        #pdb.set_trace()        
        tprob_full = fullprob[:,:,:,:,c]*p_gamma_theta
        tprob = np.sum(tprob_full)
        #pdb.set_trace()
        #print c, theta, tprob, p_gamma_theta_cluster.sum()
        theta_prob += np.log(tprob)
        #print theta, tprob, theta_prob
    #if ((field_eps <= 0.) & (tprob < 1e-50)): 
    #    theta_prob = -np.infty
    #else:        
    #    theta_prob += np.log(tprob)
    print theta, theta_prob, p_gamma_theta_cluster.sum()
    return theta_prob


#for c in range(nfiles):
#        tprob_full = fullprob[:,:,:,:]*p_gamma_theta
#        tprob = np.sum(tprob_full)
#        if ((field_eps <= 0.) & (tprob < 1e-50)): tprob = 1e-50
#        theta_prob += np.log(tprob)
#    print theta, theta_prob
#    return theta_prob



#initialization for emcee

nwalkers = 8
ndim = 4
nburn = 50
nsteps = 50


initial = [np.array([np.random.uniform(alpha_grid.min(), alpha_grid.max()), np.random.uniform(age_grid.min(), age_grid.max()), np.random.uniform(av_grid.min(), av_grid.max()), np.random.uniform(av_sig_grid.min(), av_sig_grid.max())]) for i in xrange(nwalkers)]







sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)#, threads=nthreads)
start = time.time()
pos,prob,state = sampler.run_mcmc(initial, nburn)
sampler.reset()
sampler.run_mcmc(np.array(pos),nsteps, rstate0=state)
duration = time.time()-start





                                     
                                     
                                     
                                   
                     
                                     
                    
