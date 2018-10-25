import numpy
import bigfile 
import math as m
from kdcount import correlate
from scipy.optimize import curve_fit
import numpy as np
import os
import pickle
import random
import h5py
import scipy
import multiprocessing 
import multiprocessing.pool
import time
global x,y,z,SM,lcut,NBINS,rho_crit,mid,width
rho_crit=2.77519737e11
lcuts=[9.5] 
BOXSIZE = 400.0
NBINS = 10
RMIN = 10**(-3)
RMAX = 10**2
N_RAND = 5e5

MASS_LIMIT = 10**9

LOGSCALE = True
WRAP = True
JK_SLICE = 0
z=10
        
def poiss(rmin,rmax):
    p=4./3*scipy.pi*(rmax**3-rmin**3)/400**3
    return p 

def construct_corr_info(COM,MASS_LIMIT, BOXSIZE = BOXSIZE):
    N = len(COM)  
    if N > 0:
#        print 'Total number of galaxies = ', N
#        print 'Done constructing info'
        return COM,N
    else:
#        print 'No galaxies'
        return np.ones(2),N
    
def correlate_info(data, NBINS = NBINS, RMIN = RMIN, RMAX = RMAX, BOXSIZE = BOXSIZE, WRAP = WRAP):
    if data is not None:
        if RMAX is None:
            RMAX = BOXSIZE
        
        if WRAP:
            wrap_length = BOXSIZE
        else:
            wrap_length = None
        
        dataset = correlate.points(data, boxsize = wrap_length)  
        
        binning = correlate.RBinning(np.logspace(np.log10(RMIN),np.log10(RMAX),NBINS+1))
        
#	RR=N**2*numpy.asarray([poiss(rbin[i],rbin[i+1]) for i in range(0,nbins)])
        DD = correlate.paircount(dataset, dataset, binning, np=16)
        DD = DD.sum1
        
#        print 'Done correlating'
        r = binning.centers
        return r, DD
    else:
        return None, None


if (z==8):
	pig='086'
if (z==9):
        pig='066'
if (z==10):
        pig='054'
if (z==7.5):
        pig='141'


galaxy=bigfile.BigFile('/nfs/nas-0-1/akbhowmi/my_galaxy_sample_bluetides_closest_comoving_distance_with_luminosity/'+pig+'/')
SM=galaxy.open('galaxy_stellar_mass')[:]
HM=galaxy.open('galaxy_host_mass')[:]
hostid=galaxy.open('galaxy_host_id')[:]
tag=galaxy.open('central_satellite_tag')[:]
galaxy_positions=galaxy.open('galaxy_center_of_mass')[:]



lcuts=numpy.array([7.5,8.0,8.5,9.0,9.5,10.0,10.5])
lcuts+=0.25

for lcut in reversed(lcuts):
    mask=SM>10**lcut
    data_s=galaxy_positions[mask]

    N=len(data_s)
    r, DD = correlate_info(data_s)
    
    binning = correlate.RBinning(np.logspace(np.log10(RMIN),np.log10(RMAX),NBINS+1))
    rbin=binning.edges
    RR=(N**2-N)*np.asarray([poiss(rbin[i],rbin[i+1]) for i in range(0,NBINS)])
    
    xi=DD/RR-1
    dxi=np.sqrt(DD)/RR
    r=binning.centers

    label_DD = './galaxy_correlation_functions/DD_z-%.1f_SM_cut_galaxy_COM.pickle'%(z)
    label_xi = './galaxy_correlation_functions/xi_z%.1f_%.1f_SM_cut_galaxy_COM.pickle'%(z,lcut)
    
    pickle.dump([r, DD], open(label_DD, 'w'))
    pickle.dump([r, xi,dxi,N], open(label_xi, 'w'))


    

