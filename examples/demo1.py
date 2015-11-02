#demo1 shows how to use GalSP to generate galaxy with point source
# and write the array into a fits file


import GalSP as galsp
import numpy as np
import pyfits
nstar=5 # how many point source is used in the simulation of each galaxy
ngrid=48 # size of grid
# shear for galaxy
gamma1=0.01
gamma2=0.03
#array to stor coordinates of point sources
xystaro=np.zeros((nstar,3))
xystar=np.zeros((nstar,3))
#generat the galaxy
xystaro=galsp.galaxygenerate(nstar,4.,1.2)
xystar=galsp.galaxyshear(gamma1,gamma2,xystaro)
galaxy1=galsp.galaxyingrid(ngrid,ngrid/2.+1,ngrid/2.+1,2.,3.,xystar)

#write to fits file
pyfits.writeto('galaxy.fits',galaxy1)
