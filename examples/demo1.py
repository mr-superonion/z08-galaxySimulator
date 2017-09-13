#demo1 shows how to use GalSP to generate galaxy with point source
# and write the array into a fits file

import GalSP as galsp
import numpy as np
import pyfits
import random




ndgal=12 # different galaxies generated for each shear
nshear=1 # number of shear
nrot=4 # galaxy are rotated by 4 times
nstar=80 # how many point source is used in the simulation of each galaxy
ngrid=32 # size of grid
radius=6. # extending radius of galaxy
gal_n=0.5 #parameter for Sersic index
rgal=1.5 # half light radius of galaxy
rpsf=1.5 # half light radius of PSF
truncr=8. # truncate ratio for PSF
m_psf=6. # parameter for Moffat PSF
# shear for galaxy
gamma1=np.zeros(7);gamma2=np.zeros(7)

gamma1=np.array([-0.00870246,0.01270005,0.01275759,0.01897737,-0.02373166,-0.01617663,0.02107615])
gamma2=np.array([-0.02102185, -0.02953631 ,0.02967116,-0.01457094,0.02965902,0.01259095,0.01496962])

#Creat array to store galaxy image
galimage=np.zeros((ndgal*ngrid,nshear*ngrid))
#array to store coordinates of point sources
xystaro=np.zeros((nstar,3))
xystar=np.zeros((nstar,3))
for ishear in range(nshear):
    for igal in range(ndgal): 
        # determine the intrinsic shape of galaxy
        if igal%4==0:
            #generat new galaxy
            xystaro=galsp.galaxygenerate(nstar,radius,gal_n,rgal)
            #add intrinsic ellipticity to galaxy
            #ellipticity for intrinsic galaxy
            gal_e=random.uniform(-0.1,0.1)
            xystaro=galsp.galaxyelli(gal_e,xystaro)
        else:
            xystaro=galsp.galaxyrot(np.pi/4,xystaro)
        # add shear
        xystar=galsp.galaxyshear(gamma1[ishear],gamma2[ishear],xystaro)
        #array to store one galaxy image
        galaxy1=np.zeros((ngrid,ngrid))
        galaxy1=galsp.galaxyingrid(ngrid,0.,0.,rpsf,rpsf,m_psf,truncr,xystar)
        galaxy1=galaxy1/sum(sum(galaxy1))
        #print galsp.calqua0(0,0,galaxy1,0,1)
        for iy in range(ngrid):
            for ix in range(ngrid):
                galimage[igal*ngrid+ix,ishear*ngrid+iy]=galaxy1[ix,iy]


#PSFimage
psfimage=np.zeros((ngrid,ngrid))
psfimage=galsp.mpsf(ngrid,rpsf,rpsf,m_psf,truncr,0,0,)

galimage=np.transpose(galimage)
psfimage=np.transpose(psfimage)

#write to fits file
pyfits.writeto('galaxy.fits',galimage)
pyfits.writeto('psf.fits',psfimage)
