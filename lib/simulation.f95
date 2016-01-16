subroutine galaxyGenerate(nstar,radius,xystar,sern,rgal)
implicit none
 real,parameter :: pi=acos(-1.)
 integer,intent(in) :: nstar
 real,intent(in) :: radius,sern,rgal
 real :: rstar,rstar2,thetastar,xg,yg,sx,sy,xo,yo,bm,sumflux
 real,external :: random_range
 real,intent(out) :: xystar(nstar,3)
 integer :: ig
xystar=0.
sx=0.
sy=0.
sumflux=0.
bm=2*sern-1./3.+4./405./sern
do ig=1,nstar
	rstar=random_range(0.,radius**2.)
	rstar=sqrt(rstar)
	thetastar=random_range(0.0,2.*pi)
	xg=rstar*cos(thetastar)
	yg=rstar*sin(thetastar)
	xystar(ig,1)=xg
	xystar(ig,2)=yg
	xystar(ig,3)=exp(-1.*bm*(rstar/rgal)**(1./sern))
	sx=sx+xg*exp(-1.*bm*(rstar/rgal)**(1./sern))
	sy=sy+yg*exp(-1.*bm*(rstar/rgal)**(1./sern))
	sumflux=sumflux+exp(-1.*bm*(rstar/rgal)**(1./sern))
end do
xo=sx/sumflux
yo=sy/sumflux
do ig=1,nstar
	xystar(ig,1)=xystar(ig,1)-xo
	xystar(ig,2)=xystar(ig,2)-yo
	xystar(ig,3)=xystar(ig,3)/sumflux	
end do

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function random_range(lowbound,upbound)
implicit none
real,intent(in) :: lowbound,upbound
real :: lenrandom
real :: random_range
real :: t
lenrandom=upbound-lowbound
call random_number(t)
random_range=lowbound+lenrandom*t
return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine galaxyRot(nstar,theta,xystarr,xystarr22)
implicit none
real,intent(out) :: xystarr22(nstar,3)
integer,intent(in) :: nstar
real,intent(in) :: theta
real,intent(in) :: xystarr(nstar,3)
integer :: ig

xystarr22=0.
do ig=1,nstar
	xystarr22(ig,1)=cos(theta)*xystarr(ig,1)-sin(theta)*xystarr(ig,2)
	xystarr22(ig,2)=sin(theta)*xystarr(ig,1)+cos(theta)*xystarr(ig,2)
	xystarr22(ig,3)=xystarr(ig,3)
end do


return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine galaxyShear(nstar,gamma1,gamma2,xystars,xystars2)
implicit none
integer,intent(in) :: nstar
real,intent(in) :: gamma1,gamma2
real,intent(in) :: xystars(nstar,3)
real,intent(out) :: xystars2(nstar,3)
integer :: ig

xystars2=0.
do ig=1,nstar
	xystars2(ig,1)=(1+gamma1)*xystars(ig,1)+gamma2*xystars(ig,2)
	xystars2(ig,2)=gamma2*xystars(ig,1)+(1-gamma1)*xystars(ig,2)
	xystars2(ig,3)=xystars(ig,3)
end do


return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine galaxyelli(nstar,eg,xystars,xystars2)
implicit none
integer,intent(in) :: nstar
real,intent(in) :: eg
real,intent(in) :: xystars(nstar,3)
real,intent(out) :: xystars2(nstar,3)
integer :: ig
real :: ratio

ratio=(1.-eg)**(1./2.)

do ig=1,nstar
	xystars2(ig,1)=xystars(ig,1)/ratio
	xystars2(ig,2)=xystars(ig,2)*ratio
	xystars2(ig,3)=xystars(ig,3)
end do


return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine galaxyinGrid(nstar,ngrid,x,y,re1,re2,mm,psf_trun,xystar,galaxy)
implicit none
integer,intent(in) :: nstar,ngrid
integer :: i,j,k,t1,t2,xi,yi,r1,r2
real,intent(in) :: x,y,re1,re2,mm,psf_trun
real,intent(in) :: xystar(nstar,3)
real,intent(out) :: galaxy(ngrid,ngrid)
real :: xp,yp,xs,ys,xo,yo
real,external :: gauss
real,external :: modffat
galaxy=0.
xo=ngrid/2.-x+1.
yo=ngrid/2.-y+1.
xi=aint(xo)
yi=aint(yo)
r1=aint(re1*psf_trun+1.)+3
r2=aint(re2*psf_trun+1.)+3
do i=1,nstar
	xs=xystar(i,1)
	ys=xystar(i,2)
	t1=anint(xystar(i,1))
	t2=anint(xystar(i,2))
	do k=yi+t2-r2,yi+t2+r2
	do j=xi+t1-r1,xi+t1+r1
		if(k<ngrid+1 .and. k>0 .and. j<ngrid+1 .and. j>0)then
			xp=xo+xs-j
			yp=yo+ys-k
			galaxy(j,k)=galaxy(j,k)+xystar(i,3)&
			&*modffat(re1,re2,mm,psf_trun,xp,yp)!gauss(rPSF1,rPSF2,xp,yp)!
		end if
	end do
	end do
end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mPSF(ngrid,re1,re2,mm,truncr,a,b,PSF)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: a,b,re1,re2,mm,truncr
real,intent(out) :: PSF(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer ::ia,ib
real,external :: modffat
PSF=0.
	do ib=-1*ngrid/2,ngrid/2-1
	do ia=-1*ngrid/2,ngrid/2-1
		PSF(ia,ib)=modffat(re1,re2,mm,truncr,ia+a,ib+b)
	end do
	end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gPSF(ngrid,rPSF1,rPSF2,a,b,PSF)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: a,b,rPSF1,rPSF2
real,intent(out) :: PSF(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer ::ia,ib
real,external :: gauss
PSF=0.
	do ib=-1*ngrid/2,ngrid/2-1
	do ia=-1*ngrid/2,ngrid/2-1
			PSF(ia,ib)=gauss(rPSF1,rPSF2,ia+a,ib+b)
	end do
	end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sPSF(ngrid,rPSF1,rPSF2,a,b,PSF)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: a,b,rPSF1,rPSF2
real,intent(out) :: PSF(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer ::ia,ib
real,external :: sinc2
PSF=0.
	do ib=-1*ngrid/2,ngrid/2-1
	do ia=-1*ngrid/2,ngrid/2-1
		PSF(ia,ib)=sinc2(rPSF1,rPSF2,ia+a,ib+b)
	end do
	end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dPSF(ngrid,beta,a,b,PSF)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: a,b,beta
real,intent(out) :: PSF(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer ::ia,ib
real,external :: flatdisk
PSF=0.
do ib=-1*ngrid/2,ngrid/2-1
do ia=-1*ngrid/2,ngrid/2-1
	PSF(ia,ib)=flatdisk(beta,ia+a,ib+b)
end do
end do
return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fPSF(ngrid,beta,a,b,PSF)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: a,b,beta
real,intent(out) :: PSF(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer ::ia,ib
real,external :: fang
PSF=0.
	do ib=-1*ngrid/2,ngrid/2-1
	do ia=-1*ngrid/2,ngrid/2-1
		PSF(ia,ib)=fang(beta,ia+a,ib+b)
	end do
	end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine poissonoise(ngrid,seed,pn)
implicit none
integer,intent(in) :: ngrid
integer,intent(inout) :: seed
real,intent(out) :: pn(ngrid,ngrid)
real,external :: gasdev
integer :: i,j
do j=1,ngrid
do i=1,ngrid
	pn(i,j)=gasdev(seed)
end do
end do
seed=seed+1

return
end subroutine


