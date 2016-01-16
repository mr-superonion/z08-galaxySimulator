subroutine getmax(image,ngrid,imax,x,y)
! this subroutine is used to determine the brightest pixel for the image
implicit none
integer :: ngrid,i,j,x,y
real :: image(ngrid,ngrid)
real :: imax

imax=0.
do j=1,ngrid
do i=1,ngrid
	if (image(i,j)>imax) then
		imax=image(i,j)
		x=i
		y=j
	end if

end do
end do

return
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calQua(ngrid,galaxy,a,b,Qua)
! previous edition of calQua0
implicit none
real,intent(in) :: galaxy(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer,intent(in) :: a,b,ngrid
integer :: ic,jc
real,intent(out) :: Qua

Qua=0.
do jc=-1*ngrid/2,ngrid/2-1
	do ic=-1*ngrid/2,ngrid/2-1
		Qua=Qua+((ic)**a)*((jc)**b)*galaxy(ic,jc)
	end do 
end do
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calQua0(ngrid,x,y,galaxy,a,b,Qua)
! this subroutine is used to calculate the moments of the image
! x and y are the centroid of the image
! a and b define the order of the moments
! Qua is the final outcome
implicit none
real,intent(in) :: galaxy(-1*ngrid/2:ngrid/2-1,-1*ngrid/2:ngrid/2-1)
integer,intent(in) :: a,b,ngrid
integer :: ic,jc
real,intent(in):: x,y
real,intent(out) :: Qua
Qua=0.
do jc=-1*ngrid/2+1,ngrid/2-1
	do ic=-1*ngrid/2+1,ngrid/2-1
		Qua=Qua+((ic-x)**a)*((jc-y)**b)*galaxy(ic,jc)
	end do
end do
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getep(ngrid,image,ep1,ep2)
! this subroutine is used to calculate the ellipticity of the image
! x and y are the centroid of the image
! ep1 and ep2 are ep*cos(theta) and ep*sin(theta)
implicit none
integer,intent(in) :: ngrid
real,intent(in) :: image(ngrid,ngrid)
real,intent(out) :: ep1,ep2
real :: Q11,Q20,Q02
integer :: i,j


call calQua0(ngrid,0.,0.,image,2,0,Q20)
call calQua0(ngrid,0.,0.,image,1,1,Q11)
call calQua0(ngrid,0.,0.,image,0,2,Q02)
ep1=(Q20-Q02)/(Q20+Q02)
ep2=2.*Q11/(Q20+Q02)

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getradius(image,ngrid,ratio,radius)
implicit none
real,parameter :: pi=acos(-1.)
integer :: x,y,m,n,rrr,maxradius
real :: imax,ep1,ep2,ep,bsa,sumim,sumcou
real :: imagecount(ngrid,ngrid)
real,intent(out) :: radius
integer,intent(in) :: ngrid
real,intent(in) :: image(ngrid,ngrid)
real,intent(in) :: ratio


!call getmax(image,ngrid,imax,x,y)
! set the origin at ngrid/2+1
x=ngrid/2+1
y=ngrid/2+1
call getflux(image,ngrid,sumim,maxradius)
! begin from HLR=sqrt(3)
radius=3.
do rrr=1,ngrid*ngrid
	radius=radius+1.
	imagecount=0.
	forall(m=1:ngrid,n=1:ngrid,(m-x)**2.+(n-y)**2.<=radius)
		imagecount(m,n)=image(m,n)
	end forall
	sumcou=sum(imagecount)
	if (sumcou>ratio*sumim) then 
		exit
	end if		
end do
radius=sqrt(radius)

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getflux(image,ngrid,flux,maxradius)
implicit none
integer :: x,y,m,n,rrr
real :: fluxr(ngrid/2-2)
real :: imagecount(ngrid,ngrid)
integer,intent(out) :: maxradius
real,intent(out) :: flux
integer,intent(in) :: ngrid
real,intent(in) :: image(ngrid,ngrid)

!call getmax(image,ngrid,imax,x,y)
x=ngrid/2+1
y=ngrid/2+1
fluxr=0.
maxradius=ngrid/2-2
do rrr=2,ngrid/2-2
	imagecount=0.
	forall(m=1:ngrid,n=1:ngrid,(m-x)**2.+(n-y)**2.<=rrr**2.)
		imagecount(m,n)=image(m,n)
	end forall
	fluxr(rrr)=sum(imagecount)
	if (fluxr(rrr)<=fluxr(rrr-1)) then
		maxradius=max(rrr-1,4)
		exit
	end if
end do
flux=fluxr(maxradius)
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine moment(data,n,ave,sdev)
! used to do statistics
implicit none
real :: sdev,data(n)
real*8 :: p,s,ave,var
integer :: i,n
s=0.
do i=1,n
   s=s+data(i)
end do
ave=s/n
var=0.
do i=1,n
   s=data(i)-ave
   p=s*s
   var=var+p
end do
var=var/(n-1)
sdev=sqrt(var)
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getshear(allplus,allminus,all12,times,shear1,shear2,error1,error2)
! this program is used in parallel computing
implicit none
 integer,intent(in) :: times
 real,intent(in) :: allplus(times)
 real,intent(in) :: allminus(times)
 real,intent(in) :: all12(times)
 real,intent(out) :: shear1,shear2,error1,error2
 real :: stdup1,stdup2,stddown,ng
 real*8 :: up1,up2,down
call moment(allplus,times,down,stddown)
call moment(allminus,times,up1,stdup1)
call moment(all12,times,up2,stdup2)
shear1=up1/down/2.
shear2=up2/down
ng=times+0.
error1=stdup1/down/sqrt(ng)/2.
error2=stdup2/down/sqrt(ng)

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine getshearF3(galaxy1,PSF1,ngrid,up1,up2,down)
implicit none
 real,parameter :: pi=acos(-1.)
 integer,intent(in) :: ngrid
 real,intent(in) :: galaxy1(ngrid,ngrid)
 real,intent(in) :: PSF1(ngrid,ngrid)
 real,intent(out) :: up1,up2,down
 real :: PSFweight(ngrid,ngrid),PSFu(ngrid,ngrid)
 real :: galaxyu(ngrid,ngrid)
 real :: rfour,rtarget,Q20,Q02,Q11,Q40,Q04,Q22,D4
 integer :: i,j
! calculate the HLR of Fourier Mode
PSFu=PSF1
call getradius(PSFu,ngrid,exp(-1.),rfour)
rfour=rfour*0.5
rtarget=1./sqrt(2.)/(rfour)*ngrid/2./pi

forall (j=1:ngrid, i=1:ngrid, PSFu(i,j)<PSFu(ngrid/2+1,ngrid/2+1)/1.e4)
	PSFu(i,j)=10.e7
end forall
!Gaussianization
call gPSF(ngrid,rfour,rfour,0.,0.,PSFweight)
PSFweight=PSFweight/PSFu
galaxyu=galaxy1*PSFweight

!calculate the moments of galaxy
call calQua(ngrid,galaxyu,2,0,Q20)
call calQua(ngrid,galaxyu,1,1,Q11)
call calQua(ngrid,galaxyu,0,2,Q02)
call calQua(ngrid,galaxyu,4,0,Q40)
call calQua(ngrid,galaxyu,2,2,Q22)
call calQua(ngrid,galaxyu,0,4,Q04)
D4   = Q40+2.*Q22+Q04
down = Q20+Q02-0.5*(rtarget*2*pi/(ngrid))**2*D4
up1  = -1.*(Q20-Q02)
up2  = -1.*Q11
return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

