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
real :: image1(ngrid,ngrid)
real,intent(out) :: ep1,ep2
integer :: xm,ym
real :: x,y
real :: Q11,Q20,Q02,Q10,Q01,Q00,imax
integer :: i,j

image1=image
call getmax(image1,ngrid,imax,xm,ym)
do j=1,ngrid
do i=1,ngrid
	if (image1(i,j)<0.5*imax) image1(i,j)=0.
end do
end do
call calQua0(ngrid,0.,0.,image1,0,0,Q00)
call calQua0(ngrid,0.,0.,image1,1,0,Q10)
call calQua0(ngrid,0.,0.,image1,0,1,Q01)
x=Q10/Q00
y=Q01/Q00
call calQua0(ngrid,x,y,image1,2,0,Q20)
call calQua0(ngrid,x,y,image1,1,1,Q11)
call calQua0(ngrid,x,y,image1,0,2,Q02)
ep1=(Q20-Q02)/(Q20+Q02)
ep2=2.*Q11/(Q20+Q02)

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getradius(image,ngrid,ratio,radius)
! this subroutine is used to determine the Half Width Half Maximum of the image(if the image has ellipticity, it determines the long axis)
! we first calculate the area of the pixel whose value is larger than 0.5*imax
! and then calculate the radius using pi*a*b=S approximately
implicit none
real,parameter :: pi=acos(-1.)
integer :: x,y,m, n
real :: imax,area,ep1,ep2,ep,bsa
real :: imagecount(ngrid,ngrid)
real,intent(out) :: radius
integer,intent(in) :: ngrid
real,intent(in) :: image(ngrid,ngrid)
real,intent(in) :: ratio


call getmax(image,ngrid,imax,x,y)
imagecount=0.
forall(m=1:ngrid,n=1:ngrid,image(m,n)>=ratio*imax)
	imagecount(m,n)=1.
end forall

area=sum(imagecount)
call getep(ngrid,image,ep1,ep2)
ep=sqrt(ep1**2.+ep2**2.)

bsa=sqrt((1-ep)/(1+ep))
radius=sqrt(area/bsa/pi)

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


subroutine getSNR(image,ngrid,SNR)
! this subroutine is used to calculate the signal to noise ratio approximately
implicit none
real,intent(in) :: image(ngrid,ngrid)
integer,intent(in) :: ngrid
real,intent(out) :: SNR
integer :: x,y,m,n,nn
real :: ratio,imax,p,ave,radius
ratio=.5
call getmax(image,ngrid,imax,x,y)
call getradius(image,ngrid,ratio,radius)
nn=0
p=0.
SNR=0.
do n=y-aint(radius+0.5),y+aint(radius+0.5)
do m=x-aint(radius+0.5),x+aint(radius+0.5)
	nn=nn+1
	SNR=SNR+image(m,n)
end do
end do
SNR=SNR/(radius)/2.

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


