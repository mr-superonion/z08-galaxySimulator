function gauss(beta1,beta2,xgs,ygs)
implicit none
 real,parameter:: pi=acos(-1.)
 real :: gauss
 real,intent(in) :: beta1,beta2
 real,intent(in) :: xgs,ygs
	gauss=1./(2.*pi*beta1*beta2)*exp((-0.5)*(xgs**2./beta1**2.+ygs**2./beta2**2.))
return
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sinc(k1,k2,x,y)
implicit none
 real,parameter :: pi=acos(-1.)
 real :: sinc
 real:: x,y,k1,k2

if (x**2+y**2==0.) then
	sinc=1
else if (x==0.) then
	sinc=sin(y/k2)*k2/y
else if (y==0.) then
	sinc=sin(x/k1)*k1/x
else
	sinc=sin(x/k1)*k1/x*sin(y/k2)*k2/y
end if

return
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function sinc2(k1,k2,x,y)
implicit none
 real,parameter :: pi=acos(-1.)
 real :: sinc2
 real:: x,y,k1,k2

if (x**2+y**2==0.) then
	sinc2=1.
else if (x==0.) then
	sinc2=(sin(y/k2))**2.*k2**2./y**2.
else if (y==0.) then
	sinc2=(sin(x/k1))**2.*k1**2./x**2.
else
	sinc2=(sin(y/k2))**2.*k2**2./y**2.*(sin(x/k1))**2.*k1**2./x**2.
end if

return
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function modffat(re1,re2,mm,truncr,x,y)
implicit none

 real :: modffat
 real,intent(in):: x,y
 real,intent(in) :: re1,re2,mm,truncr
 real :: r2

r2=x**2./re1**2.+y**2./re2**2.
if (r2<truncr**2.) then
	r2=r2*(2.**(1./(mm-1.))-1.)
	modffat=(1.+r2)**(-1.*mm)
else
	modffat=0.
end if

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function kaiser(alph,kk,x)
implicit none
 real,parameter :: pi=acos(-1.)
 real :: alph,kk,x,Kaiser
 
if (x**2<kk**2) then
	kaiser=sqrt(2./pi)*sinh(alph*sqrt(kk**2-x**2))/sqrt(kk**2-x**2)
else
	kaiser=0.
end if

return
end function	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function flatdisk(beta,x,y)
implicit none
 real,intent(in) :: beta,x,y
 real :: r2,flatdisk
flatdisk=0.
r2=x**2.+y**2.
if(r2<=beta**2.) flatdisk=1.

return
end function	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function fang(beta,x,y)
implicit none
 real,intent(in) :: beta,x,y
 real :: fang
fang=0.
if(abs(x)<=beta/2 .and. abs(y)<=beta/2) fang=1.

return
end function	

