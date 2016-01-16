subroutine getSNR(image,ngrid,SNR)
! this subroutine is used to calculate the signal to noise ratio approximately
implicit none
real,intent(in) :: image(ngrid,ngrid)
integer,intent(in) :: ngrid
real,intent(out) :: SNR
integer :: x,y,m,n,nn
real :: ratio,imax,radius
ratio=0.5
!call getmax(image,ngrid,imax,x,y)
x=ngrid/2+1
y=ngrid/2+1
call getradius(image,ngrid,ratio,radius)
SNR=0.
do n=y-aint(radius),y+aint(radius)
do m=x-aint(radius),x+aint(radius)
	nn=nn+1
	SNR=SNR+image(m,n)
end do
end do
SNR=SNR/(radius)/2.

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine trim_power(ngrid,powerin,powerout)
implicit none

integer,intent(in) :: ngrid
real,intent(in) :: powerin(ngrid,ngrid)
real,intent(out) :: powerout(ngrid,ngrid)

integer :: i,j
real :: temp

temp=0.
do i=1,ngrid
	temp=temp+powerin(i,1)+powerin(i,ngrid)
end do	
do i=2,ngrid-1
	temp=temp+powerin(1,i)+powerin(ngrid,i)
end do	
	
! Average 
temp=temp/(4.*ngrid-4.)
powerout=powerin-temp

	
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine snellim(galaxyu,galaxyo,ngrid)
 integer,intent(in) :: ngrid
 real,intent(in) :: galaxyu(ngrid,ngrid)
 real,intent(out) :: galaxyo(ngrid,ngrid)
 integer ::i,j
 real :: x,y
 real :: avernoise,avernum
avernoise=0.
avernum=0.
x=ngrid/2.+1.
y=ngrid/2.+1.
do j=1,ngrid
do i=1,ngrid
	if ((i-x)**2.+(j-y)**2.>=(ngrid/2.-1.5)**2.) then
		avernoise=avernoise+galaxyu(i,j)
		avernum=avernum+1.
	end if
end do
end do
avernoise=avernoise/avernum
do j=1,ngrid
do i=1,ngrid
		galaxyo(i,j)=galaxyu(i,j)-avernoise
end do
end do

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
