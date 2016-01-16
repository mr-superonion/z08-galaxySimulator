subroutine transform(image,power,ngrid)
implicit none
  integer,intent(in) :: ngrid
  real,intent(in) :: image(ngrid,ngrid)
  complex,intent(out) :: power(ngrid,ngrid)
  integer :: i,j
  real,parameter :: pi=acos(-1.)
  real wsave(4*ngrid+15)

call cffti ( ngrid, wsave )

do j=1,ngrid
  do i=1,ngrid
  	power(i,j)=cmplx(image(i,j),0.)
    power(i,j)=power(i,j)*(-1)**(i-1+j-1)
  end do
end do

call cfftf_2d(ngrid,ngrid,power,wsave)


return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transformff(ngrid,image,imageo)

integer,intent(in) :: ngrid
real,intent(in) :: image(ngrid,ngrid)
real,intent(out) :: imageo(ngrid,ngrid)
complex :: imagec(ngrid,ngrid)
call transform(image,imagec,ngrid)
imageo=abs(imagec)**2.

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transformfi(ngrid,image,imageout)
implicit none
  integer,intent(in) :: ngrid
  real,intent(in) :: image(ngrid,ngrid)
  real,intent(out) :: imageout(ngrid,ngrid)
  complex :: power(ngrid,ngrid)
  real :: image2(ngrid,ngrid)
  complex :: power2(ngrid,ngrid)
  real :: filter(ngrid,ngrid)
  complex :: filterp(ngrid,ngrid)
  integer :: i,j
  real:: wsave1(4*ngrid+15)
!call gPSF(ngrid,2.,2.,0.,0.,filter)
call cffti ( ngrid, wsave1 )

do j=1,ngrid
do i=1,ngrid
  	power(i,j)=cmplx(image(i,j),0.)
	!filterp(i,j)=cmplx(filter(i,j),0.)
end do
end do

call cfftf_2d(ngrid,ngrid,power,wsave1)
!call cfftf_2d(ngrid,ngrid,filterp,wsave1)


do j=1,ngrid
do i=1,ngrid
  	image2(i,j)=abs(power(i,j))**2.!*abs(filterp(i,j))
  	power2(i,j)=cmplx(image2(i,j),0.)
  	power2(i,j)=power2(i,j)*(-1)**(i-1+j-1)
end do
end do

call cfftb_2d(ngrid,ngrid,power2,wsave1)
imageout=real(power2)/ngrid**4.

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpolation(image,imageb,ngrid)
implicit none
  integer,intent(in) :: ngrid
  real,intent(in) :: image(ngrid,ngrid)
  real,intent(out) :: imageb(4*ngrid,4*ngrid)
  complex :: power(4*ngrid,4*ngrid)
  integer :: i,j
  real,parameter :: pi=acos(-1.)
  real :: Q10
  real:: wsave1(16*ngrid+15)
 
call cffti ( 4*ngrid, wsave1 )

imageb=0.
do j=1,ngrid
do i=1,ngrid
  	imageb(4*i-3,4*j-3)=image(i,j)
end do
end do

do j=1,4*ngrid
do i=1,4*ngrid
  	power(i,j)=cmplx(imageb(i,j),0.00)
end do
end do

call cfftf_2d(4*ngrid,4*ngrid,power,wsave1)
do j=1,4*ngrid
do i=1,4*ngrid
	!if((abs(i)>ngrid/2.and.abs(i)<7*ngrid/2).or.(abs(j)>ngrid/2.and.abs(j)<7*ngrid/2)) power(i,j)=0.
end do
end do

call cfftb_2d(4*ngrid,4*ngrid,power,wsave1)
imageb=real(power)
!call draw(imageb,4*ngrid)
!pause
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine multitransform00(image1,image2,image3,image4,image,ngrid)
implicit none
  integer :: i,j
  real:: ii,jj
  integer,intent(in) :: ngrid
  real,intent(in) :: image1(ngrid,ngrid)
  complex*16 :: power1(ngrid,ngrid)
  real,intent(in) :: image2(ngrid,ngrid)
  complex*16 :: power2(ngrid,ngrid)
  real,intent(in) :: image3(ngrid,ngrid)
  complex*16 :: power3(ngrid,ngrid)
  real,intent(in) :: image4(ngrid,ngrid)
  complex*16 :: power4(ngrid,ngrid)
  complex*16,intent(out) :: image(ngrid,ngrid)
  real wsave(4*ngrid+15)
  real,parameter :: pi=acos(-1.)
  complex*16 :: ef1,ef2,ef3,ef4
  real :: kk1,kk2,kk3,kk4
  real :: a1,b1,a2,b2,a3,b3,a4,b4
  common/shift/ a1,b1,a2,b2,a3,b3,a4,b4

  


call cffti(ngrid, wsave)
do j=1,ngrid
do i=1,ngrid
	power1(i,j)=cmplx(image1(i,j),0.)
	power2(i,j)=cmplx(image2(i,j),0.)
	power3(i,j)=cmplx(image3(i,j),0.)
	power4(i,j)=cmplx(image4(i,j),0.)
	power1(i,j)=power1(i,j)*(-1)**(i-1+j-1)
    power2(i,j)=power2(i,j)*(-1)**(i-1+j-1)
    power3(i,j)=power3(i,j)*(-1)**(i-1+j-1)
    power4(i,j)=power4(i,j)*(-1)**(i-1+j-1) 
end do
end do

call cfftf_2d(ngrid,ngrid,power1,wsave)
call cfftf_2d(ngrid,ngrid,power2,wsave)
call cfftf_2d(ngrid,ngrid,power3,wsave)
call cfftf_2d(ngrid,ngrid,power4,wsave)


do j=1,ngrid
do i=1,ngrid
	ii=i-ngrid/2-1.
	jj=j-ngrid/2-1.
	kk1=ii*a1+jj*b1
	kk2=ii*a2+jj*b2
	kk3=ii*a3+jj*b3
	kk4=ii*a4+jj*b4
	ef1=cmplx(cos(-2.*pi*kk1/ngrid),sin(-2.*pi*kk1/ngrid))
	ef2=cmplx(cos(-2.*pi*kk2/ngrid),sin(-2.*pi*kk2/ngrid))
	ef3=cmplx(cos(-2.*pi*kk3/ngrid),sin(-2.*pi*kk3/ngrid))
	ef4=cmplx(cos(-2.*pi*kk4/ngrid),sin(-2.*pi*kk4/ngrid))

	power1(i,j)=power1(i,j)*ef1
	power2(i,j)=power2(i,j)*ef2
	power3(i,j)=power3(i,j)*ef3
	power4(i,j)=power4(i,j)*ef4

end do
end do

do j=1,ngrid
do i=1,ngrid
	image(i,j)=power1(i,j)+power2(i,j)+power3(i,j)+power4(i,j)
end do
end do


return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine multitransform4(nexp,ngrid,shift,imagein,imageout,det)
implicit none
  integer :: i,j,k
  real:: ii,jj
  integer,intent(in) :: nexp
  integer,intent(in) :: ngrid
  real,intent(in) :: shift(2,nexp)
  real,intent(in) :: imagein(ngrid,ngrid,nexp)
  complex :: power1(ngrid,ngrid)
  complex :: power2(ngrid,ngrid)
  complex :: power3(ngrid,ngrid)
  complex :: power4(ngrid,ngrid)
  complex*16,intent(out) :: imageout(2*ngrid,2*ngrid)
  real wsave(4*ngrid+15)
  real,parameter :: pi=acos(-1.)
  complex*16 :: ef1,ef2,ef3,ef4
  real :: kk1,kk2,kk3,kk4
  complex*16,intent(in) :: det(4,4)
  


call cffti(ngrid, wsave)
do j=1,ngrid
do i=1,ngrid
	power1(i,j)=cmplx(imagein(i,j,1),0.)
	power2(i,j)=cmplx(imagein(i,j,2),0.)
	power3(i,j)=cmplx(imagein(i,j,3),0.)
	power4(i,j)=cmplx(imagein(i,j,4),0.)
end do
end do


call cfftf_2d(ngrid,ngrid,power1,wsave)
call cfftf_2d(ngrid,ngrid,power2,wsave)
call cfftf_2d(ngrid,ngrid,power3,wsave)
call cfftf_2d(ngrid,ngrid,power4,wsave)


do j=1,ngrid
do i=1,ngrid
	ii=i-1.
	jj=j-1.
	kk1=ii*shift(1,1)+jj*shift(2,1)
	kk2=ii*shift(1,2)+jj*shift(2,2)
	kk3=ii*shift(1,3)+jj*shift(2,3)
	kk4=ii*shift(1,4)+jj*shift(2,4)
	ef1=cmplx(cos(-2.*pi*kk1/ngrid),sin(-2.*pi*kk1/ngrid))
	ef2=cmplx(cos(-2.*pi*kk2/ngrid),sin(-2.*pi*kk2/ngrid))
	ef3=cmplx(cos(-2.*pi*kk3/ngrid),sin(-2.*pi*kk3/ngrid))
	ef4=cmplx(cos(-2.*pi*kk4/ngrid),sin(-2.*pi*kk4/ngrid))
	power1(i,j)=power1(i,j)*ef1
	power2(i,j)=power2(i,j)*ef2
	power3(i,j)=power3(i,j)*ef3
	power4(i,j)=power4(i,j)*ef4
end do
end do


do j=1,ngrid
do i=1,ngrid
	imageout(i+ngrid,j+ngrid)=det(1,1)*power1(i,j)+det(1,2)*power2(i,j)+det(1,3)*power3(i,j)+det(1,4)*power4(i,j)
	imageout(i,j+ngrid)=det(2,1)*power1(i,j)+det(2,2)*power2(i,j)+det(2,3)*power3(i,j)+det(2,4)*power4(i,j)
	imageout(i+ngrid,j)=det(3,1)*power1(i,j)+det(3,2)*power2(i,j)+det(3,3)*power3(i,j)+det(3,4)*power4(i,j)
	imageout(i,j)=det(4,1)*power1(i,j)+det(4,2)*power2(i,j)+det(4,3)*power3(i,j)+det(4,4)*power4(i,j)
end do
end do

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine multitransformR(image1,image2,image3,image4,image,det1,det2,det3,det4,ngrid)
implicit none
  integer :: i,j
  real:: ii,jj
  integer,intent(in) :: ngrid
  real,intent(in) :: image1(ngrid,ngrid)
  complex*16 :: power1(ngrid,ngrid)
  real,intent(in) :: image2(ngrid,ngrid)
  complex*16 :: power2(ngrid,ngrid)
  real,intent(in) :: image3(ngrid,ngrid)
  complex*16 :: power3(ngrid,ngrid)
  real,intent(in) :: image4(ngrid,ngrid)
  complex*16 :: power4(ngrid,ngrid)
  complex*16,intent(out) :: image(2*ngrid,2*ngrid)
  real wsave(4*ngrid+15)
  real,parameter :: pi=acos(-1.)
  complex*16 :: ef1,ef2,ef3,ef4
  real :: kk1,kk2,kk3,kk4
  real :: a1,b1,a2,b2,a3,b3,a4,b4
  common/shift/ a1,b1,a2,b2,a3,b3,a4,b4
  complex*16,intent(in) :: det1(4,4),det2(4,4),det3(4,4),det4(4,4)
  


call cffti(ngrid, wsave)
do j=1,ngrid
do i=1,ngrid
	power1(i,j)=cmplx(image1(i,j),0.)
	power2(i,j)=cmplx(image2(i,j),0.)
	power3(i,j)=cmplx(image3(i,j),0.)
	power4(i,j)=cmplx(image4(i,j),0.)
end do
end do

call cfftf_2d(ngrid,ngrid,power1,wsave)
call cfftf_2d(ngrid,ngrid,power2,wsave)
call cfftf_2d(ngrid,ngrid,power3,wsave)
call cfftf_2d(ngrid,ngrid,power4,wsave)

do j=1,ngrid
do i=1,ngrid
	ii=i-1.
	jj=j-1.
	kk1=ii*a1+jj*b1
	kk2=ii*a2+jj*b2
	kk3=ii*a3+jj*b3
	kk4=ii*a4+jj*b4
	ef1=cmplx(cos(-2.*pi*kk1/ngrid),sin(-2.*pi*kk1/ngrid))
	ef2=cmplx(cos(-2.*pi*kk2/ngrid),sin(-2.*pi*kk2/ngrid))
	ef3=cmplx(cos(-2.*pi*kk3/ngrid),sin(-2.*pi*kk3/ngrid))
	ef4=cmplx(cos(-2.*pi*kk4/ngrid),sin(-2.*pi*kk4/ngrid))
	power1(i,j)=power1(i,j)*ef1
	power2(i,j)=power2(i,j)*ef2
	power3(i,j)=power3(i,j)*ef3
	power4(i,j)=power4(i,j)*ef4
end do
end do

do j=1,ngrid
do i=1,ngrid
	image(i+ngrid,j+ngrid)=det1(1,1)*power1(i,j)+det1(1,2)*power2(i,j)+det1(1,3)*power3(i,j)+det1(1,4)*power4(i,j)
	image(i,j+ngrid)=det1(2,1)*power1(i,j)+det1(2,2)*power2(i,j)+det1(2,3)*power3(i,j)+det1(2,4)*power4(i,j)
	image(i+ngrid,j)=det1(3,1)*power1(i,j)+det1(3,2)*power2(i,j)+det1(3,3)*power3(i,j)+det1(3,4)*power4(i,j)
	image(i,j)=det1(4,1)*power1(i,j)+det1(4,2)*power2(i,j)+det1(4,3)*power3(i,j)+det1(4,4)*power4(i,j)
	!if(i==1) then
	!	image(i+ngrid,j+ngrid)=det2(1,1)*power1(i,j)+det2(1,2)*power2(i,j)+det2(1,3)*power3(i,j)+det2(1,4)*power4(i,j)
	!	image(i,j+ngrid)=det2(2,1)*power1(i,j)+det2(2,2)*power2(i,j)+det2(2,3)*power3(i,j)+det2(2,4)*power4(i,j)
	!	image(i+ngrid,j)=det2(3,1)*power1(i,j)+det2(3,2)*power2(i,j)+det2(3,3)*power3(i,j)+det2(3,4)*power4(i,j)
	!	image(i,j)=det2(4,1)*power1(i,j)+det2(4,2)*power2(i,j)+det2(4,3)*power3(i,j)+det2(4,4)*power4(i,2)
	!end if
	!if(j==1) then
	!	image(i+ngrid,j+ngrid)=det3(1,1)*power1(i,j)+det3(1,2)*power2(i,j)+det3(1,3)*power3(i,j)+det3(1,4)*power4(i,j)
	!	image(i,j+ngrid)=det3(2,1)*power1(i,j)+det3(2,2)*power2(i,j)+det3(2,3)*power3(i,j)+det3(2,4)*power4(i,j)
	!	image(i+ngrid,j)=det3(3,1)*power1(i,j)+det3(3,2)*power2(i,j)+det3(3,3)*power3(i,j)+det3(3,4)*power4(i,j)
	!	image(i,j)=det3(4,1)*power1(i,j)+det3(4,2)*power2(i,j)+det3(4,3)*power3(i,j)+det3(4,4)*power4(i,j)
	!end if
	!if(j==1 .and. i==1) then
	!	image(i+ngrid,j+ngrid)=det4(1,1)*power1(i,j)+det4(1,2)*power2(i,j)+det4(1,3)*power3(i,j)+det4(1,4)*power4(i,j)
	!	image(i,j+ngrid)=det4(2,1)*power1(i,j)+det4(2,2)*power2(i,j)+det4(2,3)*power3(i,j)+det4(2,4)*power4(i,j)
	!	image(i+ngrid,j)=det4(3,1)*power1(i,j)+det4(3,2)*power2(i,j)+det4(3,3)*power3(i,j)+det4(3,4)*power4(i,j)
	!	image(i,j)=det4(4,1)*power1(i,j)+det4(4,2)*power2(i,j)+det4(4,3)*power3(i,j)+det4(4,4)*power4(i,j)
	!end if
		
end do
end do


return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine multitransform5(image1,image2,image3,image4,image5,image,det,ngrid)
implicit none
  integer :: i,j
  real:: ii,jj
  integer,intent(in) :: ngrid
  real,intent(in) :: image1(ngrid,ngrid)
  complex :: power1(ngrid,ngrid)
  real,intent(in) :: image2(ngrid,ngrid)
  complex :: power2(ngrid,ngrid)
  real,intent(in) :: image3(ngrid,ngrid)
  complex :: power3(ngrid,ngrid)
  real,intent(in) :: image4(ngrid,ngrid)
  complex :: power4(ngrid,ngrid)
  real,intent(in) :: image5(ngrid,ngrid)
  complex :: power5(ngrid,ngrid)
  complex,intent(out) :: image(ngrid,ngrid)
  real wsave(4*ngrid+15)
  real,parameter :: pi=acos(-1.)
  complex :: ef1,ef2,ef3,ef4,ef5
  real :: kk1,kk2,kk3,kk4,kk5
  real :: a1,b1,a2,b2,a3,b3,a4,b4,a5,b5
  common/shift/ a1,b1,a2,b2,a3,b3,a4,b4,a5,b5
  complex,intent(in) :: det(5,5)
  


call cffti(ngrid, wsave)
do j=1,ngrid
do i=1,ngrid
	power1(i,j)=cmplx(image1(i,j),0.)
	power2(i,j)=cmplx(image2(i,j),0.)
	power3(i,j)=cmplx(image3(i,j),0.)
	power4(i,j)=cmplx(image4(i,j),0.)
	power5(i,j)=cmplx(image5(i,j),0.)
	power1(i,j)=power1(i,j)*(-1)**(i-1+j-1)
    power2(i,j)=power2(i,j)*(-1)**(i-1+j-1)
    power3(i,j)=power3(i,j)*(-1)**(i-1+j-1)
    power4(i,j)=power4(i,j)*(-1)**(i-1+j-1) 
    power5(i,j)=power5(i,j)*(-1)**(i-1+j-1) 
end do
end do

call cfftf_2d(ngrid,ngrid,power1,wsave)
call cfftf_2d(ngrid,ngrid,power2,wsave)
call cfftf_2d(ngrid,ngrid,power3,wsave)
call cfftf_2d(ngrid,ngrid,power4,wsave)
call cfftf_2d(ngrid,ngrid,power5,wsave)

do j=1,ngrid
do i=1,ngrid
	ii=i-ngrid/2-1.
	jj=j-ngrid/2-1.
	kk1=ii*a1+jj*b1
	kk2=ii*a2+jj*b2
	kk3=ii*a3+jj*b3
	kk4=ii*a4+jj*b4
	kk5=ii*a5+jj*b5
	ef1=cmplx(cos(-2.*pi*kk1/ngrid),sin(-2.*pi*kk1/ngrid))
	ef2=cmplx(cos(-2.*pi*kk2/ngrid),sin(-2.*pi*kk2/ngrid))
	ef3=cmplx(cos(-2.*pi*kk3/ngrid),sin(-2.*pi*kk3/ngrid))
	ef4=cmplx(cos(-2.*pi*kk4/ngrid),sin(-2.*pi*kk4/ngrid))
	ef5=cmplx(cos(-2.*pi*kk5/ngrid),sin(-2.*pi*kk5/ngrid))
	power1(i,j)=power1(i,j)*ef1
	power2(i,j)=power2(i,j)*ef2
	power3(i,j)=power3(i,j)*ef3
	power4(i,j)=power4(i,j)*ef4
	power5(i,j)=power5(i,j)*ef5
end do
end do

do j=1,ngrid
do i=1,ngrid
	image(i,j)=det(1,1)*power1(i,j)+det(1,2)*power2(i,j)+det(1,3)*power3(i,j)+det(1,4)*power4(i,j)+det(1,5)*power5(i,j)
end do
end do
!call draw(abs(image),ngrid)
!pause

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multitransform5rc(image1,image2,image3,image4,image5,image,ngrid)
implicit none
  integer :: i,j
  real:: ii,jj
  integer,intent(in) :: ngrid
  real,intent(in) :: image1(ngrid,ngrid)
  complex*16 :: power1(ngrid,ngrid)
  real,intent(in) :: image2(ngrid,ngrid)
  complex*16 :: power2(ngrid,ngrid)
  real,intent(in) :: image3(ngrid,ngrid)
  complex*16 :: power3(ngrid,ngrid)
  real,intent(in) :: image4(ngrid,ngrid)
  complex*16 :: power4(ngrid,ngrid)
  real,intent(in) :: image5(ngrid,ngrid)
  complex*16 :: power5(ngrid,ngrid)
  complex*16,intent(out) :: image(ngrid,ngrid)
  real wsave(4*ngrid+15)
  real,parameter :: pi=acos(-1.)
  complex*16 :: ef1,ef2,ef3,ef4,ef5
  real :: kk1,kk2,kk3,kk4,kk5
  real :: a1,b1,a2,b2,a3,b3,a4,b4,a5,b5
  common/shift/ a1,b1,a2,b2,a3,b3,a4,b4,a5,b5
  


call cffti(ngrid, wsave)
do j=1,ngrid
do i=1,ngrid
	power1(i,j)=cmplx(image1(i,j),0.)
	power2(i,j)=cmplx(image2(i,j),0.)
	power3(i,j)=cmplx(image3(i,j),0.)
	power4(i,j)=cmplx(image4(i,j),0.)
	power5(i,j)=cmplx(image5(i,j),0.)
	power1(i,j)=power1(i,j)*(-1)**(i-1+j-1)
    power2(i,j)=power2(i,j)*(-1)**(i-1+j-1)
    power3(i,j)=power3(i,j)*(-1)**(i-1+j-1)
    power4(i,j)=power4(i,j)*(-1)**(i-1+j-1) 
    power5(i,j)=power5(i,j)*(-1)**(i-1+j-1) 
end do
end do

call cfftf_2d(ngrid,ngrid,power1,wsave)
call cfftf_2d(ngrid,ngrid,power2,wsave)
call cfftf_2d(ngrid,ngrid,power3,wsave)
call cfftf_2d(ngrid,ngrid,power4,wsave)
call cfftf_2d(ngrid,ngrid,power5,wsave)

do j=1,ngrid
do i=1,ngrid
	ii=i-ngrid/2-1.
	jj=j-ngrid/2-1.
	kk1=ii*a1+jj*b1
	kk2=ii*a2+jj*b2
	kk3=ii*a3+jj*b3
	kk4=ii*a4+jj*b4
	kk5=ii*a5+jj*b5
	ef1=cmplx(cos(-2.*pi*kk1/ngrid),sin(-2.*pi*kk1/ngrid))
	ef2=cmplx(cos(-2.*pi*kk2/ngrid),sin(-2.*pi*kk2/ngrid))
	ef3=cmplx(cos(-2.*pi*kk3/ngrid),sin(-2.*pi*kk3/ngrid))
	ef4=cmplx(cos(-2.*pi*kk4/ngrid),sin(-2.*pi*kk4/ngrid))
	ef5=cmplx(cos(-2.*pi*kk5/ngrid),sin(-2.*pi*kk5/ngrid))
	power1(i,j)=power1(i,j)*ef1
	power2(i,j)=power2(i,j)*ef2
	power3(i,j)=power3(i,j)*ef3
	power4(i,j)=power4(i,j)*ef4
	power5(i,j)=power5(i,j)*ef5
end do
end do

do j=1,ngrid
do i=1,ngrid
	image(i,j)=power1(i,j)+power2(i,j)+power3(i,j)+power4(i,j)+power5(i,j)
end do
end do
!call draw(abs(image),ngrid)
!pause

return
end subroutine


