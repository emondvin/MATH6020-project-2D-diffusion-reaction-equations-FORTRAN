program project_part_2
!----------------------------------------------------------------------------------------
!Semilinear elliptic solver of the 2d diffusion reaction problem
!Compile using gfortran -fdefault-real-8 "thisfile.f90"
!
!Author: Vincent Emond
!----------------------------------------------------------------------------------------
implicit none

!delcarations for the program
integer :: i,stopcrit,nit
integer, parameter :: g=100, n=g*g
 !number of points along x (or y) axis of grid = g, n= number of grid points

real :: k1,k2,err1,flux,diff
real, dimension(n) :: a,b,c,x,xold
 !grid spacing h=1/g=1/number of grid points along axis
real,parameter :: h=1.0/real(g)

!used to determine how long the program runs
integer :: tstart,tfinish,clock_rate
real :: elapsed_time

!set values of k1,k2
k1=100
k2=1

!set error tolerance
err1=1e-12

!initialize number of iterations to 0, will be updated
nit=0

!initialize the grid
x=0.

stopcrit=0
call system_clock(tstart)

open (unit=31,file="P2_sol_conv.dat")
open (unit=20,file="update_p2.dat")

do while (stopcrit==0)
 !update number of iterations
 nit = nit+1

 !use this to later compare k-1 and kth iterations
 xold = x

!Here we calculate the values of the grid points which depend on boundary conditions
!Assign values of a,b,c which are used with quadratic equation

 do i=1,n

  !inner points
   if (i>g .and. mod(i,g)/=0 .and. mod(i,g)/=1 .and. i<n-g) then
    a(i)=-4.
    b(i)=x(i+1)+x(i-1)+x(i+g)+x(i-g)-4*k2-h*h*k1
    c(i)=(x(i+1)+x(i-1)+x(i+g)+x(i-g))*k2
    x(i)=(-b(i)-sqrt(b(i)**2-4*a(i)*c(i)))/(2*a(i))
   end if
  
  !i=1, bottom left bound
    a(1)= -6
    b(1)=x(2)+x(1+g)-2*k2-h*h*k1
    c(1)=(x(2)+x(1+g))*k2
    x(1)=(-b(1)-sqrt(b(1)**2-4*a(1)*c(1)))/(2*a(1))

  !top left bound i=g
    a(g)= -4
    b(g)=2+x(g-1)+x(g+g)-4*k2-h*h*k1
    c(g)=(x(g+g)+x(g-1)+2)*k2
    x(g)=(-b(g)-sqrt(b(g)**2-4*a(g)*c(g)))/(2*a(g))

  !bottom right i=n-g+1
    a(n-g+1)= -4
    b(n-g+1)=2+x(n-g+2)+x(n-g+1-g)-4*k2-h*h*k1
    c(n-g+1)=(x(n-g+2)+x(n-g+1-g)+2)*k2
    x(n-g+1)=(-b(n-g+1)-sqrt(b(n-g+1)**2-4*a(n-g+1)*c(n-g+1)))/(2*a(n-g+1))

  !top right i=n
    a(n)= -6
    b(n)=4+x(n-1)+x(n-g)-6*k2-h*h*k1
    c(n)=(x(n-1)+x(n-g)+4)*k2
    x(n)=(-b(n)-sqrt(b(n)**2-4*a(n)*c(n)))/(2*a(n))

  !top row B.C.
   if (mod(i,g)==0 .and. i/=g .and. i/=n) then
    a(i)=-5
    b(i)=2+x(i-1)+x(i+g)+x(i-g)-5*k2-h*h*k1
    c(i)=(x(i-1)+x(i+g)+x(i-g)+2)*k2
    x(i)=(-b(i)-sqrt(b(i)**2-4*a(i)*c(i)))/(2*a(i))
   endif

  !bottom row B.C.
   if (mod(i,g)==1 .and. i/=n-g+1 .and. i/=1) then
    a(i)=-3
    b(i)=x(i+1)+x(i+g)+x(i-g)-3*k2-h*h*k1
    c(i)=(x(i+1)+x(i+g)+x(i-g))*k2
    x(i)=(-b(i)-sqrt(b(i)**2-4*a(i)*c(i)))/(2*a(i))
   endif

  !left boundary
   if (i<g .and. i/=1) then
    a(i)=-3
    b(i)=x(i+1)+x(i+g)+x(i-1)-3*k2-h*h*k1
    c(i)=(x(i+1)+x(i+g)+x(i-1))*k2
    x(i)=(-b(i)-sqrt(b(i)**2-4*a(i)*c(i)))/(2*a(i))
   endif

  !right boundary
   if (i>n-g+1 .and. i/=n) then
    a(i)= -5
    b(i)=2+x(i-1)+x(i+1)+x(i-g)-5*k2-h*h*k1
    c(i)=(x(i+1)+x(i-1)+x(i-g)+2)*k2
    x(i)=(-b(i)-sqrt(b(i)**2-4*a(i)*c(i)))/(2*a(i))
   endif
  end do

write(31,*) nit,x(i)
write(20,*) nit,diff


 if (norm2(xold-x)<err1*norm2(x)) stopcrit=-10

 diff = norm2(xold-x)

 !uncomment next line to see convergence in terminal
 !write(*,'(A30,E14.7)') 'error: ',diff

end do

close(31)
close(20)

call system_clock(tfinish, clock_rate)
elapsed_time = float(tfinish-tstart) / float(clock_rate)

write(*,'(A30,E14.7)') 'CPU time in seconds: ',elapsed_time

write(*,'(A30,I6,X)') 'number of iterations: ',nit

!calculate and write flux
do i=1,n
  if (mod(i,g)==1 .or. i>n-g) then
   flux=flux+4*(1-x(i))
  end if
enddo

 write(*,'(A30,E14.7)') 'flux: ',flux

  !Write flux to data file for further use 
    open (unit=30,file="P2flux.dat",status='unknown',access='append')
    write(30,*) n, flux
    close(30)


call output(n,g,x)


! Write minimum concentration value to data file for further use
    open (unit=35,file="P2min1.dat",status='unknown',access='append')
    write(35,*) n, x(1)
    close(35)

end program

!subroutine converts the single array of dimension n to an x,y grid
subroutine output(n,g,x)

implicit none

integer,intent(in) :: n,g
 real,dimension(n),intent(in)           :: x
 real,dimension(n)                      :: w, y
 real                                   :: h
 integer                                :: i,j

h=1.0/g

!compute x coords of grid points x
do j=1,g
 do i=(j-1)*g+ 1,j*g
  w(i)=h/2.+(j-1)*h
 enddo
enddo

!compute y coords
do i =1,n
 y(i)=h/2.+h*(mod(i-1,g))
enddo

!write coords and solution vector to dat file
open(unit=30,file="tk2.dat")

do i=1,n
 write(30,*) w(i), y(i), x(i)
 if (mod(i,g)==0) then
    write(30,*)
 endif
enddo

close(30)

endsubroutine
