program project_pt3
implicit none

!----------------------------------------------------------------------------------------------------
! Author: Vincent Emond
! For use in part 3 of the project outlined in MATH*6020
! A program which uses the solution from part 1 with some modifications in order/
!  /to calculate a time dependent solution to the problem outlined in the project.
!
! Compile using the following command and flag: gfortran -fdefault-real-8 "thisfile.f90/
!   /solver.f90 output.f90"
!  where thisfile.f90 is this program. solver.90 is the chosen solver/
!  /ie. BiCGSTAB. And output.f90 is a subroutine which creates multiple/
!  / plottable data files which can be stitched together into a gif animation.
!----------------------------------------------------------------------------------------------------



INTERFACE
subroutine solveLSDIAG(n,ndiag,ioff,M,dx,b,nit,err1,err2,stopcrit)
!---------------------------------------------------------------------------------------------------
! input:  n       problem size
!         ndiag:  number of diagonals
!         ioff:   offsets (distance of sub diagonals to main diagonal)
!         M:      matrix values
!         dx:     initial guess for iteration (will be overwritten with result)
!         b:      the righ hand side of the linear system
!         nit:    maximum number of iterations to be carried out (will be overwritten)
!         err1:   tolerance for 1st stopping criterion (will be overwritten)
!         err2:   tolerance for 2nd stopping criterion (will be overwritten)
! output: sol:    solution of Ax=b
!         nit:    number of iterations taken
!         err1:   computed value for 1st stopping criterion
!         err2:   computed value for 2nd stopping criterion
!         stopcrit:  integer value that contains information which stopping criterions became active
!----------------------------------------------------------------------------------------------------

implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: m
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: b
 real,dimension(n),intent(inout)   :: dx
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit
end subroutine
END INTERFACE

!---------------------------------------------------------------------------------------------------
!declarations for the program begins
! n: Number of points spanning the area of the gxg grid
! g: Number of grid points along each axis
! A: Matrix from part 1
! M: Modification to A matrix
! E: The matrix used in solving the problem outlined in the project
! v: capital phi in equation
! x: Ax=b, its the x
! dx: change in x with each iteraiton
! xold: used in updating the solution with each iteration
! rhs: right hand side of the equation Ax=b
! k,k1,k2: the values of k parameters in the system being solved.
! t: time
! nit: number of iterations in the solver subroutine
! nitr: number of iterations for which the do loop runs. times the solver is used.
! ctr: number of times an output file is created
!----------------------------------------------------------------------------------------------------

 integer :: n,ndiag,g
 parameter (g=50, n=g*g, ndiag=5)
 real,dimension(n,ndiag) ::A,M,E
 integer,dimension(ndiag) :: ioff
 real,dimension(n) :: rhs,sol,res,y,x,defect,v,vd,dx,xold,f,df
 integer :: nit,stopcrit,i,nitr,ctr
 real :: err1,err2,err3,flux,t,dt,k1,k2,h,tol,k

 integer :: tstart,tfinish,clock_rate
 real :: elapsed_time

 ! initialise: set tolerances, max no iterations
 ! and initial guess
 !----------------------------------------------
 err1=1e-10;  err2=1e-10;  nit=n*100

 k=1000.
 k1=100.
 k2=1
 h=1/real(g)

 dt=h/1000.0
 t=0.0

 dx=1.
 x=1.
 xold=1.

 v=1.

 tol=1e-12
 err3=1

call system_clock(tstart)

 !(1) setup a test case: prescribe solution sol
 ! and generate a matrix A and rhs, such that A*sol=rhs
 !-----------------------------------------------------

!set matrix
 call genDIAG(n,ndiag,ioff,A,rhs,sol,g)


rhs=rhs/h**2
A=A/h**2


ctr=0
nitr=0

do while (stopcrit==0)

!set matrix E
E=A
E(:,3)=E(:,3)-((k1*k2)/(k2+x)**2)

M=0
M(:,3)=1.
M=M-(dt/2)*E

!set v (capital phi)
call amuxd(n,x,v,A,ndiag,ioff)
v=dt*(v-(k1*x)/(k2+x)-rhs)


!call solver
err1=1e-10
err2=1e-10

call solveLSDIAG(n,ndiag,ioff,M,dx,v,nit,err1,err2,stopcrit)

!update soln
xold=x
x=x+dx
t=t+dt
stopcrit=0




if (mod(nitr,100)==0) then
 call output(g,n,ctr,x,t)
 ctr=ctr+1
! write(*,*)norm2(dx)
endif
nitr=nitr+1



 ! print iteration
if(mod(nitr,250)==0) write(*,*) "iteration = " , nitr, "  with norm(dx) = " , norm2(dx), " @ t = ", t

if (norm2(dx)<tol) stopcrit=-10
enddo

call system_clock(tfinish, clock_rate)
elapsed_time = float(tfinish-tstart) / float(clock_rate)

write(*,*)ctr
write(*,*)nitr
write(*,*)elapsed_time



end program



subroutine genDIAG(n,ndiag,ioff,A,rhs,x,g)
implicit none
!----------------------------------------------------------
! A subroutine which creates a matrix A based on the outlined problem and its boundary conditions.
! Author: Vincent Emond
! For MATH*6020
!-----------------------------------------------------------

 integer, intent(in) :: n,ndiag,g
 integer,intent(out),dimension(ndiag) :: ioff
 real, dimension(n),intent(out) :: rhs 
 real, dimension(n,ndiag),intent(out) :: a
 real, dimension(n),intent(out) :: x
 integer :: i
 real :: h,k
 real,dimension(n) :: p, E, W, Nt, S

!define grid spacing h, 1/g where g is number of grid points
!define k to be used in calculations

   h=1.0/real(g)
   k=0 


!generate rhs vector b
!initialize values to 0
   do i=1,n
    rhs(i)=0
   enddo

!calcuate values of b for different boundary conditions
    do i=1,n
        if (mod(i,g)==0) then
            rhs(i)=-2
        endif
        if (i>n-g) then
           rhs(i)=-2
        endif
        if (i==n) then
           rhs(i)=-4
        endif
    enddo

!main diagonal
!initialize for no boundary conditions
    do i=1,n
     p(i)=-4-k*h**2
    enddo

!compute main diag based on grid position

    do i=1,n
    
    !top row in grid
     if (mod(i,g)==0) then
      p(i)=p(i)-1
     endif

    !left boundary in grid
     if (i.LE.g) then
      p(i)=p(i)+1
     endif

    !right boundary in grid
     if(i.GT.n-g) then
      p(i)=p(i)-1
     endif

    !bottom boundary in grid
     if(mod(i,g)==1) then
      p(i)=p(i)+1
     endif
    enddo

!Off diagonals
!initialize North and South = 0 and East and West = 1
   do i=1,n
     E(i)=1
     W(i)=1
     Nt(i)=0
     S(i)=0
   enddo

!compute N,S, E and W are taken care of by ioff
    do i=1,n
      if (mod(i,g).NE.0) then
       Nt(i)=1
       S(i+1)=1
      endif
    enddo

!Set up matrix
!off diagonal locations
    ioff=(/ -g,-1,0,1,g /)

!assign diagonals
    a(:,1) = W;
    a(:,2) = S;
    a(:,3) = p;
    a(:,4) = Nt;
    a(:,5) = E

end subroutine




subroutine amuxd (n,x,y,M,ndiag,ioff)
!-----------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA) 
!        f90/f95 version of the sparskit f77 subroutine
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
! ioff   = integer array of length idiag, containing the offsets of the
!          diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=A*x
!
!-----------------------------------------------------------------------
implicit none
  integer, intent(in)::  n, ndiag
  integer, intent(in),dimension(ndiag) :: ioff
  real, dimension(n), intent(in) :: x
  real, dimension(n,ndiag), intent(in) :: M
  real, dimension(n), intent(out) :: y
  integer :: j, io, i1, i2, i

      y=0.

      do j=1, ndiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do i=i1,i2
           y(i) = y(i)+m(i,j)*x(i+io)
         enddo
      enddo
end subroutine amuxd

