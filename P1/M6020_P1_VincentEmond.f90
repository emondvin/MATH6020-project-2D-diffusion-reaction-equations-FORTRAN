program project_one
!--------------------------------------------------------------------------------------------------
!A testbed for 2D diffusion reaction problems with no time dependence which can be
! computed using different iterative solvers.
!For MATH6020
!Compile using gfortran -fdefault-real-8 "thisfile.f90" "solver.f90" 
!
!Author: Vincent Emond
!--------------------------------------------------------------------------------------------------

implicit none

INTERFACE
 subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit,g)
 !---------------------------------------------------------------------------------------------------
 ! input:  n       problem size
 !         ndiag:  number of diagonals
 !         ioff:   offsets (distance of sub diagonals to main diagonal)
 !         A:      matrix values
 !         sol:    initial guess for iteration (will be overwritten with result)
 !         rhs:    the righ hand side of the linear system
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

  integer, intent(in)               :: n,ndiag,g
  real,dimension(n,ndiag),intent(in):: a
  integer, dimension(ndiag),intent(in)::ioff
  real,dimension(n),intent(in)      :: rhs
  real,dimension(n),intent(inout)   :: sol
  real,intent(inout)                :: err1,err2
  integer,intent(inout)             :: nit
  integer,intent(out)               :: stopcrit
 
 end subroutine
END INTERFACE

!---------------------------------------------------------
!declarations for the actual program begin program begin
!---------------------------------------------------------

 integer :: n,ndiag,g
 parameter (g=100, n=g*g, ndiag=5)
 real,dimension(n,ndiag) ::A
 integer,dimension(ndiag) :: ioff
 real,dimension(n) :: rhs,sol,res,y,x,defect
 integer :: nit,stopcrit,i,j
 real :: err1,err2,flux

 integer :: tstart,tfinish,clock_rate
 real :: elapsed_time
 
 !g is the number of grid points along each axis and n is the total number of grid points


 ! initialise: set tolerances, max number of iterations
 ! and initial guess
 !----------------------------------------------
 err1=1e-12;  err2=1e-12;  nit=n*100;  x=0.

 !(1) set up a test case: prescribe solution sol
 ! and generate a matrix A and rhs, such that A*sol=rhs
 !-----------------------------------------------------
 call genDIAG(n,ndiag,ioff,A,rhs,sol,g)


 !(2) call linear solver
 !-----------------------------------------------------------
 call system_clock(tstart)
 call solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit,g)
 call system_clock(tfinish, clock_rate)
 elapsed_time = float(tfinish-tstart) / float(clock_rate)

 !(3) report results: 
 !-------------------------------
 call amuxd (n,x,y,a,ndiag,ioff,g)
 res=y-rhs    ! compute residual
 defect=sol-x ! compute defect
 write(*,'(A30,I6,X,E14.7,X,E14.7,XI6)') 'done: it/err1/err2/stopcrit',nit,err1,err2,stopcrit
 write(*,'(A30,2(E14.7,X))')'max defect, max residual: ',maxval(abs(defect)),maxval(abs(res))
 write(*,'(A30,E14.7)') 'CPU time in seconds: ',elapsed_time

!Solution output
!-------------------------------
call output(n,g,sol)

!calculate the flux of the quadrant
!------------------------------------
flux=0

do i=1,n
if (mod(i,g)==1 .or. i>n-g) then
 flux=flux + 4*(1-sol(i))
end if
end do
 write(*,'(A30,E14.7)') 'flux: ',flux


  !open (unit=30,file="P1flux.dat",status='unknown',access='append')
  ! write(30,*) n, flux
  ! close(30)

   !Write minimum concentration
   !open (unit=35,file="P1minconc.dat",status='unknown',access='append')
   !	write(35,*) n, x(1)
   !close(35)


end program

!-------------------
!Begin subroutines
!-------------------
subroutine genDIAG(n,ndiag,ioff,A,rhs,x,g)
implicit none
 integer, intent(in) :: g,n,ndiag
 integer,intent(out),dimension(ndiag) :: ioff
 real, dimension(n),intent(out) :: rhs 
 real, dimension(n,ndiag),intent(out) :: a
 real, dimension(n),intent(out) :: x
 integer :: i
 real :: h,k
 real,dimension(n) :: p, E, W, Nt, S

!define grid spacing h, 1/g where g is number of grid points
!define k to be used in calculations
!---------------------------------------
  h=1.0/real(g)
  k=1
!--------------------------------------


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

subroutine amuxd (n,x,y,diag,idiag,ioff,g)
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
  integer, intent(in)::  n, idiag, g
  integer, intent(in),dimension(idiag) :: ioff
  real, dimension(n), intent(in) :: x
  real, dimension(n,idiag), intent(in) :: diag
  real, dimension(n), intent(out) :: y
  integer :: j, io, i1, i2, i

      y=0.

      do j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do i=i1,i2
           y(i) = y(i)+diag(i,j)*x(i+io)
         enddo
      enddo
end subroutine amuxd



subroutine output(n,g,sol)

implicit none

integer,intent(in) :: n,g
 real,dimension(n),intent(in)           :: sol
 real,dimension(n)                      :: w, y
 real                                   :: h
 integer                                :: i,j

h=1/real(g)

!compute x coords of grid points u
do j=1,g
 do i=(j-1)*g+ 1,j*g
  w(i)=h/2.+(j-1)*h
 enddo
enddo

!compute y coords
do i =1,n
 y(i)=h/2.+h*(mod(real(i)-1.,real(g)))
enddo

!write coords and solution vector to dat file
open(unit=30,file="P1.dat")

do i=1,n
 write(30,*) w(i), y(i), sol(i)
 if (mod(i,g)==0) then
    write(30,*)
 endif
enddo

close(30)

endsubroutine
