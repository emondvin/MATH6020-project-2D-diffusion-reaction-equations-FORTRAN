subroutine solveLSDIAG(n,ndiag,ioff,M,dx,b,nit,err1,err2,stopcrit)
!--------------------------------------------------------------------------------------------------
! Vincent Emond for MATH6020
!--------------------------------------------------------------------------------------------------
! uses the BiCGSTAB method to solve Ax=b in diagonal format
!--------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------
implicit none
 integer, intent(in)               :: n,ndiag
 real,dimension(n,ndiag),intent(in):: m
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: b
 real,dimension(n),intent(inout)   :: dx
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit

real,dimension(n) :: r,z,p,q,rtil,s,t,v 
 real :: rho, beta,rhoold,alfa,omega

 integer           :: maxit,j,i
 real              :: tol1,tol2,normx,normA,normb,norms

 !create two textfiles to write to. These will contain convergence conditions (x, r) with iteration number
open(10, file='Br.txt')
open(11, file='Bx.txt')

 ! initialization: set tolerances, max number of iterations
 ! and rough estimates of matrix and rhs norms for stopping criteria

 maxit=nit; tol1=err1; tol2=err2
 normA=maxval(abs(m)); normb=maxval(abs(b))
 norms=maxval(abs(s))
 rho=0
r=b-r

 call amuxd(n,dx,r,m,ndiag,ioff)
 r=b-r
 rtil=r
 nit=0; stopcrit=0
 do while(stopcrit==0)
   nit=nit+1; rhoold=rho
   normx=sqrt(sum(dx*dx))/n


   rho=dot_product(rtil,r)
   if (rho==0) then
	exit
   endif
   if (nit==1) then
      p=r
   else
      beta=(rho/rhoold)*(alfa/omega)
      p=r+beta*(p-omega*v)
   endif
p=p
call amuxd(n,p,v,m,ndiag,ioff)
alfa=rho/(dot_product(rtil,v))
s=r-alfa*v

   call amuxd(n,s,t,m,ndiag,ioff)
   omega=dot_product(t,s)/dot_product(t,t)
   dx=dx+alfa*p+omega*s
   r=s-omega*t

  !check norm of s
  if (sqrt(sum(s*s))<tol1) then
	exit
  end if

 ! test for convergence
 !---------------------
   err1=sqrt(sum(r*r))/n
   err2=sqrt(sum((alfa*p+omega*s)*(alfa*p+omega*s)))/n
   if (nit>maxit) stopcrit=-1
   if (err1<tol1*(norma*normx+normb)) stopcrit=stopcrit-10
   if (err2<tol2*normx) stopcrit=stopcrit-100

   ! uncomment the next line to monitor convergence progress
   !write(*,'(I6,4(E14.7,X))') nit,err1,err2,maxval(dx),minval(dx)
!Write the results of convergence conditions to two textfiles defined at the beginning of the program
   write(10,*) nit, err1
   write(11,*) nit, err2
 enddo

close(10)
close(11)
end subroutine


