subroutine solveLSDIAG(n,ndiag,ioff,A,sol,rhs,nit,err1,err2,stopcrit)
!--------------------------------------------------------------------------------------------------
! 
!--------------------------------------------------------------------------------------------------
! uses the conjugate gradient method to solve Ax=b in diagonal format
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
 real,dimension(n,ndiag),intent(in):: a
 integer, dimension(ndiag),intent(in)::ioff
 real,dimension(n),intent(in)      :: rhs
 real,dimension(n),intent(inout)   :: sol
 real,intent(inout)                :: err1,err2
 integer,intent(inout)             :: nit
 integer,intent(out)               :: stopcrit

 real,dimension(n) :: r,z,p,q
 real :: rho, beta,rhoold,alfa 

 integer           :: maxit,j,i
 real              :: tol1,tol2,normx,normA,normb
 open(10, file='CGr.txt')
 open(11, file='CGx.txt')
 ! initialization: set tolerances, max number of iterations
 ! and rough estimates of matrix and rhs norms for stopping criteria
 maxit=nit; tol1=err1; tol2=err2
 normA=maxval(abs(a)); normb=maxval(abs(rhs))


 rho=0.

 call amuxd(n,sol,r,a,ndiag,ioff)
 r=rhs-r

 nit=0; stopcrit=0
 do while(stopcrit==0)
   nit=nit+1; rhoold=rho
   normx=sqrt(sum(sol*sol))/n

   
   z=r
   rho=dot_product(r,z)
   if (nit==1) then
      p=z
   else
      beta=rho/rhoold
      p=z+beta*p
   endif
   call amuxd(n,p,q,a,ndiag,ioff)
   alfa=rho/dot_product(p,q)
   sol=sol+alfa*p
   r=r-alfa*q

     
   ! test for convergence
   !---------------------
   err1=sqrt(sum(r*r))/n
   err2=abs(alfa)*sqrt(sum(p*p))/n
   if (nit>maxit) stopcrit=-1
   if (err1<tol1*(norma*normx+normb)) stopcrit=stopcrit-10
   if (err2<tol2*normx) stopcrit=stopcrit-100

   ! uncomment the next line to monitor convergence progress
   write(*,'(I6,4(E14.7,X))') nit,err1,err2,maxval(sol),minval(sol)
   write(10, *) nit, err1
   write(11,*) nit, err2
 enddo

close(10)
close(11)
end subroutine




