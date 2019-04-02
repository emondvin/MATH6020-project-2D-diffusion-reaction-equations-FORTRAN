subroutine output(g,n,ctr,dx,t)

implicit none
  integer,intent(in) :: ctr
  integer,intent(in) :: g,n
  real,dimension(n),intent(in) :: dx
  real,intent(in) :: t
  integer ::i,j,p
  character(len=20) :: str
  character(len=9) :: suf
  external str

! call suffix(ctr, suf)
 open(13,file='proj.'//trim(str(ctr)),action='write')
 write(13,*) '# ',t
 do i=1,g
 do j=1,g
   p=(i-1)*g+j
   write(13, '(E14.4)') dx(p)
 enddo
  write(13,*)
  enddo
  close(13)

end subroutine output

character(len=20) function str(kk)
integer,intent(in) :: kk
write(str,*) kk
str=adjustl(str)
end function str
