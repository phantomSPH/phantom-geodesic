program test
 use init, only: initialise
 use metric, only: get_sourceterms, get_metric, mass1
 use step, only: step_leapfrog, step_1, step_heuns
 use utils_gr, only: dot_product_gr
 implicit none
 !integer, parameter :: ndim = 3
 real, dimension(3):: x,v,fterm
 real, dimension(0:3,0:3) :: gcov, gcon
 real :: dt, sqrtg, v4(0:3), time, en_init
 integer :: nsteps, i
 !real(8) :: energy
 !real(8) :: angmom(ndim)
 nsteps = 2*1000000

 !call sanity_checks
 !stop "sanity checks"

 print*,'START'

 call initialise(x, v)
 !Check 'causality condition thingamajig'
 call get_metric(x,gcov,gcon,sqrtg)
 v4(0) = 1.
 v4(1:3) = v(1:3)
 if (dot_product_gr(v4,v4,gcov)>0.) then
  STOP "Bad initial conditions!"
 endif

 write(1,*) x
 write(2,*) v
 time = 0.
 en_init = (1. - 2*mass1/sqrt(dot_product(x,x)))/sqrt(-dot_product_gr(v4,v4,gcov))
 write(3,*) time, (1. - 2*mass1/sqrt(dot_product(x,x)))/sqrt(-dot_product_gr(v4,v4,gcov)) - en_init


 dt = 1.e-3
 do i =1,nsteps
  time = time + dt
  call constraints(x,v)
  call get_sourceterms(x,v,fterm)
  !call step_leapfrog(x,v,fterm,dt)
  call step_heuns(x,v,fterm,dt)
  !call step_1(x,v,fterm,dt)
  if (mod(i,1000)==0) then
   print*,i
   write(1,*) x
   write(2,*) v
   v4(1:3) = v(1:3)
   write(3,*) time,(1. - 2*mass1/sqrt(dot_product(x,x)))/sqrt(-dot_product_gr(v4,v4,gcov)) - en_init
  endif

 enddo


end program test

subroutine sanity_checks
 use cons2prim, only: get_p_from_v, get_v_from_p
 use metric, only: get_metric
 use utils_gr, only: dot_product_gr
 real, dimension(3) :: x, v, pmom
 real, dimension(0:3,0:3):: gcov,gcon,gg,identity
 real :: sqrtg, v4(0:3), U0, umu(0:3), r,vx
 print*,"SANITY CHECKS:"
 !print*,"Initial position and velocity:"
 !r = 2.0+1.e-1
 !print*,"Enter r and v:"
 !read*,r,vx
 r = 3
 x = (/r,0.,0./)
 !x = (/2.9,0.,0./)
 vx= 0.0
 v = (/0.0,0.7,0.0/)
 print*,"r: ",r
 print*,'x:',x
 print*,"v:",v
 print*,""

 print*,'Covariant and contravariant metrics'
 call get_metric(x,gcov,gcon,sqrtg)

 print*,"gcov"
 do i=0,3
  print*,gcov(i,:)
 enddo

 print*,"gcon"
 do i=0,3
  print*,gcon(i,:)
 enddo


 gg = 0.
 gg = matmul(gcov,gcon)
 print*,'[gcov][gcon] (matrix multiplication)'
 do i=0,3
  print*,gg(i,:)
 enddo


 print*,''
 print*,'Four velocity stuff:'
 v4(0) = 1.
 v4(1:3) = v(:)
 identity=0
 do i=0,3
  do j=0,3
   if (i==j) identity(i,j)=1.
  enddo
 enddo
 identity(0,0)=-1.

 print*,"sqrt(gijvivj): ",sqrt(dot_product_gr(v(1:3),v(1:3),gcov(1:3,1:3)))
 print*,"dot_product_gr(v4,v4,gcov): ",dot_product_gr(v4,v4,gcov)
 U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
 print*,'U0: ',U0
 umu= U0*v4
 print*,'umu^2: ',dot_product_gr(umu,umu,gcov)


 print*,""
 print*,"Testing cons2prim"
 print*,'x:',x
 print*,"v:",v
 call get_p_from_v(pmom,v,x)
 print*,"call get_p_from_v(pmom,v,x)"
 print*,'pmom',pmom
 call get_v_from_p(pmom,v,x)
 print*,"call get_v_from_p(pmom,v,x)"
 print*,'x:',x
 print*,"v:",v

end subroutine

subroutine constraints(x,v)
 use metric, only: get_metric
 use utils_gr, only: dot_product_gr
 use cons2prim, only: get_p_from_v,get_v_from_p
 real, dimension(1:3), intent(in) :: x,v
 real, dimension(0:3) :: v4
 real, dimension(0:3,0:3) :: gcov, gcon, gg
 real :: sqrtg, trace, trace_err, offdiag_err, v_from_p(1:3), pmom(1:3)
 integer :: i,j,k

 v4(0) = 1.
 v4(1:3) = v(:)
 call get_metric(x,gcov,gcon,sqrtg)

 if (dot_product_gr(v4,v4,gcov)>0.) print*,"WARNING Causality check: dot_product_gr(v4,v4,gcov)>0"
 gg = 0.
 gg = matmul(gcov,gcon)
 trace = 0.
 do i=0,3
  trace = trace + gg(i,i)
 enddo

 trace_err = 1.e-13
 offdiag_err = 1.e-15
 if (trace-4.>trace_err) then
  print*,""
  print*, "WARNING: trace(gdown*gup)>",trace_err,"trace = ", trace
  print*, "gdown*gup /= Identity"
  do i=0,3
   print*,gg(i,:)
  enddo
 endif
 do i=0,3
  do j=0,3
   if (i/=j) then
    if (gg(i,j)>offdiag_err) then
     print*,""
     print*,"WARNING: off diagonals of gdown*gup >",offdiag_err
     print*, "gdown*gup /= Identity"
     do k=0,3
      print*,gg(k,:)
     enddo
    endif
   endif
  enddo
 enddo

 vel_tol = 1.e-15
 call get_p_from_v(pmom,v,x)
 call get_v_from_p(pmom,v_from_p,x)
 if (maxval(abs(v-v_from_p))>vel_tol) then
  print*,""
  print*,"WARNING: get_p_from_v -> get_v_from_p didn't work: |v_bef-v_aft|>",vel_tol
  print*,"v_bef-v_aft =",v-v_from_p
 endif
end subroutine constraints
