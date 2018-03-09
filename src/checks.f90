module checks
 implicit none
contains
 !----------------------------------------------------------------
 !+
 !  Subroutine to check the metric and 4velocity conditions. To be
 !  called in between desired timesteps.
 !+
 !----------------------------------------------------------------
subroutine check(x,v,passed)
 use metric_tools, only: get_metric
 use testmetric, only: test_metric_i
 real, intent(in) :: x(1:3), v(1:3)
 integer :: ntests,npass
 logical :: passed
 npass=0
 call test_metric_i(x,v,ntests,npass)
 passed = (npass==1)
end subroutine check


!THIS IS OBSOLETE NOW. USE test suite instead.
 !----------------------------------------------------------------
 !+
 !  Subroutine to do simple gr checks/tests (to be called on its own)
 !+
 !----------------------------------------------------------------
subroutine sanity_checks
 use cons2prim, only: get_p_from_v, get_v_from_p
 use metric_tools, only: get_metric
 use utils_gr, only: dot_product_gr
 real, dimension(3) :: x, v, pmom
 real, dimension(0:3,0:3):: gcov,gcon,gg,identity
 real :: sqrtg, v4(0:3), U0, umu(0:3), r,vx
 integer :: i,j
 print*,"SANITY CHECKS:"
 !print*,"Initial position and velocity:"
 !r = 2.0+1.e-1
 !print*,"Enter r and v:"
 !read*,r,vx
 r = 25.
 x = (/r,1.5,0./)
 !x = (/2.9,0.,0./)
 vx= 0.0
 v = (/0.,0.,0./)
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
end subroutine sanity_checks

end module checks
