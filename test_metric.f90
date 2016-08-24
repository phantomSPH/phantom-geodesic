module testmetric
   implicit none
contains
!----------------------------------------------------------------
!+
!  Subroutine to test the metric for different values of x and v
!  It should be used in the test suite.
!+
!----------------------------------------------------------------
   subroutine test_metric(ntests,npass)
      use metric, only: metric_type
      use utils_gr, only: get_metric3plus1
      integer, intent(inout) :: ntests,npass
      real :: x(1:3), v(1:3),alpha,beta(1:3),gij(1:3,1:3),bigv(1:3)
      integer :: i,ierr
      ierr = 0

      write(*,'(/,a)') '--> testing metric'
      write(*,'(a,/)') '    metric type = '//trim(metric_type)

      do i =0,16
         x = (/2.0+1*10**(-real(i)),0.,0./)
         bigv = (/i*0.,0.,0./)
         call get_metric3plus1(x,alpha,beta,gij)
         v = bigv*alpha-beta
         call test_metric_i(x,v,ntests,npass)
      enddo

      write(*,'(/,a,/)') '<-- metric test complete'

   end subroutine test_metric


!----------------------------------------------------------------
!+
!  Subroutine containing the actual tests of the metric
!+
!----------------------------------------------------------------
   subroutine test_metric_i(x,v,ntests,npass)
      use testutils, only: checkvalbuf
      use metric, only: get_metric
      use utils_gr, only: dot_product_gr
      integer, intent(inout) :: ntests,npass
      real, intent(in) :: x(1:3), v(1:3)
      real, dimension(0:3,0:3) :: gcov,gcon,gg
      real :: sqrtg, v4(0:3), sum
      real, parameter :: tol=1.e-15
      integer :: i,j, nerrors,ncheck,n_error
      real :: errmax

      ntests=ntests+1
      nerrors = 0

      call get_metric(x,gcov,gcon,sqrtg)
      gg = 0.
      gg = matmul(gcov,gcon)
      sum = 0
      do j=0,3
         do i=0,3
            sum = sum + gg(i,j)
         enddo
      enddo

      n_error = 0
      call checkvalbuf(sum,4.,tol,'[F]: gddgUU ',n_error,ncheck,errmax)
      nerrors = nerrors+n_error
      if (n_error>0) then
         print*, 'gdown*gup /= Identity'
         do i=0,3
            write(*,*) gg(i,:)
         enddo
      endif

      v4=(/1.,v/)
      call checkvalbuf(int(sign(1.,dot_product_gr(v4,v4,gcov))),-1,0,'[F]: &
                           &sign of dot_product_gr(v4,v4,gcov))',nerrors,ncheck)
      if (dot_product_gr(v4,v4,gcov) > 0.) then
         nerrors = nerrors + 1
         print*,'Warning: Bad combination of position and velocity... &
                &dot_product_gr(v4,v4,gcov)=',dot_product_gr(v4,v4,gcov),' > 0'
      endif

      if (nerrors/=0) then
         print*,'-- Metric test failed with:'
         print*,'    x =',x
         print*,'    v =',v
         print*,''
      else
         npass  = npass + 1
      endif

   end subroutine test_metric_i

end module testmetric
