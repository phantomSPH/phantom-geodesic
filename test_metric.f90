module testmetric
   implicit none
contains
   
   subroutine test_metric(ntests,npass)
      use metric, only: metric_type
      use utils_gr, only: get_metric3plus1
      integer, intent(inout) :: ntests,npass
      real :: x(1:3), v(1:3),alpha,beta(1:3),gij(1:3,1:3),bigv(1:3)
      integer :: i,ierr
      ierr = 0
      write(*,"(/,a)")  '--> testing metric'
      write(*,"(2a,/)") '    Metric type = ', metric_type
      x = (/0.,0.,0./)
      call get_metric3plus1(x,alpha,beta,gij)
      !  v = (/0.5,0.,0./)
      do i =0,0
         bigv = (/i*0.1,0.,0./)
         v = bigv*alpha-beta
         call metric_test(x,v,ierr)
         ntests=ntests+1
         if (ierr/=0) then
            write(*,"(/,a)") " Metric test failed",i
         else
            npass  = npass + 1
         endif
      enddo 
      write(*,"(/,a,/)") '<-- metric test complete'
   end subroutine test_metric
   
   subroutine metric_test(x,v,ierr)
      use testutils, only: checkvalbuf
      use metric, only: get_metric
      use utils_gr, only: dot_product_gr
      integer, intent(out) :: ierr
      real, intent(in) :: x(1:3), v(1:3)
      real, dimension(0:3,0:3) :: gcov,gcon,gg
      real :: sqrtg, trace, v4(0:3)
      real, parameter :: trace_err = 1.e-13, offdiag_err = 1.e-15
      integer :: i,j, n_errors
      
      n_errors = 0
      
      call get_metric(x,gcov,gcon,sqrtg)
      gg = 0.
      gg = matmul(gcov,gcon)
      trace = 0.
      do i=0,3
         trace = trace + gg(i,i)
      enddo
      ! call checkvalbuf(trace,4.,trace_err,n_error,'trace = 4')
      if (trace-4.>trace_err) then
         n_errors = n_errors + 1
         write(*,"(/,a)") " WARNING: trace(gdown*gup)>",trace_err,"trace = ", trace
      endif
      
      do i=0,3
         do j=0,3
            if (i/=j) then
               if (gg(i,j)>offdiag_err) then
                  n_errors = n_errors + 1
                  write(*,"(/,a)") " WARNING: off diagonals of gdown*gup >",offdiag_err
               endif
            endif
         enddo
      enddo
      
      if (n_errors>0) then
         print*, "gdown*gup /= Identity"
         do i=0,3
            print*,gg(i,:)
         enddo
      endif
      
      v4=(/1.,v/)
      if (dot_product_gr(v4,v4,gcov)>0.) then
         n_errors = n_errors + 1
         write(*,"(/,a)") " WARNING: Bad combination of position and velocity."
         print*,"dot_product_gr(v4,v4,gcov)=,",dot_product_gr(v4,v4,gcov)," > 0"
         print*,"x =",x
         print*,"v =",v
      endif
      
      ierr = n_errors
   end subroutine metric_test
   
end module testmetric