module testcons2prim
   implicit none
contains
   subroutine test_cons2prim(ntests,npass)
      use eos, only: get_u
      use utils_gr, only: get_metric3plus1
      integer, intent(inout) :: ntests,npass
      integer :: i,j
      real :: x(1:3),v(1:3),dens,u,p,alpha,beta(1:3),gij(1:3,1:3),bigv(1:3)
      
      write(*,"(/,a,/)") '--> testing conservative2primitive solver'
      x = (/10.,0.,0./)
      v = (/-0.5,0.,0./)
      dens = 10.
      
      call get_metric3plus1(x,alpha,beta,gij)
      
      ! p = 40./3.
      ! p = 0.
      
      ! x=(/-4.6187683284457659E-003,0.0000000000000000,0.0000000000000000/)
      ! v=(/0,0,0/)
      ! dens=10.017642316311893
      ! p=13.356856421749205
      ! u=2.0000000000000018
      ! call test_cons2prim_i(x,v,dens,u,p,ntests,npass) 
      do i=0,9
         bigv = (/i*0.1,0.,0./)
         v = bigV*alpha-beta
         do j=0,1
            p = j*.1 
            call get_u(u,p,dens)
            call test_cons2prim_i(x,v,dens,u,p,ntests,npass)
         enddo
      enddo
      
      write(*,"(/,a,/)") '<-- conservative2primitive test complete'
      
   end subroutine test_cons2prim
   
   subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
      use cons2prim, only: conservative2primitive,primitive2conservative,error_to_string
      use testutils, only: checkval
      use checks, only: check
      use testmetric, only: metric_test
      real, intent(in) :: x(1:3),v(1:3),dens,u,p
      integer, intent(inout) :: ntests,npass
      real :: rho,pmom(1:3),en
      real :: v_out(1:3),dens_out,u_out,p_out
      real, parameter :: tol = 4.e-12
      integer :: n_errors, ierr, n_error,metricfail
      integer, save :: i=0
      i = i+1
      n_errors = 0
      ! Used for initial guess in conservative2primitive
      v_out    = v
      dens_out = dens
      u_out    = u
      p_out    = p
      ntests = ntests + 1
      print*,""
      print*,"----------------------------------------------------------------"
      print*,"(cons2prim) TEST NUMBER:",i
      
      ! call check(x,v)
      call metric_test(x,v,metricfail)
      if (metricfail/=0) then
         print*,"Metric test failed so cons2prim may also fail..."
      endif
      
      ! print*,"IN :",v,dens,u,p
      call primitive2conservative(x,v,dens,u,P,rho,pmom,en)
      ! print*,' rho = ',rho,' pmom = ',pmom,' en = ',en
      ! pmom =  -1.8969934994526608E-002
      ! rho =   10.017642319400228     
      ! en   =   3.0000000000000000
      call conservative2primitive(x,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr)
      ! print*,"OUT:",v_out,dens_out,u_out,p_out
      if (ierr /= 0) print*,'ERROR: '//trim(error_to_string(ierr))
      
      ! print*,error_to_string((/1,0,1,0/))
      
      ! write(*,"(/,a,3es14.6,a,es14.6)") " v =",v," p =",p
      call checkval(ierr,0,0,n_error,'ierr = 0 for convergence')
      n_errors = n_errors + n_error
      call checkval(3,v_out,v,tol,n_error,'v_out = v')
      n_errors = n_errors + n_error
      call checkval(dens_out,dens,tol,n_error,'dens_out = dens')
      n_errors = n_errors + n_error
      call checkval(u_out,u,tol,n_error,'u_out = u')
      n_errors = n_errors + n_error
      call checkval(p_out,p,tol,n_error,'p_out = p')
      n_errors = n_errors + n_error
      
      if (n_errors>0) then
         print*,"====================="
         print*,"cons2prim test failed"
         print*,"IN :",v,dens,u,p
         print*,"OUT:",v_out,dens_out,u_out,p_out
         print*,"====================="
      else
         npass = npass + 1
      endif
   end subroutine test_cons2prim_i
   
end module testcons2prim