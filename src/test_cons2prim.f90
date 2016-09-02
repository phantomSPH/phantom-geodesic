module testcons2prim
   implicit none
contains
   subroutine test_cons2prim(ntests,npass)
      use eos, only: get_u
      use utils_gr, only: get_metric3plus1
      integer, intent(inout) :: ntests,npass
      integer :: i,j
      real :: x(1:3),v(1:3),dens,u,p,alpha,beta(1:3),gij(1:3,1:3),bigv(1:3)

      write(*,'(/,a,/)') '--> testing conservative2primitive solver'

      x = (/10.,0.,0./)
      dens = 10.
      call get_metric3plus1(x,alpha,beta,gij)
      do i=0,5
         bigv = (/i*0.1,0.,0./)
         v = bigV*alpha-beta
         do j=0,10
            p = j*.1
            call get_u(u,p,dens)
            call test_cons2prim_i(x,v,dens,u,p,ntests,npass)
         enddo
      enddo

      write(*,'(/,a,/)') '<-- conservative2primitive test complete'

   end subroutine test_cons2prim

   subroutine test_cons2prim_i(x,v,dens,u,p,ntests,npass)
      use cons2prim, only: conservative2primitive,primitive2conservative,error_to_string
      use testutils, only: checkval,checkvalbuf
      use checks, only: check
      use testmetric, only: test_metric_i
      real, intent(in) :: x(1:3),v(1:3),dens,u,p
      integer, intent(inout) :: ntests,npass
      real :: rho,pmom(1:3),en
      real :: v_out(1:3),dens_out,u_out,p_out
      real, parameter :: tol = 4.e-12
      integer :: nerrors, ierr,metricpass, j
      integer :: ncheck,dummy
      real :: errmax

      ntests = ntests + 1
      nerrors = 0

      ! Used for initial guess in conservative2primitive
      v_out    = v
      dens_out = dens
      u_out    = u
      p_out    = p

      metricpass = 0
      call test_metric_i(x,v,dummy,metricpass)
      if (metricpass==0) then
         print*,'Warning: Metric test failed so cons2prim may also fail...'
      endif

      call primitive2conservative(x,v,dens,u,P,rho,pmom,en)
      call conservative2primitive(x,v_out,dens_out,u_out,p_out,rho,pmom,en,ierr)

      ! if (ierr /= 0) print*,'ERROR: '//trim(error_to_string(ierr))
      ! print*,error_to_string((/1,0,1,0/))

      ! call checkval(ierr,0,0,n_error,'ierr = 0 for convergence')
      call checkvalbuf(ierr,0,0,'[F]: ierr (convergence)',nerrors,ncheck)
      ! nerrors = nerrors + n_error

      ! call checkval(3,v_out,v,tol,n_error,'v_out = v')
      do j=1,3
         call checkvalbuf(v_out(j),v(j),tol,'[F]: v_out',nerrors,ncheck,errmax)
         ! nerrors = nerrors + n_error
      enddo

      ! call checkval(dens_out,dens,tol,n_error,'dens_out = dens')
      call checkvalbuf(dens_out,dens,tol,'[F]: dens_out',nerrors,ncheck,errmax)
      ! nerrors = nerrors + n_error

      ! call checkval(u_out,u,tol,n_error,'u_out = u')
      call checkvalbuf(u_out,u,tol,'[F]: u_out',nerrors,ncheck,errmax)
      ! nerrors = nerrors + n_error

      ! call checkval(p_out,p,tol,n_error,'p_out = p')
      call checkvalbuf(p_out,p,tol,'[F]: p_out',nerrors,ncheck,errmax)
      ! nerrors = nerrors + n_error

      if (nerrors/=0) then
         print*,'-- cons2prim test failed with'
         print*,'  - IN:'
         print*,'     x    =',x
         print*,'     v    =',v
         print*,'     dens =',dens
         print*,'     u    =',u
         print*,'     p    =',p
         print*,'  - OUT:'
         print*,'     v    =',v_out
         print*,'     dens =',dens_out
         print*,'     u    =',u_out
         print*,'     p    =',p_out
         print*,''
      else
         npass = npass + 1
      endif
   end subroutine test_cons2prim_i

end module testcons2prim
