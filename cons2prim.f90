module cons2prim
   implicit none

   public :: conservative2primitive, primitive2conservative, get_p_from_v, get_v_from_p
   public :: error_to_string
   integer, parameter :: ierr_notconverged = 1

   private

contains

   elemental function error_to_string(ierr) result(string)
      integer, intent(in) :: ierr
      character(len=20) :: string

      select case(ierr)
      case(ierr_notconverged)
         string = 'not converged'
      case default
         string = 'unknown error'
      end select
   end function error_to_string

   !----------------------------------------------------------------
   !+
   !  Construct conserved variables from the primitive variables
   !  primitive variables are (v^i,d,u,P); v i
   !  conserved variables are (rho,pmom_i,en)
   !+
   !----------------------------------------------------------------
   subroutine primitive2conservative(x,v,dens,u,P,rho,pmom,en)
      use utils_gr, only: dot_product_gr
      use metric, only: get_metric
      real, intent(in)  :: x(1:3)
      real, intent(in) :: dens,v(1:3),u,P
      real, intent(out)  :: rho,pmom(1:3),en
      real, dimension(0:3,0:3) :: gcov, gcon
      real :: sqrtg, enth, v3d, v3_DOT_v3, U0, v4U(0:3)
      integer :: i, mu

      v4U(0) = 1.
      v4U(1:3) = v(:)

      enth = 1.+ u + P/dens

      call get_metric(x,gcov,gcon,sqrtg)
      U0  = 1./sqrt(-dot_product_gr(v4U,v4U,gcov))
      rho = sqrtg*dens*U0
      do i=1,3
         v3d = dot_product(gcov(i,:),v4U(:))
         pmom(i) = U0*enth*v3d
      enddo

      v3_DOT_v3 = 0.
      do mu=1,4
         do i=1,3
            v3_DOT_v3 = v3_DOT_v3 + gcov(i,mu)*v(mu)*v(i)
         enddo
      enddo
      en = U0*(enth*v3_DOT_v3 - (1.+u)*dot_product_gr(v4U,v4U,gcov))

   end subroutine primitive2conservative

   subroutine conservative2primitive(x,v,dens,u,P,rho,pmom,en,ierr)
      use utils_gr, only: dot_product_gr, get_metric3plus1
      use metric, only: get_metric
      use eos, only: get_enthalpy, get_u, gam
      real, intent(in)  :: x(1:3)
      real, intent(inout) :: v(1:3),dens,u,P
      real, intent(in)  :: rho,pmom(1:3),en
      integer, intent(out) :: ierr
      real, dimension(0:3,0:3) :: gcov,gcon
      real :: gij(1:3,1:3)
      real :: sqrtg,enth,lorentz_LEO,pmom2,alpha,beta(1:3),beta2,enth_old,v3d(1:3)
      real :: f,df
      integer :: niter, i
      real, parameter :: tol = 1.e-10
      integer, parameter :: nitermax = 100
      logical :: converged
      ierr = 0

      call get_metric(x,gcov,gcon,sqrtg)
      call get_metric3plus1(x,alpha,beta,gij)
      beta2 = dot_product_gr(beta,beta,gcon(1:3,1:3))
      pmom2 = dot_product_gr(pmom,pmom,gcon(1:3,1:3))

      call get_enthalpy(enth,dens,p)

      niter = 0
      converged = .false.
      do while (.not.converged)
         enth_old = enth

         lorentz_LEO = sqrt(1.+pmom2/enth_old**2)
         dens = rho*alpha/(sqrtg*lorentz_LEO)
         p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product_gr(pmom,beta,gcon(1:3,1:3))),0.)
         call get_enthalpy(enth,dens,p)

         f = enth-enth_old
         df= -1.+(gam/(gam-1.))/alpha*(1.-pmom2*p/(enth_old**3*lorentz_LEO**2*dens))
         enth = enth_old - f/df

         niter = niter + 1
         converged = (abs(enth-enth_old)/enth < tol)
         if (niter > nitermax) then
            write(*,"(/,a)") " WARNING: reached max number of iterations"
            print*,"dens,v,u,p:",dens,v,u,p
            exit
         endif
      enddo

      if (.not.converged) then
         print*,'enthold,enth,rel_err=',enth_old,enth,abs(enth-enth_old)/enth
         ierr = ierr_notconverged
         return
      endif

      lorentz_LEO = sqrt(1.+pmom2/enth**2)
      dens = rho*alpha/(sqrtg*lorentz_LEO)
      p = max(rho/sqrtg*(enth*lorentz_LEO*alpha-en-dot_product_gr(pmom,beta,gcon(1:3,1:3))),0.)
      v3d = alpha*pmom/(enth*lorentz_LEO)-beta
      do i=1,3
         v(i) = dot_product(gcon(1:3,i),v3d(1:3)) ! Raise index from down to up
      enddo
      call get_u(u,P,dens)

   end subroutine conservative2primitive

   subroutine get_v_from_p(pmom,v,x)
      ! Conservative to primitive solver for the dust case (Pressure=0)
      use metric, only: get_metric
      use utils_gr, only: dot_product_gr
      real, intent(in) :: pmom(1:3), x(1:3)
      real, intent(out) :: v(1:3)
      real, dimension(0:3,0:3) :: gcov, gcon
      real :: beta(1:3), alpha, sqrtg, pmom2, beta2,gamma,v_down(1:3)
      integer :: i

      call get_metric(x,gcov,gcon,sqrtg)
      beta  = gcov(0,1:3)
      beta2 = dot_product_gr(beta,beta,gcon(1:3,1:3))
      alpha = sqrt(beta2 - gcov(0,0))
      pmom2 = dot_product_gr(pmom,pmom,gcon(1:3,1:3))
      gamma = sqrt(1+pmom2)
      v_down = alpha*pmom/(sqrt(1+pmom2)) - beta

      !v_down = pmom/gamma

      do i=1,3
         v(i) = dot_product(gcon(1:3,i),v_down(1:3))
      enddo


      !print*,'alpha,pmom,beta',alpha,pmom,beta
      !print*,v
      !stop 'printing v'
      ! en = 0. ! ???
      ! P  = 0.
      ! rho = 0. ! ???
      ! call conservative2primitive(x,v,dens,u,P,rho,pmom,en)


   end subroutine get_v_from_p

   subroutine get_p_from_v(pmom,v,x)
      use metric, only: get_metric
      real, intent(in) :: v(1:3), x(1:3)
      real, intent(out) :: pmom(1:3)
      real :: rho, en, dens, u, P

      dens = 1. ! ???
      u = 0.
      P = 0.
      call primitive2conservative(x,v,dens,u,P,rho,pmom,en)


   end subroutine get_p_from_v

end module cons2prim
