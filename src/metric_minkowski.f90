module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Minkowski'

 public :: get_metric, get_metric_derivs, metric_type

 real, parameter, public :: mass1 = 0.  ! mass of central object
 real, parameter, public :: a     = 0.  ! spin of central object
 real, parameter, public :: rs    = 2.*mass1

 private

contains

 !----------------------------------------------------------------
 !+
 !  Compute the metric tensor in both covariant (gcov) and
 !  contravariant (gcon) form
 !+
 !----------------------------------------------------------------
 pure subroutine get_metric(x,gcov,gcon,sqrtg)
  real,    intent(in)  :: x(3)
  real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
  gcov = 0.
  gcon = 0.
  if (metric_type .eq. 'Minkowski') then
   gcov(0,0) = -1.
   gcov(1,1) = 1.
   gcov(2,2) = 1.
   gcov(3,3) = 1.
   sqrtg=1.
   gcon = gcov
  endif
 end subroutine get_metric

 pure subroutine get_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
    real,    intent(in)  :: position(3)
    real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
   !  real,    intent(out) :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
    dgcovdx = 0.
    dgcovdy = 0.
    dgcovdz = 0.
   !  dgcondx = 0.
   !  dgcondy = 0.
   !  dgcondz = 0.
 end subroutine
end module metric
