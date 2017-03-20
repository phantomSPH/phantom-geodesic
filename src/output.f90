module output
implicit none
contains

!----------------------------------------------------------------
!+
!  Write output files containing positions of all particles at a particular time.
!  (This is useful for plotting in splash, so output is always in cartesian coordinates.)
!+
!----------------------------------------------------------------
subroutine write_out(time,xall,vall,np)
   use utils_gr, only: dot_product_gr
   use metric, only: spherical2cartesian
   use metric_tools, only: coordinate_sys
   real, intent(in) :: time
   integer, intent(in) :: np
   real, dimension(1:3,np), intent(in) :: xall, vall
   ! real, dimension(:), intent(in) :: energy_init_np,angmom_init_np
   real :: x(1:3), v(1:3)
   integer :: i
   integer, save :: ifile = -1
   character(len=40) :: filename

   ifile = ifile+1
   write(filename,"(a,i5.5,a)") 'output_',ifile,'.dat'
   open(unit=50, file=filename, status='replace')
   write(50,*) time
   write(50,"(6a26)") 'x','y','z'
   do i=1,np
      x = xall(:,i)
      if (coordinate_sys == 'Spherical') call spherical2cartesian(xall(:,i),x)
      write(50,"(3e26.16)") x
   enddo
end subroutine write_out

subroutine write_ev(time,energy,angmom)
   use metric, only: metric_type
   real, intent(in) :: time, energy, angmom
   integer, save :: i=0

   if(metric_type=='Schwarzschild') then
   if (i==0) then
      open(unit=55, file='ev.dat',status='replace')
      write(55,*) '# Time, Energy, Angular momentum'
   else
      open(unit=55, file='ev.dat',position='append')
   endif
   i = i+1
   write(55,*) time, energy, angmom
   close(unit=55)
   endif
end subroutine write_ev

subroutine write_xyz(time,xall,np)
   use metric_tools, only: coordinate_sys
   use metric, only: spherical2cartesian, cartesian2spherical
   real, intent(in) :: time
   integer, intent(in) :: np
   real, dimension(1:3,np) :: xall
   integer, save :: j=0
   logical, parameter :: write_cartesian = .true.
   real :: x(3)

   if (j==0) then
      open(unit=66, file='positions.dat',status='replace')
      write(66,*) '# Number of particles (n)'
      write(66,*) np
      write(66,*) '# First Column = Time'
      write(66,*) '# Subsequent columns =  x1(1),x2(1),x3(1),.....,x1(n),x2(n),x3(n).'
   else
      open(unit=66, file='positions.dat',position='append')
   endif
   j = j+1

   ! Write to positions.dat file in cartesian
   if (write_cartesian) then
      if (coordinate_sys == 'Cartesian') then
         write(66,*) time, xall(1:3,:)
      else if (coordinate_sys == 'Spherical') then
         x = xall(:,1)
         call spherical2cartesian(xall(:,1),x)
         write(66,*) time, x
      else
          STOP "Please pick a coordinate system that I can write to file in"
      endif

   ! Write to positions.dat file in spherical
   else if (.not. write_cartesian) then
      if (coordinate_sys == 'Spherical') then
         write(66,*) time, xall(1:3,:)
      else if (coordinate_sys == 'Cartesian') then
         x = xall(:,1)
         call cartesian2spherical(xall(:,1),x)
         write(66,*) time, x
      else
          STOP "Please pick a coordinate system that I can write to file in"
      endif

   endif
   close(unit=66)
end subroutine write_xyz

subroutine write_vxyz(time,vall,np)
   real, intent(in) :: time
   integer, intent(in) :: np
   real, dimension(3,np) :: vall
   integer, save :: j

   if (j==0) then
      open(unit=70, file='velocities.dat',status='replace')
      write(70,*) '# Number of particles'
      write(70,*) np
      write(70,*) '# First Column = Time'
      write(70,*) '# Subsequent columns =  vx(1),vx(2),...vx(n),vy(1),vy(2),...vy(n),vz(1),vz(2),...,vz(n)'
   else
      open(unit=70, file='velocities.dat',position='append')
   endif
   j = j+1
   write(70,*) time, vall(1:3,:)
   close(unit=70)
end subroutine write_vxyz
end module
