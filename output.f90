module output
   implicit none
contains
   subroutine write_out(time,xall,vall,np)
      use utils_gr, only: dot_product_gr
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
      write(50,"(a12,6a26)") 'Particle #','x','y','z','vx','vy','vz'
      do i=1,np
         x = xall(:,i)
         v = vall(:,i)
         write(50,"(I12,6e26.16)") i,x,v
      enddo
   end subroutine write_out

   subroutine write_ev(time,energy,angmom)
      real, intent(in) :: time, energy, angmom
      integer, save :: i=0

      if (i==0) then
         open(unit=55, file='output.ev',status='replace')
         write(55,*) '# Time, Energy, Angular momentum'
      else
         open(unit=55, file='output.ev',position='append')
      endif
      i = i+1
      write(55,*) time, energy, angmom
      close(unit=55)
   end subroutine write_ev

   subroutine write_xyz(time,xall,np)
      real, intent(in) :: time
      integer, intent(in) :: np
      real, dimension(1:3,np) :: xall
      integer, save :: j=0

      if (j==0) then
         open(unit=66, file='positions.dat',status='replace')
         write(66,*) '# Number of particles (n)'
         write(66,*) np
         write(60,*) '# First Column = Time'
         write(60,*) '# Subsequent columns =  x(1),x(2),...x(n),y(1),y(2),...y(n),z(1),z(2),...,z(n)'
      else
         open(unit=66, file='positions.dat',position='append')
      endif
      j = j+1
      write(66,*) time, xall(1,:),xall(2,:),xall(3,:)
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
      write(70,*) time, vall(1,:),vall(2,:),vall(3,:)
      close(unit=70)
   end subroutine write_vxyz
end module
