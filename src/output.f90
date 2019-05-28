module output
 implicit none

 logical :: write_cartesian = .true.

contains

!----------------------------------------------------------------
!+
!  Write output files containing positions of all particles at a particular time.
!  (This is useful for plotting in splash, so output is always in cartesian coordinates.)
!+
!----------------------------------------------------------------
subroutine write_out(time,xall,vall,np)
 use utils_gr,     only: dot_product_gr
 use metric,       only: spherical2cartesian
 use metric_tools, only: coordinate_sys
 real,    intent(in) :: time
 integer, intent(in) :: np
 real,    intent(in) :: xall(3,np),vall(3,np)
 integer, parameter :: iu = 50
 integer, save      :: ifile = -1
 character(len=40)  :: filename
 real    :: x(3),v(3)
 integer :: i

 ifile = ifile+1
 write(filename,"(a,i5.5,a)") 'output_',ifile,'.dat'
 open(unit=iu, file=filename, status='replace')
 write(iu,*) '# grtest dumpfile'
 write(iu,*) time
 write(iu,*) '#',np
 write(iu,"('#',6a26)") 'x','y','z','vx','vy','vz'
 do i=1,np
    x = xall(:,i)
    if (coordinate_sys == 'Spherical') call spherical2cartesian(xall(:,i),x,vall(:,i),v)
    write(iu,"(6e26.16)") x,v
 enddo
 close(iu)
end subroutine write_out

!----------------------------------------------------------------
!+
!  Write one output file containing the time evolution of total energy and  total angular momentum.
!+
!----------------------------------------------------------------
subroutine write_ev(time,energy,angmom)
 real, intent(in) :: time, energy, angmom
 integer, save :: i=0
 integer, parameter :: iu = 55

 if (i==0) then
    open(unit=iu, file='ev.dat',status='replace')
    write(iu,*) '# Time, Energy, Angular momentum'
    i = i+1
 else
    open(unit=iu, file='ev.dat',position='append')
 endif

 write(iu,*) time, energy, angmom
 close(iu)

end subroutine write_ev

!----------------------------------------------------------------
!+
!  Write one output file containing the time evolution of position for each particle in columns.
!+
!----------------------------------------------------------------
subroutine write_xyz(time,xall,np)
 use metric_tools, only: coordinate_sys
 use metric,       only: metric_type, spherical2cartesian, cartesian2spherical
 real,    intent(in) :: time
 integer, intent(in) :: np
 real, dimension(1:3,np), intent(in) :: xall
 real, dimension(1:3,np) :: x
 character(len=6) :: nstring
 integer, parameter :: iu = 66
 logical, save :: first = .true.
 integer :: i,i1,i2
 character(len=11) :: ihead1,ihead2,ihead3,column_labels(np*3)
 character(len=30) :: fmt

 if (first) then
    open(unit=iu, file='positions.dat',status='replace')
    write(iu,*) '# Number of particles (n)'
    write(iu,*) np

    !-- Create an array of strings for the column labels
    do i=1,np
      select case(coordinate_sys)
      case('Cartesian')
         write(ihead1,'("x",i0)') i
         write(ihead2,'("y",i0)') i
         write(ihead3,'("z",i0)') i
      case('Spherical')
         write(ihead1,'("r",i0)') i
         write(ihead2,'("theta",i0)') i
         write(ihead3,'("phi",i0)') i
      end select
      i1 = (i-1)*3+1
      i2 = i1 + 2
      column_labels(i1:i2) = [ihead1,ihead2,ihead3]
    enddo

    !-- Create the format string for column labels
    write(nstring,'(i0)') np*3
    fmt = '("# Time",'//trim(nstring)//'a30)'

    !-- Write column labels to the file
    write(iu,fmt) column_labels(:)

    !-- Write some other info to file
    write(iu,*) '# Coordinate system: ',coordinate_sys,'. Metric: ',metric_type,'. Written to file in cartesian:',write_cartesian

    first = .false.
 else
    open(unit=iu, file='positions.dat',position='append')
 endif

 ! Write to positions.dat file in cartesian
 if (write_cartesian) then
    if (coordinate_sys == 'Cartesian') then
       write(iu,*) time, xall(1:3,:)
    else if (coordinate_sys == 'Spherical') then
       do i = 1,np
          call spherical2cartesian(xall(:,i),x(:,i))
       enddo
       write(iu,*) time, x(1:3,:)
    else
       STOP "Please pick a coordinate system that I can write to file in"
    endif

    ! Write to positions.dat file in spherical
 else if (.not. write_cartesian) then
    if (coordinate_sys == 'Spherical') then
       write(iu,*) time, xall(1:3,:)
    else if (coordinate_sys == 'Cartesian') then
       do i=1,np
          call cartesian2spherical(xall(:,i),x(:,i))
       enddo
       write(iu,*) time, x(1:3,:)
    else
       STOP "Please pick a coordinate system that I can write to file in"
    endif

 endif
 close(iu)
end subroutine write_xyz

!----------------------------------------------------------------
!+
!  Write one output file containing the time evolution of velocity for each particle in columns.
!+
!----------------------------------------------------------------
subroutine write_vxyz(time,vall,np)
 use metric_tools, only: coordinate_sys
 use metric,       only: metric_type
 integer, intent(in) :: np
 real,    intent(in) :: time, vall(3,np)
 logical, save       :: first = .true.
 integer, parameter  :: iu = 70
 integer :: i,i1,i2
 character(len=11) :: ihead1,ihead2,ihead3,column_labels(np*3)
 character(len=30) :: fmt
 character(len=6)  :: nstring

 if (first) then
    open(unit=iu, file='velocities.dat',status='replace')
    write(iu,*) '# Number of particles'
    write(iu,*) np

    !-- Create an array of strings for the column labels
    do i=1,np
      write(ihead1,'("vx",i0)') i
      write(ihead2,'("vy",i0)') i
      write(ihead3,'("vz",i0)') i
      i1 = (i-1)*3+1
      i2 = i1 + 2
      column_labels(i1:i2) = [ihead1,ihead2,ihead3]
    enddo

    !-- Create the format string for column labels
    write(nstring,'(i0)') np*3
    fmt = '("# Time",'//trim(nstring)//'a30)'

    !-- Write column labels to the file
    write(iu,fmt) column_labels(:)

    !-- Write some other info to the file
    write(iu,*) '# Coordinate system: ',coordinate_sys,' Metric: ',metric_type
    first = .false.
 else
    open(unit=iu, file='velocities.dat',position='append')
 endif
 write(iu,*) time, vall(1:3,:)
 close(iu)
end subroutine write_vxyz
end module
