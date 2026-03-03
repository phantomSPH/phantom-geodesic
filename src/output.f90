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
subroutine write_out(time,xall,vall,np,mall)
 use utils_gr,     only: dot_product_gr
 use metric,       only: spherical2cartesian
 use metric_tools, only: coordinate_sys
 real,    intent(in) :: time
 integer, intent(in) :: np
 real,    intent(in) :: xall(3,np),vall(3,np),mall(np)
 integer, parameter :: iu = 50
 integer, save      :: ifile = -1
 character(len=40)  :: filename
 real    :: x(3),v(3)
 integer :: i
 ifile = ifile+1
 write(filename,"(a,i5.5,a)") 'output_',ifile,'.dat'
 open(unit=iu, file=filename, status='replace')
 write(iu,"('#',a)") ' grtest dumpfile'
 write(iu,*) time
 write(iu,"('#',i20)") np
 write(iu,"('#',7a26)") 'x','y','z','vx','vy','vz','m'
 do i=1,np
    x = xall(:,i)
    v = vall(:,i)
    if (coordinate_sys == 'Spherical') call spherical2cartesian(xall(:,i),x,vall(:,i),v)
    write(iu,"(7e26.16)") x,v,mall(i)
 enddo
 close(iu)
end subroutine write_out

!----------------------------------------------------------------
!+
!  Write one output file containing the time evolution of total energy and  total angular momentum.
!+
!----------------------------------------------------------------
subroutine write_ev(time,energy,energy_init,angmom,angmom_init)
 real, intent(in) :: time,energy,energy_init,angmom,angmom_init
 integer, save :: i=0
 integer, parameter :: iu = 55
 character(len=120)   :: ev_file
 character(len=120)   :: file_input
 real :: diff_en, diff_angmom
 real :: d_en, d_ang
 
 call  get_environment_variable("ev_file", file_input)
 
 if (LEN_TRIM(file_input) > 0) then
    ev_file=trim(file_input)
 else
    ev_file=trim('ev.dat')
 endif

 diff_en = energy - energy_init
 diff_angmom = angmom - angmom_init

 d_en = abs(diff_en) / energy_init
 d_ang = abs(diff_angmom) / angmom_init

 if (i==0) then
    open(unit=iu, file=ev_file,status='replace')
    write(iu,*) '# Time, Energy, dE, Angular momentum, dL'
    i = i+1
 else
    open(unit=iu, file=ev_file,position='append')
 endif

 write(iu,*) time, energy, d_en,  angmom, d_ang
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
 character(len=6) :: nstring
 integer, parameter :: iu2=33
 logical, save :: first2 = .true.
 integer :: i,i3,i4
 character(len=11) :: ihead12,ihead22,ihead32,column_labels2(np*3)
 character(len=30) :: fmt
 character(len=120)   :: positions_file,positions_all
 CHARACTER(len=120)   :: file_input

 call  get_environment_variable("positions_file", file_input)

 if (len_trim(file_input) > 0) THEN
    positions_file=trim(file_input)
 else
    positions_file=trim('positions.dat')
 endif

 write(positions_all,"(a,i3.3,a)") 'positions.dat'
 ! write the file with all the data points.
 if (first2) then
    open(unit=iu2, file=positions_file,status='replace')
    write(iu2,*) '# Number of particles'
    write(iu2,*) np

    !-- Create an array of strings for the column labels
    do i=1,np
      select case(coordinate_sys)
       case('Cartesian')
          write(ihead12,'("x",i0)') i
          write(ihead22,'("y",i0)') i
          write(ihead32,'("z",i0)') i
       case('Spherical')
          write(ihead12,'("r",i0)') i
          write(ihead22,'("theta",i0)') i
          write(ihead32,'("phi",i0)') i
       end select
       i3 = (i-1)*3+1
       i4 = i3 + 2
       column_labels2(i3:i4) = [ihead12,ihead22,ihead32]

    enddo
    !-- Create the format string for column labels
    write(nstring,'(i0)') np*3
    fmt = '("# ",'//trim(nstring)//'a30," Time")'

    !-- Write column labels to the file
    write(iu2,fmt) column_labels2(:)

    !-- Write some other info to file
    write(iu2,*) '# Coordinate system: ',coordinate_sys,'. Metric: ',metric_type,'. Written to file in cartesian:',write_cartesian
    first2 = .false.
 else
    open(unit=iu2, file=positions_file,position='append')
 endif

 write(iu2,*) xall(1:3,:),time
 close(iu2)

end subroutine write_xyz

!----------------------------------------------------------------
!+
!  Write one output file containing the time evolution of velocity for each particle in columns.
!+
!----------------------------------------------------------------
subroutine write_vxyz(time,vall,np,force)
 use metric_tools, only: coordinate_sys
 use metric,       only: metric_type
 integer, intent(in) :: np
 real,    intent(in) :: time, vall(3,np), force(3,np)
 logical, save       :: first = .true.
 integer, parameter  :: iu = 70
 integer :: i,i1,i2
 character(len=11) :: ihead1,ihead2,ihead3,column_labels(np*3)
 character(len=6)  :: nstring
 character(len=120)   :: velocities_file
 CHARACTER(len=120)   :: file_input

 call  get_environment_variable("velocities_file", file_input)

 IF (LEN_TRIM(file_input) > 0) THEN
        velocities_file=trim(file_input)
 else
        velocities_file=trim('velocities.dat')
 endif

 if (first) then
    open(unit=iu, file=velocities_file,status='replace')
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
    write(iu,'(A)', advance='no') '# Time'        ! no newline yet
    write(iu,'(1X,'//trim(nstring)//'A30)') column_labels
    write(iu,*)                                     

    !-- Write some other info to the file
    write(iu,*) '# Coordinate system: ',coordinate_sys,' Metric: ',metric_type
    first = .false.
 else
    open(unit=iu, file=velocities_file,position='append')
 endif
 write(iu,*) time, vall(1:3,:)
 close(iu)
end subroutine write_vxyz
end module
