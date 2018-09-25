module infile
use options,      only:dt,tmax,dnout_ev,dtout,write_pos_vel
use step,         only:steptype
use metric_tools, only:coordinate_sys

implicit none

contains

subroutine init_infile(filename)
 character(len=*), intent(in) :: filename
 integer :: ierr
 logical :: iexist

 inquire(file=filename,exist=iexist)
 if (iexist) call read_paramsfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_paramsfile(filename)
    print*,' Edit '//trim(filename)//' and rerun code'
    stop
 endif

end subroutine init_infile

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_paramsfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file'

 write(iunit,'(/,"#",30("-"),a,30("-"))') ' Timestepping '
 call write_inopt(steptype,'steptype',"(1 = Leapfrog  |  2 = RK2  |  3 = Euler  |  4 = Heun's  |  5 = L&R05)",iunit)
 call write_inopt(tmax,'tmax','simulation end time',iunit)
 call write_inopt(dt  ,'dt'  ,'timestep (fixed)'   ,iunit)

 write(iunit,'(/,"#",30("-"),a,30("-"))') ' Output '
 call write_inopt(dtout,'dtout',"output between dump files (-ve don't write dumps)",iunit)
 call write_inopt(write_pos_vel,'write_pos_vel',"write positions and velocities to one file (logical)",iunit)
 call write_inopt(dnout_ev,'dnout_ev',"frequency of writing to ev file (number of steps)",iunit)

 write(iunit,'(/,"#",30("-"),a,30("-"))') ' Coordinates'
 call write_inopt(coordinate_sys,'coordinate_sys',"'Cartesian' or 'Spherical'",iunit)

 close(iunit)

end subroutine write_paramsfile

subroutine read_paramsfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(steptype,'steptype',db,min=1,errcount=nerr)
 call read_inopt(tmax,'tmax',db,min=0.,errcount=nerr)
 call read_inopt(dt,'dt',db,min=0.,errcount=nerr)
 call read_inopt(dtout,'dtout',db,errcount=nerr)
 call read_inopt(write_pos_vel,'write_pos_vel',db,errcount=nerr)
 call read_inopt(dnout_ev,'dnout_ev',db,min=0,errcount=nerr)
 call read_inopt(coordinate_sys,'coordinate_sys',db,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif

end subroutine read_paramsfile

end module infile
