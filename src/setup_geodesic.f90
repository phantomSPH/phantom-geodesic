module setup
 use set_geodesic, only:setgeodesic,gtypelist,ngtypes,print_geodesic_choices,iprec
 use prompting,    only:prompt
 implicit none
 integer :: gtype
contains

!--- Setup up a single test particle

subroutine setpart(xall,vall,np,mall)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, intent(inout), dimension(:)   :: mall
 !integer :: gtype
 character(len=120)      :: filename
 integer                 :: ierr
 logical                 :: iexist

 call print_geodesic_choices
 gtype = iprec
 filename = 'grtest'//'.params'                                ! moddump should really know about the output file prefix...
 inquire(file=filename,exist=iexist)
 if (iexist) call read_option(filename,np,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_option(filename,np)
    print*,' Edit '//trim(filename)//' and rerun'
    stop
 endif

 if (gtype == 11) then
      call get_binary_mass(np, mall)
 endif

 allocate(xall(3,np),vall(3,np))
 call setgeodesic(xall,vall,mall,np,gtype)

end subroutine setpart

! Next this subroutine reads the mass of the particles so that
! they can be used for calculating force between the particles
subroutine get_binary_mass(np, mall)
 integer, intent(in) :: np
 real, allocatable,dimension(:),intent(inout) :: mall
 character(len=120)      :: filename_binary
 integer                 :: ierr
 logical                 :: iexist
 filename_binary = 'binary'//'.params'

 allocate(mall(np))
 inquire(file=filename_binary,exist=iexist)
 if (iexist) call read_binary_option(filename_binary,np,mall,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_binary_option(filename_binary,np,mall)
    print*,' Edit with mass of stars '//trim(filename_binary)//' and rerun'
    stop
 endif


end subroutine get_binary_mass
!
!---Read/write option file--------------------------------------------------
!
subroutine write_option(filename,np)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i,np

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# tde setup file'

 call write_inopt(gtype,'gtype','option to use',iunit)
 call write_inopt(np,'np','number of particles',iunit)

 close(iunit)

end subroutine write_option

subroutine read_option(filename,np,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer,          intent(out) :: np
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(gtype,'gtype',db,min=1,errcount=nerr)
 call read_inopt(np,'np',db,min=1,errcount=nerr)
 print*,gtype,"gtype in read option",np,"np in read option"
 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif

end subroutine read_option

subroutine write_binary_option(filename_binary,np,mall)
  use infile_utils, only:write_inopt
  character(len=*), intent(in) :: filename_binary
  integer, intent(in) :: np
  real, dimension(np), intent(inout) :: mall
  integer, parameter :: iunit = 23
  integer :: i

  mall(:) = 0.

  print "(a)",' writing setup options file '//trim(filename_binary)
  open(unit=iunit,file=filename_binary,status='replace',form='formatted')
  write(iunit,"(a)") '# binary setup file'

  do i = 1, np
     call write_inopt(mall(i), 'm'//char(i+48), 'mass', iunit)
  end do
  close(iunit)
end subroutine write_binary_option

subroutine read_binary_option(filename_binary,np,mall,ierr)
  use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
  character(len=*), intent(in)  :: filename_binary
  integer,          intent(out) :: ierr
  integer,          intent(in)  :: np
  real, dimension(np), intent(inout) :: mall
  integer, parameter :: iunit = 24
  integer :: nerr
  type(inopts), allocatable :: db(:)
  integer :: i

  print "(a)",'reading setup options from '//trim(filename_binary)
  nerr = 0
  ierr = 0
  call open_db_from_file(db,filename_binary,iunit,ierr)
  do i = 1, np
    call read_inopt(mall(i), 'm'//char(i+48), db, min=0., errcount=nerr)
 end do
  call close_db(db)
  if (nerr > 0) then
      print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
      ierr = nerr
  endif


end subroutine read_binary_option
end module setup
