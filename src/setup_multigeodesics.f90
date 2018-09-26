module setup
 use set_geodesic, only:setgeodesic

 implicit none

 public :: setpart

 private

 integer :: npart,gtype
 real    :: r1,dr

contains

!--- Setup up multiple test particles

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real :: r0
 integer :: i

 call init_setupfile('multigeodesic.setup')

 np = npart

 allocate(xall(3,np),vall(3,np))

 do i=1,np
    r0 = r1+(i-1)*dr
    call setgeodesic(xall(:,i),vall(:,i),gtype,r0)
 enddo

end subroutine setpart


subroutine init_setupfile(filename)
 use prompting,    only:prompt
 use set_geodesic, only:iprec,ngtypes,print_geodesic_choices
 character(len=*), intent(in) :: filename
 integer :: ierr
 logical :: iexist

 !
 ! Defaults
 !
 npart = 10
 r1    = 4.
 dr    = 4.
 gtype = iprec

 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then

    call prompt(' Enter number of particles',npart,0)
    call prompt(' Enter starting r1',r1,0.)
    call prompt(' Enter starting dr',dr,0.)
    call print_geodesic_choices()
    call prompt(' Enter geodesic choice:',gtype,1,ngtypes)

    call write_setupfile(filename)
    print*,' Edit '//trim(filename)//' and rerun code'
    stop
 endif

end subroutine init_setupfile

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use set_geodesic, only:ngtypes,gtypelist
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# multigeodesics setup file'

 call write_inopt(npart,'npart','number of particles',iunit)
 call write_inopt(r1,'r1','starting radius',iunit)
 call write_inopt(dr,'dr','spacing of particles',iunit)
 call write_inopt(gtype,'gtype','geodesic type',iunit)

 write(iunit,"(/,10('-'),a,10('-'))") ' Geodesic choices'
 do i=1,ngtypes
    write(iunit,"('#',i10,' : ',a)") i,gtypelist(i)
 enddo

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
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

 call read_inopt(npart,'npart',db,min=1,errcount=nerr)
 call read_inopt(r1,'r1',db,min=0.,errcount=nerr)
 call read_inopt(dr,'dr',db,min=0.,errcount=nerr)
 call read_inopt(gtype,'gtype',db,min=1,errcount=nerr)

 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif

end subroutine read_setupfile

end module setup
