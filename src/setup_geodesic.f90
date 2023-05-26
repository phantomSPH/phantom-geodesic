module setup
 use set_geodesic, only:setgeodesic,gtypelist,ngtypes,print_geodesic_choices,iprec
 use prompting,    only:prompt
 implicit none
 integer :: gtype
contains

!--- Setup up a single test particle

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 !integer :: gtype
 character(len=120)      :: filename
 integer                 :: ierr
 logical                 :: iexist
 np = 1
 allocate(xall(3,np),vall(3,np))

 call print_geodesic_choices
 gtype = iprec
 print*,gtype,"gtype initial"
 ! call prompt(' Enter geodesic choice:',gtype)
 filename = 'grtest'//'.params'                                ! moddump should really know about the output file prefix...
 inquire(file=filename,exist=iexist)
 if (iexist) call read_option(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_option(filename)
    print*,' Edit '//trim(filename)//' and rerun'
    stop
 endif
 print*,gtype,"gtype from file"
 call setgeodesic(xall(:,1),vall(:,1),gtype)
end subroutine setpart
!
!---Read/write option file--------------------------------------------------
!
subroutine write_option(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# tde setup file'

 call write_inopt(gtype,'gtype','option to use',iunit)

 close(iunit)

end subroutine write_option

subroutine read_option(filename,ierr)
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

 call read_inopt(gtype,'gtype',db,min=1,errcount=nerr)
 print*,gtype,"gtype in read option"
 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif

end subroutine read_option

end module setup
