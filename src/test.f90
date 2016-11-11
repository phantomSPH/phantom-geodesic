program test
   use testmetric, only: test_metric
   use testcons2prim, only: test_cons2prim
   integer :: ntests, npass
   ntests=0
   npass=0

   write(*,"(a)") "Running tests...."
   call test_metric(ntests,npass)
   call test_cons2prim(ntests,npass)

   write(*,'(/,a,i10)') 'Number of tests = ',ntests
   write(*,'(a,i10)')   'Number passed   = ',npass
   write(*,'(a,i10,/)') 'Number failed   = ',ntests-npass

   if (ntests==npass) then
      write(*,'(a)') "-----> PASSED <-----"
   else
      write(*,'(a)') "-----> FAILED <-----"
   endif
end program test
