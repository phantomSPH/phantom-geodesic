program test
   use testmetric, only: test_metric
   use testcons2prim, only: test_cons2prim
   integer :: ntests, npass
   ntests=0
   npass=0
   
   write(*,"(a)") "Running tests...."
   call test_metric(ntests,npass)
   call test_cons2prim(ntests,npass)
   
   print*,'Number of tests = ',ntests
   print*,'Number passed   = ',npass
   print*,""
   if (ntests==npass) then
      print*,"-----> PASSED <-----"
   else 
      print*,"-----> FAILED <-----"
   endif
end program test