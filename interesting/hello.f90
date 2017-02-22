program main
!
!first fortran 90 program
!gfortran hell.f90 && ./a.out
!
write(*,*) 'Hello'
call timestamp()

end program main

subroutine timestamp()
implicit none

        character(8) date
        character(5) time

        call date_and_time(date,time)
        write(*,'(a8,2x,a10)') date,time
end subroutine

