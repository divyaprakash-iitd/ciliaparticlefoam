program check_mod
    use soft_cilia
    use iso_c_binding, only: C_INT, c_double
    implicit none

    integer(C_INT)  :: n
    real(c_double) :: h
    call say_hello()
    h = 0.1d0
    !call generateellipse(n)
    print *, "n = ", n

    call generatecilia(n,h) 

end program check_mod
