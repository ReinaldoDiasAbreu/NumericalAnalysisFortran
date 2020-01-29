! #####################################################
! title: Euler's method
! creation: 25 Jan 2020
! modification: 27 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Euler method for solving ordinary differential equations
! #####################################################

program Euler
    implicit none
    integer :: n, i
    double precision :: h, a, alf, b
    double precision, allocatable, dimension(:) :: x, y
    double precision, allocatable, dimension(:,:) :: M
    character(len=100) :: l, act="write", formato = '(2(F10.3,1X))'

    print *, "Enter domain range: "
    read (*,*) a, b
    print *, "Enter the step: "
    read(*,*) h
    print *, "Enter the initial value of Y'(x): "
    read(*,*) alf 

    n = idint((b-a) / h)

    allocate(x(n))
    allocate(y(n))
    allocate(M(n,2))

    x(1) = a
    y(1) = alf
    do i=1, n 
        y(i+1) = y(i) + h * func(x(i), y(i))
        x(i+1) = a + i * h
    end do

    print *, "Output file name: "
    read(*,*) l

    M(:,1) = x
    M(:,2) = y

    do i = 1, n
        write(*,*) M(i,1), " - ", M(i,2)
    end do

    call open_file(10, l, act)
    call write_file(10, formato, M, n)
    call close_file(10)

    deallocate(x)
    deallocate(y)
    deallocate(M)
    
    contains

    real*8 function func(x, y) result(r)
        implicit none
        double precision, intent(in) :: x, y
        ! Differential equation with Y'(x) isolated
        r = cos(x)      
    end function func
    
end program Euler

! Auxiliary functions

subroutine open_file(ref, arq, act)
    implicit none
    integer :: ref
    character(len=100), intent(in) :: arq, act
    open(unit=ref, file=arq, status='UNKNOWN', action=act)
end subroutine open_file

subroutine close_file(ref)
    implicit none
    integer :: ref
    close(unit=ref)
end subroutine close_file

subroutine read_file(ref, formato, matriz, lines)
    implicit none
    integer :: i, ref, lines
    character(len=100), intent(in) :: formato
    real*8 :: matriz(10,2)
    do i=1, lines
        read(unit=ref, fmt=formato) matriz(i,1), matriz(i,2)
        print *, '|', matriz(i,1), '|', matriz(i,2), '|'
    end do
end subroutine read_file

subroutine write_file(ref, formato, matriz, lines)
    implicit none
    integer :: i, ref, lines
    character(len=100), intent(in) :: formato
    real*8 :: matriz(10,2)
    do i=1, lines
        write(unit=ref, fmt=formato) matriz(i, 1), matriz(i, 2)
    end do
end subroutine write_file





