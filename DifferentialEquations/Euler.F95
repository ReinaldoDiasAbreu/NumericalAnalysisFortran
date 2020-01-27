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
    logical :: l

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

    print *, "Display result? (T/F) "
    read(*,*) l
    print *,""
    if(l .EQV. .TRUE.) then
        M(:,1) = x
        M(:,2) = y
        do i = 1, n
            print *, M(i,:)
        end do  
    end if

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




