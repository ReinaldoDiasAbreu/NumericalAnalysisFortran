! #####################################################
! title: Simpson method (3/8 Rule)
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Simpson method for numerical integration
! #####################################################

program Simpson
    implicit none
    integer :: n
    real*8 :: a, b

    write(*,*) "Integration limits (a,b): "
    read(*,*) a, b

    write(*,*) "Degree of integration: "
    read(*,*) n

    call SimpsonMethod(a,b,n)

    
end program Simpson

subroutine SimpsonMethod(a, b, n)
    ! a, b - Integration limits
    ! n    - Degree of integration
    implicit none
    integer :: i, n
    real*8 :: area, h, h2, a, b, x(n)
    real*8 :: x1, x2, x3, x4, A1, A2, A3, A4
   
    h = (b-a)/n
    h2 = (b-a)/(n-1)
    x(1) = a
    do i = 2, n
        x(i) = x(i-1) + h2
    end do
    area = 0.00
    do i=1, n
        x1 = x(i)
        x2 = x1 + h/3
        x3 = x2 + h/3
        x4 = x(i+1)
        A1 = 1.0/8.0
        A2 = 3.0/8.0
        A3 = 3.0/8.0
        A4 = 1.0/8.0 
        area = area + ((A1*func(x1) + A2*func(x2) + A3*func(x3) + A4*func(x4)) * h)
    end do
    
    print *, area
    
    contains

    real*8 function func(x)
        implicit none
        real*8 :: x, y 
        !f(x) function to integrate 
        y = 2*x
        return
    end function func

end subroutine SimpsonMethod

