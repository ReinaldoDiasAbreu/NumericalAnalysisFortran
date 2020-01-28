! #####################################################
! title: Trapezoid Rule
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Trapezoidal rule for numerical integration
! #####################################################

program Trapezoidal
    implicit none
    real*8 :: h, a, b

    write(*,*) "Integration limits (a,b): "
    read(*,*) a, b

    write(*,*) "Domain range: "
    read(*,*) h

    call TrapezoidalRule(a, b, h)

    
end program Trapezoidal

subroutine TrapezoidalRule(a, b, h)
    ! a, b - Integration limits
    ! h    - domain range
    implicit none
    integer :: ntrap
    real*8 :: area, h, a, b, jump

    ntrap = idint((b-a)/h)
    jump = a
    ! An interpolating polynomial was used for the function images
    area = func(a)
    do while(jump < b)
        area = area + (2*func(jump))
        jump = jump + h
    end do
    area = area + func(b)

    print *, "Solution: ", (area * (h/2))

    contains

    real*8 function func(x)
        implicit none
        real*8 :: x, y 
        !f(x) function to integrate 
        y = 2*x
        return
    end function func

end subroutine TrapezoidalRule

