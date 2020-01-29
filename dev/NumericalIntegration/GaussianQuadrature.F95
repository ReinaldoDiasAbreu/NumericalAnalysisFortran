! #####################################################
! title: Gaussian Quadrature
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Gaussian quadrature for numerical integration
! #####################################################

program GaussianQuadrature
    implicit none
    integer :: n
    real*8 :: a, b

    write(*,*) "Integration limits (a,b): "
    read(*,*) a, b

    write(*,*) "Numbers Intervals: "
    read(*,*) n

    call Gaussian_Quadrature(a,b,n)

    
end program GaussianQuadrature

subroutine Gaussian_Quadrature(a, b, n)
    ! a, b - Integration limits
    ! n    - Degree of integration
    implicit none
    integer :: i, n
    real*8 :: area, h, h2, a, b, x(n)
    real*8 :: w1, w2, w3, t1, t2, t3, alpha, bet, x1, x2, x3

    h = (b-a)/n
    h2 = (b-a)/(n-1)
    x(1) = a
    do i = 2, n
        x(i) = x(i-1) + h2
    end do
    
    w1 = 5.0/9.0
    t1 = -(dsqrt(dble(3.0/5.0)))
    w2 = 8.0/9.0
    t2 = 0.0
    w3=w1
    t3 = -(t1)
    area = 0.0
    
    do i=1, n
        alpha = (x(i+1) - x(i))/2.0
        bet = (x(i+1) +x(i))/2.0
        x1 = alpha*t1+bet;
        x2 = alpha*t2+bet;
        x3 = alpha*t3+bet;
        
        area = area+ ((w1*func(x1) + w2*func(x2) + w3*func(x3)) * (h/2.0))
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

end subroutine Gaussian_Quadrature

