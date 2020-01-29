
program main
    use numerical_analysis
    implicit none
    real*8, dimension(10):: x, y
    real*8 :: a, b, h, r
    integer :: t

    x = (/1,2,3,4,5,6,7,8,9,10/)
    y = (/2,4,6,8,10,12,14,16,18,20/)
    a = 1.0
    b = 6.0
    h = 0.01
    t = 10
    r = DiscreetTrapezoidal(x,y, a, b, h, t)

    print *, r
    
end program main


