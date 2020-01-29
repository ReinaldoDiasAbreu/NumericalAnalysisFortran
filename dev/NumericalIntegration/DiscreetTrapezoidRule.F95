! #####################################################
! title: Discreet Trapezoid Tule
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Trapezoidal rule for numerical integration
! #####################################################

program Trapezoidal
    implicit none
    integer :: n
    real*8, allocatable, dimension(:) :: x, y
    real*8 :: h, a, b

    write(*,*) "Vector Size: "
    read(*,*) n

    write(*,*) "Integration limits (a,b): "
    read(*,*) a,b

    write(*,*) "Domain range: "
    read(*,*) h

    allocate(x(n))
    allocate(y(n))

    write(*,*) "Domain vector: "
    read(*,*) x

    write(*,*) "Image vector: "
    read(*,*) y

    call TrapezoidalRule(x, y, a, b, h, n)

    deallocate(x)
    deallocate(y)
    
end program Trapezoidal

subroutine TrapezoidalRule(x, y, a, b,  h, n)
    ! x, y - domain and function image
    ! a, b - Integration limits
    ! h    - domain range
    ! n    - vector size
    implicit none
    integer :: n, ntrap
    real*8, dimension(n) :: x, y
    real*8 :: area, h, a, b, jump

    ntrap = idint((b-a)/h)
    jump = a
    ! An interpolating polynomial was used for the function images
    area = Newton_Polynomial(x, y, n, a)
    do while(jump < b)
        area = area + (2*Newton_Polynomial(x, y, n, jump))
        jump = jump + h
    end do
    area = area + Newton_Polynomial(x, y, n, b)

    print *, "Solution: ", (area * (h/2))

    contains

    real*8 function Lagrange_Polynomial(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        double precision, dimension(n) :: x, y
        double precision :: ent
        double precision, allocatable, dimension(:)  :: L, inter
        real*8 :: num

        allocate(L(n))
        allocate(inter(n))

        do i=1, n
            L(i) = 1
            do j = 1, n 
                if (j /= i) then
                    L(i) = L(i) * (ent - x(j)) / (x(i)- x(j)) 
                end if
            end do
        end do

        inter = L * y

        num = sum(inter)
        return

        deallocate(L)
        deallocate(inter)

    end function Lagrange_Polynomial

    real*8 function Newton_Polynomial(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        double precision, dimension(n) :: x, y
        double precision :: ent, r
        double precision, allocatable, dimension(:)  :: inter
        real*8 :: num
    
        allocate(inter(n))
    
        do i=1, n
            inter(i) = y(i)
        end do
    
        do j=1, n-1
            do i=n,j+1,-1
                inter(i) = (inter(i)-inter(i-1))/(x(i)-x(i-j))
            end do
        end do
            
        r = inter(n)
        do i=n-1,1,-1
                r = r * (ent - x(i)) + inter(i)
        end do

        num = r
        return
    
        deallocate(inter)
    
    end function Newton_Polynomial

    real*8 function GregoryNewton_Polynomial(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        double precision, dimension(n) :: x, y
        double precision :: ent, r, u
        double precision, allocatable, dimension(:)  :: inter
        real*8 :: num
    
        allocate(inter(n))
    
        do i = 1, n
            inter(i) = y(i)
        end do
    
        do j = 1, n-1
            do i = n, j+1, -1
                inter(i) = (inter(i) - inter(i-1))
            end do
        end do
    
        u = (ent - x(1)) / (x(2) - x(1))
        r = inter(n)
    
        do i = n-1, 1, -1
            r = r * (u - i + 1) / i + inter(i)
        end do
    
        num = r
        return
    
        deallocate(inter)
    
    end function GregoryNewton_Polynomial

end subroutine TrapezoidalRule

