! #####################################################
! title: Newton Interpolator Polynomial
! creation: 25 Jan 2020
! modification: 27 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Newton Interpolator Polynomial
! #####################################################

program Newton
    implicit none
    integer :: n
    double precision, allocatable, dimension(:) :: x, y
    double precision :: ent

    write(*,*) "Vector Size: "
    read(*,*) n

    write(*,*) "Value to Interpolate: "
    read(*,*) ent

    allocate(x(n))
    allocate(y(n))

    write(*,*) "Domain vector: "
    read(*,*) x

    write(*,*) "Image vector: "
    read(*,*) y

    call Newton_Polynomial(x, y, n, ent)

    deallocate(x)
    deallocate(y)
    
end program Newton

subroutine Newton_Polynomial(x, y, n, ent)
    ! x, y - Domain and Image
    ! n    - Vectors Size
    ! ent  - Value to Interpolate
    implicit none
    integer :: n, i, j
    double precision, dimension(n) :: x, y
    double precision :: ent, r
    double precision, allocatable, dimension(:)  :: inter

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

    write(*,*) "Solution: ", r

    deallocate(inter)

end subroutine Newton_Polynomial

