! #####################################################
! title: Lagrange Interpolator Polynomial
! creation: 25 Jan 2020
! modification: 27 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Exemplo
! #####################################################

program Lagrange
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

    call Lagrange_Polynomial(x, y, n, ent)

    deallocate(x)
    deallocate(y)
    
end program Lagrange

subroutine Lagrange_Polynomial(x, y, n, ent)
    ! x, y - Domain and Image
    ! n    - Vectors Size
    ! ent  - Value to Interpolate
    implicit none
    integer :: n, i, j
    double precision, dimension(n) :: x, y
    double precision :: ent
    double precision, allocatable, dimension(:)  :: L, inter

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

    write(*,*) "Solution: ", sum(inter)

    deallocate(L)
    deallocate(inter)

end subroutine Lagrange_Polynomial

