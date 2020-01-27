! #####################################################
! title: Tridiagonal Systems
! creation: 25 Jan 2020
! modification: 27 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Algorithm for Tridiagonal Systems
! #####################################################

program Tridiagonal
    implicit none
    integer          :: order
    double precision, allocatable, dimension(:) :: a, b, c, d

    write(*, *) "Order: "
    read(*,*) order

    allocate(a(order))
    allocate(b(order))
    allocate(c(order))
    allocate(d(order))

    write(*,*) "First vector: "
    read(*,*) a
    write(*,*) "Second vector: "
    read(*,*) b
    write(*,*) "Third vector: "
    read(*,*) c
    write(*,*) "Independent values: "
    read(*,*) d

    call Tridiagonal_System(a, b, c, d, order)

    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(d)


end program Tridiagonal


subroutine Tridiagonal_System(a, b, c, d, order)
    !a, b ,c  - system vectors 
    !d        - independent values
    !order    - system order
    implicit none
    integer          :: order, i
    double precision, dimension(order) :: a, b, c, d
    double precision, allocatable, dimension(:) :: x

    allocate(x(order))
    x = 0
    c(1) = c(1) / b(1)  
    do i=2, order-1  
        c(i) = c(i) / (b(i) - a(i) * c(i-1))  
    end do 
 
    d(1)=d(1)/b(1)  
    do i=2, order
        d(i) = (d(i) - a(i) * d(i-1)) / (b(i) - a(i) * c(i-1))  
    end do
 
    x(order) = d(order)  
    do i=order-1, 1, -1
        x(i) = d(i) - c(i) * x(i+1)  
    end do
    
    write(*,*) "Solution: ", x

    deallocate(x)

end subroutine Tridiagonal_System

