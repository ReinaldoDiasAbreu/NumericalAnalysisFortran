! #####################################################
! title: Gauss-Seidel
! creation: 25 Jan 2020
! modification: 25 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Gauss-Seidel Iterative Method with Initial Kick
! #####################################################

program Exec
    implicit none
    integer :: i, j, m, inte
    double precision  :: Matrix(3,3), V(3), X(3), tol

    write(*,*) "Order (m): "
    read(*,*) m

    write(*,*) "Matrix Elements: "
    do i=1,m
        read (*,*) (Matrix(i,j), j = 1, m)
    end do

    write(*,*) "Independent values: "
    read (*,*) V

    write(*,*) "First kick: "
    read (*,*) X

    write(*,*) "Tolerance: "
    read (*,*) tol

    write(*,*) "Maximum interactions: "
    read (*,*) inte

    call Gauss_Seidel(Matrix, m, V, X, tol, inte)
    
end program Exec

subroutine Gauss_Seidel(Matrix, m, V, X, tol, inte)
    ! Matrix - System Matrix
    ! m      - Matrix order
    ! V      - Independent values
    ! X      - First kick
    ! tol    - Tolerance
    ! inte   - Maximum interactions

    implicit none
    integer :: m, inte, p, i, j
    double precision :: Matrix(m,m), V(m), X(m), Xant(m), tol, aux, d

    aux = 0
    Xant = X

    do p = 1, inte
        do i = 1, m
            do j = 1, m 
                if (i /= j) then
                    aux = aux + Matrix(i,j) * Xant(j)
                end if
            end do
            Xant(i) = ((V(i) - aux) / (Matrix(i,i)))
            aux = 0
        end do
        d = abs(maxval(Xant) - maxval(X)) / abs(maxval(Xant))
        if (d <= tol) then
            exit
        end if
    end do

    write(*,*) X
end

subroutine Print_Matrix(Matrix, V, m, n)
    implicit none
    integer :: m, n, i
    real :: Matrix(m,n), V(m)

    do i = 1, m
        print *, Matrix(i,:), "|", V(i)
    end do

end
