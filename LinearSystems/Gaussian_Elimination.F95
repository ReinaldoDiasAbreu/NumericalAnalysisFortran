! #####################################################
! title: Gaussian elimination
! creation: 25 Jan 2020
! modification: 25 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Gaussian elimination with partial pivoting of the lines
! #####################################################

program GaussianElimination
    implicit none
    integer :: i, j, m
    double precision  :: Matrix(3,3), V(3)

    write(*,*) "Order (m): "
    read(*,*) m

    write(*,*) "Matrix Elements: "
    do i=1,m
        read (*,*) (Matrix(i,j), j = 1, m)
    end do

    write(*,*) "Independent values: "
    read (*,*) V

    call gaussian(Matrix, m, V)
    

end program GaussianElimination

subroutine gaussian(Matrix, m, V)
    ! Matrix - System Matrix
    ! m      - Matrix Order
    ! V      - independent values

    implicit none
    integer :: m, i, j, p, linha_pivo
    double precision  :: V(m), Matrix(m,m), pivo, aux(m), auxv, mult
    
    do p = 1, m-1
        pivo = Matrix(p,p)

        do i = p+1, m
            if (abs(Matrix(i,p)) > abs(pivo)) then
                linha_pivo = i
                pivo = Matrix(i,p)
                ! permutation of matrix lines
                aux = Matrix(p,:)
                Matrix(p,:) = Matrix(linha_pivo,:)
                Matrix(linha_pivo,:) = aux
                ! permutation of independent values
                auxv = V(p)
                V(p) = V(linha_pivo)
                V(linha_pivo) = auxv
            end if
        end do

        do i = p+1, m
            mult = Matrix(i,p) / pivo
            do j = 1, m
                Matrix(i,j) = Matrix(i,j) - (mult * Matrix(p,j))
            end do
            V(i) = V(i) - mult * V(p)
        end do

    end do

    ! Substitutions of x values
    V(m) = V(m) / Matrix(m,m)
    do i = m-1, 1, -1
        do j = m, i+1, -1
            auxv = auxv + Matrix(i,j)*V(j)
        end do
        V(i) = (V(i) - auxv) / Matrix(i,i)
        auxv = 0
    end do

    write(*,*) "Solution: ", V

end

