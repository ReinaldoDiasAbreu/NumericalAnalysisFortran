! #####################################################
! title: Exemple
! creation: 25 Jan 2020
! modification: 27 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Exemplo
! #####################################################

program Exemple
    implicit none
   
    
end program Exemple

subroutine Print_Matrix(Matrix, V, m, n)
    implicit none
    integer :: m, n, i
    double precision :: Matrix(m,n), V(m)

    do i = 1, m
        print *, Matrix(i,:), "|", V(i)
    end do
end



