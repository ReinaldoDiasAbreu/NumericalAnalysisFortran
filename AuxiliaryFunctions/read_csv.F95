! #####################################################
! title: Read CSV
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Read CSV files
! #####################################################

program CSV
    implicit none
    integer :: ref, lines
    character(len=100) :: arq, act, formato
    real*8 :: matriz(10,2)

    ref = 10
    arq = "texto.csv"
    act = "read"
    formato = '(2(F10.3,1X))'
    lines = 10

    call open_file(ref, arq, act)
    call read_file(ref, formato, matriz, lines)
    call close_file(ref)
    act = "write"
    call open_file(ref, arq, act)
    call write_file(ref, formato, matriz, lines)
    call close_file(ref)

end program CSV

subroutine open_file(ref, arq, act)
    implicit none
    integer :: ref
    character(len=100), intent(in) :: arq, act
    open(unit=ref, file=arq, status='UNKNOWN', action=act)
end subroutine open_file

subroutine close_file(ref)
    implicit none
    integer :: ref
    close(unit=ref)
end subroutine close_file

subroutine read_file(ref, formato, matriz, lines)
    implicit none
    integer :: i, ref, lines
    character(len=100), intent(in) :: formato
    real*8 :: matriz(10,2)
    do i=1, lines
        read(unit=ref, fmt=formato) matriz(i,1), matriz(i,2)
        print *, '|', matriz(i,1), '|', matriz(i,2), '|'
    end do
end subroutine read_file

subroutine write_file(ref, formato, matriz, lines)
    implicit none
    integer :: i, ref, lines
    character(len=100), intent(in) :: formato
    real*8 :: matriz(10,2)
    do i=1, lines
        write(unit=ref, fmt=formato) matriz(i, 1), matriz(i, 2)
    end do
end subroutine write_file
