! Exemple file

program main
    use numerical_analysis
    implicit none
    integer :: i, n
    real*8 :: a, b, alf1, alf2, h
    real*8, allocatable, dimension(:,:) :: solution

    ! Solving system of differential equations
    ! NOTE: The isolated differential equations are 
    !       already written in the function RungeKuttaSist()
    ! Y1'(x,y) = (1.5*y1)-0.5*y1*y2
    ! Y2'(x,y1,y2) = (1.5*y1)-0.5*y1*y2

    ! Defining interval
    a = 0
    b = 20
    ! Setting initial values
    alf1 = 5
    alf2 = 2
    ! Setting step
    h = 0.05
    ! NOTE : With the step value we can calculate the output value, according to documentation
    n = idint((b - a) / h)
    allocate(solution(n+1,3))
    ! Running
    solution = RungeKuttaSist(a, b, alf1, alf2, h)
    !Printing solution
    do i=1, n 
        write(*,*) solution(i,1), "|", solution(i,2), "|", solution(i,3)
    end do
    
end program main


