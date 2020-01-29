! #####################################################
! title: Numerical_Analysis
! creation: 28 Jan 2020
! modification: 28 Jan 2020
! author: Reinaldo Junio Dias de Abreu
! contributors: @reinaldo,
! description: Numerical Analysis is a module that brings 
!               together functions for numerical analysis
! #####################################################

module numerical_analysis

    contains

    ! ------- Interpolation -------

    function Lagrange(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        real*8, dimension(n) :: x, y
        real*8, allocatable, dimension(:) :: L, inter
        real*8 :: ent, Lagrange

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

        Lagrange = sum(inter)

    end function Lagrange

    function Newton(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        real*8, dimension(n) :: x, y
        real*8, allocatable, dimension(:)  :: inter
        real*8 :: Newton, ent

        allocate(inter(n))
        
        do i=1, n
            inter(i) = y(i)
        end do
    
        do j=1, n-1
            do i=n,j+1,-1
                inter(i) = (inter(i)-inter(i-1))/(x(i)-x(i-j))
            end do
        end do

        Newton = inter(n)
        do i=n-1,1,-1
            Newton = Newton * (ent - x(i)) + inter(i)
        end do

        deallocate(inter)

    end function Newton
    
    function GregoryNewton(x, y, n, ent)
        ! x, y - Domain and Image
        ! n    - Vectors Size
        ! ent  - Value to Interpolate
        implicit none
        integer :: n, i, j
        double precision, dimension(n) :: x, y
        double precision :: ent, u, GregoryNewton
        double precision, allocatable, dimension(:)  :: inter

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
        GregoryNewton = inter(n)

        do i = n-1, 1, -1
            GregoryNewton = GregoryNewton * (u - i + 1) / i + inter(i)
        end do

        deallocate(inter)

    end function GregoryNewton

    ! ------- Linear Systems -------

    function GaussianElimination(Matrix, m, V)
        ! Matrix - System Matrix
        ! m      - Matrix Order
        ! V      - independent values
        implicit none
        integer                 :: m, i, j, p, linha_pivo
        real*8, dimension(m,m)  :: Matrix
        real*8, dimension(m)    :: V
        real*8                  :: pivo, auxv, mult
        real*8, allocatable, dimension(:) :: aux
        real*8, dimension(m) :: GaussianElimination

        allocate(aux(m))

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

        deallocate(aux)

        GaussianElimination = V

    end function GaussianElimination

    function Tridiagonal(a, b, c, d, order)
        ! a, b ,c  - system vectors 
        ! d        - independent values
        ! order    - system order
        implicit none
        integer          :: order, i
        real*8, dimension(order) :: a, b, c, d
        real*8, allocatable, dimension(:) :: x
        real*8, dimension(order) :: Tridiagonal

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
        
        Tridiagonal = x

        deallocate(x)

    end function Tridiagonal

    ! ------- Differential Equations -------

    function Euler(a, b, h, alf)
        ! a, b - domain range
        ! h    - step
        ! alf  - initial value of Y'(x,y)
        implicit none
        integer :: n, i
        real*8 :: h, a, alf, b
        real*8, allocatable, dimension(:) :: x, y
        real*8, allocatable, dimension(:,:) :: Euler
        
        n = idint((b - a) / h) + 2

        allocate(x(n))
        allocate(y(n))
        allocate(Euler(n,2))

        x(1) = a
        y(1) = alf
        do i=1, n
            y(i+1) = y(i) + h * func(x(i), y(i))
            x(i+1) = a + i * h
        end do

        Euler(:,1) = x
        Euler(:,2) = y

        deallocate(x)
        deallocate(y)

        contains

        real*8 function func(x, y) result(r)
            implicit none
            double precision, intent(in) :: x, y
            ! Differential equation with Y'(x,y) isolated
            r = cos(x)
        end function func

    end function Euler

    function RungeKutta2(a, b, h, alf)
        ! a, b - domain range
        ! h    - step
        ! alf  - initial value of Y'(x,y)
        implicit none
        integer :: n, i
        real*8 :: h, a, alf, b, K1, K2
        real*8 , allocatable, dimension(:) :: x, y, U
        real*8 , allocatable, dimension(:,:) :: RungeKutta2

        n = idint((b - a) / h) + 2

        allocate(x(n))
        allocate(y(n))
        allocate(U(n))
        allocate(RungeKutta2(n, 2))

        x(1) = a
        y(1) = alf

        do i = 1, n-1
            k1 = func(x(i),y(i))
            x(i+1) = a + i*h
            U(i+1) = y(i) + h*k1
            k2 = func(x(i+1), U(i+1))
            y(i+1) = y(i) + (h/2)*(k1+k2)
        end do

        RungeKutta2(:,1) = x
        RungeKutta2(:,2) = y

        deallocate(x)
        deallocate(y)
        deallocate(U)
        
        contains

        real*8 function func(x, y) result(r)
            implicit none
            double precision, intent(in) :: x, y
            ! Differential equation with Y'(x,y) isolated
            r = cos(x)      
        end function func
    end function RungeKutta2

    function RungeKutta4(a, b, h, alf)
        ! a, b - domain range
        ! h    - step
        ! alf  - initial value of Y'(x,y)
        implicit none
        integer :: n, i
        real*8 :: h, a, alf, b, K1, K2, K3, K4
        real*8, allocatable, dimension(:) :: x, y, U
        real*8, allocatable, dimension(:,:) :: RungeKutta4
      
        n = idint((b - a) / h) + 2

        allocate(x(n))
        allocate(y(n))
        allocate(U(n))
        allocate(RungeKutta4(n, 2))

        x(1) = a
        y(1) = alf

        do i=1, n-1
            k1 = func(x(i), y(i))
            U(i+1) = y(i) + h*k1
            x(i+1) = a + i*h
            k2 = func(x(i)+h/2, y(i)+h*k1/2)
            k3 = func(x(i)+h/2, y(i)+h*k2/2)
            k4 = func(x(i)+h, y(i) + h*k3)
            y(i+1) = y(i) + (h/6)*(k1+2*k2+2*k3+k4)
        end do

        RungeKutta4(:,1) = x
        RungeKutta4(:,2) = y

        deallocate(x)
        deallocate(y)
        deallocate(U)
        
        contains

        real*8 function func(x, y) result(r)
            implicit none
            double precision, intent(in) :: x, y
            ! Differential equation with Y'(x,y) isolated
            ! Exemple: Y'(x) = cos(x,y)
            r = cos(x)
        end function func
    end function RungeKutta4

    function RungeKuttaSist(a, b, alf1, alf2, h)
        ! a, b - domain range
        ! h    - step
        ! alf1, alf2  - initial value of Y1'(x) e Y1'(x)
        implicit none
        integer :: i, n
        real*8 :: h, a, b, alf1, alf2
        real*8 :: xt, yt1, yt2
        real*8 :: k11, k12, k21, k22, k31, k32, k41, k42
        real*8, allocatable, dimension(:) :: x, y, y1, y2
        real*8, allocatable, dimension(:,:) :: RungeKuttaSist
      
        n = idint((b - a) / h) + 1

        allocate(x(n))
        allocate(y(n))
        allocate(y1(n))
        allocate(y2(n))
        allocate(RungeKuttaSist(n, 3))

        x(1) = a
        y1(1) = alf1
        y2(1) = alf2   
        
        do i=1, n
            xt = x(i)
            yt1 = y1(i)
            yt2 = y2(i)
            k11 = func1(xt, yt1, yt2)
            k12 = func2(xt, yt1, yt2)

            xt=x(i)+(h/2)
            yt1=y1(i)+(h*k11)/2
            yt2=y2(i)+(h*k12)/2
            k21 = func1(xt,yt1,yt2)
            k22 = func2(xt,yt1,yt2)
            xt=x(i)+(h/2)
            yt1=y1(i)+(h*k21)/2
            yt2=y2(i)+(h*k22)/2
            k31 = func1(xt, yt1, yt2)
            k32 = func2(xt, yt1, yt2) 
            xt = x(i)+h
            yt1=y1(i)+(h*k31)
            yt2=y2(i)+(h*k32)
            k41 = func1(xt, yt1, yt2)
            k42 = func2(xt, yt1, yt2)
            
            x(i+1) = a+(i*h)
            
            y1(i+1) = y1(i) + (h/6)*(k11+2*(k21+k31)+k41)
            y2(i+1) = y2(i) + (h/6)*(k12+2*(k22+k32)+k42)
            
        end do

        RungeKuttaSist(:,1) = x
        RungeKuttaSist(:,2) = y1
        RungeKuttaSist(:,3) = y2

        contains

        real*8 function func1(x, y1, y2) result(r)
            implicit none
            double precision, intent(in) :: x, y1, y2
            ! Differential equation with Y'(x,y) isolated
            ! Exemple: Y'(x,y) = (1.5*y1)-0.5*y1*y2
            r = (1.5*y1)-0.5*y1*y2
        end function func1

        real*8 function func2(x, y1, y2) result(r)
            implicit none
            double precision, intent(in) :: x, y1, y2
            ! Differential equation with Y'(x,y1,y2) isolated
            ! Exemple: Y'(x,y1,y2) = (1.5*y1)-0.5*y1*y2
            r = (-0.5)*y2+(y2*y1)
        end function func2

    end function RungeKuttaSist

    ! ------- Numerical Integration -------

    function TrapezoidalRule(a, b, h)
        ! a, b - Integration limits
        ! h    - domain range / precision
        implicit none
        integer :: ntrap
        real*8 :: area, h, a, b, jump, TrapezoidalRule
    
        ntrap = idint((b-a)/h)
        jump = a
        ! An interpolating polynomial was used for the function images
        area = func(a)
        do while(jump < b)
            area = area + (2*func(jump))
            jump = jump + h
        end do
        area = area + func(b)
    
        TrapezoidalRule = (area * (h/2))
    
        contains
    
        function func(x)
            implicit none
            real*8 :: x, func
            !f(x) function to integrate 
            func = 2*x
        end function func
    
    end function TrapezoidalRule

    function Simpson(a, b, n)
        ! a, b - Integration limits
        ! n    - Degree of integration
        implicit none
        integer :: i, n
        real*8 :: Simpson, h, h2, a, b, x(n)
        real*8 :: x1, x2, x3, x4, A1, A2, A3, A4
    
        h = (b-a)/n
        h2 = (b-a)/(n-1)
        x(1) = a
        do i = 2, n
            x(i) = x(i-1) + h2
        end do
        Simpson = 0.00
        do i=1, n
            x1 = x(i)
            x2 = x1 + h/3
            x3 = x2 + h/3
            x4 = x(i+1)
            A1 = 1.0/8.0
            A2 = 3.0/8.0
            A3 = 3.0/8.0
            A4 = 1.0/8.0 
            Simpson = Simpson + ((A1*func(x1) + A2*func(x2) + A3*func(x3) + A4*func(x4)) * h)
        end do
        
        contains

        function func(x)
            implicit none
            real*8 :: x, func
            !f(x) function to integrate 
            func = 2*x
        end function func

    end function Simpson

    function GaussianQuadrature(a, b, n)
        ! a, b - Integration limits
        ! n    - Degree of integration
        implicit none
        integer :: i, n
        real*8 :: GaussianQuadrature, h, h2, a, b, x(n)
        real*8 :: w1, w2, w3, t1, t2, t3, alpha, bet, x1, x2, x3

        h = (b-a)/n
        h2 = (b-a)/(n-1)
        x(1) = a
        do i = 2, n
            x(i) = x(i-1) + h2
        end do
        
        w1 = 5.0/9.0
        t1 = -(dsqrt(dble(3.0/5.0)))
        w2 = 8.0/9.0
        t2 = 0.0
        w3=w1
        t3 = -(t1)
        GaussianQuadrature = 0.0
        
        do i=1, n
            alpha = (x(i+1) - x(i))/2.0
            bet = (x(i+1) +x(i))/2.0
            x1 = alpha*t1+bet;
            x2 = alpha*t2+bet;
            x3 = alpha*t3+bet;
            
            GaussianQuadrature = GaussianQuadrature + ((w1*func(x1) + w2*func(x2) + w3*func(x3)) * (h/2.0))
        end do
        
        contains

        function func(x)
            implicit none
            real*8 :: x, func
            !f(x) function to integrate 
            func = 2*x
        end function func

    end function GaussianQuadrature

    function DiscreetTrapezoidal(x, y, a, b,  h, n)
        ! x, y - domain and function image
        ! a, b - Integration limits
        ! h    - domain range
        ! n    - vector size
        implicit none
        integer :: n, ntrap
        real*8, dimension(n) :: x, y
        real*8 :: h, a, b, jump, DiscreetTrapezoidal
    
        ntrap = idint((b-a)/h)
        jump = a
        ! An interpolating polynomial was used for the function images
        DiscreetTrapezoidal = GregoryNewton(x, y, n, a)
        do while(jump < b)
            DiscreetTrapezoidal = DiscreetTrapezoidal + (2 * GregoryNewton(x, y, n, jump))
            jump = jump + h
        end do
        DiscreetTrapezoidal = DiscreetTrapezoidal + GregoryNewton(x, y, n, b)
        DiscreetTrapezoidal = (DiscreetTrapezoidal * (h/2))
    
    end function DiscreetTrapezoidal

end module numerical_analysis

