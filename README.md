# Numerical Analysis Fortran
Repository designed to maintain fortran scripts for numerical analysis.

The file numerical_analysis.F95, contains the following functions:

Some functions must be edited to include the equations necessary for the analysis, such as the differential equations that must be included in the functions before executing them.

- Linear Systems

    - Gaussian Elimination

        ```
            GaussianElimination(Matrix, m, V)

        ! Matrix - System Matrix (real*8, dimension(m,m))
        ! m      - Matrix Order (integer)
        ! V      - independent values (real*8, dimension(m))
        ! Returns a real*8 vector of size m
        ```

    - Tridiagonal Systems

        ```
            Tridiagonal(a, b, c, d, order)

        ! a, b ,c  - system vectors (real*8, dimension(order))
        ! d        - independent values (real*8, dimension(order))
        ! order    - system order (integer)
        ! Returns a real*8 vector of size order
        ```


- Interpolation Methods

    - Newton Interpolator Polynomial

        ```
            Newton(x, y, n, ent)

        ! x, y - Domain and Image (real*8, dimension(n))
        ! n    - Vectors Size (integer)
        ! ent  - Value to Interpolate (real*8)
        ! Returns real*8 value
        ```

    - Lagrange Interpolator Polynomial

        ```
            Lagrange(x, y, n, ent)

        ! x, y - Domain and Image (real*8, dimension(n))
        ! n    - Vectors Size (integer)
        ! ent  - Value to Interpolate (real*8)
        ! Returns real*8 value
        ```

    - Gregory Newton Interpolator Polynomial

        ```
            GregoryNewton(x, y, n, ent)

        ! x, y - Domain and Image (real*8, dimension(n))
        ! n    - Vectors Size (integer)
        ! ent  - Value to Interpolate (real*8)
        ! Returns real*8 value
        ```

- Differential Equations
    
    - Euler Method 

        ```
            Euler(a, b, h, alf)

        ! a, b - domain range (real*8)
        ! h    - step (real*8)
        ! alf  - initial value of Y'(x,y) (real*8)
        ! Returns an array [ real*8 array(n+2,2) ]
        ```
    - Second-order Runge Kutta


        ```
            RungeKutta2(a, b, h, alf)

        ! a, b - domain range (real*8)
        ! h    - step (real*8)
        ! alf  - initial value of Y'(x,y) (real*8)
        ! Returns an array [ real*8 array(n+2,2) ]
        ```
    
    - Fourth-order Runge Kutta

        ```
            RungeKutta4(a, b, h, alf)

        ! a, b - domain range (real*8)
        ! h    - step (real*8)
        ! alf  - initial value of Y'(x,y) (real*8)
        ! Returns an array [ real*8 array(n+2,2) ]
        ```
    
    - Runge Kutta method for systems of differential equations

        ```
            RungeKuttaSist(a, b, alf1, alf2, h)
            
        ! a, b - domain range (real*8)
        ! h    - step (real*8)
        ! alf1, alf2  - initial value of Y1'(x) e Y1'(x) (real*8)
        ! Returns an array [ real*8 array(n+1,3) ]
        ```

- Numerical Integration

    - Trapezoid Rule

        ```
            TrapezoidalRule(a, b, h)

        ! a, b - Integration limits (real*8)
        ! h    - domain range / precision (real*8)
        ! Returns real*8 value
        ```

    - Discreet Trapezoid Rule

        ```
            DiscreetTrapezoidal(x, y, a, b,  h, n)

        ! x, y - domain and function image (real*8, dimension(n) )
        ! a, b - Integration limits (real*8)
        ! h    - domain range (real*8)
        ! n    - vector size (integer)
        ! Returns real*8 value
        ```

    - Gaussian Quadrature

        ```
            GaussianQuadrature(a, b, n)

        ! a, b - Integration limits (real*8)
        ! n    - Degree of integration (integer)
        ! Returns real*8 value
        ```

    - Simpson method (3/8 Rule)

        ```
            Simpson(a, b, n)

        ! a, b - Integration limits (real*8)
        ! n    - Degree of integration (integer)
        ! Returns real*8 value
        ```

To compile run the following command:

```
gfortran main.F95 numerical_analysis.F95
```

The main.F95 file, for example, is the file that makes use of the functions.
