Module data_pb
    USE numerics

    IMPLICIT NONE

    contains 

    FUNCTION compute_u(x,y,t) result(u_val)
        REAL(dp), INTENT(IN) :: x, y, t
        REAL(dp) :: u_val

        u_val = (x*y)**2 + 2 * t * x * y + sin(2 * pi * x * y)
    END FUNCTION compute_u

    FUNCTION f(x,y) result(f_val)
        REAL(dp), INTENT(IN) :: x, y
        REAL(dp) :: f_val, T
        T = 1.0_dp

        f_val = (x*y)**2 + 2 * T *x * y - 4 &
              + sin(2*pi*x*y) * (1.d0 + 8*pi**2 - 16*pi**4*(x*y)**2) &
              + 32*pi**3 * (x*y) * cos(2*pi*x*y) 
    END FUNCTION f

End module data_pb