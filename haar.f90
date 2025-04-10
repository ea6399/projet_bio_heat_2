module haar_functions

    USE numerics

    implicit none

contains

    recursive integer function factorial(n) result(res)
        integer, intent(in) :: n
        if (n <= 1) then
            res = 1
        else
            res = n * factorial(n - 1)
        end if
    end function factorial

    integer function find_j(i, J) result(petit_j_result)
        integer, intent(in) :: i, J
        integer :: petit_j, M
        M = 2**J
        if (i > 2*M) then
            print*, "Problem: i > 2*M"
            stop
        end if
        petit_j = 0
        do while (2**petit_j < i)
            petit_j = petit_j + 1
        end do
        petit_j = petit_j - 1
        petit_j_result = 2**petit_j
    end function find_j

    real(dp) function haar(i, x, a, b, J) result(res)
        integer, intent(in) :: i, J
        real(dp), intent(in) :: x, a, b
        integer :: m, k
        res = 0.0_dp
        if (i == 1) then
            if (x >= a .and. x <= b) res = 1.0_dp
        else
            m = find_j(i, J)
            k = i - m - 1
            if (x >= a + real(k, dp)*(b - a)/real(m, dp) .and. &
               x <  a + (real(k, dp) + 0.5_dp)*(b - a)/real(m, dp)) then
                res = 1.0_dp
            else if (x >= a + (real(k, dp) + 0.5_dp)*(b - a)/real(m, dp) .and. &
                   x <= a + (real(k, dp) + 1.0_dp)*(b - a)/real(m, dp)) then
                res = -1.0_dp
            end if
        end if
    end function haar

    real(dp) function P_haar(x, i, n, a, b, J) result(res)
        real(dp), intent(in) :: x, a, b
        integer, intent(in) :: i, n, J
        real(dp) :: xi_1, xi_2, xi_3
        integer :: m, k
        real(dp), parameter :: EPSILON = 1.0e-10_dp

        res = 0.0_dp
        m = find_j(i, J)
        k = i - m - 1

        if (abs(x - a) < EPSILON) return

        if (i == 1) then
            if (x <= a) then
                res = 0.0_dp
            else if (x > a .and. x <= b) then
                res = (x - a)**n / real(factorial(n), dp)
            else
                select case(n)
                case(1)
                    res = (b - a)
                case(2)
                    res = 0.5_dp * (b - a)**2 + (b - a)*(x - b)
                case(3)
                    res = (b - a)**3 / 6.0_dp + 0.5_dp*(b - a)**2*(x - b) + 0.5_dp*(x - b)**2*(b - a)
                case(4)
                    res = (b - a)**3 /6.0_dp*(x - b) + 0.25_dp*(b - a)**2*(x - b)**2 + &
                          (x - b)**3/6.0_dp*(b - a) + (b - a)**4 /24.0_dp
                end select
            end if
        else
            xi_1 = a + real(k, dp)*(b - a)/real(m, dp)
            xi_2 = a + (real(k, dp) + 0.5_dp)*(b - a)/real(m, dp)
            xi_3 = a + (real(k, dp) + 1.0_dp)*(b - a)/real(m, dp)

            if (x <= xi_1) then
                res = 0.0_dp
            else if (x > xi_1 .and. x <= xi_2) then
                res = (x - xi_1)**n / real(factorial(n), dp)
            else if (x >= xi_2 .and. x <= xi_3) then
                res = ( (x - xi_1)**n - 2.0_dp*(x - xi_2)**n ) / real(factorial(n), dp)
            else if (x > xi_3 .and. x <= b) then
                res = ( (x - xi_1)**n - 2.0_dp*(x - xi_2)**n + (x - xi_3)**n ) / real(factorial(n), dp)
            end if
        end if
    end function P_haar

end module haar_functions