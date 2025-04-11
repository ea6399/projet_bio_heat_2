PROGRAM main

    USE numerics
    USE haar_functions
    USE matrix
    USE data_pb

    IMPLICIT NONE
    ! Déclaration des variables locales
    integer :: i, j, r, s


                    ! ------------------------------ !
                    !          Début du code         !
                    ! ------------------------------ !

        

        ! This is the main program
        !Initialisation des variables
        CALL initialize_globals()
        
        write(*, '(/,T25,A,/)') " Toutes les variables globales ont été initialisées avec succès !"

        ! Affichage de points de collocation
        ! write(*, '(/,T5,A,/)') "Points de collocation :"
        DO i = 1, N
            write(*, '(1X, F10.5, 1X, F10.5)') x_c(i), y_c(i)
        END DO

        ! Initialisation des ondelettes de Haar
        DO j = 1, N
            DO i = 1, N
                hx(i,j) = haar(i, x_c(j), x_min, x_max, N)
                hy(i,j) = haar(i, y_c(j), y_min, y_max, N)
            END DO
        END DO

        ! Affichage des ondelettes de Haar
        write(*, '(/,T5,A,/)') "Fonctions Haar calcule :"
        call print_matrix(hx)
        call print_matrix(hy)

        ! Affichage des intégrales premières
        write(*, '(/,T5,A,/)') "Intégrales premières calcules :"
        DO j = 1, N
            DO i = 1, N
                p1x(i,j) = P_haar(x_c(j), i, 1, x_min, x_max, N)
                p1y(i,j) = P_haar(y_c(j), i, 1, y_min, y_max, N)
            END DO
        END DO
        call print_matrix(p1x)
        call print_matrix(p1y)

        ! Affichage des intégrales secondes
        write(*, '(/,T5,A,/)') "Intégrales secondes calculees :"
        DO j = 1, N
            DO i = 1, N
                    p2x(i,j) = P_haar(x_c(j), i, 2, x_min, x_max, N)
                    p2y(i,j) = P_haar(y_c(j), i, 2, y_min, y_max, N)
            END DO
        END DO
        call print_matrix(p2x)
        call print_matrix(p2y)



                              ! -------------------------------- !
                              !   Calcul des inconnus (aij)      !
                              ! -------------------------------- !

        ! Initialisation de V
        DO s = 1, N
            DO r = 1, N
                v_old(r,s,0) = initialize_v(x_c(r), y_c(s))
            END DO
        END DO

        write(*, '(/,T25,A,/)') "V initialisé."

        ! Calcul des coefficients a_ij

        

        WRITE(*, '(/,T25,A,/)') "Calcul des inconnus a_ij calculés"

        ! Résolution du système A * a = b




        contains

        function initialize_v(x,y) result(val)
            use numerics
            implicit none
            real(dp) :: val
            real(dp), intent(in) :: x, y
            
            val = sin(pi*x) * sin(pi*y)

        end function initialize_v

        function v_boundary_x0(x, y, t) result(val)
            real(dp), intent(in) :: x, y, t
            real(dp) :: val
            ! Dirichlet,  x = 0 donc g_1(x,t)
            val = 0.0_dp  
        end function v_boundary_x0

        function v_boundary_0y(x, y, t) result(val)
            real(dp), intent(in) :: x, y, t
            real(dp) :: val
            ! x = 0 donc g_1(y,t)
            val = 0.0_dp  
        end function v_boundary_0y

        function compute_v(x, y, t) result(v_val)
            use numerics
            implicit none
            real(dp), intent(in) :: x, y, t  ! Position et temps actuels
            real(dp) :: v_val                ! Valeur de v(x,y,t)
            
            ! Variables locales
            integer :: i, j
            real(dp) :: delta_x, delta_y, dx, dy
            real(dp) :: term1, term2, term3, term4, term5, wavelet_part
            
            !------------------------------------------
            ! 1. Calcul des pas spatiaux deltax et deltay
            !------------------------------------------
            ! dx = (x_max - x_min)/N
            ! dy = (y_max - y_min)/N
            delta_x = dx
            delta_y = dy

            
            !------------------------------------------
            ! 2. Évaluation des termes aux limites 
            !------------------------------------------
            
            ! Terme 1 : v(x,0,t) + y*(v(x,Δy,t) - v(x,0,t))/Δy
            term1 = v_boundary_x0(x, 0.0_dp, t) + y * (v_boundary_x0(x, delta_y, t) - v_boundary_x0(x, 0.0_dp, t)) / delta_y
            
            ! Terme 2 : v(0,y,t) - v(0,0,t)
            term2 = v_boundary_0y(0.0_dp, y, t) - v_boundary_0y(0.0_dp, 0.0_dp, t)
            
            ! Terme 3 : x*(v(Δx,y,t) - v(0,y,t))/Δx - x*(v(Δx,0,t) - v(0,0,t))/Δx
            term3 = x * (v_boundary_0y(delta_x, y, t) - v_boundary_0y(0.0_dp, y, t)) / delta_x &
                - x * (v_boundary_0y(delta_x, 0.0_dp, t) - v_boundary_0y(0.0_dp, 0.0_dp, t)) / delta_x
            
            ! Terme 4 : -y*(v(0,Δy,t) - v(0,0,t))/Δy
            term4 = -y * (v_boundary_0y(0.0_dp, delta_y, t) - v_boundary_0y(0.0_dp, 0.0_dp, t)) / delta_y
            
            ! Terme 5 : xy*(v(Δx,Δy,t) - v(Δx,0,t) - v(0,Δy,t) + v(0,0,t))/(ΔxΔy)
            term5 = x * y * (v_boundary_0y(delta_x, delta_y, t) - v_boundary_0y(delta_x, 0.0_dp, t) &
                        - v_boundary_0y(0.0_dp, delta_y, t) + v_boundary_0y(0.0_dp, 0.0_dp, t)) / (delta_x * delta_y)
            
            !---------------------
            ! 3. Partie ondelette :
            !---------------------
            wavelet_part = 0.0_dp
            do i = 1, N
                do j = 1, N
                    wavelet_part = wavelet_part + a_ij(i,j) * p2x(i, j) * p2y(i, j)
                end do
            end do
            
            !----------------
            ! 4. Somme finale 
            !----------------
            v_val = term1 + term2 + term3 + term4 + term5 + wavelet_part
            
        end function compute_v 

END PROGRAM main