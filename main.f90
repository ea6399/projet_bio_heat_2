PROGRAM main

    USE numerics
    USE haar_functions
    USE matrix
    USE data_pb

    IMPLICIT NONE
    ! Déclaration des variables locales
    integer :: i, j


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

        WRITE(*, '(/,T25,A,/)') "Calcul des inconnus a_ij :"

        ! do r = 1, N
        !     do s = 1, N
        !         idx = (r-1)*N + s
        !         do i = 1, N
        !             do j = 1, N
        !                 A_mat(idx, (i-1)*N + j) = -epsilon * dt * h(i, r) * h(j, s)
        !             end do
        !         end do
        !         A_mat(idx, idx) = A_mat(idx, idx) + 1.0_dp  ! Terme diagonal
        !         B_vec(idx) = v_old(r, s)  
        !     end do
        ! end do




        


























        contains

        function v_boundary_x0(x, y, t) result(val)
            real(dp), intent(in) :: x, y, t
            real(dp) :: val
            ! y = 0 donc h_1(x,t)
            val = 0.0_dp  
        end function v_boundary_x0

        function v_boundary_0y(x, y, t) result(val)
            real(dp), intent(in) :: x, y, t
            real(dp) :: val
            ! x = 0 donc g_1(y,t)
            val = 0.0_dp  
        end function v_boundary_0y

        function closest_index(x, grid) result(idx)
            real(dp), intent(in) :: x
            real(dp), intent(in) :: grid(:)
            integer :: idx
            idx = minloc(abs(grid - x), dim=1)
        end function closest_index

        function compute_v(x, y, t) result(v_val)
            use numerics
            implicit none
            real(dp), intent(in) :: x, y, t  ! Position et temps actuels
            real(dp) :: v_val                ! Valeur de v(x,y,t)
            
            ! Variables locales
            integer :: i, j, idx_x, idx_y
            real(dp) :: delta_x, delta_y, dx, dy
            real(dp) :: term1, term2, term3, term4, term5, wavelet_part
            
            !------------------------------------------
            ! 1. Calcul des pas spatiaux Δx et Δy
            !------------------------------------------
            ! dx = (x_max - x_min)/N
            ! dy = (y_max - y_min)/N
            delta_x = dx
            delta_y = dy
            
            !------------------------------------------
            ! 2. Évaluation des termes aux limites (exemples)
            !------------------------------------------
            ! Supposons que v(x,0,t), v(0,y,t), etc. sont stockés dans des tableaux
            ! ou calculés via des fonctions externes (à adapter selon votre cas)
            
            ! Indices des points voisins pour les différences finies
            idx_x = min(floor(x/dx) + 1, N-1)
            idx_y = min(floor(y/dy) + 1, N-1)
            
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
            
            !------------------------------------------
            ! 3. Partie ondelette : Σ a_ij p2_i(x) p2_j(y)
            !------------------------------------------
            wavelet_part = 0.0_dp
            do i = 1, N
                do j = 1, N
                    wavelet_part = wavelet_part + a_ij(i,j) * p2x(i, closest_index(x, x_c)) * p2y(j, closest_index(y, y_c))
                end do
            end do
            
            !----------------
            ! 4. Somme finale 
            !----------------
            v_val = term1 + term2 + term3 + term4 + term5 + wavelet_part
            
        end function compute_v 

END PROGRAM main