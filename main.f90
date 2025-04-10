PROGRAM main

    USE numerics
    USE haar_functions
    USE matrix

    IMPLICIT NONE
    ! Déclaration des variables locales
    integer :: i, j


                    ! ------------------------------ !
                    !          Début du code         !
                    ! ------------------------------ !

        

        ! This is the main program
        !Initialisation des variables
        CALL initialize_globals()
        
        write(*, '(/,T5,A,/)') " Toutes les variables globales ont été initialisées avec succès !"

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


        

      



      

    

END PROGRAM main