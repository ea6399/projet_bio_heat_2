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
              h(i,j) = haar(i, x_c(j), x_min, x_max, N)
          END DO
      END DO

      ! Affichage des ondelettes de Haar
      write(*, '(/,T5,A,/)') "Fonctions Haar :"
      call print_matrix(h)

      

    

END PROGRAM main