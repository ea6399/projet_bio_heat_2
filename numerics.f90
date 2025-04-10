MODULE numerics

      ! This module contains the global variables declarations

      IMPLICIT NONE
      !------------------------------------------
      ! Paramètres de résolution et discrétisation
      !------------------------------------------
      integer, parameter :: dp = kind(1.d0)  ! Double précision
      integer :: M = 16                      ! Niveau de résolution des ondelettes (2^M points)
      integer :: N = 2 * M                   ! Nombre total de points de collocation (N = 2*M)
      real(dp) :: dt = 0.01d0                ! Pas de temps Δt
      real(dp) :: epsilon = 1.0d0            ! Coefficient de diffusion ε

      !------------------------------------------
      ! Domaine spatial et temporel
      !------------------------------------------
      real(dp) :: x_min = 0.0d0, x_max = 1.0d0  ! Domaine en x
      real(dp) :: y_min = 0.0d0, y_max = 1.0d0  ! Domaine en y
      real(dp) :: T_total = 1.0d0               ! Temps final

      !------------------------------------------
      ! Points de collocation et maillages
      !------------------------------------------
      real(dp), allocatable :: x_c(:), y_c(:)   ! Points de collocation (x_r = (r-0.5)/N)
      real(dp), allocatable :: t_grid(:)        ! Grille temporelle

      !------------------------------------------
      ! Fonctions et coefficients des ondelettes de Haar
      !------------------------------------------
      real(dp), allocatable :: h(:,:)         ! Fonctions Haar h_i(x_c) [i, x]
      real(dp), allocatable :: p1(:,:)        ! Intégrales premières p_{i,1}(x_c)
      real(dp), allocatable :: p2(:,:)        ! Intégrales secondes p_{i,2}(x_c)
      real(dp), allocatable :: C(:)           ! Constantes d'intégration C_i = 1/(4m²)

      real(dp), allocatable :: a_ij(:,:)      ! Coefficients HWCM1 pour v_{xxyy} [i,j]
      real(dp), allocatable :: b_ij(:,:)      ! Coefficients HWCM2 pour u_{xxyy} [i,j]

      !------------------------------------------
      ! Conditions aux limites et source exacte
      !------------------------------------------
      real(dp), allocatable :: phi(:,:)       ! Condition initiale φ(x,y)
      real(dp), allocatable :: psi(:,:)       ! Condition finale ψ(x,y)
      real(dp), allocatable :: f_exact(:,:)   ! Source exacte f(x,y)

      !------------------------------------------
      ! Paramètres de bruit et randomisation
      !------------------------------------------
      real(dp) :: sigma = 0.0d0              ! Intensité du bruit σ (%)
      integer :: seed = 123456789            ! Graine pour le générateur aléatoire

      !------------------------------------------
      ! Matrices système et vecteurs de travail
      !------------------------------------------
      real(dp), allocatable :: A_mat(:,:)    ! Matrice du système linéaire
      real(dp), allocatable :: B_vec(:)      ! Vecteur second membre
      integer, allocatable :: ipiv(:)        ! Pivots pour la résolution (LAPACK)


END MODULE numerics