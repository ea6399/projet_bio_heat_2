MODULE numerics

      ! This module contains the global variables declarations

      IMPLICIT NONE
      !------------------------------------------
      ! Paramètres de résolution et discrétisation
      !------------------------------------------
      integer, parameter :: dp = kind(1.d0)  ! Double précision
      integer, parameter :: M = 4                      ! Niveau de résolution des ondelettes (2^M points)
      integer, parameter :: N = 2 ** M                  ! Nombre total de points de collocation (N = 2*M)
      real(dp) :: dt = 0.01d0                                   ! Pas de temps Δt
      real(dp), parameter :: epsilon = 1.0d0            ! Coefficient de diffusion ε

      !------------------------------------------
      ! Domaine spatial et temporel
      !------------------------------------------
      real(dp), parameter :: x_min = 0.0d0, x_max = 1.0d0  ! Domaine en x
      real(dp), parameter :: y_min = 0.0d0, y_max = 1.0d0  ! Domaine en y
      real(dp), parameter :: T_total = 1.0d0               ! Temps final

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



      contains


       subroutine initialize_globals()
        implicit none
        integer :: i, r
        real(dp) :: dx, dy
        integer :: num_points

        ! Initialisation
        num_points = ceiling(T_total/dt) + 1


        ! Allocation des tableaux
        allocate(x_c(N), y_c(N), t_grid(num_points))
        allocate(h(N, N), p1(N, N), p2(N, N), C(N))
        allocate(a_ij(N, N), b_ij(N, N))
        allocate(phi(N, N), psi(N, N), f_exact(N, N))
        allocate(A_mat(N*N, N*N), B_vec(N*N), ipiv(N*N))

        ! Initialisation des points de collocation
        dx = (x_max - x_min)/N
        dy = (y_max - y_min)/N
        do r = 1, N
            x_c(r) = x_min + (r - 0.5d0)*dx
            y_c(r) = y_min + (r - 0.5d0)*dy
        end do

        ! Initialisation de la grille temporelle
        do i = 0, num_points - 1
            t_grid(i+1) = i*dt
        end do

    end subroutine initialize_globals


END MODULE numerics