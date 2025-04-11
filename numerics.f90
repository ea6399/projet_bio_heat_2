MODULE numerics

      ! This module contains the global variables declarations

      IMPLICIT NONE
      !------------------------------------------
      ! Paramètres de résolution et discrétisation
      !------------------------------------------
      integer, parameter :: dp = kind(1.d0)  ! Double précision
      integer, parameter :: M = 2**1                   ! Niveau de résolution des ondelettes (M^J points)
      integer, parameter :: N = 2 * M                  ! Nombre total de points de collocation (N = 2*M)
      real(dp) :: dt = 0.01d0                          ! Pas de temps Δt
      real(dp), parameter :: epsilon = 1.0d0           ! Coefficient de diffusion ε
      real(dp), parameter :: pi = atan(1.0d0) * 4.0d0   ! Constante π

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
      real(dp), allocatable :: hx(:,:)         ! Fonctions Haar h_i(x_c) [i, x]
      real(dp), allocatable :: hy(:,:)         ! Fonctions Haar h_i(y_c) [i, y]
      real(dp), allocatable :: p1x(:,:)       ! Intégrales premières p_{i,1}(x_c)
      real(dp), allocatable :: p1y(:,:)       ! Intégrales premières p_{i,1}(y_c) 
      real(dp), allocatable :: p2x(:,:)       ! Intégrales secondes p_{i,2}(x_c)
      real(dp), allocatable :: p2y(:,:)       ! Intégrales secondes p_{i,2}(y_c) 
      real(dp), allocatable :: u(:,:,:)       ! Inconnu u(x,y,t) [i,j,k]
      real(dp), allocatable :: v(:,:,:)       ! Inconnu v(x,y,t) [i,j,k]
      real(dp), allocatable :: v_old(:,:,:)     ! Inconnu v(x,y,t-Δt) [i,j]
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
      integer :: idx



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
        allocate(hx(N, N), hy(N,N), p1x(N, N),p1y(N,N), p2x(N, N),p2y(N,N), C(N))
        allocate(u(N, N, num_points), v(N, N, num_points), v_old(N, N, num_points))
        allocate(a_ij(N,N), b_ij(N,N))
        allocate(phi(N, N), psi(N, N), f_exact(N, N))
        allocate(A_mat(N*N, N*N), B_vec(N*N), ipiv(N*N))

        x_c = 0.0_dp
        y_c = 0.0_dp
        t_grid = 0.0_dp
        hx = 0.0_dp
        hy = 0.0_dp
        p1x = 0.0_dp
        p1y = 0.0_dp
        p2x = 0.0_dp
        p2y = 0.0_dp
        u = 0.0_dp
        v = 0.0_dp
        v_old = 0.0_dp
        a_ij = 0.0_dp
        b_ij = 0.0_dp
        phi = 0.0_dp
        psi = 0.0_dp
        f_exact = 0.0_dp
        A_mat = 0.0_dp
        B_vec = 0.0_dp
        ipiv = 0
        C = 0.0_dp

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