Module matrix
    USE numerics

    Implicit None

    Contains

    ! Subroutine pour afficher une matrice 2D
    Subroutine print_matrix(mat)
        Implicit None
        Real(dp), Intent(in) :: mat(:,:)  ! Matrice à afficher
        Integer :: rows, cols ! Dimensions de la matrice
        Integer :: i, j

        rows = size(mat, 1)
        cols = size(mat, 2)

        do i = 1, rows
        write(*, '(*(f8.2))') (mat(i, j), j = 1, cols)
        end do

        print *  ! Saut de ligne final
    end subroutine print_matrix

endmodule matrix