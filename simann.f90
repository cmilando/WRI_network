!--------------------------------------------------------------------------
! Author: CWM
! Date started: 6.21.2024
! Purpose: Subroutines for simulated annealing
!
! NOTES: 
! - This only works for S as an integer, and with no other confounders
! - All subroutines need to be in this file
! - For doing in parallel, you may need to add THREADPRIVATE somewhere
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Subroutine for swapping discrete neighbors
! 
! Note: This does not need to change for a new application
!--------------------------------------------------------------------------
subroutine get_neighbor(x0, x1, l)
    
    implicit none

    integer :: l                ! length of the index vector
    real    :: p0, p1           ! random numbers, [0, 1)
    integer :: p0_int, p1_int   ! converted to fractions of l
    integer :: x0(l), x1(l)     ! the old and new index choices

    ! Initialize x1 to be x0
    x1 = x0
    
    ! -- SWAP STEP -----------------------
    ! now swap two indices
    ! random_number is [0, 1)
    ! FORTRAN arrays start at 1
    ! because this will give you [0,1] you need to adjust l by 1
    call random_number(p0)
    p0_int = floor(p0 * real(l)) + 1

 10 call random_number(p1)
    p1_int = floor(p1 * real(l)) + 1
    ! these cannot be the same so, try again if so
    ! Note: This is how you ensure you have all different indices
    if (x0(p1_int) .eq. x0(p0_int) ) go to 10
    ! -- END SWAP ------------------------
    
    ! Note: if you wanted to introduce other options
    ! for joint annealing, here is where you would do it
    ! and the steps themselves would have different probabilities
    ! which likely would be adjusted on the fly
    ! so that the closer you got to the end, the more likely
    ! swaps were and the less likely you were to add in or remove neighbors
    
    ! Note2: if you do this though, you need to update the SCORE to have
    ! a penalty for the additional stations
    
    ! -- UPDATE THE NEIGHBORS -- 
    x1(p0_int) = x0(p1_int)
    x1(p1_int) = x0(p0_int)

end subroutine get_neighbor

!--------------------------------------------------------------------------
! quantile function
!--------------------------------------------------------------------------
subroutine get_metric(x, nrow, ncol, q, xout)
  implicit none
  
  integer :: nrow, ncol
  real(kind = 8) :: x(nrow, ncol)
  real(kind = 8) :: q
  real(kind = 8) :: xout(ncol)

  real(kind = 8) :: col(nrow)
  integer :: i, j, k, idx_low, idx_high
  real(kind = 8) :: pos, frac, temp
  
  do j = 1, ncol
  
     ! copy column j into a working array
     do i = 1, nrow
        col(i) = x(i, j)
     end do

     ! sort the column in ascending order (simple insertion sort for clarity)
     do i = 2, nrow
        temp = col(i)
        k = i - 1
        do while (k >= 1 .and. col(k) > temp)
           col(k+1) = col(k)
           k = k - 1
        end do
        col(k+1) = temp
     end do

     ! compute quantile position (R type 7: (n-1)*q + 1)
     pos = (nrow - 1) * q + 1
     idx_low = int(floor(pos))
     idx_high = int(ceiling(pos))
     frac = pos - idx_low

     if (idx_high > nrow) idx_high = nrow

     if (idx_low == idx_high) then
        xout(j) = col(idx_low)
     else
        xout(j) =  col(idx_low) + frac * (col(idx_high) - col(idx_low))
     end if
  end do
  
end subroutine get_metric

!--------------------------------------------------------------------------
! rmse
!--------------------------------------------------------------------------
subroutine get_rmse(a, b, n, rmse)

  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: a(n), b(n)
  double precision, intent(out) :: rmse

  integer :: i
  double precision :: diff, sumsq

  sumsq = 0.0d0
  do i = 1, n
     diff = a(i) - b(i)
     sumsq = sumsq + diff * diff
  end do

  rmse = sqrt ( sumsq / dble(n) )
  
end subroutine get_rmse

!--------------------------------------------------------------------------
! OLS SVD
! 
! The reason I chose this over other options is because I suspect that 
! columns of X will be collinear
!--------------------------------------------------------------------------

subroutine ols_svd(X, Y, n, p, rcond, beta, rank, info)
  !-----------------------------------------------------------------------
  ! Solve the least squares problem:  minimize ||X*beta - Y||_2
  !
  ! Equivalent to:  beta = MASS::ginv(t(X) %*% X) %*% t(X) %*% Y
  ! but computed via LAPACK's DGELSD routine, which uses SVD internally.
  !
  ! This approach is more numerically stable and can handle
  ! rank-deficient X matrices by using a tolerance (rcond).
  !
  ! Inputs:
  !   X(n, p)   : design matrix (n observations x p predictors)
  !   Y(n)      : response vector (length n)
  !   n         : number of rows
  !   p         : number of columns
  !   rcond     : cutoff for small singular values (e.g., 1e-12)
  !
  ! Outputs:
  !   beta(p)   : estimated coefficients
  !   rank      : effective rank of X (after applying rcond)
  !   info      : 0 if successful; >0 means convergence failure
  !-----------------------------------------------------------------------
  implicit none
  integer, intent(in)    :: n, p
  real(kind=8), intent(in)  :: X(n, p), Y(n), rcond
  real(kind=8), intent(out) :: beta(p)
  integer, intent(out)   :: rank, info

  ! Local variables
  integer :: m, nrhs, lda, ldb, lwork, liwork
  real(kind=8), allocatable :: A(:,:), B(:,:), S(:), work(:)
  integer, allocatable :: iwork(:)

  ! Declare LAPACK routine
  external dgelsd

  ! Set dimensions
  m    = n          ! Number of rows in X
  nrhs = 1          ! Number of RHS vectors (Y is a single column)
  lda  = max(1, m)  ! Leading dimension of A
  ldb  = max(1, max(m, p))  ! Leading dimension of B (must be >= max(m,p))

  ! Allocate working arrays
  allocate(A(lda, p))
  allocate(B(ldb, nrhs))
  allocate(S(min(m, p)))   ! Singular values

  ! Copy inputs because DGELSD overwrites A and B
  A = X
  B(1:m,1) = Y
  if (ldb > m) B(m+1:ldb,1) = 0.0d0   ! Pad with zeros if necessary

  ! Workspace query (first call with lwork = -1 and liwork = -1)
  ! to get the optimal work array sizes for DGELSD.
  lwork  = -1
  liwork = -1
  allocate(work(1))
  allocate(iwork(1))

  call dgelsd(m, p, nrhs, A, lda, B, ldb, S, rcond, rank, &
     &        work, lwork, iwork, liwork, info)

  ! Optimal sizes are returned in work(1) and iwork(1)
  lwork  = int(work(1))
  liwork = iwork(1)
  
  ! write(*,*) info
  ! write(*,*) lwork
  ! write(*,*) liwork

  deallocate(work)
  deallocate(iwork)
  allocate(work(lwork))
  allocate(iwork(liwork))

  ! Actual solve step:
  ! DGELSD overwrites B with the solution in its first p entries.
  call dgelsd(m, p, nrhs, A, lda, B, ldb, S, rcond, rank, & 
   &          work, lwork, iwork, liwork, info)
  
  ! write(*,*) B(1:p, 1)
  
  ! Extract beta
  ! if (info == 0) then
  if (info >= 0) then
     beta = B(1:p, 1)
  else
     beta = 0.0d0
  end if
  
  ! write(*,*) beta
  
  ! Clean up
  deallocate(A, B, S, work, iwork)
  
end subroutine ols_svd

!--------------------------------------------------------------------------
! Cost function for SCORE
! 
! this is what needs to change for a new application
!--------------------------------------------------------------------------
subroutine get_score(S, df_wide, magic_n, nsites, ndays, q, score_cols, & 
 &                   X_matrix, n_predictors, ID_vector, Y_matrix, lambda_vec,  & 
 &                   SCORE_z1, SCORE_z2, SCORE)

    implicit none
    
    ! Inputs / Outputs
    integer, intent(in)        :: magic_n                  ! current value of n subset
    integer, intent(in)        :: nsites                   ! total number of sites
    integer, intent(in)        :: ndays                    ! total number of days
    
    integer, intent(in)        :: score_cols               ! number of columns in the score matrix, so 3 (med, min, max)
    real(kind=8), intent(in)   :: q(score_cols)
    integer                    :: S(nsites)                ! S, a specific iteration of the magic indicator
    
    real(kind=8), intent(in)   :: df_wide(nsites, ndays)
    
    integer, intent(in)        :: n_predictors             ! total number of predictors in X

    real(kind=8), intent(in)   :: X_matrix(nsites * ndays, n_predictors)
    integer, intent(in)        :: ID_vector(nsites * ndays)
    real(kind=8), intent(in)   :: Y_matrix(nsites * ndays)
    real(kind=8), intent(in)   :: lambda_vec(2)
    
    ! Target
    real(kind=8)  :: SCORE                                 ! Score to be calculated
    real(kind=8)  :: SCORE_z1                              ! Component of z1
    real(kind=8)  :: SCORE_z2                              ! Component of z2
    
    ! Local variables
    integer       :: i, s_count
    integer       :: S_ones(magic_n)
    
    real(kind=8)  :: df_sub(magic_n, ndays)
    real(kind=8)  :: metric_matrix(ndays, score_cols)
    real(kind=8)  :: best_matrix(ndays, score_cols)
    real(kind=8)  :: tmp_xout(ndays)
    real(kind=8)  :: tmp_rmse
    
    integer       :: total_rows
    integer       :: m_sub_rows(magic_n * ndays)
    real(kind=8)  :: X_sub(magic_n * ndays, n_predictors)
    real(kind=8)  :: Y_sub(magic_n * ndays)
    real(kind=8)  :: tmp_beta(n_predictors)
    real(kind=8)  :: rcond
    integer       :: rank, info
    real(kind=8)  :: pred(nsites * ndays)
    
    ! Initialize
    SCORE    = 0.0
    SCORE_z1 = 0.0
    SCORE_z2 = 0.0

    s_count = 0
    total_rows = nsites * ndays
    
    ! Loop through S to find indices where S(i) == 1
    ! keep this
    do i = 1, nsites
        if (S(i) == 1) then
            s_count = s_count + 1
            S_ones(s_count) = i
        end if
    end do
    
    ! -----------------
    ! --- PART 1 ------
    ! -----------------
    
    if (lambda_vec(1) > 0.0) then
      
      ! filter the data frame
      ! df_sub <- df_wide[S_ones, ]
      do i = 1, magic_n
          df_sub(i, :) = df_wide(S_ones(i), :)
      end do
      
      ! get just the points for this sample
      ! subroutine get_metric(x, nrow, ncol, q, xout)
      ! yy <- get_metrics(df_sub)
      do i = 1, score_cols
      
          ! best
          tmp_xout = 0.0
          call get_metric(df_wide, nsites, ndays, q(i), tmp_xout)
          best_matrix(:, i) = tmp_xout
          
          ! subset
          tmp_xout = 0.0
          call get_metric(df_sub, magic_n, ndays, q(i), tmp_xout)
          metric_matrix(:, i) = tmp_xout
          
      end do
      
      ! get day-wise errors
      ! zz = vector("numeric", ncol(df_best))
      ! for(i in 1:(ncol(df_best))) {
      !  a = unlist(unname(as.vector(df_best[, i])))
      !  b = unlist(unname(as.vector(yy[, i])))
      !  zz[i] = get_rmse(a, b)
      ! }
      ! get_rmse(a, b, n, rmse)
      do i = 1, score_cols
          tmp_rmse = 0.0
          call get_rmse(best_matrix(:, i), metric_matrix(:, i), ndays, tmp_rmse)
          SCORE_z1 = SCORE_z1 + tmp_rmse
      end do
      
    else 
      
      SCORE_z1 = 0.0
    
    end if
    
    ! -----------------
    ! --- PART 2 ------
    ! -----------------
    
    if (lambda_vec(2) > 0.0) then
    
      ! m_sub_rows <- which(ID_vector %in% S_ones)
      ! Y_sub <- matrix(Y_matrix[m_sub_rows, ], ncol = 1)
      ! X_matrix_sub <- X_matrix[m_sub_rows, ]
      s_count = 0
      do i = 1, total_rows
          if (any(S_ones == ID_vector(i))) then
              s_count = s_count + 1
              m_sub_rows(s_count) = ID_vector(i)
              Y_sub(s_count) = Y_matrix(i)
              X_sub(s_count, :) = X_matrix(i, :)
          end if
      end do
      
      !# classic OLM invervsion lets go
      !# β = (XTX)−1XTy
      !beta_vector <- MASS::ginv(t(X_matrix_sub) %*% X_matrix_sub) %*% 
      !  t(X_matrix_sub) %*% Y_sub
      !ols_svd(X, Y, n, p, rcond, beta, rank, info)
      tmp_beta = 0.0
      rcond = 1e-6
      rank = 1
      info = 1
      
      ! *** NEED TO FIX THIS BECAUSE I TOOK A SHORTCUT FOR THE INFO PARAMETER
      call ols_svd(X_sub, Y_sub, magic_n * ndays, n_predictors, rcond,  &
       &           tmp_beta, rank, info)
      ! ********
       
      ! write(*,*) tmp_beta
    
      !# now predict on the full set
      !pred_sub <- X_matrix %*% beta_vector
      pred = MATMUL(X_matrix, tmp_beta)
      
      !# and, you guessed it get_rmse
      !z2 <- get_rmse(Y_matrix, pred_sub)
      SCORE_z2 = 0.0
      call get_rmse(Y_matrix, pred, total_rows, SCORE_z2)
    
    else
    
      SCORE_z2 = 0.0
    
    end if
    
    ! -----------------
    ! --- PART 3 ------
    ! -----------------
    
    SCORE = lambda_vec(1) * SCORE_z1 + lambda_vec(2) * SCORE_z2
    
end subroutine get_score

! ------------------------------------------------------------------------------------
! Simulated Annealing algorithm
!
! informed by sci-kits sko SA.py
!
! This only needs to change in a few places for a new application:
! - The inputs and Outputs
! - and then the two places where get_score is calculated
! - and the subroutine call obvi
! ------------------------------------------------------------------------------------
subroutine simann(S, df_wide, magic_n, nsites, ndays, q, score_cols, & 
      & X_matrix, n_predictors, ID_vector, Y_matrix, lambda_vec, &
      & SCORE_z1, SCORE_z2, SCORE, cooling_rate, verbose)
    
    implicit none
    
    ! Inputs / Outputs that need to change for new applications
    integer                    :: magic_n     ! current value of n subset
    integer, intent(in)        :: nsites      ! total number of sites
    integer, intent(in)        :: ndays       ! total number of days
    integer, intent(in)        :: score_cols  ! number of columns in the score matrix, so 3 (med, min, max)
    real(kind=8), intent(in)   :: q(score_cols)
    integer                    :: S(nsites)   ! S, a specific iteration of the magic indicator
    real(kind=8), intent(in)   :: df_wide(nsites, ndays)
    
    
    integer, intent(in)        :: n_predictors             ! total number of predictors in X
    real(kind=8), intent(in)   :: X_matrix(nsites * ndays, n_predictors)
    integer, intent(in)        :: ID_vector(nsites * ndays)
    real(kind=8), intent(in)   :: Y_matrix(nsites * ndays)
    real(kind=8), intent(in)   :: lambda_vec(2)
    
    ! Inputs / Outputs that don't need to change for new applications
    integer              :: Sprev(nsites)                !  , a prevoious iteration of the magic indicator
    integer              :: Sbest(nsites)                !  , the best

    real(kind=8)         :: SCORE                        !
    real(kind=8)         :: SCORE_z1                              ! Component of z1
    real(kind=8)         :: SCORE_z2                              ! Component of z2
    real(kind=8)         :: SCOREprev                    !  
    real(kind=8)         :: SCOREbest                    !  
    real(kind=8)         :: SCOREbest1, SCOREbest2       ! generations of SCOREbest
    
    real(kind=8)         :: SCORE_z1best                    !
    real(kind=8)         :: SCORE_z2best                    !
    
    real(kind=8)         :: df             ! comparison of SCOREs
    real(kind=8)         :: rv             ! random variable to check

    integer :: verbose
    real(kind=8) :: rel_tol, abs_tol

    ! Simulated annealing parmeters
    real(kind = 8) :: Temp_curr  ! 
    real(kind = 8) :: Temp_max   ! 
    real(kind = 8) :: Temp_min   ! 
    real(kind = 8) :: cooling_rate 
    real(kind = 8) :: curr_rel_tol
    real(kind = 8) :: abs_score_diff
    integer        :: LoC
    integer        :: stay_counter
    integer        :: max_stay_counter
    integer        :: cycle_i
    integer        :: iter_i  ! long of chain, may_stay_counter

    ! --------------------------------------------------------
    ! Initialize
    call random_seed()
    rel_tol          = 1.0e-2
    abs_tol          = 1.0e-3
    LoC              = 500
    max_stay_counter = 500
    stay_counter     = 0
    cycle_i          = 0
    SCOREbest        = 0.
    SCORE_z1best     = 0.
    SCORE_z2best     = 0.
    Temp_max         = 500000.
    Temp_min         = 0.0000000001
    Temp_curr        = Temp_max

    ! **************
    ! INTIAL VALUES
    
    ! //////////////////////////////////////////////
    ! >>>>>>> this call changes with get_score <<<<<<<<
    call get_score(S, df_wide, magic_n, nsites, ndays, q, score_cols, & 
    &  X_matrix, n_predictors, ID_vector, Y_matrix, lambda_vec,  & 
    &                   SCORE_z1, SCORE_z2,SCORE)
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! //////////////////////////////////////////////
    
    Sprev   = S
    Sbest   = S
    SCOREprev  = SCORE
    SCOREbest  = SCORE
    SCOREbest1 = SCORE
    SCOREbest2 = SCORE
    ! **************

    if(verbose .eq. 1) write(*, *) "Started"

    ! Simulated Annealing Loop
    ! converted to FORTRAN from python in scikit-opt on 10.4.2023 by CWM
    do 
        
        ! Update SCOREbest1 and SCOREbest2
        SCOREbest2 = SCOREbest1
        SCOREbest1 = SCORE

        ! -----------------
        ! ** PHASE 1 **
        ! -----------------
        ! Loop for all L in a specific chain
        do iter_i = 1, LoC

            !if(verbose .eq. 1 .and. mod(iter_i, 100) .eq. 0) write(*, *) iter_i
            
            ! get a new x. so S becomes a new version of Sprev
            ! so x_current = Sprev
            ! and x_new = S
            call get_neighbor(Sprev, S, nsites)

            ! calculate a new y. so this updates the value of SCORE
            ! similarly, y_current = SCOREprev
            ! and y_new = SCORE
            
            ! //////////////////////////////////////////////
            ! >>>>>>> this call changes with get_score <<<<<<<<
            call get_score(S, df_wide, magic_n, nsites, ndays, q, score_cols, & 
            &  X_matrix, n_predictors, ID_vector, Y_matrix, lambda_vec,  & 
            &                   SCORE_z1, SCORE_z2,SCORE)
            ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            ! //////////////////////////////////////////////
            
            ! ** METROPOLIS **
            df = SCORE - SCOREprev
            call random_number(rv)
            if (df < 0.0 .or. exp(-df / Temp_curr) > rv) then
                Sprev  = S             ! x_current = x_new
                SCOREprev = SCORE      ! y_current = y_new
                if (SCORE < SCOREbest) then
                    if(verbose .eq. 1) write(*, *) "best updated", cycle_i, SCORE
                    Sbest  = S
                    SCOREbest = SCORE
                    SCORE_z1best = SCORE_z1
                    SCORE_z2best = SCORE_z2
                end if
            end if
            
        end do

        ! -----------------
        ! ** PHASE 2 **
        ! -----------------
        ! Update cycle
        cycle_i = cycle_i + 1
        
        ! Cool down
        Temp_curr = Temp_curr * cooling_rate

        ! Get counter conditions
        ! so you have [SCOREbest2, SCOREbest1, and SCOREbest]
        ! SCOREbest was just updated, and so you should have the best past 3
        ! so now, check stay counter
        ! Need to check that this works correctly ....
        abs_score_diff = abs(SCOREbest1 - SCOREbest2)
        curr_rel_tol = rel_tol * max(abs(SCOREbest1), abs(SCOREbest2))
        if(abs_score_diff <= max(curr_rel_tol, abs_tol)) then
            stay_counter = stay_counter + 1
        else 
            stay_counter = 0
        end if

        ! check two stop conditions 
        ! Condition !: temperature   
         if (Temp_curr .lt. Temp_min) then
            if(verbose .eq. 1) write(*,*) "Cooled to final temperature"
            if(verbose .eq. 1) write(*,*) "Stay counter: ", stay_counter
            if(verbose .eq. 1) write(*,*) "abs score diff", abs_score_diff
            if(verbose .eq. 1) write(*,*) "curr_rel_tol", curr_rel_tol
            if(verbose .eq. 1) write(*,*) "abs_tol", abs_tol
            go to 99
        end if

        ! Condition 2: stability counter
        if (stay_counter .gt. max_stay_counter) then
            if(verbose .eq. 1) write(*,*) "Stayed unchanged in ", stay_counter, " iterations"
            go to 99
        end if

    end do

    ! Update the outputs
99  S = Sbest
    SCORE = SCOREbest
    SCORE_z1 = SCORE_z1best
    SCORE_z2 = SCORE_z2best
    if(verbose .eq. 1) write(*,*) "N. cycles = ", cycle_i

end subroutine simann


