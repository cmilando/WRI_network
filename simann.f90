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
    
    ! -- UPDATE THE NEIGHBORS -- 
    x1(p0_int) = x0(p1_int)
    x1(p1_int) = x0(p0_int)

end

!--------------------------------------------------------------------------
! Cost function for SCORE
! 
! this is what needs to change for a new application
!--------------------------------------------------------------------------
subroutine get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)

    implicit none
    
    ! Inputs / Outputs
    integer, intent(in)        :: np                          ! total number of site-point pairs
    integer, intent(in)        :: magic_n                     ! total number of site-point pairs
    integer, intent(in)        :: nsites                      ! total number of sites
    integer, intent(in)        :: S(nsites)                   ! S, a specific iteration of the magic indicator
    real(kind=8), intent(in)   :: site_pairs_sub(np, 6)       ! X, outcome data (columns: site_id, pop_id, dist, population, metric1, score)
    real(kind=8), intent(in)   :: row_lookup(nsites, 3)       ! row lookup (columns: site_id, start, end)
    real(kind=8), intent(in)   :: penalty
    
    real(kind=8)  :: SCORE                       ! Score to be calculated
    
    ! Local variables
    integer                   :: i, j, start_idx, end_idx, s_count
    integer, allocatable      :: pop_id_list(:), S_ones(:)
    real(kind=8), allocatable :: score_list(:)
    integer                   :: freq, idx
    real(kind=8)              :: score2

    SCORE = 0.0

    ! Allocate arrays to store intermediate data
    allocate(pop_id_list(np))
    allocate(score_list(np))
    allocate(S_ones(magic_n))

    s_count = 0

    ! Loop through S to find indices where S(i) == 1
    do i = 1, nsites
        if (S(i) == 1) then
            s_count = s_count + 1
            S_ones(s_count) = i
        end if
    end do

    ! Loop through selected sites (where S = 1) and gather corresponding site pairs
    idx = 0

    do i = 1, magic_n
        ! ok this is where its starting to get weird
        start_idx = int(row_lookup(S_ones(i), 2))
        end_idx   = int(row_lookup(S_ones(i), 3))
        !write(*, *) start_idx, end_idx

        do j = start_idx, end_idx
            idx = idx + 1
            pop_id_list(idx) = int(site_pairs_sub(j, 2))
            score_list(idx) = site_pairs_sub(j, 6)
        end do
    end do

    ! Calculate frequency of each population id
    do i = 1, idx
        freq = count(pop_id_list(1:idx) == pop_id_list(i))

        ! Adjust score based on frequency and penalty
        score2 = score_list(i) * 1.0 / (penalty * freq)
        
        ! create a running sum
        SCORE = SCORE + score2
    end do

    SCORE = SCORE * (-1.0)

    !write(*, *) 'SCORE: ', SCORE

    ! Deallocate arrays
    deallocate(pop_id_list)
    deallocate(score_list)
    deallocate(S_ones)

end


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
subroutine simann(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE, cooling_rate, verbose)
    
    !  Temp_max, Temp_min, cooling_rate

    implicit none
    
    ! Inputs / Outputs that need to change for new applications
    integer, intent(in)        :: np                     ! total number of site-point pairs
    integer, intent(in)        :: magic_n                ! total number of site-point pairs
    integer, intent(in)        :: nsites                 ! total number of sites

    real(kind=8), intent(in)   :: site_pairs_sub(np, 6)  ! X, outcome data (columns: site_id, pop_id, dist, population, metric1, score)
    real(kind=8), intent(in)   :: row_lookup(nsites, 3)  ! row lookup (columns: site_id, start, end)
    real(kind=8), intent(in)   :: penalty
    
    ! Inputs / Outputs that don't need to change for new applications
    integer              :: S(nsites)              ! S, a specific iteration of the magic indicator
    integer              :: Sprev(nsites)                !  , a prevoious iteration of the magic indicator
    integer              :: Sbest(nsites)                !  , the best

    real(kind=8)         :: SCORE                        !  
    real(kind=8)         :: SCOREprev                    !  
    real(kind=8)         :: SCOREbest                    !  
    real(kind=8)         :: SCOREbest1, SCOREbest2       ! generations of SCOREbest

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
    Temp_max         = 500000.
    Temp_min         = 0.0000000001
    Temp_curr        = Temp_max

    ! **************
    ! INTIAL VALUES
    
    ! //////////////////////////////////////////////
    ! >>>>>>> this call will need to change <<<<<<<<
    call get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)
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
            ! >>>>>>> this call will need to change <<<<<<<<
            call get_score(S, magic_n, np, nsites, site_pairs_sub, row_lookup, penalty, SCORE)
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
    if(verbose .eq. 1) write(*,*) "N. cycles = ", cycle_i

END