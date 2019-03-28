! Compile file with
! R CMD SHLIB tree_climb.f90


! update_trait randomly updates a trait according to an Ornstein-Uhlenbeck
! process with mean 0.
!  trait: the value to be updated (value of variable changes)
!  deltat: time period between old and updated trait
!  sigma: random drift speed
!  theta: OU process parameter. High theta means values decay quickly towards 0
!  rands: just for pre-allocated memory. Mostly kept for sake of
!         minimizing changes to old code.

subroutine update_trait(trait, deltat, sigma, theta, rands)
  double precision, intent(INOUT) :: trait
  double precision, intent(IN) :: deltat, sigma, theta
  double precision, dimension(2) :: rands
  double precision, parameter :: pi = 3.1415926535879
  double precision :: rand_var
  
  ! Generate N(0,1) random_var by the Box-Muller tranform.
  ! We could get a second normal random number
  ! for little effort, but we don't bother to do so.
  call random_number(rands)
  rand_var=sqrt(-2.0*log(rands(1)))*cos(2.0*pi*rands(2))

  if (theta==0) then
    ! For theta=0, OU process becomes just Brownian motion
    trait = trait + sigma * sqrt(deltat) * rand_var
 else
    ! Ornstein-Uhlenbeck process value x at time t
    ! given value was x' at (earlier) time t' is normally distributed
    ! with mean x' exp(-theta(t-t')) and variance
    ! sigma^2/(2 theta)(1-exp(-2 theta(t-t')))
    ! Here x' - old 'trait' value, t-t' = deltat.
    trait = trait*exp(-theta*deltat) + &
            sigma*sqrt((1-exp(-2*theta*deltat))/(2*theta))*rand_var
 endif
end subroutine update_trait

! n.edge = n + n.events - 2
! res = tree.climb(1, 1, matrix(0, n.edge, 2), numeric(n.edge))

! arguments: n, n_events, leaves, changes, changer, nodes, times, time, a,&
! edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, ws

recursive subroutine tree_climb(n, n_events, leaves, changes, changer, nodes, times, time, a,&
edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, theta, ws)

  implicit none
  integer, intent(IN) :: n, n_events, time, n_samples, n_leaves, ml
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: leaves, changer, nodes
  integer, dimension(n_samples), intent(IN) :: se

  integer :: a, place, newtime, tmpws, ws, ows
  integer, dimension(n + n_events - 2, 2) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: t_el

  double precision :: time_in, tmptrait, trait, otrait, sumtime, rand_var
  double precision, dimension(2) :: rands
  double precision, dimension(n + n_events - 2) :: edge_length, edge_trait
  double precision, dimension(ml, n_samples) :: samples

  !open(12, file="test.txt", status="old", position="append", action="write")
! write(12,*) time, ws, trait
!  write(12,*) samples(1,1), samples(442,15), samples(442,1), "start"
  ows = ws
  otrait = trait
  time_in = 0.0
  place = changer(time)

  if (ws <= n_samples) then
    do while (time == se(ws))
      call update_trait(trait, t_el(ws)-time_in, sigma, theta, rands)
      !call random_number(rands)
      !trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
      !write(12,*) place, ws, trait, "straight_in"
      samples(place, ws) = trait
                !write(12,*) samples(place,ws)
      time_in = t_el(ws)
      ws = ws + 1
      if (ws > n_samples) exit      
    enddo
  endif
  
  newtime = time + 1
  sumtime = times(newtime) - time_in
  time_in = 0.0
  
  if (newtime < n_events) then
    do while (place .ne. changer(newtime))
      if (changer(newtime) < place) place = place + changes(newtime)
      
      time_in = 0.0
  
      if (ws <= n_samples) then
        do while (newtime == se(ws))
          call update_trait(trait, t_el(ws)-time_in, sigma, theta, rands)
          !call random_number(rands)
          !trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
          !write(12,*) place, ws, trait, "along"
          samples(place, ws) = trait
                    !write(12,*) samples(place,ws)
          time_in = t_el(ws)
          ws = ws + 1
          if (ws > n_samples) exit
        enddo
      endif
      
      newtime = newtime + 1
      sumtime = sumtime + times(newtime) - time_in
      if (newtime == n_events) exit
    enddo
  endif
  
  edge_length(a) = sum(times((time + 1):newtime)) ! sumtime (doesn't include time to samples)
  call update_trait(trait, sumtime, sigma, theta, rands)
  !call random_number(rands)
  !trait = trait + sigma * sqrt(-2.0*sumtime*log(rands(1)))*cos(2.0*pi*rands(2))
  edge_trait(a) = trait
    
  if (newtime < n_events) then
    edge(a,:) = (/ nodes(time), nodes(newtime) /)

    if (changes(newtime) == 1) then
      a = a + 1
      tmpws = ws
      tmptrait = trait
      ! arguments: n, n_events, leaves, changes, changer, nodes, times, time, a,&
      ! edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, ws
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newtime, a, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, theta, tmpws)
    endif
  else
    edge(a,:) = (/ nodes(time), place /)
  endif
  
  a = a + 1

  ws = ows
  trait = otrait
  time_in = 0.0
  place = changer(time) + 1

  if (ws <= n_samples) then
    do while (time == se(ws))
      call update_trait(trait, t_el(ws)-time_in, sigma, theta, rands)
      !call random_number(rands)
      !trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
      !write(12,*) place, ws, trait, "straight_in"
      samples(place, ws) = trait
                !write(12,*) samples(place,ws)
      time_in = t_el(ws)
      ws = ws + 1
      if (ws > n_samples) exit
    enddo
  endif
  
  newtime = time + 1
  sumtime = times(newtime)
  if (newtime < n_events) then
    do while (place .ne. changer(newtime))
      if (changer(newtime) < place) place = place + changes(newtime)
      time_in = 0.0

      if (ws <= n_samples) then  
        do while (newtime == se(ws))
          call update_trait(trait, t_el(ws)-time_in, sigma, theta, rands)
          !call random_number(rands)
          !trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
          !write(12,*) place, ws, trait, "along"
          samples(place, ws) = trait
          !write(12,*) samples(place,ws)
          time_in = t_el(ws)
          ws = ws + 1
          if (ws > n_samples) exit          
        enddo
      endif

      newtime = newtime + 1
      sumtime = sumtime + times(newtime) - time_in
      if (newtime == n_events) exit
    enddo
  endif

  edge_length(a) = sum(times((time + 1):newtime)) ! sumtime (doesn't include time to samples)
  call update_trait(trait, sumtime, sigma, theta, rands)
  !call random_number(rands)
  !trait = trait + sigma * sqrt(-2.0*sumtime*log(rands(1)))*cos(2.0*pi*rands(2))
  edge_trait(a) = trait

  if (newtime < n_events) then
    edge(a,:) = (/ nodes(time), nodes(newtime) /)
    if (changes(newtime) == 1) then
      a = a + 1
      tmpws = ws
      tmptrait = trait
      ! arguments: n, n_events, leaves, changes, changer, nodes, times, time, a,&
      ! edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, ws
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newtime, a, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, theta, tmpws)
    endif
  else
    edge(a,:) = (/ nodes(time), place /) 
  endif

!  write(12,*) samples(1,1), samples(442,15), samples(442,1), "end"
!  close(12)

end subroutine

