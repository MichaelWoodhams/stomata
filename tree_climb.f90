
! rand_normal: return a normally distributed random number with mean 0, variance 1.
! This uses the Box-Muller tranform. We could get a second normal random number
! for little effort, but we don't bother to do so.
! 'rands' is just pre-allocated memory. It changes, but value isn't expected to be
! used outside this function.
! Result is returned in rand_var. This should be a function, but I'm modifying someone
! else's code and haven't used Fortran since Fortran 77, and haven't figured out
! how to define the return type of the function when I'm not inside a 'program' block.
! Yay for cargo cult programming!

subroutine rand_normal(rands, rand_var)
  double precision, dimension(2) :: rands
  double precision, intent(OUT) :: rand_var
  double precision, parameter :: pi = 3.1415926535879

  call random_number(rands)
  rand_var=sqrt(-2.0*log(rands(1)))*cos(2.0*pi*rands(2))
end subroutine rand_normal

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
      call rand_normal(rands,rand_var)
      trait = trait + sigma * sqrt(t_el(ws) - time_in) * rand_var
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
          call rand_normal(rands,rand_var)
          trait = trait + sigma * sqrt(t_el(ws) - time_in) * rand_var
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
  call rand_normal(rands,rand_var)
  trait = trait + sigma * sqrt(sumtime) * rand_var
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
      call rand_normal(rands,rand_var)
      trait = trait + sigma * sqrt(t_el(ws) - time_in) * rand_var
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
          call rand_normal(rands,rand_var)
          trait = trait + sigma * sqrt(t_el(ws) - time_in) * rand_var
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
  call rand_normal(rands,rand_var)
  trait = trait + sigma * sqrt(sumtime) * rand_var
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

