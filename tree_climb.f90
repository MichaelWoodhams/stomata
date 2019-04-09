! Compile file with
! R CMD SHLIB tree_climb.f90

! Glossary of variables:

! "Events" construct: an ordered sequence of event, each of which is
! a speciation or extinction at a particular node on the tree.
!   event 0 is the root of the tree, always a speciation.
!   even n_events is the end of time, where extant tree tips live.
!
! n_events, const integer: number of events
! changes, const integer(n_events-1): changes(event-1) = +1 if speciation,
!      -1 if extinction event. (No record for final event, as that is
!      just the end of time.)
! changer, const integer(n_events): The 'place' at which extinction/speciation
!      occurs.
! leaves, const integer(n_events): The number of leaves present at
!      immediately before each event. Is cumulative sum of 'changes'
! nodes, const integer(n_events): node identifier number for the node
!      which is where this event happens. (Final leaves have node
!      numbers 1...n. These are not represented in nodes array.)
! times, const double(n_events): clock time prior to each event
!      i.e. times(i) is interval between event i-1 and event i.
! event, integer: index into events. The root event for which each call
!      to tree_climb is processing the descendents of.
! newevent, integer: index into events, advanced as we traverse each
!      edge.
!
! "Edge" construct: a set of tree edges with properties.
! n, const integer: number of tips on tree.
!      n is used to derive the number of edges, which is n+n_events-2
! edge, integer(n+n_events-2,2): parent/child node numbers for each edge.
!     We fill this in as we go.
! edge_length, double(n+n_events-2): length of edge
! edge_trait, double(n+n_events-2): trait at end of edge
! thisEdge, integer: indexes edge arrays. Starts at 1 and increments everytime a
!     new edge is encountered
!
! "Sample" construct. We have several sample times, at which we want
! to know all trait values of leaves which exist at that time.
!
! n_samples, const integer: number of times at which sampling occurs
! ml, const integer: maximum number of leaves at any sampling time.
! samples, double(ml, n_samples): samples(place,sample_number) = trait value
!      at given sample time, for edge 'place' from left at that slice through
!      tree.
! se, const integer(n_samples): event number immediately preceeding sample time
! t_el, const double(n_samples): t_el(nextSample) is clock time between most recent
!      event (number se(nextSample)) and sample time.
! nextSample, integer: the next sampling time we will meet, i.e. earliest sampling
!    time we have not already passed.
! onextSample, integer: stores old nextSample
! tempnextSample, integer: used to prevent nextSample being overwritten.
!
! Others:
! 
! place, integer: indexes how far across the tree we are. E.g. at event 'event',
!     'place'=1 means we are on the leftmost branch extant at 'event'.
! trait, double: the current value of the randomly drifting trait
! sigma, double: variance/time of trait
! theta, double: OU parameter, strength of attraction of 0
! time_in, double: How far along current edge we are in 'clock' time?
! sumtime, double: ? cumulative time along this branch?



! n_leaves, const integer: ???


! Redundant variables: I think these can all be removed.
! leaves, const integer(n_events): number of extant leaves immediately prior(??)
!    to indexed event. (Is cumulative sum of 'changes')
! rands, double(2): just storage for two random numbers.
! n is used only in dimensions of edge arrays. Could replace by n_edges.


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

! If there are sample times between current event and next event,
! advance clock to last of those sample times, and record
! the trait values on the current edge.
! Else do nothing.
! nextSample, time_in can change (increase).

subroutine advance_to_sample_time(nextSample, n_samples, ml, event, se, trait, t_el, &
  time_in, sigma, theta, rands, samples, place)

  implicit none
  integer, intent(INOUT) :: nextSample
  integer, intent(IN) :: n_samples, ml, event, place, se(n_samples)
  double precision, intent(INOUT) :: trait, time_in, rands(2), samples(ml,n_samples)
  double precision, intent(IN) :: sigma, theta, t_el(n_samples)
  
  do while (nextSample <= n_samples .and. event == se(nextSample))
     call update_trait(trait, t_el(nextSample)-time_in, sigma, theta, rands)
     samples(place, nextSample) = trait
     time_in = t_el(nextSample)
     nextSample = nextSample + 1
  enddo
end subroutine advance_to_sample_time


! n.edge = n + n.events - 2
! res = tree.climb(1, 1, matrix(0, n.edge, 2), numeric(n.edge))

! arguments: n, n_events, leaves, changes, changer, nodes, times, event, thisEdge,&
! edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, baseNextSample

recursive subroutine tree_climb(n, n_events, leaves, changes, changer, nodes, times, event, thisEdge,&
edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, baseTrait, sigma, theta, baseNextSample)

  implicit none
  integer, intent(IN) :: n, n_events, event, n_samples, n_leaves, ml
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: leaves, changer, nodes
  integer, dimension(n_samples), intent(IN) :: se

  integer :: thisEdge, place, newevent, baseNextSample, tmpnextSample, nextSample
  integer, dimension(n + n_events - 2, 2) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: t_el

  double precision :: time_in, tmptrait, trait, baseTrait, sumtime, rand_var
  double precision, dimension(2) :: rands
  double precision, dimension(n + n_events - 2) :: edge_length, edge_trait
  double precision, dimension(ml, n_samples) :: samples

  nextSample = baseNextSample
  trait = baseTrait
  time_in = 0.0
  place = changer(event)

  call advance_to_sample_time(nextSample, n_samples, ml, event, se, trait, t_el, &
                         time_in, sigma, theta, rands, samples, place)

  newevent = event + 1
  sumtime = times(newevent) - time_in
  time_in = 0.0
  
  do while (place .ne. changer(newevent) .and. newevent < n_events )
     if (changer(newevent) < place) place = place + changes(newevent)
     
     time_in = 0.0
     
     call advance_to_sample_time(nextSample, n_samples, ml, newevent, se, trait, t_el, &
          time_in, sigma, theta, rands, samples, place)
     
     newevent = newevent + 1
     sumtime = sumtime + times(newevent) - time_in
  enddo
  
  edge_length(thisEdge) = sum(times((event + 1):newevent)) ! sumtime (doesn't include time to samples)
  call update_trait(trait, sumtime, sigma, theta, rands)
  edge_trait(thisEdge) = trait
    
  if (newevent == n_events) then
    edge(thisEdge,:) = (/ nodes(event), place /)
  else 
    edge(thisEdge,:) = (/ nodes(event), nodes(newevent) /)

    if (changes(newevent) == 1) then
      thisEdge = thisEdge + 1
      tmpnextSample = nextSample
      tmptrait = trait
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newevent, thisEdge, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, theta, tmpnextSample)
    endif
  endif
  
  thisEdge = thisEdge + 1

  nextSample = baseNextSample
  trait = baseTrait
  time_in = 0.0
  place = changer(event) + 1

  call advance_to_sample_time(nextSample, n_samples, ml, event, se, trait, t_el, &
                         time_in, sigma, theta, rands, samples, place)
  
  newevent = event + 1
  sumtime = times(newevent)
  time_in = 0.0 
  do while (place .ne. changer(newevent) .and. newevent < n_events)
     if (changer(newevent) < place) place = place + changes(newevent)
     time_in = 0.0
     
     call advance_to_sample_time(nextSample, n_samples, ml, newevent, se, trait, t_el, &
          time_in, sigma, theta, rands, samples, place)
     
     newevent = newevent + 1
     sumtime = sumtime + times(newevent) - time_in
  enddo

  edge_length(thisEdge) = sum(times((event + 1):newevent)) ! sumtime (doesn't include time to samples)
  call update_trait(trait, sumtime, sigma, theta, rands)
  edge_trait(thisEdge) = trait

  if (newevent == n_events) then
    edge(thisEdge,:) = (/ nodes(event), place /) 
  else
    edge(thisEdge,:) = (/ nodes(event), nodes(newevent) /)
    if (changes(newevent) == 1) then
      thisEdge = thisEdge + 1
      tmpnextSample = nextSample
      tmptrait = trait
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newevent, thisEdge, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, theta, tmpnextSample)
    endif
  endif

end subroutine

