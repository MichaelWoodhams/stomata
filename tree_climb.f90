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
! thisEdge, integer: indexes edge arrays. Starts at 1 and increments everytime
!     a new edge is encountered
!
! "Sample" construct. We have several sample times, at which we want
! to know all trait values of leaves which exist at that time.
!
! n_samples, const integer: number of times at which sampling occurs
! maxLeaves, const integer: maximum number of leaves at any sampling time.
! samples, double(maxLeaves, n_samples): samples(place,sample_number) = 
!      trait value at given sample time, for edge 'place' from left at that 
!      slice through tree.
! preSampleEvent, const integer(n_samples): event number immediately preceeding
!      sample time.
! t_el, const double(n_samples): t_el(nextSample) is clock time between most 
!      recent event (number preSampleEvent(nextSample)) and sample time.
! nextSample, integer: the next sampling time we will meet, i.e. earliest 
!    sampling time we have not already passed.
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
! n is used only in dimensions of edge arrays. Could replace by n_edges.


! update_trait randomly updates a trait according to an Ornstein-Uhlenbeck
! process with mean 0.
!  trait: the value to be updated (value of variable changes)
!  deltat: time period between old and updated trait
!  sigma: random drift speed
!  theta: OU process parameter. High theta means values decay quickly towards 0

subroutine update_trait(trait, deltat, sigma, theta)
  double precision, intent(INOUT) :: trait
  double precision, intent(IN) :: deltat, sigma, theta
  double precision, parameter :: pi = 3.1415926535879
  double precision :: rand_var, rand1, rand2
  
  ! Generate N(0,1) random_var by the Box-Muller tranform.
  ! We could get a second normal random number
  ! for little effort, but we don't bother to do so.
  call random_number(rand1)
  call random_number(rand2)
  rand_var=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)

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

subroutine advance_to_sample_time(nextSample, n_samples, maxLeaves, event, &
  preSampleEvent, trait, t_el, time_in, sigma, theta, samples, place)

  implicit none
  integer, intent(INOUT) :: nextSample
  integer, intent(IN) :: n_samples, maxLeaves, event, place, preSampleEvent(n_samples)
  double precision, intent(INOUT) :: trait, time_in, samples(maxLeaves,n_samples)
  double precision, intent(IN) :: sigma, theta, t_el(n_samples)
  
  do while (nextSample <= n_samples .and. event == preSampleEvent(nextSample))
     call update_trait(trait, t_el(nextSample)-time_in, sigma, theta)
     samples(place, nextSample) = trait
     time_in = t_el(nextSample)
     nextSample = nextSample + 1
  enddo
end subroutine advance_to_sample_time

! We are at event 'event', which is a speciation event.
! We are going to process one of the two subtrees, which one
! is determined by value of 'place'. (For place =  place-of-'event',
! we do left subtree, if place = ditto + 1, right subtree.)
recursive subroutine do_subtree(place, n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, t_el, samples, trait, sigma, &
     theta, nextSample)

    implicit none
  integer, intent(IN) :: n, n_events, event, n_samples, maxLeaves
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: changer, nodes
  integer, dimension(n_samples), intent(IN) :: preSampleEvent

! place, nextSample change, but returned value not used by calling function.
  integer, intent(INOUT) :: nextSample, thisEdge, place 
  integer, dimension(n + n_events - 2, 2), intent(INOUT) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: t_el

! trait changes, but returned value not used by calling function.
  double precision, intent(INOUT) :: trait
  double precision, dimension(n + n_events - 2), intent(INOUT) :: edge_length, edge_trait
  double precision, dimension(maxLeaves, n_samples), intent(INOUT) :: samples

  integer :: newevent, tmpnextSample
  double precision :: time_in, tmptrait, sumtime

  ! On entry:
  ! 'event' is a speciation event, and we are processing the 'place'th
  !   edge, which is one of the two edges from this event.
  !   At start of this edge, we have 'trait' value. 
  
  ! How far have we progressed along branch (in 'trait' value)?
  time_in = 0.0
  ! Closely associated with this is 'sumtime' which is how much
  ! trait evolution time we've passed, but not updated 'trait' for.
  
  ! Usually changes nothing (if no sample times between 'event' and
  ! 'event+1'), but if there are, updates 'trait' to time of last of
  ! these sample times, updates preSampleEvent to point to latest of
  ! these sample times, and time_in to the time diff between event
  ! and last sample time (where trait was updated to.)
  call advance_to_sample_time(nextSample, n_samples, maxLeaves, event, &
       preSampleEvent, trait, t_el, time_in, sigma, theta, samples, place)

  newevent = event + 1
  ! In stepping to next event, sumtime is how much time has passed which
  ! has not been used to update 'trait'
  ! This version of program is being bug-compatible with old version.
  ! This is an instance of that. Old program was inconsistent.
  if (place == changer(event)) then
     sumtime = times(newevent) - time_in
  else
     sumtime = times(newevent)
  endif
  time_in = 0.0 ! no more 'hidden' accounted-for time, as sumtime accounts for it.

  ! do while (new event is not occuring on this edge)
  !      and (haven't reached tip of tree)
  do while (place .ne. changer(newevent) .and. newevent < n_events )
     ! As 'place' counts how far we are from the left of the tree,
     ! if an event happened to the left of us, 'place' changes
     if (changer(newevent) < place) place = place + changes(newevent)
     
     time_in = 0.0

     ! I believe there is a bug here: time being advanced should make use
     ! of sumtime. I'll come back to fix this later.
     call advance_to_sample_time(nextSample, n_samples, maxLeaves, newevent, &
          preSampleEvent, trait, t_el, time_in, sigma, theta, samples, place)
     
     newevent = newevent + 1
     sumtime = sumtime + times(newevent) - time_in
  enddo
  ! At this point, one of three things has happened:
  ! * We've reached the end of the tree (newevent == n_events)
  ! * We've reached a speciation event on this branch
  ! * We've reached an extinction event on this branch
  
  ! Update edge info: edge length, trait value at end.
  edge_length(thisEdge) = sum(times((event + 1):newevent)) 
  call update_trait(trait, sumtime, sigma, theta)
  edge_trait(thisEdge) = trait
    
  if (newevent == n_events) then
     ! 'reached end of tree' case. Terminating node is numbered 'place'.
    edge(thisEdge,:) = (/ nodes(event), place /)
  else
    ! two 'reached event' cases: terminating node number is nodes(newevent)
    edge(thisEdge,:) = (/ nodes(event), nodes(newevent) /)

    ! If 'reached extinction event' case, nothing more to be done. Otherwise:
    if (changes(newevent) == 1) then
       ! Found a speciation event at end of this edge.
       thisEdge = thisEdge + 1 ! Ready to fill in next edge
       ! TO DO: are next two needed? Could be avoided by
       ! making local copy in tree_climb?
      tmpnextSample = nextSample ! temp save to not trash nextSample.
      tmptrait = trait ! Ditto trait
      ! Starting from speciation event 'newevent', process both branches:
      call tree_climb(n, n_events, changes, changer, nodes, times, &
           newevent, thisEdge, edge, edge_length, edge_trait, n_samples, &
           preSampleEvent, maxLeaves, t_el, samples, tmptrait, &
           sigma, theta, tmpnextSample)
    endif
  endif

end subroutine do_subtree

! n.edge = n + n.events - 2
! res = tree.climb(1, 1, matrix(0, n.edge, 2), numeric(n.edge))

! See variable glossary for explanation of parameters.
recursive subroutine tree_climb(n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, t_el, samples, baseTrait, sigma, &
     theta, baseNextSample)

  implicit none
  integer, intent(IN) :: n, n_events, event, n_samples, maxLeaves
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: changer, nodes
  integer, dimension(n_samples), intent(IN) :: preSampleEvent

  integer :: thisEdge, place, newevent, baseNextSample, tmpnextSample, nextSample
  integer, dimension(n + n_events - 2, 2), intent(INOUT) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: t_el

  double precision :: time_in, tmptrait, trait, baseTrait, sumtime
  double precision, dimension(n + n_events - 2) :: edge_length, edge_trait
  double precision, dimension(maxLeaves, n_samples) :: samples

  ! On entry to this routine:
  ! 'event' is a speciation event.
  ! 'place' is the place of the edge leading up to this event.
  ! 'baseTrait' is the trait value at this event
  ! 'baseNextSample' is the index of next sampling time after this event
  
  nextSample = baseNextSample
  trait = baseTrait
  place = changer(event)
  ! As this is speciation node at 'place', it has two child edges:
  ! at 'place' and 'place+1'. We call do_subtree on each of them.
  
  call do_subtree(place, n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, t_el, samples, trait, sigma, &
     theta, nextSample)

  thisEdge = thisEdge + 1 ! moving to a new edge
  nextSample = baseNextSample
  trait = baseTrait
  place = changer(event) + 1 ! Do the place+1 edge this time.

  call do_subtree(place, n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, t_el, samples, trait, sigma, &
     theta, nextSample)
end subroutine

