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
!      occurs. Value of changer(n_events) is never used.
! nodes, const integer(n_events): node identifier number for the node
!      which is where this event happens. (Final leaves have node
!      numbers 1...n. These are not represented in nodes array.)
! times, const double(n_events): clock time prior to each event
!      i.e. times(i) is interval between event i-1 and event i.
!      value of times(1) is never used.
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
! time_E_to_S, const double(n_samples): time_E_to_S(nextSample) is clock time between most 
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

! Lots of this never changes and gets passed through many calls.
! Had I the patience, I'd learn about modules and store all this
! in what passes for a 'common' block in modern Fortran.

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

subroutine advance_trait(trait, timeBehind, nextSample, samples, event, &
     place, finalUpdate, n_events, times, n_samples, maxLeaves, &
     preSampleEvent, time_E_to_S, sigma, theta, debug)

  implicit none
  double precision, intent(INOUT) :: trait, timeBehind, samples(maxLeaves,n_samples)
  integer, intent(INOUT) :: nextSample
  integer, intent(IN) :: event, place
  logical, intent(IN) :: finalUpdate, debug

  ! And the invariant stuff passed everywhere
  
  integer, intent(IN) :: n_events, n_samples, maxLeaves, preSampleEvent(n_samples)
  double precision, intent(IN) :: times(n_events), time_E_to_S(n_samples), sigma, theta
  ! Local variables
  double precision :: time_in

  time_in=0
  if (nextSample <= n_samples .and. event == preSampleEvent(nextSample)+1) then
     ! Stepping forward in time to 'event' will take us past (at least) one
     ! sample time. So we need to advance 'trait' to this new time,
     ! record value in 'samples', and set timeLeft appropriately.
     ! Do loop is necessary in case one step include several sample times.
     do while (nextSample <= n_samples .and. event == preSampleEvent(nextSample)+1)
        if (debug) then
           print '(A,F3.1,A,I1,A,I1,A,I1)',&
                "Advancing time by ",timeBehind+time_E_to_S(nextSample)-time_in,&
                " for sample ",nextSample," after event=",event-1," place=",place               
                
        endif
        call update_trait(trait, timeBehind+time_E_to_S(nextSample)-time_in, sigma, theta)
        timeBehind=0
        ! time_in is how far into this time step we've advanced the value of
        ! 'trait'. 
        time_in = time_E_to_S(nextSample) 
        samples(place, nextSample) = trait
        nextSample = nextSample + 1
     enddo
     timeBehind = times(event)-time_in
  else
     ! No samples in this step
     timeBehind = timeBehind+times(event)
  end if
  if (finalUpdate) then
     if (debug) then
        print '(A,F3.1,A,I1,A,I1)', "Advancing time by ",timeBehind,&
             " for end of branch at event=",event," place=",place
     endif
     call update_trait(trait, timeBehind, sigma, theta)
     timeBehind = 0
  end if
end subroutine advance_trait

! We are at event 'event', which is a speciation event.
! We are going to process one of the two subtrees, which one
! is determined by value of 'place'. (For place =  place-of-'event',
! we do left subtree, if place = ditto + 1, right subtree.)
recursive subroutine do_subtree(place, n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, time_E_to_S, samples, trait, sigma, &
     theta, nextSample, debug)

  implicit none
  logical, intent(IN) :: debug
  integer, intent(IN) :: n, n_events, event, n_samples, maxLeaves
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: changer, nodes
  integer, dimension(n_samples), intent(IN) :: preSampleEvent

! place, nextSample change, but returned value not used by calling function.
  integer, intent(INOUT) :: nextSample, thisEdge, place 
  integer, dimension(n + n_events - 2, 2), intent(INOUT) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: time_E_to_S

! trait changes, but returned value not used by calling function.
  double precision, intent(INOUT) :: trait
  double precision, dimension(n + n_events - 2), intent(INOUT) :: edge_length, edge_trait
  double precision, dimension(maxLeaves, n_samples), intent(INOUT) :: samples

  integer :: newevent, tmpnextSample
  double precision :: timeBehind, tmptrait

  ! On entry:
  ! 'event' is a speciation event, and we are processing the 'place'th
  !   edge, which is one of the two edges from this event.
  !   At start of this edge, we have 'trait' value. 
  ! We are going to step forward through time, by advancing event pointer
  ! 'newevent'. 'timeBehind' records how much time has passed which has not
  ! yet been used in altering 'trait' value.
  
  if (debug) then
     print '(A,I1,A,I1,A,I1,A,I1)', "Enter do_subtree, event=",event," (node=",nodes(event),"), place=",place
  endif

  timeBehind=0
  newevent = event+1
  ! Step forward through events until either newevent is a
  ! speciation/extinction on this branch, or is end of tree
  do while (place .ne. changer(newevent) .and. newevent .ne. n_events)
     call advance_trait(trait, timeBehind, nextSample, samples, newevent, &
          place, .false., n_events, times, n_samples, maxLeaves, &
          preSampleEvent, time_E_to_S, sigma, theta, debug)
     if (changer(newevent) < place) place = place + changes(newevent)
     newevent = newevent+1
  end do
  ! Now newevent is on this branch, or is end of tree.
  ! Update trait with any remaining time:
  call advance_trait(trait, timeBehind, nextSample, samples, newevent, &
       place, .true., n_events, times, n_samples, maxLeaves, preSampleEvent, &
       time_E_to_S, sigma, theta, debug)

  ! Update edge info: edge length, trait value at end.
  edge_length(thisEdge) = sum(times((event + 1):newevent)) 
  edge_trait(thisEdge) = trait
    
  if (newevent == n_events) then
     ! 'reached end of tree' case. Terminating node is numbered 'place'.
     edge(thisEdge,:) = (/ nodes(event), place /)
     thisEdge = thisEdge+1
     if (debug) then
        print '(A,I1,A,I1)', "Reached end of tree at event=",newevent,&
              ", node and place=",place
     endif
  else
     ! two 'reached event' cases: extinction or speciation.
     ! terminating node number is nodes(newevent)
    edge(thisEdge,:) = (/ nodes(event), nodes(newevent) /)
    thisEdge = thisEdge + 1
    if (changes(newevent) == 1) then
       if (debug) then
          print '(A,I1,A,I1,A,I1)', "Speciation at event=",newevent,&
               " (node=",nodes(newevent),"), place=",place
       endif
       ! Found a speciation event at end of this edge.
       ! TO DO: are next two needed? Could be avoided by
       ! making local copy in tree_climb?
      tmpnextSample = nextSample ! temp save to not trash nextSample.
      tmptrait = trait ! Ditto trait
      ! Starting from speciation event 'newevent', process both branches:
      call tree_climb(n, n_events, changes, changer, nodes, times, &
           newevent, thisEdge, edge, edge_length, edge_trait, n_samples, &
           preSampleEvent, maxLeaves, time_E_to_S, samples, tmptrait, &
           sigma, theta, tmpnextSample, debug)
   else
      ! extinction event. Nothing need doing except debug output.
      if (debug) then
         print '(A,I1,A,I1,A,I1)', "Extinction at event=",newevent,&
              " (node=",nodes(newevent),"), place=",place
      endif
   endif
  endif
  if (debug) then
     print '(A,I1,A,I1,A,I1)', "Exiting do_subtree, event=",event,&
          " (node=",nodes(event),"), place=",place
  endif

end subroutine do_subtree

! n.edge = n + n.events - 2
! res = tree.climb(1, 1, matrix(0, n.edge, 2), numeric(n.edge))

! See variable glossary for explanation of parameters.
recursive subroutine tree_climb(n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, time_E_to_S, samples, baseTrait, sigma, &
     theta, baseNextSample, debug)

  implicit none
  logical, intent(IN) :: debug
  integer, intent(IN) :: n, n_events, event, n_samples, maxLeaves
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: changer, nodes
  integer, dimension(n_samples), intent(IN) :: preSampleEvent

  integer :: thisEdge, place, newevent, baseNextSample, tmpnextSample, nextSample
  integer, dimension(n + n_events - 2, 2), intent(INOUT) :: edge

  double precision, intent(IN) :: sigma, theta
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: time_E_to_S

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
     preSampleEvent, maxLeaves, time_E_to_S, samples, trait, sigma, &
     theta, nextSample, debug)

  nextSample = baseNextSample
  trait = baseTrait
  place = changer(event) + 1 ! Do the place+1 edge this time.

  call do_subtree(place, n, n_events, changes, changer, nodes, &
     times, event, thisEdge, edge, edge_length, edge_trait, n_samples, &
     preSampleEvent, maxLeaves, time_E_to_S, samples, trait, sigma, &
     theta, nextSample, debug)
end subroutine

