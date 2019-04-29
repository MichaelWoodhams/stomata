! This is simple call around tree_climb so that I can gdb it.

program debug
  implicit none
  integer :: n ! number of live tips of tree
  integer :: n_events ! number of speciation/extinction events.
  integer, dimension(:), allocatable :: leaves ! len n_events
  integer, dimension(:), allocatable :: changes ! len n_events-1
  integer, dimension(:), allocatable :: changer ! len n_events
  integer, dimension(:), allocatable :: nodes ! len n_events
  double precision, dimension(:), allocatable :: times ! len n_events
  integer :: time
  integer :: thisEdge
  ! Initialize with nonsense values. I'd use NaN, but Fortran doesn't make this easy.
  integer, dimension(:,:), allocatable :: edge ! dim (n+n_events-2,2)
  double precision, dimension(:), allocatable :: edge_length ! len n+n_events-2
  double precision, dimension(:), allocatable :: edge_trait ! len n+n_events-2
  integer :: n_samples ! =16 normally. Number of times at which to sample
  integer, dimension(:), allocatable :: preSampleEvent ! len n_samples
  integer :: maxLeaves ! max number of leaves at any time. Typically equal or slightly bigger than n
  double precision, dimension(:), allocatable :: time_E_to_S ! len n_samples
  double precision, dimension(:,:),allocatable :: samples ! dim(maxLeaves,n_samples)
  double precision :: trait 
  double precision :: sigma 
  double precision :: theta 
  integer :: nextSample 
  integer :: i, j
  integer, allocatable :: seed(:)
  character (len=20) :: string
  logical :: debugFlag

  ! Seed RNG. Will later reuse n for something else.
  call random_seed(size = n)
  allocate(seed(n))
  seed = [(i,i=1,n)]
  call random_seed(put=seed)

  string='R_output.txt'
  if (command_argument_count() > 0) then
     call get_command_argument(1, string)
  endif
  ! For debugging, provide the optional file name argument, then another arg.
  debugFlag = (command_argument_count() > 1)
  
  ! A lack of error checking here.
  open(10,file=string)
  read(10,*) n
  read(10,*) n_events
  read(10,*) n_samples  
  read(10,*) sigma
  read(10,*) theta

  read(10,*) string
  if (string /= "changes") then
    stop "changes not found"
  end if
  allocate(changes(n_events-1))
  do i=1,n_events-1
     read(10,*) changes(i)
  end do

  read(10,*) string
  if (string /= "changer") then
    stop "changer not found"
  end if
  allocate(changer(n_events))
  do i=1,n_events
     read(10,*) changer(i)
  end do

  read(10,*) string
  if (string /= "nodes") then
    stop "nodes not found"
  end if
  allocate(nodes(n_events))
  do i=1,n_events
     read(10,*) nodes(i)
  end do
  
  read(10,*) string
  if (string /= "preSampleEvent") then
    stop "preSampleEvent not found"
  end if
  allocate(preSampleEvent(n_samples))
  do i=1,n_samples
     read(10,*) preSampleEvent(i)
  end do

  read(10,*) string
  if (string /= "times") then
    stop "times not found"
  end if
  allocate(times(n_events))
  do i=1,n_events
     read(10,*) times(i)
  end do

  read(10,*) string
  if (string /= "time_E_to_S") then
    stop "time_E_to_S not found"
  end if
  allocate(time_E_to_S(n_samples))
  do i=1,n_samples
     read(10,*) time_E_to_S(i)
  end do

  close(10) 

  maxLeaves=0
  allocate(leaves(n_events))
  leaves(1)=1
  do i=2,n_events
     leaves(i)=leaves(i-1)+changes(i-1)
     maxLeaves=max(maxLeaves,leaves(i))
  end do

  allocate(edge(n+n_events-2,2))
  allocate(edge_length(n+n_events-2))
  allocate(edge_trait(n+n_events-2))
  allocate(samples(maxLeaves,n_samples))
  ! Initialize to nonsense values
  edge=-1
  edge_length=-1
  edge_trait=-100
  samples=-100

  time=1
  thisEdge=1
  trait=0.0
  nextSample=1

  call tree_climb(n, n_events, changes, changer, nodes, times, time, &
       thisEdge, edge, edge_length, edge_trait, n_samples, preSampleEvent, &
       maxLeaves, time_E_to_S, samples, trait, sigma, theta, nextSample, debugFlag)

  open(20,file="f95out.txt",action='write')
  write(20,*) "edge"
  do i=1,n+n_events-2
     write(20,"(2I6)") edge(i,1), edge(i,2)
  end do
  write(20,*) "edge_length"
  do i=1,n+n_events-2
     write(20,"(F10.7)") edge_length(i)
  end do
  write(20,*) "edge_trait"
  do i=1,n+n_events-2
     write(20,"(F10.5)") edge_trait(i)
  end do
  write(20,*) "samples"
  do i=1,maxLeaves
     do j=1,n_samples
        write(20,"(F10.4)",advance='no') samples(i,j)
     end do
     write(20,*)
  end do
  close(20)

end program


