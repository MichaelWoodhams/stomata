! This is simple call around tree_climb so that I can gdb it.

program debug
  implicit none
  integer, parameter :: n=5 ! number of live tips of tree
  integer, parameter :: n_events=17 ! number of speciation/extinction events.
  integer, dimension(n_events) :: leaves = [1,2,3,2,3,4,5,6,5,4,5,4,5,6,5,4,5]
  integer, dimension(n_events-1) :: changes = [1,1,-1,1,1,1,1,-1,-1,1,-1,1,1,-1,-1,1]
  integer, dimension(n_events) :: changer = [1,2,1,1,3,4,5,2,1,3,4,1,2,3,1,1,5]
  integer, dimension(n_events) :: nodes = [12,13,11,14,15,16,17,10,9,18,8,19,20,7,6,21,0]
  double precision, dimension(n_events) :: times = &
       [0.236686751,0.288092867,0.349368826,0.355765762,0.182190466,&
        0.087765656,0.190416062,0.110437159,0.016545403,0.038370095,&
        0.089007730,0.005424049,0.090686383,0.002863372,0.060970358,&
        0.121784159,0.107278600]
  integer :: time=1
  integer :: a=1
  ! Initialize with nonsense values. I'd use NaN, but Fortran doesn't make this easy.
  integer, dimension(n + n_events - 2, 2) :: edge=-1
  double precision, dimension(n + n_events - 2) :: edge_length=-1 
  double precision, dimension(n + n_events - 2) :: edge_trait=-10000
  integer, parameter :: n_samples=16
  integer, dimension(n_samples) :: se=[3,3,3,3,4,4,4,5,6,6,6,6,7,10,10,12]
  integer :: n_leaves = -1000  ! test my theory this is never used.
  integer, parameter :: ml = 6
  double precision, dimension(n_samples) :: t_el = &
       [0.19681548,0.30955564,0.33210367,0.35465170,0.04398200, &
        0.08907807,0.13417413,0.01962776,0.06715029,0.08969832, &
        0.11224635,0.15734242,0.10221455,0.02705402,0.07215008,0.04536240]
  double precision, dimension(ml, n_samples) :: samples=-10000
  double precision :: trait = 0.0
  double precision :: sigma = 1.0
  double precision :: theta = 0.0
  integer :: ws = 1
 
  call tree_climb(n, n_events, leaves, changes, changer, nodes, times, &
       time, a, edge, edge_length, edge_trait, n_samples, se, n_leaves, &
       ml, t_el, samples, trait, sigma, theta, ws)

end program


