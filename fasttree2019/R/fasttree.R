fasttree = function(n, lambda, mu, frac = 1, traits = FALSE, sigma = 1, alpha=0,
                     treeDepth = 1, sampleTimes = 1, seed=1L, traitseed=NULL, startTrait=0)
{
  #  print('start')
  if (any(sampleTimes <= 0 | sampleTimes > treeDepth)) stop('All values of sampleTimes must be between 0 and treeDepth')

  # start off doing things in reverse
  # (start with maxleaf extant species, create species at rate mu, destroy at rate lambda)
  # events = -1 with prob l/(l+mu), 1 with prob mu/(l+mu)

  # calculate mean of distribution to predict number of likely events before
  # all maxleaf nodes are removed, then run enough events that we hit zero species

  maxleaf = round(n / frac)
  m = (lambda - mu)/(mu + lambda) # average rate of species decrease (in reverse)

  sample.events = integer(length(sampleTimes))
  set.seed(seed)
  if (is.null(traitseed)) { traitseed <- seed }

  while (any(sample.events <= 2))
  {
    changes = NULL
    # track events
    # run slightly more than the expected number of events, should usually be enough
    # but will repeat if we haven't hit zero yet
    # (number chosen to roughly optimise run time for 1700 leaves, selfishly,
    # but should be pretty fast for all cases)


    while (sum(changes) > -maxleaf)
      changes = rbind(changes, sample(c(-1, 1), size = ceiling(1.2*maxleaf/m), prob = c(lambda, mu), replace = TRUE))
    # track number of leaves initially (maxleaf), then after each event
    leaves = maxleaf + c(0, cumsum(changes))
    # find out where we hit zero leaves and crop just before that
    n.events = min(which(leaves == 0)) - 1

    # now turn things the right way round, excluding the first change in number of leaves
    changes = -changes[(n.events-1):1]
    leaves = leaves[n.events:1]
    # generate times till next event based on current number of leaves
    times = rexp(n = n.events, rate = leaves * (mu + lambda))
    # for each event, randomly allocate a leaf for the event to happen to
    changer = ceiling(runif(n.events) * leaves)
    # allocate edge IDs using protocol for phylo class: first IDs go to tips
    # (i.e. 1:(maxleaf+L) below) then remaining IDs to internal nodes
    # note that the first maxleaf IDs go to extant tips
    nodes = 0*leaves
    ii = which(changes == -1) # events at which an extinction happens
    L = length(ii) # number of extinct tips
    nodes[ii] = maxleaf + (L:1) # number these (most recent has smallest number)
    ii = which(changes == 1) # events at which speciation happens
    nodes[ii] = maxleaf + L + 1:length(ii) # number these
    n.edge = maxleaf + n.events - 2 # overall number of edges

    # Scaling of times. sum(times[-1]) is length of tree before scaling.
    times = times*treeDepth/sum(times[-1]) # now first split is time 0, leaves at time treeDepth
    curr.time = cumsum(times[-1]) # time from first split to each event
    #tree.time = sum(times[-1]) # time from first split to current day
    #sampleTimes2 = sampleTimes * tree.time # use sampleTimes as proportion of way through the tree
    n.samples = length(sampleTimes)


    ## New version
    # sample.events = sapply(sampleTimes, function(x) max(which(curr.time <= x))+1L) # work out which event each sample is after
    # n.leaves = leaves[sample.events] # number of leaves at each sampled point
    # max.leaves = max(n.leaves) # max number of leaves at any sampled point
    # t.el = sampleTimes # for sampleTimes before the first event, just set as the time elapsed
    # ii = which(sample.events > -Inf)
    # t.el[ii] = t.el[ii] - curr.time[sample.events[ii] - 1] # otherwise, how far after the last event each sample point is

    ## Old version
    sample.events = sapply(sampleTimes, function(x) which.min(curr.time < x)) # work out which event each sample is after
  }
  n.leaves = leaves[sample.events + 1] # number of leaves at each sampled point
  max.leaves = max(n.leaves) # max number of leaves at any sampled point
  t.el = sampleTimes - curr.time[sample.events - 1] # how far after the last event each sample point is

  # Rescale sigma for tree length.
  #scaled.sigma <- sigma*tree.time/

  # use precompiled Fortran code to convert tree data to
  # phylo form, using a recursive function to calculate edges between nodes
  # also carries along trait information along edges, as well as at sampled points

  #  print('pre fortran')

  #  recover()

  # parameters come from main function and remain unchanged unless comments say otherwise
  # Modified for rewritten Fortran.
  test = .Fortran('tree_climb',
                  n = as.integer(n),
                  n_events = as.integer(n.events),
                  changes = as.integer(changes),
                  changer = as.integer(changer),
                  nodes = as.integer(nodes),
                  times = as.double(times),
                  event = 1L,
                  thisEdge = 1L,
                  edge = matrix(0L, n + n.events - 2, 2),
                  edge_length = numeric(n + n.events - 2),
                  edge_trait = numeric(n + n.events - 2),
                  n_samples = as.integer(n.samples),
                  se = as.integer(sample.events),
                  ml = as.integer(max.leaves),
                  t_el = as.double(t.el),
                  samples = matrix(as.double(0.0), max.leaves, n.samples),
                  trait = startTrait,
                  sigma = as.double(sigma),
                  alpha = as.double(alpha),
                  ws = 1L,
                  rngseed = traitseed,
                  debug = FALSE
  )

  if (length(test$edge_length) == 0) recover()

  # print('post test')

  res = list(edge = test$edge,
             edge.length = test$edge_length,
             edge.trait = test$edge_trait,
             tip.label = paste('t', 1:((maxleaf + n.events)/2), sep=''),
             root.edge = times[1],
             Nnode = ((maxleaf + n.events)/2) - 1
  )

  #  print('post res')0

  class(res) = 'phylo'

  # return tree and extra outputs if traits = TRUE
  if (!traits)
    res
  else
    list(tree = res,  # output tree
         samples = test$samples, # max.leaves x n.samples matrix of traits taken at sampleTimes (padded with zeros)
         n.leaves = n.leaves, # number of leaves at each of the sample points
         sampleTimes = sampleTimes) # times of each sample
}
