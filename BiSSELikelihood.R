# BiSSE Likelihood estimation
# Sebastian Dunn
# BIOINF702 2019 Sem1

# Refs:
  # - Wayne P. Maddison, Peter E. Midford, Sarah P. Otto, Estimating a Binary 
  #     Character's Effect on Speciation and Extinction, Systematic Biology, 
  #     Volume 56, Issue 5, October 2007, Pages 701–710, https://doi.org/10.1080/10635150701607033
  # - Diversitree R documentation on make.bisse function.,
  #     https://www.rdocumentation.org/packages/diversitree/versions/0.9-11/topics/make.bisse

# First load packages we'll need
library(ape)
library(diversitree)
library(deSolve)

# Set up simulated example to act as default values. Set seed to make this repeatable
set.seed(4)
exampleRates <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01) # lambda0, lambda1, mu0, mu1, q01, q10
exampleTree <- tree.bisse(exampleRates, max.t=30, x0=0)
exampleTipStates = exampleTree$tip.state

# The BiSSE Likelihood function. 
# Args:
#     - tree = an APE Phylo object, assumed to be bifurcating and ultrametric
#     - rates = a vector of length 6 with these values in order: lambda0, lambda1, mu0, mu1, q01, q10
#     - tipStates = a vector of length Ntip(tree) with a binary state, 1 or 0, for each tip. Each state has the tip label as a name
#
# Returns:
#     - An estimation of the likelihood for these data under this model. 
BiSSELikelihood <- function(rates = exampleRates, tree = exampleTree, tipStates = exampleTipStates) {
  names(rates) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10") # names the rates for ease of referencing
  
  # Set up arrays to hold the DN0(t) and DN1(t) probabilities calculated at each internal node. 
  # These represent the probability of the clade observed downstream of the node being created,  
  # given the time, parameters, and starting in state 0 or 1 at that node. 
  dN0values = rep(NA, Nnode(tree) + Ntip(tree))
  dN1values = rep(NA, Nnode(tree) + Ntip(tree))
  
  # Set up parallel arrays for calculated extinction rates.
  # These represent the chance at each node of the lineage going extinct before present day,
  # given state 0 or 1 at the node
  E0values = rep(NA, Nnode(tree) + Ntip(tree))
  E1values = rep(NA, Nnode(tree) + Ntip(tree))
  
  # Keep track of height, which we need for the integration over the branch
  nodeHeights = rep(0.0, Nnode(tree) + Ntip(tree))
  
  # NB: I know the previous 5 vectors could be a single matrix, but this
  # is easier for me to conceptualise at this point
  
  # We now integrate along each branch, storing the probabilities of 0 and 1 at
  # each internal node, and using those as the precondition for the branch.
  # We're climbing backwards from the tips to the root by iterating in postorder over the edges.
  for (edgeIndex in postorder(tree)) {
    rootWards = tree$edge[edgeIndex, 1] # index of rootWards node
    tipWards = tree$edge[edgeIndex, 2] # index of tipWards node
    edgeLength = tree$edge.length[edgeIndex]
    
    tipWardsHeight = 0
    rootWardsHeight = 0
    preProb0 = 0 # probability of 0 at the start of this branch
    preProb1 = 0 # probability of 1 at start of this branch
    E0 = 0 # prob of extinction for the lineage at the tipWards node, if in state 0
    E1 = 0 # prob of extinction for the lineage at the tipWards node, if in state 1
    
    if (tipWards <= Ntip(tree)) {
      # if this is a tip, preProbs are just whether the tip is in that state or not
      tipState = tipStates[tree$tip.label[tipWards]]
      preProb0 = ifelse(tipState == 0, 1, 0)
      preProb1 = ifelse(tipState == 1, 1, 0)
      # NB: E0 = E1 = 0 if this is a tip
      
      # Store the height for rootWards. Tip is at 0. Overwriting these values from other
      # branches shouldn't matter as the tree is assumed to be ultrametric. 
      rootWardsHeight = edgeLength
      nodeHeights[rootWards] = rootWardsHeight
    } else {
      # The indexed values should have been calculated already, based on postorder traversal
      # The preconditions for an internal branch are the values at its tipWards node. 
      preProb0 = dN0values[tipWards]
      preProb1 = dN1values[tipWards]
      E0 = E0values[tipWards]
      E1 = E1values[tipWards]
      
      # Calculate and record height
      tipWardsHeight = nodeHeights[tipWards]
      rootWardsHeight = tipWardsHeight + edgeLength
      nodeHeights[rootWards] = rootWardsHeight
    }

    # Integrate over the branch, calculating dN0(t) and dN1(t) at rootWard end of the branch.
    # This involves solving 4 parallel differential equations, laid out in diffEquations below, 
    # for each step defined in branchSlices. These steps move from the tipWards height to the
    # rootWards height. I'm naively using 200 steps on each branch, though I could optimise this
    # based on branch lengths. 
    # I define the integration method "lsoda" so the same one can be used by our external BiSSE function for comparison.
    branchSlices = seq(from = tipWardsHeight, to = rootWardsHeight, length.out = 200)
    integrationOutput = ode(y = c(preProb0, preProb1, E0, E1), times = branchSlices, 
                            func = diffEquations, parms = rates, method = "lsoda")
    
    # Grab the values of dNO(t) and dN1(t) from the last slice, at the end of the branch. 
    # I use length(integrationOutput)/5 because length is the total number of cells, not rows, and each row has 5 values. 
    valuesAtEndOfBranch = integrationOutput[length(integrationOutput)/5, ]
    
    # Save these numbers, or combine them with other incoming branches if we can.
    if (is.na(dN0values[rootWards])) {
      # This is the first of two branches that need to combine into this node. 
      # Just save the current values for now.
      dN0values[rootWards] <- valuesAtEndOfBranch[2]
      dN1values[rootWards] <- valuesAtEndOfBranch[3]
      E0values[rootWards] <- valuesAtEndOfBranch[4]
      E1values[rootWards] <- valuesAtEndOfBranch[5]
    } else {
      # This is the second branch, therefore we combine the values at the node as in Maddison et al (2007)
      dN0values[rootWards] <- dN0values[rootWards] * valuesAtEndOfBranch[2] * rates[1] # rates[1] is lambda0
      dN1values[rootWards] <- dN1values[rootWards] * valuesAtEndOfBranch[3] * rates[2] # rates[2] is lambda1
      E0values[rootWards] <- (E0values[rootWards] + valuesAtEndOfBranch[4]) / 2 # average the extinction rates when combining at a node
      E1values[rootWards] <- (E1values[rootWards] + valuesAtEndOfBranch[5]) / 2
    }
  } 
  
  q01 = rates["q01"]
  q10 = rates["q10"]
  
  rootIndex = Ntip(tree) + 1 # root is always at Ntips + 1 in R. 
  dN0 = dN0values[rootIndex]
  dN1 = dN1values[rootIndex]
  
  # calculate likelihood at the root, combining 0 and 1 states based on the equilibrium frequency of their transition rates
  lik = (q01 / (q01 + q10)) * dN1 + (q10 / (q01 + q10)) * dN0
  lnL = unname(log(lik))
  
  txt0 = paste("lambda0:", rates["lambda0"], "mu0:", rates["mu0"], "q01:", rates["q01"], sep="\t")
  txt1 = paste("lambda1:", rates["lambda1"], "mu1:", rates["mu1"], "q10:", rates["q10"], sep="\t")
  txtLnL = paste("\tlnL:", lnL)
  cat(paste(txt0, txt1, txtLnL, sep = "\n"), "\n\n")
  
  return(lnL)
}

# The 4 dependent differential equations for the values of dN0(t), dN1(t), extinction0 and extinction1.
# This function has the right arguments to be given to the ode solver in deSolve.
# Args:
#     - time = the timestep we are at in the integration. This is given by ODE and taken from branchSlices above
#     - y = the values that we calculate and that change at each timestep. 
#     - parms = constants needed at every time step in this calculation. For us it is the rates vector of lambda0, lambda1, mu0, mu1, q01, q10
#
# Return:
#     - A list of vectors, where the first element contains the calculated values that mirror the y input. 
#       This will become the y input for the next timestep.
diffEquations <- function(time, y, parms) {
  preProb0 = y[1]
  preProb1 = y[2]
  E0 = y[3]
  E1 = y[4]
  
  # See Maddison et al (2007) equation 3a
  dDN0.dt = -1 * (parms["lambda0"] + parms["mu0"] + parms["q01"]) * preProb0 +
            parms["q01"] * preProb1 + 2 * parms["lambda0"] * E0 * preProb0
  
  # See Maddison et al (2007) equation 3b
  dDN1.dt = -1 * (parms["lambda1"] + parms["mu1"] + parms["q10"]) * preProb1 +
            parms["q10"] * preProb0 + 2 * parms["lambda1"] * E1 * preProb1
  
  # See Maddison et al (2007) equation 7a
  dE0.dt = parms["mu0"] - (parms["lambda0"] + parms["mu0"] + parms["q01"]) * E0 +
           parms["q01"] * E1 + parms["lambda0"] * E0 * E0
  
  # See Maddison et al (2007) equation 7b
  dE1.dt = parms["mu1"] - (parms["lambda1"] + parms["mu1"] + parms["q10"]) * E1 +
           parms["q10"] * E0 + parms["lambda1"] * E1 * E1
  
  return(list(c(dDN0.dt, dDN1.dt, dE0.dt, dE1.dt)))
}

# Uncomment to run with example data on reload
# BiSSELikelihood()
