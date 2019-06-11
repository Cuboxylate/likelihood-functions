# Set up simulated example to act as default values. Set seed to make this repeatable
set.seed(5)
bisseParams <- c(0.25, 0.2, 0.1, 0.1, 0.01, 0.01) # lambda0, lambda1, mu0, mu1, q01, q10
bisseTree = tree.bisse(pars = bisseParams, max.taxa = 50)
bisseTipStates = bisseTree$tip.state

# And now some Birth-Death data
set.seed(1)
bdRates = c(0.2, 0.05)
bdTree = tree.bd(pars = exampleRates, max.taxa = 50)