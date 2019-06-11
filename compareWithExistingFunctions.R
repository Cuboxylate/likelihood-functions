# This script compares my BiSSE and Birth-Death likelihood functions to
# those available in the Diversitree package.
# Please see https://www.rdocumentation.org/packages/diversitree/versions/0.9-11

# TO RUN: 
#     - Please change assignmentLocation below to the location of my assignment folder
#     - If you don't have Diversitree installed, uncomment the next line

# install.packages("diversitree")
assignmentLocation = "~/Uni/BIOINF702/likelihood-functions/"

library("ape")
library("diversitree")

setwd(assignmentLocation)
source(paste0(assignmentLocation, "BiSSELikelihood.R"), echo=FALSE)
source(paste0(assignmentLocation, "birthDeathLikelihood.R"), echo=FALSE)
source(paste0(assignmentLocation, "createExampleData.R"), echo=FALSE)

## BiSSE Likelihood

# Find the likelihood using my function
lnL1 = BiSSELikelihood(rates = bisseParams, tree = bisseTree, tipStates = bisseTipStates)

# Find the likelihood using Diversitree
bisseLik2 <- make.bisse(bisseTree, bisseTipStates, control = list("backend" = "deSolve")) # construct the likelihood functin
lnL2 = bisseLik2(bisseParams, root = ROOT.EQUI) # combine at the root by equilibrium frequencies

## Birth-Death Likelihood

# Find likelihood using my function
lnL3 = birthDeathLikelihood(rates = bdRates, tree = bdTree)

# Find likelihood using Diversitree function
bdLik2 <- make.bd(bdTree)
lnL4 = bdLik2(bdRates)

text1 = paste("My BiSSE lnL:", lnL1, sep = "\t\t")
text2 = paste("Diversitree BiSSE lnL:", lnL2, sep = "\t")
text3 = paste("My BD lnL:", lnL3, sep = "\t\t")
text4 = paste("Diversitree BD lnL:", lnL4, sep = "\t")
cat(paste(text1, text2, text3, text4, sep = "\n"), "\n")