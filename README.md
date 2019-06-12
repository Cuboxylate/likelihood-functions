# Likelihood function implementations for BIOINF702 Assignments 3 and 4.

Find my implementations of the BiSSE Log-likeihood and Birth-Death Log-likelihood functions under BiSSELikelihood.R and birthDeathLikelihood.R.

To compare my outputs with the published functions in Densitree, run compareWithExistingFunctions.R. NOTE: this script will need you to set a variable - the path on your computer to this repository. 

To create plots using real data, run generateLikelihoodPlots.R. This will create two pdfs in the working directory. This script assumes BiSSELikelihood.R and birthDeathLikelihood.R are loaded, and the working directory is already set to this repository.

I implemented some likelihood optimisation for a few nested models in maximumLikelihoodEstimation.R.

See Assignment3_discussion.pdf and Assignment4_discussion.pdf for written assignment answers, as well as examples of likelihood plots from generateLikelihoodPlots.R and ML outputs from maximumLikelihoodEstimation.R.
