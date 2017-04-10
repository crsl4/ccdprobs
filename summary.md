# Notebook BISTRO
Bret Larget, Claudia Solis-Lemus (2016)

## To do now
- acabar de checar q el orden este bien en test
- checar las alphas y lambdas ahora, no deben haber cosas raras esta vez

- grid of branch lengths on likelihood when we stop, to try to understand what is happening, first we need to fix the NR problem


- Write up manuscript, and figure out simulation study for small datasets:
  - Rerun many datasets of increasing size with the current state of bistro (fixed tree, so run mrbayes first): Edit bistroOneRep and bistroAllRep: run with fixed topology and without
  - Create scripts to analyze output files: ESS, correct bl and p in posterior interval, MAP tree = true tree (or average PP for the true tree); do plots

## Later:
- how often do we have the bad alpha cases in the bl?
- add prints to see if the weird big alphas are linked to L(0,0): do we want a minimum variance? weird negative var! add checks for the vc matrix: tiny v, and v<0; why alpha==1; add a check for positive definite (determinant)=> what to do? make the correlation coefficient = 1 or -1 depending on the sign of s: sigma1*sigma2
r s
s t
- test with weightMean option
- do a rmd file to summarize result better, what is the problem with distance weights!? what is the problem with bl? (use Rstudio): check conclusions below. first rerun to see things are fixed now with alpha problem
*** 2 problems: long BL, true tree sampled only once in whales
- sim16 1500: 1.23%: long BL, some true BL outside of densities
- sim16 1500 indep: 0.65%
- sim16 fix tree: 0.51%
- sim16 500: 0.15%
- sim cats dogs 1500: 0.5%, long BL, but plots well centered, most sampled tree is different from true tree
- sim cats dogs fixed tree: 1.16%, long BL, but plots well centered
cats mean tree = true tree
- sim whales 1500: 0.14%: true tree only sampled once, it is not even the one with w=0.8: comparison plot, most sampled tree is different from true tree; mean tree tiny difference from true tree
- sim whales fixed tree: 0.1%, long BL, bl 54.53=20!!


## Current problems:
- individual MLE cambia del joint MLE:
branch lengths before joint MLE: 0.050790 0.028995
Initial gradient:  -602.892870 -1312.026804
Initial hessian:
-8979.288156   948.334247
  948.334247   862.319441
branch lengths after joint MLE:  0.050529 46.345215
Gradient: -0.000000 -0.000000
Hessian:
-14671.025746      0.000000
     0.000000      0.000000
But we still have weird hessians, which lead to weird cov (non positive matrices)
- When finding the individual MLE, the first derivative is 0, but when using the same individual MLE in the joint likelihood, then the gradient is far from zero, why?
When we optimized t1, we had t2 at something else!
Do we care?
Also, see these BL:
MLE BL: 0.005544
1st derivative: -0.000000
2nd derivative: -118928.778204
MLE BL: 0.046195
1st derivative: 0.000000
2nd derivative: -28237.935574
branch lengths before joint MLE: 0.005544 0.046195
Initial gradient: 257.686841 900.635353
Initial hessian:
-152586.908387   -1237.719476
  -1237.719476  -48395.345476
branch lengths after joint MLE: 0.007418 0.077757
Gradient: 0.000000 0.000000
Hessian:
-103429.536925    -363.135084
   -363.135084  -16626.654090
   -------------------------
   MLE BL: 0.025094
   1st derivative: -0.000000
   2nd derivative: -50764.555167
   MLE BL: 0.032610
   1st derivative: -0.000000
   2nd derivative: -39232.473999
branch lengths before joint MLE: 0.025094 0.032610
Initial gradient: -1036.520865 -1394.911543
Initial hessian:
-4077.546622  1122.497826
 1122.497826   999.863442
branch lengths after joint MLE:  0.024296 48.201786
Gradient: -0.000000 -0.000000
Hessian:
-23281.606224      0.000000
     0.000000      0.000000

## Check with Bret
- problem 1: individual MLE does not match joint MLE, do we care?
- problem 2: if we use the individual MLE, we get hessians with very big numbers => cov with very small numbers, and unstable => positive definite problems
- problem 3: if we use the joint MLE, the NR is a nightmare, and we don't fix the hessian problems => hessian close to zero, flat lik, problems!
- problem 4: if we use the individual MLE, we are not in the joint MLE, so the fisher info does not behave well:
branch lengths before joint MLE: 0.050160 0.028767
Gradient:  -609.960221 -1331.223637
Hessian:
-9282.820332   979.369411
  979.369411   886.695211




## Jordan
- Create bootstrap sample of trees for different datasets:
  - Compare bootstrap frequencies of trees and clades to MrBayes PP (with plots):
    - artiodactyl
    - cats dogs
    - primates
    - whales
    - 024,027 (or more) in datasets
  - Compare Bret's mean tree to Megan's mean tree, as well as with most likely tree in MrBayes
  - Later: weight trees by distance to "center tree", and compare weights to MrBayes PP
  - Compare performance of bistro with different scales

## Known problems
- Problem with long branches: only use joint BL density if there is a long branch (ESS is better with independent BL if there is not a long BL)
- Biased and spread proposal densities for Q in some causes
- Bootstrap sample of trees does not provide a good estimate of clade distribution:
  - Idea: Use a "center tree" and distance to this tree as weights

## Theoretical things for the future
- Prove crude SE is close (smaller than) real SE with the formula
- Study correlation matrices: For one tree, keep gamma parameters, and get a correlation matrix used to sample and compare to the one in mrbayes and the sample correlation from bistro. Keep in mind: 3 matrices: the one from the gamma that we sampled from (the observed inverse info matrix), the sample correlation from bistro and the sample correlation from mrbayes (which we are using as theoretical one, but we don’t know). One idea is to shrink the covariance matrix from the observed information matrix, in a controlled way towards a identity: divide the off diagonals by another scale. Add a column that prints the root, and save the gamma parameters to study the actual correlation that we use (I think we don’t need this anymore because we are satisfied with the idea that long branch causes correlation to neighbor branches, so we only need to study this correlation)
- **Mean trees** Why minimizing geodesic to get the frechet mean instead of bret algorithm: `argmin_a \sum_t \sum_s (t(s)-a(s))^2` (s=split, t=list of input trees)?
Bret algorithm: give a score to each split (see below), sort the splits, and
input split into tree if they are compatible. The score of a split is the square of the parent edge of the split, and the best score for the subtree. **Question:** How are the Frechet mean and the Bret mean different?
- Sequential IS?


## Performance improvements
- Test if first branches to estimate are poor in comparison to deep branches: to test without randomize, we need to call sort canonical before sampling branch lengths
- Add proposal to mcmc to loop over internal nodes, and slide the two neighboring edges
- Use Paul Lewis measure of information from sample as well as ESS

## Code improvements

### Minor fixes
- Use max number of cores if specified bigger number? Give warning!
- fatal error: boost/dynamic_bitset.hpp: No such file or directory; sometimes Boost, sometimes boost
- Alignment summary print number of taxa and sequence length
- Some functions need the root to have two children, some need to have three children, we need to change the functions so that they all work either way
- Get rid of old unused functions in bistro
- change back --no-reweight to use raw bootstrap counts



### Big fixes
- Bootstrap multithreads
- Add a bool to Edge:: calculate to only do logl, instead of dlogl and ddlogl always
- Include summarizeBL.jl and summarizeOut.r in bistro to get understandable output files for BL and rates
- Having an option to read a list of trees to use as input for ccdprobs, instead of reading the sequences: get a sample of fastME


### Known errors
- Bret error: fix error here: A new large data set took a long time to get through the MCMC. It then did the bootstrap. When it stopped, it printed out to the .topCounts file, but there were no end lines. It then froze without actually doing any importance sampling. Something is wrong at the stage of going from bootstrap to importance sampling with this data set.


### Threads problem
- Get rid of logwt to get ESS from multimap
- Try thread guard from book; try page 28 as well; google mem_fn wrapper
- Number of cores -1 by convention?
- New valgrind error

## Notes:
ln -s Documents/phylo/projects/present/CladeCondProb/ccdprobs/Bistro/Code/bistro/bistro /usr/local/bin/bistro

-----------------------
# URS project
Jordan Vonderwell (Fall 2016-Spring 2017)

## Objectives:
Presentations at Undergraduate symposium (end of spring), final research report (end of spring)

## Project global goals:
- Find out if parsimony weights with bootstrap really cover the posterior of trees.
- Determine the parsimony scale that is better, we wish to make this data-specific (maybe depending on the range of parsimony score?).
- Parallelize boostrap?
- Crude SE match empirical standard errors; compare time with MrBayes

## Contract:
- tasks for which the student is responsible
- number of hours/week the student is expected to work: 6
- aspects of a research process (see list p. 2) the student will be learning
- number of credits – 2 or 3 – appropriate for the student’s time invested
- If a student is expected to work with you 4-6 hours per week, 2 credits is appropriate. For 7-12 hours a week, 3 credits are appropriate. If you find during the semester that the student is averaging more or fewer hours than expected, please contact the URS Director and the credit load will be adjusted accordingly.


## Background:
- Tree thinking book
- Phylogenetic Handbook
- Bayesian introduction: small code example of MCMC
- Importance sampling: small code example
- Explain big picture of our project
- Continuous time MC
- Bootstrap
- Parsimony
- Example with the coins
- Later: Github account to fork the repository (whenever ready to make changes to code)

### Project
- Describe code, how to compile, how to run, explain parameters and output
- Git clone the repository (cannot make changes to code at this point)
- Run a small example to explain input and output
- Run a small example in MrBayes to explain input and output
- Perl script to compare these two outputs? What do we try to compare?
- Design simulations/real-life dataset examples to compare MrBayes and bistro

### Written report
Goal: identify a proposal density for the posterior distribution of trees

- introduction of importance sampling
- known difficulties of importance sampling: bias (recall R tests). Basically, we need good proposals
- introduction of phylogenetic trees
  - what is a phylogenetic tree
  - how do we estimate it with bayesian inference (mcmc)
  - downside of mcmc (slow)
  - importance sampling for trees: we need a good proposal density for topologies that matches the posterior probabilites (PP) from mrbayes. Some ideas:
    - bootstrap sample of trees (frequency does not match PP, bias)
    - weight based on distance to mean tree (frequency matches PP?)


-------------------------
## Previously in BranchLengths/Code/test
- change exp case to only use mean, not var
- we do not do calculate(qmatrix) for the edges after setting them, I think
- study 3 taxa case with ugly bl, why ESS so low? compare mvnormal and gammabeta, and study NR, are we failing to reach the MLE?
- revisit the debug.txt case: MLE not reached comparing C++ and R, do we reach it now?
- after satsified with new code, create simulations2, and redo exerything from simulations1. simulations2: mvnormal, after bug fixed
- run for real life data, only after fixing bug
better do arty3, birds4, arty6, cats-dogs12


### later:
- why gammabeta is much worse than mvnormal for 4 taxa real life data? do we have this analysis in R?
- save seed per replicate that failed to rerun (maybe does not work)
- small example for bret: is vs mcmc,why independent sampler does work for us? centered correctly, not enough


### much later:
- different types of exceptions to distinguish
