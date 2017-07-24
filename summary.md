# Notebook BISTRO
Bret Larget, Claudia Solis-Lemus (2016)

## To do now


- check paper that does mcmc and IS (someone in cimat told me)
- read steve derivation: exponential families
- wait for steve derivation for a better way to compute alpha and lambda



## Check with Bret
- bistro and mb mean trees are very different (see for example 059. mb mean tree = map, but != bistro mean tree and map). Bistro more balanced, is there a bias towards balance with bootstrap? So, the problem is not with spread, but with center. Could consensus tree of bootstrap be better than mean tree?
- example for bistro manuscript
- manuscript for dirichlet in BA: many things in red missing: how to cite this unpublished work in the main bistro paper? or do we want to start with the dirichlet one?

## To knitr Rmd files:
```R
library(knitr)
pandoc("input.Rmd", format = "latex")
```


## Jordan
### Final steps:
1. Clean up and add comments to the changes you made to files in `Bistro/Code/bistro`. Also, put the boost folder and everything that you needed to compile bistro in Ubuntu.
2. Put the scripts that you created to parse the output in the `Bistro/Scripts` folder with a README file on how to use them
3. Fork the ccdprobs repository on Github (you need a Github account for this), and replace your local version of ccdprobs with the cloned version from your repository:
`git clone https://github.com/YOUR_USERNAME/ccdprobs.git`
Make sure to backup your Bistro/Code/bistro folder, and all the files you added to `Bistro/Scripts`
4. Replace your files in the ccdprobs file, and do:
```
git add .
git commit -m "message explaining your changes"
git push
```

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
- BL sampling looks good for now, but we could still improve it by studying how many cases fall into the half normal or alpha==1 cases. Specially, we need to study the alpha==1 cases that are joint gamma, but num small, so we are forcing alpha to be >1
- weights by distance to mean tree seem to work fine, except in the case of very short branches (but this is the same problem with all bayesian approaches)

## Theoretical things for the future
- Prove crude SE is close (smaller than) real SE with the formula
- Study correlation matrices: For one tree, keep gamma parameters, and get a correlation matrix used to sample and compare to the one in mrbayes and the sample correlation from bistro. Keep in mind: 3 matrices: the one from the gamma that we sampled from (the observed inverse info matrix), the sample correlation from bistro and the sample correlation from mrbayes (which we are using as theoretical one, but we don’t know). One idea is to shrink the covariance matrix from the observed information matrix, in a controlled way towards a identity: divide the off diagonals by another scale. Add a column that prints the root, and save the gamma parameters to study the actual correlation that we use (I think we don’t need this anymore because we are satisfied with the idea that long branch causes correlation to neighbor branches, so we only need to study this correlation)
- **Mean trees** Why minimizing geodesic to get the frechet mean instead of bret algorithm: `argmin_a \sum_t \sum_s (t(s)-a(s))^2` (s=split, t=list of input trees)?
Bret algorithm: give a score to each split (see below), sort the splits, and
input split into tree if they are compatible. The score of a split is the square of the parent edge of the split, and the best score for the subtree. **Question:** How are the Frechet mean and the Bret mean different?
- Sequential IS?
- distribution on the norm of multivariate normal: iid case, with covariance? plot distance vs logl in multivariate normal case; and compare to what we get from the bootstrap trees: compare the percentiles of bootstrap with the ones with multivariate normal
  - sample of bootstrap trees->mean tree-> distances from mean tree-> plot logl vs distances: difference between 95 percentile logl and the logl of the mean
  - sample of multivariate normal (not independent)->distances from mean->plot logl vs distances: difference between 95 percentile logl to true mean logl
  - problem: which covariance in the multinormal? maybe the difference in logl does not depend on the standard deviation
  - we do this to calculate the scale in the weights by distance
  - we can get quantiles of logl-max logl of mrbayes runs, and compare to the quantiles we get from normal (or bivariate normal): `dnorm(x)/dnorm(0)`: compare to different datasets
  - we want to justify the use of the scale in weight



## Performance improvements
- Test if first branches to estimate are poor in comparison to deep branches: to test without randomize, we need to call sort canonical before sampling branch lengths
- Add proposal to mcmc to loop over internal nodes, and slide the two neighboring edges
- Use Paul Lewis measure of information from sample as well as ESS
- If MLE is close to zero, then pick three points (0.001,0.002,0.003) and get a parabola. Then extrapolate the parabola to the identify the range in which the maxlogl(0) and maxlogl(0)-3 live. This is a range in x. Take 5 points in this range, evaluate the logl again, and compute the mean to then estimate alpha for the pareto: this is instead of what we are doing for the truncated normal

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
