# The Generalized Mutually Best Pairs condition

 This repository provides Matlab code to verify the Generalized Mutually Best Pairs (GMBP) condition. The GMBP condition is useful to identify in the classical school choice problem whether the student-proposing deferred acceptance algorithm (DA) is efficient and whether there is a unique allocation that respects priorities. More information about the condition can be found in \[[Cantillon, Chen, Pereyra, 2022](https://arxiv.org/abs/2212.02881)\]. The code was tested on Matlab R2021b and should be compatible with older versions.

 ## File structure
 * UtilitiesPriorities.m
    + generate student preferences and school priorities with chosen parameters for the simulations
    + verify for the given market of student preferences and school priorities, whether the market satisfies GMBP, the sequential Mutually Best Pairs \[[Salonen and Salonen, 2018](https://www.sciencedirect.com/science/article/pii/S016548961730001X)\]
    + produce graphs visualizing the simulated results
* Algorithm
  + gmbp.m \[Code to verify GMBP\]
  + mbp.m  \[Code to verify MBP\]
  + dastudent.m  \[Code to generate matching produced by the student-proposing DA algorithm\]
  + TTC.m \[Code to generate the matching produced by the Top-trading cycle (TTC) mechanism\]
  + pareto.m \[Code to verify if the outcome produced by the student-proposing DA is efficient\]
