# Polytopic underapproximation of stochastic reach-avoid set

This directory contains code to recreate the examples from

  Abraham P. Vinod and Meeko M. K. Oishi, "Scalable Underapproximative Verification of Stochastic LTI Systems Using Convexity and Compactness", ACM Hybrid Systems: Computation and Control, 2018, 1--10.

# Requirements
    - MATLAB:
        - Statistics and Machine Learning (for qscmvnv.m)
	- Optimization toolbox (for fmincon)
        - Global optimization toolbox (for pattersearch) 
    - CVX:
        - Convex optimization problem parser
        - (Optional) Gurobi
    - MPT3:
        - Polyhedral computations and plotting

# Instructions

- Obtain Figure 5 and 6 and the relevant computation times (Table 2) using main_CWH.m

## Contact details

* Abraham P. Vinod ([aby.vinod@gmail.com](mailto:aby.vinod@gmail.com))
* Meeko M. K. Oishi ([oishi@unm.edu](mailto:oishi@unm.edu))
