# puncturing-optimization
Source code for article
Svistunov German und  Vasiliy Usatyuk. 
Interleaved Cyclic Group Decomposition Length Adaptive MET QC-LDPC Codes.
IEEE International Black Sea Conference on Communications and Networking 3-6 June 2019. Sochi, Russia
 
It construct interleaver for improving perfomance under 5G eMBB lenght adaptation for short length using k-step recover puncture pattern optimization.

In reposit contained two methods:
1. Baseline, Syndrom based (one of trivial approach, exist more efficient MacKay syndrom entropy approach), find_puncture_pattern.m;
2. k-step recover puncture pattern optimization based on paper Jeongseok Ha, Jaehong Kim, D. Klinc and S. W. McLaughlin, "Rate-compatible punctured low-density parity-check codes with short block lengths," in IEEE Transactions on Information Theory, vol. 52, no. 2, pp. 728-738, Feb. 2006.


To undestand how it work run example.m
[P, Rates, Pairs]=GroupAndSort2(parity_check_matrix,rate_without_puncturing,step_in_rate)
P - Variable nodes to puncturing.
Rates - rate after punturing nodes in P.
Maximal chain of punturing depend from random initial nodes, run several times until get maximal rate.
