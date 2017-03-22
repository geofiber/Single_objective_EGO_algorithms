# The standard and parallel efficient global optimization (EGO) algorithms

The standard efficient global optimization (EGO) [1] is a widely used surrogate-based optimization algorithm for expensive single-objective optimiztion. The EGO algorithm starts by building an initial Kriging model (aka. Gaussion process model) using some initial design points which are often produced by a experiment design method, such as Latin Hypercube Sampling (LHS) method. Then, in each iteration, the point with the highest expected improvement (EI) value is selected  by using a traditional optimization algorithm, such as genetic algorithm (GA). The selected point is evaluated using the real expensive objective function and used to update the Kriging model. In such a way, the EI criterion guides the search toward the optimum of the real problem.

The pseudo EI criterion is an extendtion of the standard EI criterion [2] which allows the EGO algorithm to select multiple design points in each iteration (cycle). Then these candidate points can be evaluated in parallel that may save some wall-clock time.

The dace toolbox [3] is used for building the Kriging models in the implementations.

Reference:

1. Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of expensive black-box functions. Journal of Global Optimization 13(4), 455-492 (1998).

2. D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for parallel, Journal of Global Optimization. doi:10.1007/s10898-016-0484-7

3. Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002. Available at: http://www2.imm.dtu.dk/~hbn/dace/.
