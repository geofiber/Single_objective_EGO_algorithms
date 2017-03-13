# The pseudo EI criterion for parallel EGO algorithm
This is the MATLAB of the pseudo EI criterion for parallel EGO algorithm described in [1]. If you use this code for publication, please cite this reference.

The pseudo EI criterion is an extendtion of the standard EI criterion [2] which allows the EGO algorithm to select multiple design points in each iteration (cycle). Then these candidate points can be evaluated in parallel that may save some wall-clock time.

The dace toolbox [2] is used for building the Kriging model in the implementation.

Reference:

1. D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for parallel, Journal of Global Optimization. doi:10.1007/s10898-016-0484-7

2. Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of expensive black-box functions. Journal of Global Optimization 13(4), 455-492 (1998).

3. Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002. Available at: http://www2.imm.dtu.dk/~hbn/dace/.
