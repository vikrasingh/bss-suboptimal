%-------------------------------------------------------------------------------- % Authors: Vikram Singh and Min Sun % % % % Copyright (2025): Vikram Singh and Min Sun % % This code is distributed under the terms of the GNU General Public License % % 2.0. % % % % Permission to use, copy, modify, and distribute this software for % % any purpose without fee is hereby granted, provided that this entire % % notice is included in all copies of any software which is or includes % % a copy or modification of this software and in all copies of the % % supporting documentation for such software. % % This software is being provided "as is" without any express or % % implied warranty. In particular, the authors do not make any % % representation or warranty of any kind concerning the merchantability % % of this software or its fitness for any particular purpose. % %-------------------------------------------------------------------------------- This package provides the testing and algorithm files to solve the following Best Subset Selection (BSS) problem

min_{b in R^p} ||y-Xb||{2}^{2} subject to ||b||{0}<=k , l<=b<=u

where l,u,b are in R^p, y in R^n, X in R^{nxp}, ||.||{2} is l-2 norm, ||.||{0} is pseudo-norm which counts the no. of non-zero entries of b. The following five suboptimal algorithms have been implemented
and compared to solve the BSS problem.
1. Discrete First Order
2. Forward Selection
3. Genetic Algorithm
4. Sequential Feature Swapping
5. Sequential Floating Feature Selection.
